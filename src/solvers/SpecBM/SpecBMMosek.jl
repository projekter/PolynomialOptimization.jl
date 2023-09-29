struct SpecBMSubsolverMosek
    task::Mosek.MSKtask
    sparsemats::Dict{Int,Vector{Int64}}
    M₁afeidx::Vector{Int64}
    M₁val::Vector{Float64}
    M₂barvaridx::Vector{Int32}
    M₂numterm::Vector{Int64}
    M₂ptrterm::Vector{Int64}
    M₂termidx::Vector{Int64}
end

function specbm_setup_primal_subsolver(::Val{:Mosek}, num_psds, r, rdims, Σr, ρ)
    task = Mosek.maketask()
    taskptr = task.task
    # Note that here, we use the macros for direct access. Often, this is more efficient as no index conversion is
    # necessary, so we do it everywhere for consistency. However, they are not part of the documented API, so we'll hope
    # that the interface does not change.
    # first num_psds vars: γⱼ; following vars: triangles corresponding to vec(Sⱼ); last variable: objective
    Mosek.@MSK_appendvars(taskptr, num_psds +1)
    Mosek.@MSK_putobjsense(taskptr, Mosek.MSK_OBJECTIVE_SENSE_MAXIMIZE.value)
    Mosek.@MSK_putcj(taskptr, num_psds, -1.)
    # γⱼ ≥ 0
    Mosek.@MSK_putvarboundsliceconst(taskptr, zero(Int32), num_psds, Mosek.MSK_BK_LO.value, 0., Inf)
    # objective variable: no bounds
    Mosek.@MSK_putvarbound(taskptr, num_psds, Mosek.MSK_BK_FR.value, -Inf, Inf)
    # Sⱼ ⪰ 0
    Mosek.@MSK_appendbarvars(taskptr, length(r), convert(Vector{Int32}, r))
    # γⱼ + tr(Sⱼ) ≤ ρ
    Mosek.@MSK_appendcons(taskptr, num_psds)
    Mosek.@MSK_putconboundsliceconst(taskptr, zero(Int32), num_psds, Mosek.MSK_BK_UP.value, -Inf, ρ)
    ur = unique(r)
    # with the same breath, also create sparse symmats for extracting svec(S)
    sparsemats = sizehint!(Dict{Int,Vector{Int64}}(), length(ur))
    cfz = Iterators.countfrom(zero(Int32))
    let dims=Vector{Int32}(undef, maximum(rdims)), range=collect(Int32(0):Int32(max(num_psds, length(dims)) -1)),
        values=ones(Float64, length(range)), traces=sizehint!(Dict{Int,Int}(), length(ur)),
        rows=range, columns=similar(dims), nzs=ones(Int64, length(dims))
        Mosek.@MSK_putaijlist64(taskptr, num_psds, range, range, values)
        for (j, rⱼ) in zip(cfz, r)
            sparseidx = get!(() -> let idx=Ref{Int64}()
                Mosek.@MSK_appendsparsesymmat(taskptr, rⱼ, rⱼ, range, range, values, idx)
                idx[]
            end, traces, rⱼ)
            Mosek.@MSK_putbaraij(taskptr, j, j, 1, Ref(sparseidx), values)
        end
        # We loop once more, which allows us to overwrite range as "rows"
        @inbounds for rⱼ in ur
            dimⱼ = packedsize(rⱼ)
            k = 1
            for col in 0:rⱼ-1
                rows[k] = col
                columns[k] = col
                values[k] = 1.
                k += 1
                for row in col+1:rⱼ-1
                    rows[k] = row
                    columns[k] = col
                    values[k] = inv(sqrt2) # triangle is doubled, we don't want this
                    k += 1
                end
            end
            @assert(k -1 == dimⱼ)
            output = Vector{Int64}(undef, dimⱼ)
            fill!(dims, Int32(rⱼ))
            Mosek.@MSK_appendsparsesymmatlist(taskptr, dimⱼ, dims, nzs, rows, columns, values, output)
            sparsemats[dimⱼ] = output
        end
    end
    # quadratic objective constraint: F x + ⟨F̄, X̄⟩ + g ∈ D
    qc = Ref{Int64}()
    conedim = num_psds + Σr +2
    Mosek.@MSK_appendrquadraticconedomain(taskptr, conedim, qc)
    Mosek.@MSK_appendafes(taskptr, conedim)
    Mosek.@MSK_putafefentry(taskptr, 0, num_psds, 1.) # put the objective in the quadratic cone already
    Mosek.@MSK_putafeg(taskptr, 1, .5) # we need to minimize the actual square, so we use the rotated quadratic cone, but we
                                       # then have to fix the second variable to 1/2.
    Mosek.@MSK_appendaccseq(taskptr, qc[], conedim, 0, C_NULL)
    # the rest is set dynamically, but we also need further caches. The aim is to do one Cholesky factorization of the
    # matrix M and then just assign this data portion verbatim. Note that we only have direct access to the lower
    # triangular part, i.e., M = L * Lᵀ, so here we must put Lᵀ in the afes.
    # And we need to look at the right part of Lᵀ (in dense storage format), which is the bottom part in L. This is bad, as
    # it is not in a contiguous storage. But the rows of the right part of Lᵀ are indeed in contiguous storage, and Mosek
    # provides a routine just to set these rows.
    M₂barvaridx = collect(Int32(0):Int32(num_psds -1))
    M₂numterm = convert(Vector{Int64}, rdims)
    M₂ptrterm = isempty(rdims) ? Int64[] : accumulate(+, Iterators.flatten((zero(Int64),
                                                                            Iterators.take(rdims, length(rdims) -1))))
    M₂termidx = collect(Iterators.flatten(sparsemats[dimⱼ] for dimⱼ in rdims))
    # We put the γ data into the afes columnwise
    M₁afeidx = collect(Int64(2):Int64(num_psds + Σr +1))
    M₁val = similar(M₁afeidx, Float64)
    return SpecBMSubsolverMosek(task, sparsemats, M₁afeidx, M₁val, M₂barvaridx, M₂numterm, M₂ptrterm, M₂termidx)
end

specbm_finalize_primal_subsolver!(data::SpecBMSubsolverMosek) = Mosek.deletetask(data.task)

#= old - assumes that we don't use barvars, but variables in the scalar psd cone
let indices=Vector{Int32}(undef, maximum(rdims)), values=ones(Float64, length(indices) +1),
    indices2=similar(indices), values2=similar(values), i=num_psds, sqrt2=sqrt(2) # zero-based
    @inbounds for (j, (rⱼ, dimⱼ)) in enumerate(zip(r, rdims))
        istart = i
        # first encode γⱼ + tr(Sⱼ) ≤ ρ, we need rⱼ +1 entries
        indices[1] = j -1 # γⱼ
        for col in 2:rⱼ+1
            indices[col] = i
            i += col
        end
        # we directly call the macro - there's no bounds check involved, so we don't even need views.
        Mosek.@MSK_putarow(task.task, Int32(j -1), Int32(rⱼ +1), indices, values)

        # next encode Sⱼ ∈ 𝕊ʳ₊, we need rⱼ(rⱼ +1)÷2 entries. We have the vectorized upper triangle; MOSEK expects the
        # scaled and vectorized lower triangle.
        indices[1:dimⱼ] = 0:dimⱼ-1
        k = 1
        firstindex = 0
        firstindexΔ = 1
        for col in 1:rⱼ
            indices2[k] = col
            values2[k] = 1.
            k += 1

            curindexΔ = col
            curindex = firstindex + curindexΔ
            for row in col+1:rⱼ
                indices2[k] = curindex # col + row * (row -1) ÷ 2 -1
                values2[k] = sqrt2
                k += 1

                curindexΔ += 1
                curindex += curindexΔ
            end

            firstindexΔ += 1
            firstindex += firstindexΔ
        end
        Mosek.@MSK_putafefentrylist(task, dimⱼ, indices, indices2, values)
    end
end=#

function mosekdump(task, fn)
    Mosek.writedata(task, "temptask.task")
    try
        newtask = Mosek.maketask() do t
            Mosek.readdata(t, "temptask.task")
            Mosek.writedata(t, fn)
        end
    finally
        rm("temptask.task")
    end
end

function specbm_primal_subsolve!(mastersolver::SpecBMMastersolverData{R}, cache::SpecBMCache{R,F,ACV,SpecBMSubsolverMosek}) where {R,F,ACV}
    # Now we have the matrix M and can in principle directly invoke Mosek using putqobj. However, this employs a sparse
    # Cholesky factorization for large matrices. In our case, the matrix M is dense and not very large, so we are better of
    # calculating the dense factorization by ourselves and then using the conic formulation. This also makes it easier to use
    # other solvers which have a similar syntax.
    Mfact = cholesky!(cache.M, RowMaximum(), tol=sqrt(eps(R)), check=false)
    data = cache.subsolver
    num_psds = length(cache.m₁)
    cfz = Iterators.countfrom(zero(Int32))
    taskptr = data.task.task
    Mosek.@MSK_putcslice(taskptr, 0, length(cache.m₁), cache.m₁)
    for (i, m₂) in zip(cfz, cache.m₂s)
        Mosek.@MSK_putbarcj(taskptr, i, length(m₂), data.sparsemats[length(m₂)], m₂)
    end
    pLᵀ = transpose(@view(Mfact.L[invperm(Mfact.p), 1:Mfact.rank]))
    for (idx, col) in zip(zero(Int32):Int32(num_psds -1), eachcol(pLᵀ))
        copyto!(data.M₁val, col)
        Mosek.@MSK_putafefcol(
            taskptr,
            idx,
            Mfact.rank,
            data.M₁afeidx,
            data.M₁val
        )
    end
    let
        pL₂ᵀ = @view(pLᵀ[:, num_psds+1:end])
        buf = gettmp(cache, size(pL₂ᵀ, 2))
        # in the afes, our row index starts with 2 - as 0 is the index of the upper bound variable (-> objective) and 1 is the
        # index of the fixed 1/2
        for (i, row) in zip(Iterators.countfrom(Int64(2)), eachrow(pL₂ᵀ))
            copyto!(buf, row)
            Mosek.@MSK_putafebarfrow(
                taskptr,
                i,
                num_psds,
                data.M₂barvaridx,
                data.M₂numterm,
                data.M₂ptrterm,
                cache.Σr,
                data.M₂termidx,
                buf
            )
        end
    end
    for row in Int64(Mfact.rank +2):Int64(size(Mfact, 2) +1)
        Mosek.@MSK_emptyafebarfrow(
            taskptr,
            row
        )
    end
    # Bug in Mosek < 10.1.11 (that's why we disable the method for these versions): The changed data will not be processed
    # correctly in the solution. In principle, a very inefficient solution is to store the task in TASK format (others are also
    # broken prior to 10.1.11) and re-load it. We don't even offer this as a workaround, but hypothetically speaking, it works.
    Mosek.@MSK_writedata(taskptr, "task.jtask")
    #Mosek.@MSK_readdata(taskptr, "task.task")
    Mosek.optimize(data.task)
    solutionsta = Mosek.getsolsta(data.task, Mosek.MSK_SOL_ITR)
    solutionsta == Mosek.MSK_SOL_STA_OPTIMAL || error("Subsolver failed with status ", solutionsta)
    Mosek.@MSK_getxxslice(taskptr, Mosek.MSK_SOL_ITR.value, 0, length(mastersolver.γstars), mastersolver.γstars)
    for (j, Sⱼ) in zip(cfz, mastersolver.Sstar_psds)
        Sⱼunscaled = PackedMatrix(LinearAlgebra.checksquare(Sⱼ), vec(Sⱼ), :L)
        Mosek.@MSK_getbarxj(taskptr, Mosek.MSK_SOL_ITR.value, j, Sⱼunscaled)
        packed_scale!(Sⱼunscaled)
    end
    return
end