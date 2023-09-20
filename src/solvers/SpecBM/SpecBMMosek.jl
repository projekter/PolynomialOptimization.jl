struct SpecBMSubsolverMosek
    task::Mosek.MSKtask
    sparsemats::Dict{Int,Vector{Int64}}
    M₁afeidx::Vector{Int64}
    M₁val::Vector{Float64}
    M₂barvaridx::Vector{Int32}
    M₂numterm::Vector{Int64}
    M₂ptrterm::Vector{Int64}
    M₂termidx::Vector{Int64}
    lupermutations::Dict{Int,Vector{Int}}
end

function specbm_setup_primal_subsolver(::Val{:Mosek}, num_psd, r, rdims, Σr, ρ)
    task = Mosek.maketask()
    # Note that here, we use the macros for direct access. Often, this is more efficient as no index conversion is
    # necessary, so we do it everywhere for consistency. However, they are not part of the documented API, so we'll hope
    # that the interface does not change.
    # first num_psd vars: γⱼ; following vars: triangles corresponding to vec(Sⱼ); last variable: objective
    Mosek.@MSK_appendvars(task.task, num_psd +1)
    Mosek.@MSK_putobjsense(task.task, Mosek.MSK_OBJECTIVE_SENSE_MINIMIZE.value)
    Mosek.@MSK_putcj(task.task, num_psd, 1.)
    # γⱼ ≥ 0
    Mosek.@MSK_putvarboundsliceconst(task.task, zero(Int32), num_psd, Mosek.MSK_BK_LO.value, 0., Inf)
    # objective variable: no bounds
    Mosek.@MSK_putvarbound(task.task, num_psd, Mosek.MSK_BK_FR.value, -Inf, Inf)
    # Sⱼ ⪰ 0
    Mosek.@MSK_appendbarvars(task.task, length(r), convert(Vector{Int32}, r))
    # γⱼ + tr(Sⱼ) ≤ ρ
    Mosek.@MSK_appendcons(task.task, num_psd)
    Mosek.@MSK_putconboundsliceconst(task.task, zero(Int32), num_psd, Mosek.MSK_BK_UP.value, -Inf, ρ)
    ur = unique(r)
    # with the same breath, also create sparse symmats for extracting vec(S)
    sparsemats = sizehint!(Dict{Int,Vector{Int64}}(), length(ur))
    cfz = Iterators.countfrom(zero(Int32))
    lupermutations = sizehint!(Dict{Int,Vector{Int}}(), length(ur))
    let dims=Vector{Int32}(undef, maximum(rdims)), range=collect(Int32(0):Int32(max(num_psd -1, length(dims)))),
        values=ones(Float64, length(range)), traces=sizehint!(Dict{Int,Int}(), length(ur)),
        rows=range, columns=similar(dims), nzs=ones(Int64, length(dims))
        Mosek.@MSK_putaijlist64(task.task, num_psd, range, range, values)
        for (j, rⱼ) in zip(cfz, r)
            sparseidx = get!(() -> let idx=Ref{Int64}()
                Mosek.@MSK_appendsparsesymmat(task.task, rⱼ, rⱼ, range, range, values, idx)
                idx[]
            end, traces, rⱼ)
            Mosek.@MSK_putbaraij(task.task, j, j, 1, Ref(sparseidx), values)
        end
        # We loop once more, which allows us to overwrite range as "rows"
        @inbounds for rⱼ in ur
            dimⱼ = packedsize(rⱼ)
            lupermutation = Vector{Int}(undef, dimⱼ)
            k = 1
            for col in 0:rⱼ-1
                curidx = col +1
                Δ = rⱼ
                for row in 0:col-1
                    rows[k] = row
                    columns[k] = col
                    values[k] = .5 # triangle is doubled, we don't want this
                    lupermutation[k] = curidx
                    curidx += (Δ -= 1)
                    k += 1
                end
                rows[k] = col
                columns[k] = col
                values[k] = 1.
                lupermutation[k] = curidx
                k += 1
            end
            @assert(k -1 == dimⱼ)
            output = Vector{Int64}(undef, dimⱼ)
            fill!(dims, Int32(rⱼ))
            # We flip rows and cols, as MOSEK wants the lower triangle
            Mosek.@MSK_appendsparsesymmatlist(task.task, dimⱼ, dims, nzs, columns, rows, values, output)
            sparsemats[dimⱼ] = output
            lupermutations[dimⱼ] = lupermutation
        end
    end
    # quadratic objective constraint: F x + ⟨F̄, X̄⟩ + g ∈ D
    qc = Ref{Int64}()
    conedim = num_psd + Σr +2
    Mosek.@MSK_appendrquadraticconedomain(task.task, conedim, qc)
    Mosek.@MSK_appendafes(task.task, conedim)
    Mosek.@MSK_putafefentry(task.task, 0, num_psd, 1.) # put the objective in the quadratic cone already
    Mosek.@MSK_putafeg(task.task, 1, .5) # we need to minimize the actual square, so we use the rotated quadratic cone, but we
                                         # then have to fix the second variable to 1/2.
    Mosek.@MSK_appendaccseq(task.task, qc[], conedim, 0, C_NULL)
    # the rest is set dynamically, but we also need further caches. The aim is to do one Cholesky factorization of the
    # matrix M and then just assign this data portion verbatim. Note that we only have direct access to the lower
    # triangular part, i.e., M = L * Lᵀ, so here we must put Lᵀ in the afes.
    # And we need to look at the right part of Lᵀ (in dense storage format), which is the bottom part in L. This is bad, as
    # it is not in a contiguous storage. But the rows of the right part of Lᵀ are indeed in contiguous storage, and Mosek
    # provides a routine just to set these rows.
    M₂barvaridx = collect(Int32(0):Int32(num_psd -1))
    M₂numterm = convert(Vector{Int64}, rdims)
    M₂ptrterm = isempty(rdims) ? Int64[] : accumulate(+, Iterators.flatten((zero(Int64),
                                                                            Iterators.take(rdims, length(rdims) -1))))
    M₂termidx = vcat((sparsemats[dimⱼ] for dimⱼ in rdims)...)
    # We put the γ data into the afes columnwise
    M₁afeidx = collect(Int64(2):Int64(num_psd + Σr +1))
    M₁val = similar(M₁afeidx, Float64)
    return SpecBMSubsolverMosek(task, sparsemats, M₁afeidx, M₁val, M₂barvaridx, M₂numterm, M₂ptrterm,
        M₂termidx, lupermutations)
end

specbm_finalize_primal_subsolver(data::SpecBMSubsolverMosek) = Mosek.deletetask(data.task)

#= old - assumes that we don't use barvars, but variables in the scalar psd cone
let indices=Vector{Int32}(undef, maximum(rdims)), values=ones(Float64, length(indices) +1),
    indices2=similar(indices), values2=similar(values), i=num_psd, sqrt2=sqrt(2) # zero-based
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

function specbm_primal_subsolve(data::SpecBMSubsolverMosek, Mfact, m₁, m₂s, Σr, γstars, Sstars, tmp)
    num_psd = length(m₁)
    cfz = Iterators.countfrom(zero(Int32))
    taskptr = data.task.task
    Mosek.@MSK_putcslice(taskptr, 0, num_psd, m₁)
    for (i, m₂) in zip(cfz, m₂s)
        Mosek.@MSK_putbarcj(taskptr, i, length(m₂), data.sparsemats[length(m₂)], m₂)
    end
    pLᵀ = transpose(@view(Mfact.L[invperm(Mfact.p), 1:Mfact.rank]))
    for (idx, col) in zip(zero(Int32):Int32(num_psd -1), eachcol(pLᵀ))
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
        pL₂ᵀ = @view(pLᵀ[:, num_psd+1:end])
        buf = @view(tmp[1:size(pL₂ᵀ, 2)])
        # in the afes, our row index starts with 2 - as 0 is the index of the upper bound variable (-> objective) and 1 is the
        # index of the fixed 1/2
        for (i, row) in zip(Iterators.countfrom(Int64(2)), eachrow(pL₂ᵀ))
            copyto!(buf, row)
            Mosek.@MSK_putafebarfrow(
                taskptr,
                i,
                num_psd,
                data.M₂barvaridx,
                data.M₂numterm,
                data.M₂ptrterm,
                Σr,
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
    #Mosek.@MSK_writedata(taskptr, "task.task")
    #Mosek.@MSK_readdata(taskptr, "task.task")
    Mosek.optimize(data.task)
    solutionstap = Ref{Int32}()
    Mosek.@MSK_getsolsta(taskptr, Mosek.MSK_SOL_ITR.value, solutionstap)
    solutionsta = Mosek.Solsta(solutionstap[])
    solutionsta == Mosek.MSK_SOL_STA_OPTIMAL || error("Subsolver failed with status ", solutionsta)
    Mosek.@MSK_getxxslice(taskptr, Mosek.MSK_SOL_ITR.value, 0, length(γstars), γstars)
    for (j, Sⱼ) in zip(cfz, Sstars)
        Mosek.@MSK_getbarxj(taskptr, Mosek.MSK_SOL_ITR.value, j, tmp)
        @inbounds copyto!(Sⱼ, @view(tmp[data.lupermutations[length(Sⱼ)]]))
    end
    return
end