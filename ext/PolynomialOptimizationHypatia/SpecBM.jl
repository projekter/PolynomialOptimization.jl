struct SpecBMSubsolverHypatia{R}
    model::Models.Model{R}
    triupos::Vector{Int}
    ρpos::UnitRange{Int}
    solver::Solvers.Solver{R}
end

function SpecBM.setup_primal_subsolver(::Val{:Hypatia}, num_psds, r, rdims, Σr, ρ::R) where {R}
    @inbounds begin
        sidedim = num_psds + Σr
        # first num_psds vars: γⱼ; following vars: triangles corresponding to vec(Sⱼ); last variable: objective
        c = Vector{R}(undef, sidedim +1)
        c[end] = one(R)
        # construct the matrix G in the following way (h - G x ∈ [Nonneg for all γ, all PSDs, Nonneg for tr-ρ-bound, quad obj])
        # [-I(length(γ)) 0                               0
        #  0             -swap triu to tril for each PSD 0
        #  I(length(γ))  extract diagonals for each PSD  0
        #  0             0                               -1
        #  0             0                               0
        #  dummy-entries for Cholesky factor             0]
        # Note we need dense storage for the Cholesky stuff due to permutations!
        colptr = Vector{Int}(undef, sidedim +2)
        rowval = Vector{Int}(undef, 2num_psds + Σr + sum(r) + 1 + sidedim^2)
        nzvals = similar(rowval, R)
        colptr[end] = length(rowval) +1
        triustart = sidedim + num_psds + 3
        triupos = FastVec{Int}(buffer=sidedim)
        # first length(γ) columns are easy
        i = 1
        @simd for j in 1:num_psds
            colptr[j] = i
            rowval[i] = j
            nzvals[i] = -one(R)
            i += 1
            rowval[i] = j + sidedim
            nzvals[i] = one(R)
            i += 1
            unsafe_push!(triupos, i)
            for tri in triustart:triustart+sidedim-1
                rowval[i] = tri
                i += 1
            end
        end
        # rest
        j = num_psds +1
        i₁ = num_psds
        i₂ = sidedim +1
        for rⱼ in r
            # triu-tril swap. Columns first, so in the input order.
            rδ = 1
            for incol in 1:rⱼ
                i₁ += rδ
                rδ += 1
                cδ = incol
                outpos = i₁
                for inrow in incol:rⱼ
                    # swap
                    colptr[j] = i
                    rowval[i] = outpos
                    nzvals[i] = -one(R)
                    i += 1
                    outpos += cδ
                    cδ += 1
                    # trace
                    if incol == inrow
                        rowval[i] = i₂
                        nzvals[i] = one(R)
                        i += 1
                    end
                    # dummy-entry struff
                    unsafe_push!(triupos, i)
                    for tri in triustart:triustart+sidedim-1
                        rowval[i] = tri
                        i += 1
                    end
                    j += 1
                end
            end
            i₂ += 1
        end
        # final one for objective
        colptr[j] = i
        rowval[i] = sidedim + num_psds +1
        nzvals[i] = -one(R)

        G = SparseMatrixCSC(sidedim + num_psds + 2 + sidedim, sidedim +1, colptr, rowval, nzvals)
        h = zeros(R, size(G, 1))
        h[sidedim+1:sidedim+num_psds] .= ρ
        h[sidedim+num_psds+2] = inv(R(2))
        cones = Vector{Cones.Cone{R}}(undef, num_psds +3)
        cones[1] = Hypatia.Cones.Nonnegative{R}(num_psds)
        cones[2:end-2] = Hypatia.Cones.PosSemidefTri{R,R}.(rdims)
        cones[end-1] = Hypatia.Cones.Nonnegative{R}(num_psds)
        cones[end] = Hypatia.Cones.EpiPerSquare{R}(sidedim +2)
        return SpecBMSubsolverHypatia{R}(
            Models.Model{R}(
                c, Matrix{R}(undef, 0, length(c)), Vector{R}(undef, 0), G, h, cones
            ),
            finish!(triupos),
            sidedim+1:sidedim+num_psds,
            Solvers.Solver{R}(
                verbose=false,
                syssolver=Solvers.SymIndefSparseSystemSolver{Float64}(),
                preprocess=false,
                reduce=false
            )
        )
    end
end

SpecBM.adjust_penalty_subsolver!(data::SpecBMSubsolverHypatia, ρ) = @inbounds data.model.h[data.ρpos] .= ρ

SpecBM.finalize_primal_subsolver!(data::SpecBMSubsolverHypatia) = nothing

function SpecBM.primal_subsolve!(mastersolver::SpecBM.MastersolverData{R},
    cache::SpecBM.Cache{R,F,ACV,SpecBMSubsolverHypatia{R}}) where {R,F,ACV}
    sidedim = size(cache.M, 1)
    Mfact = cholesky!(cache.M, RowMaximum(), tol=sqrt(eps(R)), check=false)
    data = cache.subsolver
    num_psds = length(cache.m₁)
    model = data.model
    copyto!(model.c, cache.m₁)
    copyto!(model.c, num_psds +1, cache.m₂, 1, length(cache.m₂))
    rmul!(@view(model.c[1:end-1]), -one(R))
    pLᵀ = transpose(@view(Mfact.L[invperm(Mfact.p), 1:Mfact.rank]))
    vals = nonzeros(model.G)
    if sidedim > Mfact.rank
        @inbounds for (j, (pos, col)) in enumerate(zip(data.triupos, eachcol(pLᵀ)))
            copyto!(vals, pos, col, 1, Mfact.rank)
            fill!(@view(vals[pos+Mfact.rank:pos+sidedim-1]), zero(R))
        end
    else
        @inbounds for (j, (pos, col)) in enumerate(zip(data.triupos, eachcol(pLᵀ)))
            copyto!(vals, pos, col, 1, sidedim)
        end
    end
    # While in principle, we'd have to negate all the entries due to Hypatia's extra minus, we just ignore this due to the
    # squaring in the cone.
    solver = data.solver
    Solvers.load(solver, model)
    Solvers.solve(solver)
    Solvers.get_status(solver) ∈ (Hypatia.Solvers.Optimal, Hypatia.Solvers.NearOptimal) ||
        error("Subsolver failed with status ", Hypatia.Solvers.get_status(solver))
    copyto!(mastersolver.γstars, 1, solver.result.x, 1, num_psds)
    i = num_psds +1
    for Sⱼ in mastersolver.Sstar_psds
        copyto!(vec(Sⱼ), 1, solver.result.x, i, length(Sⱼ))
        i += length(Sⱼ)
    end
    return
end