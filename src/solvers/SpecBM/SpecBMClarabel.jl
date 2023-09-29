import Clarabel

struct SpecBMSubsolverClarabel{R}
    solver::Clarabel.Solver{R}
    P::SparseMatrixCSC{R,Int}
    q::Vector{R}
    A::SparseMatrixCSC{R,Int}
    b::Vector{R}
    cones::Vector{Clarabel.SupportedCone}
    settings::Clarabel.Settings
end

function specbm_setup_primal_subsolver(::Val{:Clarabel}, num_psds, r, rdims, Σr, ρ::R) where {R}
    solver = Clarabel.Solver()
    sidedim = num_psds + Σr
    P = SparseMatrixCSC(sidedim, sidedim, accumulate(+, Iterators.flatten((1, sidedim:-1:1))),
        collect(Iterators.flatten(j:sidedim for j in 1:sidedim)), zeros(R, packedsize(sidedim)))
    q = Vector{R}(undef, sidedim)
    # construct the matrix A in the following way (A x + s = b, s = [Nonneg for γ, all PSDs, Nonneg for tr-ρ-bound]):
    # [-I(length(γ)) 0
    #  0             -swap triu to tril for each PSD
    #  I(length(γ))  extract diagonals for each PSD]
    colptr = Vector{Int}(undef, sidedim +1)
    rowval = Vector{Int}(undef, 2num_psds + Σr + sum(r))
    nzvals = similar(rowval, R)
    colptr[end] = length(rowval) +1
    # first length(γ) columns are easy
    i = 1
    @inbounds @simd for j in 1:num_psds
        colptr[j] = i
        rowval[i] = j
        nzvals[i] = -one(R)
        i += 1
        rowval[i] = j + sidedim
        nzvals[i] = one(R)
        i += 1
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
                j += 1
            end
        end
        i₂ += 1
    end
    A = SparseMatrixCSC(sidedim + num_psds, sidedim, colptr, rowval, nzvals)
    b = zeros(R, size(A, 1))
    b[end-num_psds+1:end] .= ρ
    cones = Vector{Clarabel.SupportedCone}(undef, num_psds +2)
    @inbounds cones[1] = Clarabel.NonnegativeConeT(num_psds)
    @inbounds cones[2:end-1] = Clarabel.PSDTriangleConeT.(r)
    @inbounds cones[end] = Clarabel.NonnegativeConeT(num_psds)
    settings = Clarabel.Settings(verbose=false)
    return SpecBMSubsolverClarabel{R}(Clarabel.Solver{R}(), P, q, A, b, cones, settings)
end

specbm_finalize_primal_subsolver!(data::SpecBMSubsolverClarabel) = nothing

function specbm_primal_subsolve!(mastersolver::SpecBMMastersolverData{R}, cache::SpecBMCache{R,F,ACV,SpecBMSubsolverClarabel{R}}) where {R,F,ACV}
    data = cache.subsolver
    num_psds = length(cache.m₁)
    copyto!(data.q, cache.m₁)
    copyto!(data.q, num_psds +1, cache.m₂, 1, length(cache.m₂))
    # while Clarabel uses upper triangles, the matrix a was set up such that we can use the lower triangle
    trttp!('L', parent(cache.M), nonzeros(data.P))
    rmul!(data.q, -one(R))
    rmul!(data.P, R(2))
    rmul_offdiags!(PackedMatrix(size(data.P, 1), nonzeros(data.P), :L), R(2))
    Clarabel.setup!(data.solver, data.P, data.q, data.A, data.b, data.cones, data.settings)
    result = Clarabel.solve!(data.solver)
    result.status == Clarabel.SOLVED || error("Subsolver failed with status ", result.status)
    copyto!(mastersolver.γstars, 1, result.x, 1, num_psds)
    i = num_psds +1
    for Sⱼ in mastersolver.Sstar_psds
        copyto!(vec(Sⱼ), 1, result.x, i, length(Sⱼ))
        i += length(Sⱼ)
    end
    return
end