"""
    poly_solutions(:mvhankel, result::Result, ϵ=1e-6, δ=1e-3; verbose=false)

Performs a multivariate Hankel decomposition of the full moment matrix that was obtained via optimization of the problem, using
the [best currently known algorithm](https://doi.org/10.1016/j.laa.2017.04.015). This method is not deterministic due to a
random sampling, hence negligible deviations are expected from run to run.
The parameter `ϵ` controls the bound below which singular values are regarded as zero.
The parameter `δ` controls up to which threshold supposedly identical numerical values for the same variables from different
cliques must match.
Note that for sparsity patterns different from no sparsity and [`SparsityCorrelative`](@ref Relaxation.SparsityCorrelative),
this method not be able to deliver results, although it might give partial results for
[`SparsityCorrelativeTerm`](@ref Relaxation.SparsityCorrelativeTerm) if some of the cliques did not have a term sparsity
pattern. Consider using the `:heuristic` method in such a case.
This function returns an iterator.

See also [`poly_optimize`](@ref).
"""
poly_solutions(::Val{:mvhankel}, result::Result{Rx,V}, ϵ::R=R(1 // 1_000_000), δ::R=R(1 // 1_000);
    verbose::Bool=false) where {Nr,Nc,Rx<:AbstractRelaxation{<:Problem{<:SimplePolynomial{<:Any,Nr,Nc}}},R<:Real,V<:Union{R,Complex{R}}}  =
    poly_solutions(Val(:mvhankel), result.relaxation, result.moments, ϵ, δ, verbose)

function poly_solutions(::Val{:mvhankel}, relaxation::AbstractRelaxation{<:Problem{<:SimplePolynomial{<:Any,Nr,Nc}}},
    moments::MomentVector{R,V,Nr,Nc}, ϵ::R, δ::R, verbose::Bool) where {Nr,Nc,R<:Real,V<:Union{R,Complex{R}}}
    @verbose_info("Preprocessing for decomposition")
    nvars = Nr + 2Nc
    deg = 2degree(relaxation)
    @assert(deg > 1) # It is in principle possible to have deg = 0 - but only if the problem itself was to optimize zero over
                     # the whole space. Let's ignore this.
    # Potentially scale the moment matrix. Note that the vector of moments does not know about complex structure, so
    # monomial_count is the correct function.
    λ = maximum(absval, @view(moments[monomial_count(nvars, deg -2)+1:monomial_count(nvars, deg -1)])) /
        maximum(absval, @view(moments[monomial_count(nvars, deg -1)+1:end]))
    if λ > ϵ
        for d in 1:deg-1
            rmul!(@view(moments[monomial_count(nvars, d -1)+1:monomial_count(nvars, d)]), λ^d)
        end
        rmul!(@view(moments[monomial_count(nvars, deg -1)+1:end]), λ^deg) # the last might be incomplete
    end
    # for each variable clique, we can perform the original decomposition algorithm
    cliques = groupings(relaxation).var_cliques
    solutions_cl = FastVec{Union{Matrix{V},Missing}}(buffer=length(cliques))
    @verbose_info("Starting solution extraction per clique")
    extraction_time = @elapsed begin
        for (i, clique) in enumerate(cliques)
            @verbose_info("Investigating clique ", clique)
            a1 = basis(relaxation, i)
            a2 = Relaxation.truncate_basis(a1, degree(relaxation) -1)
            if !isreal(relaxation)
                # this is the transpose of what we'd get with moment_matrix, but this is not important. Since the monomials
                # in the matrix will be multiplied by variables (which are the un-conjugated ones), we must make sure that
                # the un-conjugated ones are of the smaller degree.
                a1 = conj(a1)
            end
            try
                result_cl = poly_solutions(Val(Symbol("mvhankel-scaled")), moments, a1, a2, clique, ϵ)
                λ > ϵ && rmul!(result_cl, inv(λ))
                unsafe_push!(solutions_cl, result_cl)
                @verbose_info("Potential solutions:\n", result_cl)
            catch e
                if e isa MonomialMissing
                    unsafe_push!(solutions_cl, missing)
                    @verbose_info("Not all monomials were present to allow for a solution extraction in this clique")
                else
                    rethrow()
                end
            end
        end
    end
    # undo the scaling
    if λ > ϵ
        λ = inv(λ)
        for d in 1:deg-1
            rmul!(@view(moments[monomial_count(nvars, d -1)+1:monomial_count(nvars, d)]), λ^d)
        end
        rmul!(@view(moments[monomial_count(nvars, deg -1)+1:end]), λ^deg) # the last might be incomplete
    end
    # now we combine all the solutions. In principle, this is Iterators.product, but we only look for compatible entries (as
    # variables in the cliques can overlap).
    @verbose_info("Found all potential individual solutions in ", extraction_time, " seconds. Building iterator.")
    return PolynomialSolutions{R,V,Nr,Nc,
                               SimpleVariable{Nr,Nc,SimplePolynomials.smallest_unsigned(Nr + 2Nc)},isone(length(cliques))}(
        cliques,
        finish!(solutions_cl),
        δ,
        verbose
    )
end

function poly_solutions(::Val{Symbol("mvhankel-scaled")}, moments::MomentVector{R,V}, a1, a2, variables, ϵ::R) where {R<:Real,V<:Union{R,Complex{R}}}
    Hd1d2 = moment_matrix(moments, Val(:throw), a1, a2)
    UrSrVrbar = svd!(Hd1d2)
    ϵ *= UrSrVrbar.S[1]
    numEntries = findfirst(<(ϵ), UrSrVrbar.S)
    if isnothing(numEntries)
        Ur, Sr, Vrbar = UrSrVrbar.U, UrSrVrbar.S, UrSrVrbar.V
    else
        numEntries -= 1
        @inbounds Ur, Sr, Vrbar = UrSrVrbar.U[:, 1:numEntries], UrSrVrbar.S[1:numEntries], UrSrVrbar.V[:, 1:numEntries]
    end
    Ms = Vector{Matrix{V}}(undef, length(variables))
    for i in 1:length(variables)
        @inbounds Ms[i] = (Ur' ./ Sr) * moment_matrix(moments, Val(:throw), a1, a2, variables[i]) * Vrbar
    end
    v = eigvecs(sum((2rand(V) - one(V)) * Mi for Mi in Ms))
    len = size(Sr, 1)
    result = Matrix{V}(undef, length(variables), len)
    for j in 1:len
        @inbounds vj = @view(v[:, j])
        for i in 1:length(variables)
            @inbounds Mv = Ms[i] * vj
            den = zero(R)
            vre = zero(R)
            if V <: Complex
                vim = zero(R)
            end
            for (Mvᵢ, vjᵢ) in zip(Mv, vj)
                vre += real(Mvᵢ) * real(vjᵢ) + imag(Mvᵢ) * imag(vjᵢ)
                if V <: Complex
                    vim += real(vjᵢ) * imag(Mvᵢ) - real(Mvᵢ) * imag(vjᵢ)
                end
                den += abs2(vjᵢ)
            end
            if V <: Complex
                @inbounds result[i, j] = V(vre, vim) / den
            else
                @inbounds result[i, j] = vre / den
            end
        end
    end
    # we delete duplicates (which might arise due to floating point inefficiencies)
    keep = fill(true, len)
    redo = false
    @inbounds for i in 1:len -1
        keep[i] || continue
        for j in i+1:len
            keep[j] || continue
            @views if norm(result[:, i] - result[:, j]) < ϵ
                keep[j] = false
                redo = true
            end
        end
    end
    if redo
        result = keepcol!(result, keep)
    end
    return result
end