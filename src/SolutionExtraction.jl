export moment_matrix, poly_solutions, poly_solution_badness, poly_all_solutions, optimality_certificate

struct MonomialMissing <: Exception end

"""
    moment_matrix(problem::Union{PolyOptProblem,AbstractSparsity}; max_deg=Inf, prefix=1)

After a problem has been optimized, this function assembles the associated moment matrix (possibly by imposing a degree bound
`max_deg`, and possibly multiplying each monomial by the monomial `prefix`, which does not add to `max_deg`).

See also [`poly_optimize`](@ref), [`sparse_optimize`](@ref).
"""
function moment_matrix(moments::Dict{M,R}, rows, cols, prefix, missing_func=(_) -> R(NaN)) where {M,R<:Real}
    if M <: MonomialComplexContainer
        return Complex{R}[
            let key = k * j * prefix, key_conj = conj(key), key_re, key_im;
                if key_conj < key
                    key_re = MonomialComplexContainer(key_conj, true)
                    key_im = MonomialComplexContainer(key_conj, false)
                    @inbounds Complex{R}(haskey(moments, key_re) ? moments[key_re] : missing_func(key_re),
                                         haskey(moments, key_im) ? -moments[key_im] : -missing_func(key_im))
                elseif key_conj == key
                    key_re = MonomialComplexContainer(key, true)
                    @inbounds Complex{R}(haskey(moments, key_re) ? moments[key_re] : missing_func(key_re), zero(R))
                else
                    key_re = MonomialComplexContainer(key, true)
                    key_im = MonomialComplexContainer(key, false)
                    @inbounds Complex{R}(haskey(moments, key_re) ? moments[key_re] : missing_func(key_re),
                                         haskey(moments, key_im) ? moments[key_im] : missing_func(key_im))
                end
            end for k in rows, j in cols]
    else
        return R[let key = k * j * prefix; haskey(moments, key) ? moments[key] : missing_func(key) end for k in rows, j in cols]
    end
end

function moment_matrix(problem::PolyOptProblem; max_deg=Inf, prefix=1)
    b = sort(problem.basis, by=degree)
    if max_deg < Inf
        upto = searchsortedfirst(b, max_deg+1, by=x -> x isa Int ? x : degree(x))
        if upto ≤ lastindex(b)
            deleteat!(b, upto:lastindex(b))
        end
    end
    if problem.complex
        return Hermitian(moment_matrix(last_moments(problem), b, conj.(b), prefix))
    else
        @assert iszero(imag(prefix))
        return Symmetric(moment_matrix(last_moments(problem), b, b, prefix))
    end
end
moment_matrix(sparse::AbstractSparsity; kwargs...) = moment_matrix(sparse_problem(sparse); kwargs...)

struct PolynomialSolutions{O,has_length}
    var_to_idx
    cliques
    solutions_cl
    δ
    verbose
end

"""
    poly_solutions(problem::Union{PolyOptProblem,AbstractSparsity}, ϵ=1e-6, δ=1e-3; verbose=false)

Performs a multivariate Hankel decomposition of the full moment matrix that was obtained via optimization of the problem, using
the [best currently known algorithm](https://doi.org/10.1016/j.laa.2017.04.015). This method is not deterministic due to a
random sampling, hence negligible deviations are expected from run to run.
The parameter `ϵ` controls the bound below which singular values are regarded as zero.
The parameter `δ` controls up to which threshold supposedly identical numerical values for the same variables from different
cliques must match.
Note that for sparsity patterns different from [`SparsityNone`](@ref) and [`SparsityCorrelative`](@ref), this method not be
able to deliver results, although it might give partial results for [`SparsityCorrelativeTerm`](@ref) if some of the cliques
did not have a term sparsity pattern. Consider using [`poly_solutions_heuristic`](@ref) in such a case.
This function returns an iterator.

See also [`poly_optimize`](@ref), [`sparse_optimize`](@ref), [`poly_solutions_heuristic`](@ref).
"""
function poly_solutions(state::AbstractSparsity, ϵ::R=1e-6, δ::R=1e-3; verbose::Bool=false) where {R<:Real}
    if !(state isa SparsityNone || state isa SparsityCorrelative)
        @warn("poly_solutions requires a certain group structure in the moment matrix. For term sparsity, this is usually not satisfied and the results are probably useless or not to be trusted. Consider using poly_solutions_heuristic instead.")
    end
    @verbose_info("Preprocessing for decomposition")
    problem = sparse_problem(state)
    moments = problem.last_moments
    isempty(moments) && error("The problem was not optimized yet.")
    deg = 2problem.degree
    λ = maximum(degree(m) == deg - 1 ? abs(v) : 0 for (m, v) in moments) /
        maximum(degree(m) == deg ? abs(v) : 0 for (m, v) in moments)
    moments_scaled = λ ≤ ϵ ? moments : typeof(moments)(m => λ^degree(m) * v for (m, v) in moments)
    # for each variable clique, we can perform the original decomposition algorithm
    _, cliques = sparse_groupings(state)
    solutions_cl = FastVec{Union{AbstractMatrix,Missing}}(buffer=length(cliques))
    if problem.gröbner_basis isa EmptyGröbnerBasis
        missing_func = (_) -> throw(MonomialMissing())
    else
        missing_func = (mon) -> let r = rem(mon, problem.gröbner_basis)
            isempty(r) ? zero(R) : λ^degree(mon) * sum(coefficient(t) *
                                                       get((_) -> throw(MonomialMissing()), moments, monomial(t)) for t in r)
        end
    end
    @verbose_info("Starting solution extraction per clique")
    extraction_time = @elapsed begin
        for variables in cliques
            @verbose_info("Investigating clique ", variables)
            a1 = problem.basis
            a2 = filter(mon -> degree(mon) < problem.degree, problem.basis)
            if problem.complex
                # this is the transpose of what we'd get with moment_matrix, but this is not important. Since the monomials
                # in the matrix will be multiplied by variables (which are the un-conjugated ones), we must make sure that
                # the un-conjugated ones are of the smaller degree.
                a1 = conj.(a1)
            end
            try
                result = poly_solutions_scaled(moments_scaled, a1, a2, variables, ϵ, missing_func)
                λ > ϵ && (result ./= λ)
                unsafe_push!(solutions_cl, result)
                @verbose_info("Potential solutions:\n", result)
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
    # now we combine all the solutions. In principle, this is Iterators.product, but we only look for compatible entries (as
    # variables in the cliques can overlap).
    @verbose_info("Found all potential individual solutions in ", extraction_time, " seconds. Building iterator.")
    return PolynomialSolutions{problem.complex ? ComplexF64 : Float64, length(cliques) == 1}(
        Dict{eltype(problem.variables),Int}((v, i) for (v, i) in zip(problem.variables, 1:length(problem.variables))),
        cliques,
        finish!(solutions_cl),
        δ,
        verbose
    )
end
poly_solutions(problem::PolyOptProblem, args...; kwargs...) = poly_solutions(SparsityNone(problem), args...; kwargs...)

function Base.iterate(iter::PolynomialSolutions)
    length(iter.cliques) == 0 && return nothing
    # we first fabricate an artificial state whose only purpose is to be iterated in the very next step
    state = ones(Int, length(iter.cliques))
    @inbounds state[1] = 0
    return iterate(iter, state)
end

function Base.iterate(iter::PolynomialSolutions{O}, state::Vector{Int}) where {O}
    verbose = iter.verbose
    # first, we should increase the current state by one
    hasnext = false
    @inbounds for j in 1:length(state)
        if !ismissing(iter.solutions_cl[j]) && state[j] < size(iter.solutions_cl[j], 2)
            state[j] += 1
            fill!(@view(state[1:j-1]), 1)
            hasnext = true
            break
        end
    end
    hasnext || return nothing
    # then, try to get the solution for this state
    solution = Vector{O}(undef, length(iter.var_to_idx))
    @label whiletrue
        fill!(solution, O(NaN))
        @inbounds for (j, (variables, subsolutions)) in enumerate(zip(iter.cliques, iter.solutions_cl))
            ismissing(subsolutions) && continue
            for (var, val) in zip(variables, @view(subsolutions[:, state[j]]))
                var_idx = iter.var_to_idx[var]
                old_val = solution[var_idx]
                if isnan(old_val)
                    solution[var_idx] = val
                elseif abs(old_val - val) > iter.δ
                    @verbose_info("Skipped a solution: Discrepancy ", abs(old_val - val), " > δ = ", iter.δ)
                    @inbounds for j in 1:length(state)
                        if !ismissing(iter.solutions_cl[j]) && state[j] < size(iter.solutions_cl[j], 2)
                            state[j] += 1
                            fill!(@view(state[1:j-1]), 1)
                            @goto whiletrue
                        end
                    end
                    return nothing
                end
            end
        end
    return solution, state
end

Base.IteratorSize(::Type{PolynomialSolutions{O,true}}) where {O} = Base.HasLength()
Base.IteratorSize(::Type{PolynomialSolutions{O,false}}) where {O} = Base.SizeUnknown()
Base.IteratorEltype(::Type{PolynomialSolutions}) = Base.HasEltype()
Base.eltype(::Type{PolynomialSolutions{O}}) where {O} = Vector{O}
function Base.length(iter::PolynomialSolutions{O,true}) where {O}
    @inbounds cl = iter.solutions_cl[1]  # true = there is only one clique
    if ismissing(cl)
        return 0
    else
        return length(iter.solutions_cl[1])
    end
end

function poly_solutions_scaled(moments::Dict{M,R}, a1, a2, variables, ϵ::R, missing_func) where {M,R<:Real}
    O = M <: MonomialComplexContainer ? Complex{R} : R
    Hd1d2 = moment_matrix(moments, a1, a2, 1, missing_func)
    UrSrVrbar = svd!(Hd1d2)
    ϵ *= UrSrVrbar.S[1]
    numEntries = findfirst(x -> x < ϵ, UrSrVrbar.S)
    if numEntries === nothing
        Ur, Sr, Vrbar = UrSrVrbar.U, UrSrVrbar.S, UrSrVrbar.V
    else
        @inbounds Ur, Sr, Vrbar = UrSrVrbar.U[:, 1:numEntries-1], UrSrVrbar.S[1:numEntries-1], UrSrVrbar.V[:, 1:numEntries-1]
    end
    Ms = Vector{Matrix{O}}(undef, length(variables))
    for i in 1:length(variables)
        @inbounds Ms[i] = (Ur' ./ Sr) * moment_matrix(moments, a1, a2, variables[i], missing_func) * Vrbar
    end
    v = eigvecs(sum((2rand(O) - one(O)) * Mi for Mi in Ms))
    len = size(Sr, 1)
    result = Matrix{O}(undef, length(variables), len)
    for j in 1:len
        @inbounds vj = @view(v[:, j])
        for i in 1:length(variables)
            @inbounds Mv = Ms[i] * vj
            if M <: MonomialComplexContainer
                @inbounds result[i, j] = Complex{R}(
                    (dot(real.(Mv), real.(vj)) + dot(imag.(Mv), imag.(vj))),
                    (dot(real.(vj), imag.(Mv)) - dot(real.(Mv), imag.(vj)))
                ) / sum(abs2.(vj))
            else
                @inbounds result[i, j] = (dot(real.(Mv), real.(vj)) + dot(imag.(Mv), imag.(vj))) / sum(abs2.(vj))
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
        @inbounds result = result[:, keep]
    end
    return result
end

"""
    poly_solution_badness(problem::Union{PolyOptProblem,AbstractSparsity}, solution)

Determines the badness of a solution by comparing the value of the objective with the value according to the last optimization
that was done on the problem, and also by checking the violation of the constraints.
The closer the return value is to zero, the better. If the return value is too large, `solution` probably has nothing to do
with the actual solution.

See also [`poly_optimize`](@ref), [`sparse_optimize`](@ref), [`poly_solutions`](@ref), [`poly_all_solutions`](@ref).
"""
function poly_solution_badness(problem::PolyOptProblem, solution::Vector)
    # check whether we can certify optimality
    any(isnan, solution) && return Inf
    solution_map = problem.variables => solution
    violation = abs(problem.objective(solution_map) - last_objective(problem))
    for constr in problem.constraints
        # those are real values, they might just not know it yet
        if constr.type == pctEqualitySimple || constr.type == pctEqualityGröbner || constr.type == pctEqualityNonneg
            new_violation = abs(constr.constraint(solution_map))
        elseif constr.type == pctNonneg
            new_violation = -real(constr.constraint(solution_map))
        elseif constr.type == pctPSD
            new_violation = -real(eigmin((solution_map,) .|> constr.constraint))
        else
            @assert(false)
        end
        new_violation > violation && (violation = new_violation)
    end
    return violation
end
poly_solution_badness(state::AbstractSparsity, args...) = poly_solution_badness(sparse_problem(state), args...)

function default_solution_method end

"""
    poly_all_solutions(problem::Union{PolyOptProblem,AbstractSparsity}, ϵ=1e-6, δ=1e-3; verbose=false,
        rel_threshold=100, abs_threshold=Inf, method::Symbol)

Obtains a vector of all the solutions to a previously optimized problem; then iterates over all of them and grades and sorts
them by their badness. Every solution of the returned vector is a tuple that first contains the optimal point and second the
badness at this point. Solutions that are `rel_threshold` times worse than the best solution or worse than `abs_threshold` will
be dropped from the result.
The currently accepted methods are
- `:mvhankel` to perform a multivariate Hankel decomposition as in [`poly_optimize`](@ref), default for dense problems and
  those with correlative sparsity only.
- `:heuristic` to perform a heuristic method as in [`poly_solutions_heuristic`](@ref), default for problems with term sparsity.

See also [`poly_optimize`](@ref), [`sparse_optimize`](@ref), [`poly_solutions`](@ref), [`poly_solution_badness`](@ref),
[`poly_solutions_heuristic`](@ref).
"""
function poly_all_solutions(state::AbstractSparsity, ϵ::R=1e-6, δ::R=1e-3; verbose::Bool=false,
    rel_threshold::Float64=100., abs_threshold::Float64=Inf, method::Symbol=default_solution_method(state, missing)) where {R<:Real}
    if method === :heuristic
        sol_itr = poly_solutions_heuristic(sparse_problem(state); verbose)
    elseif method === :mvhankel
        sol_itr = poly_solutions(state, ϵ, δ; verbose)
    else
        error("Unknown solution extraction method: $method")
    end
    solutions_fv = FastVec{Tuple{Vector{sparse_problem(state).complex ? ComplexF64 : Float64},Float64}}()
    Base.haslength(sol_itr) && sizehint!(solutions_fv, length(sol_itr))
    for solution in sol_itr
        push!(solutions_fv, (solution, poly_solution_badness(state, solution)))
    end
    solutions = finish!(solutions_fv)
    isempty(solutions) && return solutions
    sort!(solutions, by=last)
    # cut off solutions that are too bad
    threshold = min(rel_threshold * solutions[1][2], abs_threshold)
    return filter(x -> x[2] <= threshold, solutions)
end
poly_all_solutions(problem::PolyOptProblem, args...; kwargs...) = poly_all_solutions(SparsityNone(problem), args...; kwargs...)

"""
    optimality_certificate(problem::PolyOptProblem, ϵ=1e-6)

This function applies the flat extension/truncation criterion to determine whether the optimality of the given problem can be
certified, in which case it returns `:Optimal`. If no such certificate is found, the function returns `:Unknown`. The criterion
is meaningless for sparse problems.
The parameter `ϵ` controls the bound below which singular values are considered to be zero, and its negative below which
eigenvalues are considered to be negative.

See also [`poly_optimize`](@ref).
"""
function optimality_certificate(problem::PolyOptProblem, ϵ::R=1e-6) where {R<:Real}
    isempty(problem.last_moments) && error("The problem was not optimized yet.")
    if problem.complex
        d₀ = maxhalfdegree(problem.objective)
        dconstr = maximum(c -> maxhalfdegree(c.constraint), problem.constraints; init=0)
        d₁ = max(2, dconstr)
        if d₀ == 1 && dconstr == 1
            # this is the one case that is excluded in the loop, but the rank-1 condition is also valid here
            rank(moment_matrix(problem, max_deg=1), atol=ϵ) == 1 && return :Optimal
        end
        n = length(problem.variables)
        for t in problem.degree:-1:max(d₀, d₁)
            Mₜ = moment_matrix(problem, max_deg=t)
            rkMₜ = rank(Mₜ, atol=ϵ)
            rkMₜ == 1 && return :Optimal
            max_deg = t - d₁
            m₁₁ = moment_matrix(problem; max_deg)
            if rkMₜ == rank(m₁₁, atol=ϵ)
                fail = false
                if n > 1
                    @inbounds for i in 1:n-1, j in i:n
                        zᵢ, zⱼ = problem.variables[i], problem.variables[j]
                        m₁₂ = moment_matrix(problem; max_deg, prefix=zᵢ)
                        m₁₃ = moment_matrix(problem; max_deg, prefix=zⱼ)
                        m₂₂ = moment_matrix(problem; max_deg, prefix=zᵢ*conj(zᵢ))
                        m₂₃ = moment_matrix(problem; max_deg, prefix=zⱼ*conj(zᵢ))
                        m₃₃ = moment_matrix(problem; max_deg, prefix=zⱼ*conj(zⱼ))
                        if min(eigvals!(Hermitian([m₁₁       m₁₂       m₁₃
                                                   conj(m₁₂) m₂₂       m₂₃
                                                   conj(m₁₃) conj(m₂₃) m₃₃]))) < -ϵ
                            fail = true
                            break
                        end
                    end
                else
                    @inbounds z = problem.variables[1]
                    m₁₂ = moment_matrix(problem; max_deg, prefix=z)
                    m₂₂ = moment_matrix(problem; max_deg, prefix=z*conj(z))
                    if min(eigvals!(Hermitian([m₁₁ m₁₂; conj(m₁₂) m₂₂]))) < -ϵ
                        fail = true
                    end
                end
                fail || return :Optimal
            end
        end
    else
        d₀ = maxhalfdegree(problem.objective)
        d₁ = max(1, maximum(c -> maxhalfdegree(c.constraint), problem.constraints; init=0))
        for t in problem.degree:-1:max(d₀, d₁)
            rkMₜ = rank(moment_matrix(problem, max_deg=t), atol=ϵ)
            if rkMₜ == 1 || rkMₜ == rank(moment_matrix(problem, max_deg=t-d₁), atol=ϵ)
                return :Optimal
            end
        end
    end
    return :Unknown
end