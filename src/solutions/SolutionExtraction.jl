export poly_solutions, poly_solution_badness, poly_all_solutions

struct MonomialMissing <: Exception end

struct PolynomialSolutions{R<:Real,V<:Union{R,Complex{R}},Nr,Nc,Var<:SimpleVariable{Nr,Nc},has_length}
    cliques::Vector{Vector{Var}}
    solutions_cl::Vector{Union{Missing,Matrix{V}}}
    δ::R
    verbose::Bool
end

function Base.show(io::IO, ::MIME"text/plain", s::PolynomialSolutions)
    println(io, "Solution iterator for a polynomial optimization problem")
    if isone(length(s.cliques))
        print(io, "Number of solutions: ", size(s.solutions_cl[1], 2))
    else
        println(io, "Number of solutions: unknown")
        println(io, "Number of cliques: ", length(s.cliques))
        for (i, (clique, clique_sol)) in enumerate(zip(s.cliques, s.solutions_cl))
            println(io, "Clique #", i, ": ", length(clique.variables), " variable")
            isone(length(clique)) || print(io, "s")
            print(io, ", ", size(clique_sol, 2), " subsolution")
            isone(size(clique_sol, 2)) || print(io, "s")
        end
        println(io, "Maximally allowed discrepancy between same variables in different cliques: ", s.δ)
        print(io, "Verbose display while assembling solutions: ", s.verbose ? "enabled" : "disabled")
    end
end

absval(x) = isnan(x) ? zero(x) : abs(x)

"""
    poly_solutions([method, ]result::Result, args...; verbose, kwargs...)

Extracts solutions from a polynomial optimization result using the method `method`. Depending on the chosen method, the result
may be an iterator or a vector. Consult the documentation of the methods for further information.
If `method` is omitted, a default method will be chosen according to the relaxation that was used for the optimization.

See also [`poly_optimize`](@ref), [`poly_all_solutions`](@ref).
"""
function poly_solutions end

poly_solutions(result::Result, args...; kwargs...) =
    poly_solutions(Relaxation.default_solution_method(result), args...; kwargs...)

poly_solutions(method::Symbol, rest...; kwrest...) = poly_solutions(Val(method), rest...; kwrest...)

function Base.iterate(iter::PolynomialSolutions)
    length(iter.cliques) == 0 && return nothing
    # we first fabricate an artificial state whose only purpose is to be iterated in the very next step
    state = ones(Int, length(iter.cliques))
    @inbounds state[1] = 0
    return iterate(iter, state)
end

function Base.iterate(iter::PolynomialSolutions{R,V,Nr,Nc}, state::Vector{Int}) where {R<:Real,V<:Union{R,Complex{R}},Nr,Nc}
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
    solution = Vector{V}(undef, Nr + Nc)
    @label whiletrue
        fill!(solution, V(NaN))
        @inbounds for (j, (variables, subsolutions)) in enumerate(zip(iter.cliques, iter.solutions_cl))
            ismissing(subsolutions) && continue
            for (var, val) in zip(variables, @view(subsolutions[:, state[j]]))
                var_idx = variable_index(var)
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

Base.IteratorSize(::Type{PolynomialSolutions{R,V,Nr,Nc,<:SimpleVariable{Nr,Nc},true} where {R<:Real,V<:Union{R,Complex{R}},Nr,Nc}}) =
    Base.HasLength()
Base.IteratorSize(::Type{PolynomialSolutions{R,V,Nr,Nc,<:SimpleVariable{Nr,Nc},false} where {R<:Real,V<:Union{R,Complex{R}},Nr,Nc}}) =
    Base.SizeUnknown()
Base.IteratorEltype(::Type{PolynomialSolutions}) = Base.HasEltype()
Base.eltype(::Type{<:PolynomialSolutions{R,V}}) where {R<:Real,V<:Union{R,Complex{R}}} = Vector{V}
function Base.length(iter::PolynomialSolutions{R,V,Nr,Nc,<:SimpleVariable{Nr,Nc},true}) where {R<:Real,V<:Union{R,Complex{R}},Nr,Nc}
    @inbounds cl = iter.solutions_cl[1]  # true = there is exactly one clique
    if ismissing(cl)
        return 0
    else
        return size(iter.solutions_cl[1], 2)
    end
end

"""
    poly_solution_badness(result::Result, solution)

Determines the badness of a solution by comparing the value of the objective with the value according to the optimization given
in `result`, and also by checking the violation of the constraints.
The closer the return value is to zero, the better. If the return value is too large, `solution` probably has nothing to do
with the actual solution.

See also [`poly_optimize`](@ref), [`poly_solutions`](@ref), [`poly_all_solutions`](@ref).
"""
function poly_solution_badness(result::Result, solution::Vector)
    # check whether we can certify optimality
    any(isnan, solution) && return Inf
    problem = poly_problem(result)
    violation = abs(problem.objective(solution) - result.objective)
    for constr in problem.constr_zero
        new_violation = abs(constr(solution))
        new_violation > violation && (violation = new_violation)
    end
    for constr in problem.constr_nonneg
        new_violation = -real(constr(solution)) # they are real, but they might not know it yet
        new_violation > violation && (violation = new_violation)
    end
    for constr in problem.constr_psd
        new_violation = -eigmin((eltype(solution) <: Complex ? Hermitian : Symmetric)((solution,) .|> constr))
        new_violation > violation && (violation = new_violation)
    end
    return violation
end

Relaxation.default_solution_method(result::Result) = Relaxation.default_solution_method(result.relaxation)

"""
    poly_all_solutions([method, ]result::Result, args...; verbose=false, rel_threshold=100,
        abs_threshold=Inf, kwargs...)

Obtains a vector of all the solutions to a previously optimized problem; then iterates over all of them and grades and sorts
them by their badness. Every solution of the returned vector is a tuple that first contains the optimal point and second the
badness at this point. Solutions that are `rel_threshold` times worse than the best solution or worse than `abs_threshold` will
be dropped from the result.

See also [`poly_optimize`](@ref), [`poly_solutions`](@ref).
"""
function poly_all_solutions(@nospecialize(method::Val), result::Result, args...; verbose::Bool=false,
    rel_threshold::Float64=100., abs_threshold::Float64=Inf, kwargs...)
    sol_itr = poly_solutions(method, result, args...; verbose, kwargs...)
    solutions_fv = FastVec{Tuple{Vector{eltype(result)},Float64}}()
    Base.haslength(sol_itr) && sizehint!(solutions_fv, length(sol_itr))
    for solution in sol_itr
        push!(solutions_fv, (solution, poly_solution_badness(result, solution)))
    end
    solutions = finish!(solutions_fv)
    isempty(solutions) && return solutions
    sort!(solutions, by=last)
    # cut off solutions that are too bad (maybe if the list grows too fast, we could do this during construction?)
    threshold = min(rel_threshold * solutions[1][2], abs_threshold)
    return filter(x -> x[2] ≤ threshold, solutions)
end

poly_all_solutions(result::Result, args...; kwargs...) =
    poly_all_solutions(Relaxation.default_solution_method(result), result, args...; kwargs...)

poly_all_solutions(method::Symbol, rest...; kwrest...) = poly_all_solutions(Val(method), rest...; kwrest...)

include("./MVHankel.jl")
include("./Heuristic.jl")