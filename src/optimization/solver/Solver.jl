module Solver

using ..SimplePolynomials, ..PolynomialOptimization, MultivariatePolynomials, LinearAlgebra, SparseArrays, Reexport
using ..SimplePolynomials: SimpleMonomialOrConj, SimpleConjMonomial, _get_I
@reexport using ..FastVector
using ..FastVector: overallocation
using ..PolynomialOptimization: @assert, @inbounds, @verbose_info, @capture, @unroll, FastKey, StackVec, Problem, sort_along!
import ..PolynomialOptimization: MomentVector
using ..SimplePolynomials.MultivariateExponents: ExponentsAll, ExponentsDegree, Unsafe, unsafe
using ..Relaxation: AbstractRelaxation, RelaxationGroupings
import LinearAlgebra: issuccess, UpperOrUnitUpperTriangular, LowerOrUnitLowerTriangular
# We re-export things that implementations of solvers (which is the only place where this module should be use'd) will most
# likely need
export
    poly_problem, Problem, MomentVector, StackVec, FastKey, sort_along!, @verbose_info, @capture, # from PolynomialOptimization
    AbstractRelaxation, RelaxationGroupings, # from Relaxation
    SimpleMonomialOrConj, SimpleConjMonomial, monomial_index, _get_I, # from SimplePolynomials
    overallocation, # from FastVector (not exported)
    poly_optimize, solver_methods, @solver_alias,
    issuccess # from LinearAlgebra

const solver_methods = Symbol[]

function default_solver_method()
    isempty(solver_methods) &&
        error("No solver method is available. Load a solver package that provides such a method (e.g., Mosek)")
    @inbounds return solver_methods[begin]
end

"""
    @solver_alias(alias, original)

Defines the solver identifier `alias` to map to the same optimization routine as `original`.
"""
macro solver_alias(alias, original)
    quote
        $Solver.poly_optimize(::Val{$(QuoteNode(alias))}, relaxation::$AbstractRelaxation, args...; kwargs...) =
            $poly_optimize(Val($(QuoteNode(original))), relaxation, args...; kwargs...)
    end
end

function poly_optimize end

"""
    poly_optimize(::Val{method}, relaxation::AbstractRelaxation,
        groupings::RelaxationGroupings; representation, verbose, kwargs...)

This is the central entry point that a solver has to implement. It has to carry out the optimization and must return a tuple of
three values:
1. An internal state that can be used later on to access the solver again to extract solutions or (if possible) reoptimization.
2. The success status as returned by the solver and understood by the [`issuccess`](@ref issuccess(::Val, ::Any))
   implementation.
3. The minimum value of the objective.

See also [`moment_setup!`](@ref), [`sos_setup!`](@ref).
"""
poly_optimize(::Val, ::AbstractRelaxation, ::RelaxationGroupings)

"""
    issuccess(::Val{method}, status)

A solver must implement this method for all of its possible methods to indicate whether a status `status` signifies success.
"""
issuccess(::Val, ::Any)

struct _Copied{X}
    data::X
end

"""
    poly_optimize(::Val{method}, oldstate, relaxation::AbstractRelaxation,
        groupings::RelaxationGroupings; representation, verbose, kwargs)

A solver that supports re-optimization of an already optimized problem with changed rotations on the DD and SDD representations
should implement this method. It is guaranteed that only the rotations change, and only in a structure-preserving way (diagonal
and nondiagonal rotations will not interchange). The return value is as documented in [`poly_optimize`](@ref), and the
`oldstate` parameter holds the first return value of the previous call to [`poly_optimize`](@ref).
"""
function poly_optimize(method::Val, ::Any, relaxation::AbstractRelaxation, groupings::RelaxationGroupings; representation,
    kwargs...)
    @info("The chosen solver does not support re-optimization. Starting from the beginning.")
    if representation isa Rerepresent
        # we can disable any constraints on not changing the representation, as everything is new anyway
        representation = Rerepresent(representation, false)
    end
    return _Copied(poly_optimize(method, relaxation, groupings; representation, kwargs...))
end


include("./Helpers.jl")
include("./Interface.jl")
include("./SOSHelpers.jl")
include("./MomentHelpers.jl")
include("./AbstractAPISolver.jl")
include("./AbstractSparseMatrixSolver.jl")
include("./PrimalMomentHelpers.jl")

end