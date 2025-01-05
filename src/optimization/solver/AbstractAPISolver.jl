export AbstractAPISolver

"""
    AbstractAPISolver{K<:Integer,T,V} <: AbstractSolver{T,V}

Superclass for a solver that requires new variables/constraints to be added via API calls. Solvers that are of this type must
implement [`append!`](@ref append!(::AbstractAPISolver{K}, ::K) where {K<:Integer}) in such a way that they directly add a
variable (moment-case) to or constraint (SOS-case) to the solver.
Concrete types that inherit from `AbstractAPISolver` must have a property `mon_to_solver::Dict{FastKey{K},T}`.
"""
abstract type AbstractAPISolver{K<:Integer,T,V} <: AbstractSolver{T,V} end

"""
    append!(solver::AbstractAPISolver{K}, key::K)

Appends at least one new variable (moment-case) or constraint (SOS-case) to the solver `state` that represents the monomial
given by `key`.
"""
Base.append!(::AbstractAPISolver{K}, ::K) where {K<:Integer}

@inline function mindex(solver::AbstractAPISolver{<:Integer,T}, monomials::SimpleMonomialOrConj{Nr,Nc}...) where {T,Nr,Nc}
    idx = monomial_index(monomials...)
    dictidx = Base.ht_keyindex(solver.mon_to_solver, FastKey(idx))
    @inbounds return (dictidx < 0 ?
        append!(solver, idx) :
        solver.mon_to_solver.vals[dictidx]
    )::T
end

"""
    MomentVector(relaxation::AbstractRelaxation, moments::Vector{<:Real},
        solver::AbstractAPISolver)

Given the moments vector as obtained from an [`AbstractAPISolver`](@ref), convert it to a [`MomentVector`](@ref). Note that this
function is not fully type-stable, as the result may be based either on a dense or sparse vector depending on the relaxation.
"""
function MomentVector(relaxation::AbstractRelaxation{<:Problem{<:SimplePolynomial{<:Any,Nr,Nc}}}, moments::Vector{V},
    solver::AbstractAPISolver{K}) where {Nr,Nc,K<:Integer,V<:Real}
    @assert(length(moments) ≥ length(solver.mon_to_solver))
    # In any case, our variables will not be in the proper order. First figure this out.
    # Dict does not preserve the insertion order. While we could use OrderedDict instead, we need to do lots of insertions and
    # lookups and only once access the elements in the insertion order; so it is probably better to do the sorting once.
    mon_pos = convert.(K, keys(solver.mon_to_solver))
    var_vals = collect(values(solver.mon_to_solver))
    sort_along!(var_vals, mon_pos)
    # Now mon_pos is in insertion order. There might still not be a 1:1 correspondence between x and mon_pos, as slack
    # variables could be present if the solver requires those.
    if length(var_vals) < length(moments)
        var_vals .+= one(Cint)
        @inbounds moments = moments[var_vals] # remove the slack
    end
    # Now we have the 1:1 correspondence, but we want the actual monomial order.
    sort_along!(mon_pos, moments)
    max_mons = mon_pos[end] # total number of monomials in ExponentsAll (which in the complex case can exceed even the length
                            # of the dense relaxation significantly - because the exponents don't know anything about complex
                            # structure)
    @assert(length(moments) ≤ max_mons)
    # Finally, x is ordered!
    if length(moments) == mon_pos # (real) dense case
        solution = moments
    elseif 3length(mon_pos) < max_mons
        solution = SparseVector(max_mons, mon_pos, moments)
    else
        solution = fill(NaN, max_mons)
        copyto!(@view(solution[mon_pos]), moments)
    end
    return MomentVector(relaxation, ExponentsAll{Nr+2Nc,K}(), solution)
end