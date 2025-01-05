export AbstractSolver, mindex

"""
    AbstractSolver{T,V<:Real}

Abstract supertype for any solver state. `T` is the type that is returned by calling [`mindex`](@ref) with such a state; `V` is
the type of the coefficients used in the solver. Using the default implementation for [`mindex`](@ref), `T` will be `UInt`.
Most likely, `V` will be `Float64`.
"""
abstract type AbstractSolver{T,V<:Real} end

"""
    mindex(::AbstractSolver{T}, monomials::SimpleMonomialOrConj...)::T

Calculates the index that the product of all monomials will have in the SDP represented by `state`.
The default implementation calculates the one-based monomial index according to a dense deglex order and returns an `UInt`.
The returned index is arbitrary as long as it is unique for the total monomial.
"""
@inline mindex(::AbstractSolver{T}, monomials::SimpleMonomialOrConj{Nr,Nc}...) where {T,Nr,Nc} =
    monomial_index(ExponentsAll{Nr+2Nc,UInt}(), monomials...)::T

include("./Support.jl")
include("./DataTypes.jl")
include("./Representations.jl")
include("./Counters.jl")

include("./MomentInterface.jl")
include("./SOSInterface.jl")
include("./Extraction.jl")