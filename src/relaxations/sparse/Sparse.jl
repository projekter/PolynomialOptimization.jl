"""
    AbstractRelaxationSparse{Prob} <: AbstractRelaxation{Prob}

An `AbstractRelaxationSparse` is a relaxation of a polynomial optimization problem that applies sparsity methods to reduce the
size of the associated problem, possibly at the expense of lowering the objective bound.
"""
abstract type AbstractRelaxationSparse{Prob<:Problem} <: AbstractRelaxation{Prob} end

basis(relaxation::AbstractRelaxationSparse) = basis(relaxation.parent)

function basis(relaxation::AbstractRelaxationSparse, i::Int)
    1 ≤ i ≤ length(relaxation.groupings.var_cliques) || throw(ArgumentError("Unknown clique index: $i"))
    return filter(Base.Fix2(SimplePolynomials.effective_variables_in, Set(relaxation.groupings.var_cliques[i])),
        basis(relaxation))
end

include("./Chordal.jl")
include("./Correlative.jl")
include("./Term.jl")
include("./CorrelativeTerm.jl")
include("./CliqueMerged.jl")