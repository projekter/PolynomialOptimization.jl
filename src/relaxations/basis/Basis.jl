"""
    AbstractRelaxationBasis{Prob} <: AbstractRelaxation{Prob}

An `AbstractRelaxationBasis` is a relaxation of a polynomial optimization problem that is built using a single basis for
everything (objective and constraints). The groupings for the individual elements will come from a degree truncation of the
same shared basis for all constituents of the problem (intersected with a parent grouping).
"""
abstract type AbstractRelaxationBasis{Prob<:Problem} <: AbstractRelaxation{Prob} end

basis(relaxation::AbstractRelaxationBasis) = relaxation.basis

function basis(relaxation::AbstractRelaxationBasis, i::Int)
    i == 1 || throw(ArgumentError("Unknown clique index: $i"))
    return relaxation.basis
end

function groupings(problem::Problem{Prob}, basis::AbstractVector{M}, degree::Integer, parent) where
    {Nr,Nc,I<:Integer,M<:SimpleMonomial{Nr,Nc,I},Prob<:SimplePolynomial{<:Any,Nr,Nc}}
    return embed(RelaxationGroupings(
        SimpleMonomialVector{Nr,Nc,I}[basis],
        [SimpleMonomialVector{Nr,Nc,I}[truncate_basis(basis, degree - maxhalfdegree(x))] for x in problem.constr_zero],
        [SimpleMonomialVector{Nr,Nc,I}[truncate_basis(basis, degree - maxhalfdegree(x))] for x in problem.constr_nonneg],
        [SimpleMonomialVector{Nr,Nc,I}[truncate_basis(basis, degree - maxhalfdegree(x))] for x in problem.constr_psd],
        [filter(∘(!, isconj), variables(problem.objective))]
    ), parent, true)
end

MultivariatePolynomials.degree(relaxation::AbstractRelaxationBasis) = relaxation.degree

function Base.show(io::IO, m::MIME"text/plain", relaxation::AbstractRelaxationBasis)
    _show(io, m, relaxation)
    print(io, "\nRelaxation degree: ", relaxation.degree)
end

struct _DummyMonomial
    degree::Int
end

MultivariatePolynomials.degree(d::_DummyMonomial) = d.degree

function truncate_basis(v::AbstractVector{M} where {M<:SimpleMonomial}, maxdeg::Integer)
    idx = searchsortedlast(v, _DummyMonomial(maxdeg), by=degree)
    if idx < firstindex(v)
        return @view(v[begin:end])
    else
        return @view(v[1:idx])
    end
end

include("./Dense.jl")
include("./Custom.jl")
include("./Newton.jl")