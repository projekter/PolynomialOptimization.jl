struct SparsityCorrelativeTerm{P<:Problem,S<:SparsityTerm{<:Integer,P}} <: AbstractRelaxationSparse{P}
    s::S

    function SparsityCorrelativeTerm(relaxation::SparsityCorrelative;
        method::Union{TermMode,<:AbstractVector{TermMode}}=TERM_MODE_BLOCK, kwargs...)
        s = SparsityTerm(relaxation; method, kwargs...)
        new{typeof(poly_problem(s)),typeof(s)}(s)
    end
end

@doc """
    SparsityCorrelativeTerm(relaxation::AbstractRelaxation; method=TERM_MODE_BLOCK, kwargs...)

Analyze both the [correlative as well as the term sparsity](http://arxiv.org/abs/2005.02828v2) of the problem.
This is the most versatile kind of sparsity analysis, combining the effects of correlative sparsity with term analysis per
clique. However, it is nothing more than first performing correlative sparsity analysis, followed by term sparsity analysis.
This constructor will take all keyword arguments and distribute them appropriately to the [`SparsityCorrelative`](@ref)
and [`SparsityTerm`](@ref) constructors.
The returned object will be a very thin wrapper around [`SparsityTerm`](@ref), with the only difference in printing;
[`SparsityCorrelativeTerm`](@ref) objects by default print the clique grouping. Note that the same can be achieved for any
relaxation (or groupings of a relaxation) if the IO parameter `bycliques` is set to `true`.

See also [`SparsityCorrelative`](@ref), [`SparsityTermBlock`](@ref), [`SparsityTermChordal`](@ref), [`TermMode`](@ref).
"""
SparsityCorrelativeTerm(relaxation::AbstractRelaxation; high_order_zero=missing,
    high_order_nonneg=missing, high_order_psd=missing, low_order_zero=missing, low_order_nonneg=missing,
    low_order_psd=missing, chordal_completion::Bool=true, verbose::Bool=false,
    method::Union{TermMode,<:AbstractVector{TermMode}}=TERM_MODE_BLOCK, termkwargs...) =
    SparsityCorrelativeTerm(SparsityCorrelative(relaxation; high_order_zero, high_order_nonneg, high_order_psd,
        low_order_zero, low_order_nonneg, low_order_psd, chordal_completion, verbose); verbose, method, termkwargs...)

"""
    SparsityCorrelativeTerm(relaxation::SparsityCorrelative; method=TERM_MODE_BLOCK, kwargs...)

This form allows to wrap an already created correlative sparsity pattern into a term sparsity pattern.

See also [`SparsityCorrelative`](@ref), [`TermMode`](@ref).
"""
SparsityCorrelativeTerm(::SparsityCorrelative) # the main documentation should be for the Problem form, so put this later

Base.getproperty(relaxation::SparsityCorrelativeTerm, f::Symbol) = getproperty(getfield(relaxation, :s), f)
Base.propertynames(relaxation::SparsityCorrelativeTerm) = propertynames(getfield(relaxation, :s))

_show(io::IO, m::MIME"text/plain", x::SparsityCorrelativeTerm) =
    _show(IOContext(io, :bycliques => true), m, x, typeof(x).name.name)

function iterate!(relaxation::SparsityCorrelativeTerm; kwargs...)
    isnothing(iterate!(getfield(relaxation, :s); kwargs...)) && return nothing
    return relaxation
end

default_solution_method(::SparsityCorrelativeTerm) = :heuristic