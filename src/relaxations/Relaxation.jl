module Relaxation

using ..IntPolynomials, .IntPolynomials.MultivariateExponents, ..PolynomialOptimization, MultivariatePolynomials,
    ..PolynomialOptimization.FastVector, SortingAlgorithms
import StatsBase, Graphs, CliqueTrees
using ..PolynomialOptimization: @assert, @capture, @inbounds, Problem, @unroll, @verbose_info, issubset_sorted
import ..PolynomialOptimization: poly_problem, iterate!

export AbstractRelaxation, basis, groupings, iterate!

"""
    AbstractRelaxation

This is the general abstract type for any kind of relaxation of a polynomial optimization problem.
Its concrete types can be used for analyzing and optimizing the problem.

See also [`poly_problem`](@ref), [`Problem`](@ref), [`poly_optimize`](@ref).
"""
abstract type AbstractRelaxation{Prob<:Problem} end

@doc raw"""
    RelaxationGroupings

Contains information about how the elements in a certain (sparse) polynomial optimization problem combine.
Groupings are contained in the fields `obj`, `zero`, `nonneg`, and `psd`:
- ``\sum_i \mathit{obj}_i^\top \sigma_i \overline{\mathit{obj}_i}`` is the SOS representation of the objective with
  ``\sigma_i \succeq 0``
- ``\sum_i \mathit{zero}_{k, i}^\top f_k`` is the prefactor for the kᵗʰ equality
  constraint with ``f_k`` a free vector
- ``\sum_i \mathit{nonneg}_{k, i}^\top \sigma_{k, i} \overline{\mathit{nonneg}_{k, i}}`` is the SOS representation of
  the prefactor of the kᵗʰ nonnegative constraint with ``\sigma_{k, i} \succeq 0``
- ``\sum_i (\mathit{psd}_{k, i}^\top \otimes \mathbb1) Z_{k, i} (\overline{\mathit{psd}_{k, i}} \otimes \mathbb1)``
  is the SOS matrix representation of the prefactor of the kᵗʰ PSD constraint with ``Z_{k, i} \succeq 0``
The field `var_cliques` contains a list of sets of variables, each corresponding to a variable clique in the total problem. In
the complex case, only the declared variables are returned, not their conjugates.
"""
struct RelaxationGroupings{Nr,Nc,I<:Integer,V<:IntVariable{Nr,Nc}}
    # These can all be different types of monomial vectors, but we don't want to have a type explosion with nested relaxation
    # groupings. Keep it dynamic.
    obj::Vector{IntMonomialVector{Nr,Nc,I}}
    zeros::Vector{Vector{IntMonomialVector{Nr,Nc,I}}}
    nonnegs::Vector{Vector{IntMonomialVector{Nr,Nc,I}}}
    psds::Vector{Vector{IntMonomialVector{Nr,Nc,I}}}
    var_cliques::Vector{Vector{V}}
end

_lensort(x) = (-length(x), x) # use as "by" argument for sort to sort with descending length, standard tie-breaker

# We don't store the association of a grouping with a clique; we just find it anew every time. While this is not too efficient,
# it also doesn't occur so often that we actually need this information.
_findclique(grouping::IntMonomialVector{Nr,Nc}, var_cliques::(Vector{Vector{V}} where V<:IntVariable)) where {Nr,Nc} =
    findlast(Base.Fix1(⊆, effective_variables(grouping)), var_cliques) # put the grouping in the smallest fitting clique

function _show_groupings(io::IO, grouping::Vector{<:IntMonomialVector}, cliques)
    if get(io, :bycliques, false)::Bool
        noc = IOContext(io, :bycliques => false)
        for i in 1:length(cliques)
            inclique = filter(let i=i; g -> _findclique(g, cliques) == i end, grouping)
            if !isempty(inclique)
                print(io, "\n> Clique #", i, ": ")
                _show_groupings(noc, inclique, cliques)
            end
        end
        return
    end
    lg = length(grouping)
    print(io, lg, " block", isone(lg) ? "" : "s")
    iszero(lg) && return
    lensorted = sort(grouping, by=_lensort)
    len = floor(Int, log10(length(@inbounds lensorted[begin]))) +1
    limit = get(io, :limit, false)::Bool
    for block in Iterators.take(lensorted, limit ? 5 : length(lensorted))
        # we must do the printing manually to avoid all the type cluttering. We can assume that a grouping is never empty.
        print(io, "\n  ", lpad(length(block), len, " "), " [")
        a, rest = Iterators.peel(block)
        show(io, "text/plain", a)
        for x in Iterators.take(rest, limit ? 10 : length(block) -1)
            print(io, ", ")
            show(io, "text/plain", x)
        end
        limit && length(block) > 10 && print(io, ", ...")
        print(io, "]")
    end
    limit && length(lensorted) > 5 && print(io, "\n  ", lpad("⋮", len, " "))
end

function Base.show(io::IO, m::MIME"text/plain", groupings::RelaxationGroupings{Nr,Nc,I}) where {Nr,Nc,I<:Integer}
    println(io, "Groupings for the relaxation of a polynomial optimization problem\nVariable cliques\n================")
    bycliques = get(io, :bycliques, false)::Bool
    if bycliques
        cliquelen = floor(Int, log10(length(groupings.var_cliques))) +1
    end
    for (i, clique) in enumerate(groupings.var_cliques)
        if bycliques
            print(io, "#", lpad(i, cliquelen, " "), ": ")
        end
        print(io, "[")
        a, rest = Iterators.peel(clique)
        show(io, "text/plain", a)
        for x in rest
            print(io, ", ")
            show(io, "text/plain", x)
        end
        println(io, "]")
    end
    print(io, "\nBlock groupings\n===============\nObjective: ")
    _show_groupings(io, groupings.obj, groupings.var_cliques)
    for (name, f) in (("Equality", :zeros), ("Nonnegative", :nonnegs), ("Semidefinite", :psds))
        block = getproperty(groupings, f)::Vector{Vector{IntMonomialVector{Nr,Nc,I}}}
        if !isempty(block)
            for (i, constr) in enumerate(block)
                print(io, "\n", name, " constraint #", i, ": ")
                _show_groupings(io, constr, groupings.var_cliques)
            end
        end
    end
end

function embed!(to::AbstractVector{X}, new::X, olds::AbstractVector{X}) where {X}
    for oldᵢ in olds
        if new ⊆ oldᵢ
            push!(to, new)
            return true
        end
    end
    temp = sort!(intersect.((new,), olds), by=_lensort)
    @inbounds for i in length(temp):-1:2
        tempᵢ = temp[i]
        for tempⱼ in @view(temp[1:i-1])
            if tempᵢ ⊆ tempⱼ
                deleteat!(temp, i)
                break
            end
        end
    end
    append!(to, temp)
    return false
end

function embed(news::AbstractVector{X}, olds::AbstractVector{X}, news_is_clean::Bool) where {X}
    complete = Threads.Atomic{Bool}(true)
    nos = Threads.nthreads()
    batch = div(length(news), nos, RoundUp)
    tos = [FastVec{X}(buffer=batch) for _ in 1:nos]
    @inbounds Threads.@threads for i in 1:nos
        items = @view(news[(i-1)*batch+1:min(i*batch, length(news))])
        completeᵢ = true
        toᵢ = tos[i]
        for new in items
            completeᵢ &= embed!(toᵢ, new, olds)
        end
        Threads.atomic_and!(complete, completeᵢ)
        if X <: AbstractVector
            for toᵢᵢ in toᵢ
                sort!(toᵢᵢ)
            end
        end
        sort!(toᵢ, by=_lensort)
    end
    to = sizehint!(X[], sum(length, tos, init=0))
    for toᵢ in tos
        append!(to, toᵢ)
    end
    result = Base._groupedunique!(sort!(to, by=_lensort, alg=TimSort)) # use TimSort as the subslices are already presorted
    news_is_clean && complete[] && return result
    # it is not guaranteed that news is completely subset-free, as it might have been constructed from different sources
    deletes = fill(false, length(result))
    @inbounds Threads.@threads for i in 2:length(result)
        resultᵢ = result[i]
        for j in 1:i-1
            if !deletes[j] && resultᵢ ⊆ result[j]
                deletes[i] = true
            end
        end
    end
    deleteat!(result, deletes)
    return result
end

function embed(new::RG, old::RG, new_is_clean::Bool) where {Nr,Nc,I<:Integer,RG<:RelaxationGroupings{Nr,Nc,I}}
    (length(new.zeros) == length(old.zeros) && length(new.nonnegs) == length(old.nonnegs) &&
        length(new.psds) == length(old.psds)) ||
        throw(ArgumentError("Cannot embed two relaxation groupings for different optimization problems"))
    newobj = embed(new.obj, old.obj, new_is_clean)
    newzeros = embed.(new.zeros, old.zeros, new_is_clean)
    newnonnegs = embed.(new.nonnegs, old.nonnegs, new_is_clean)
    newpsds = embed.(new.psds, old.psds, new_is_clean)
    newcliques = embed(new.var_cliques, old.var_cliques, new_is_clean)
    return RG(newobj, newzeros, newnonnegs, newpsds, newcliques)
end
embed(new::RelaxationGroupings, ::Nothing, ::Bool) = new

Base.:(==)(g₁::G, g₂::G) where {G<:RelaxationGroupings} = g₁.obj == g₂.obj && g₁.zeros == g₂.zeros &&
    g₁.nonnegs == g₂.nonnegs && g₁.psds == g₂.psds && g₁.var_cliques == g₂.var_cliques

"""
    poly_problem(relaxation::AbstractRelaxation)

Returns the original problem associated with a relaxation.
"""
poly_problem(relaxation::AbstractRelaxation) = relaxation.problem

"""
    basis(relaxation::AbstractRelaxation[, clique::Int]) -> IntMonomialVector

Constructs the basis that is associated with a given polynomial relaxation. If `clique` is given, only the monomials that are
relevant for the given clique must be returned.
"""
function basis end

"""
    groupings(relaxation::AbstractRelaxation) -> RelaxationGroupings

Analyze the current state and return the bases and cliques as indicated by its relaxation in a [`RelaxationGroupings`](@ref)
struct.
"""
groupings(relaxation::AbstractRelaxation) = relaxation.groupings
groupings(::Nothing) = nothing

"""
    iterate!(relaxation::AbstractRelaxation)

Some sparse polynomial optimization relaxations allow to iterate their sparsity, which will lead to a more dense representation
and might give better bounds at the expense of a more costly optimization. Return `nothing` if the iterations converged
(`state` did not change any more), else return the new state. Note that `state` will be modified.
"""
iterate!(::AbstractRelaxation) = nothing

"""
    Relaxation.XXX(problem::Problem[, degree]; kwargs...)

This is a convenience wrapper for `Relaxation.XXX(Relaxation.Dense(problem, degree))` that works for any
[`AbstractRelaxation`](@ref) `XXX`.
`degree` is the degree of the Lasserre relaxation, which must be larger or equal to the halfdegree of all polynomials that are
involved. If `degree` is omitted, the minimum required degree will be used.
Specifying a degree larger than the minimal only makes sense if there are inequality or PSD constraints present, else it
needlessly complicates calculations without any benefit.

The keyword arguments will be passed on to the constructor of `XXX`.
"""
function (r::Type{<:AbstractRelaxation})(problem::Problem, args...; kwargs...)
    r <: Dense && throw(MethodError(r, (problem, args...)))
    # to avoid infinite recursion if the arguments did not match
    return r(Dense(problem, args...); kwargs...)
end

function _show(io::IO, m::MIME"text/plain", x::AbstractRelaxation, name=typeof(x).name.name)
    groups = groupings(x)
    # we don't want to print the fully parameterized type type
    print(io, "Relaxation.", name, " of a polynomial optimization problem")
    bycliques = get(io, :bycliques, false)::Bool
    if bycliques
        cliquesizes_psd = [Dict{Int,Int}() for _ in 1:length(groups.var_cliques)]
        cliquesizes_lin = [Dict{Int,Int}() for _ in 1:length(groups.var_cliques)]
        for grouping in groups.obj
            cliquesize = cliquesizes_psd[_findclique(grouping, groups.var_cliques)]
            l = length(grouping)
            cliquesize[l] = get!(cliquesize, l, 0) +1
        end
        for constr in groups.zeros, grouping in constr
            cliquesize = cliquesizes_lin[_findclique(grouping, groups.var_cliques)]
            l = length(grouping)
            cliquesize[l] = get!(cliquesize, l, 0) +1
        end
        for constr in groups.nonnegs, grouping in constr
            cliquesize = cliquesizes_psd[_findclique(grouping, groups.var_cliques)]
            l = length(grouping)
            cliquesize[l] = get!(cliquesize, l, 0) +1
        end
        for (constr, constr_mat) in zip(groups.psds, poly_problem(x).constr_psd), grouping in constr
            cliquesize = cliquesizes_psd[_findclique(grouping, groups.var_cliques)]
            l = length(grouping) * size(constr_mat, 1)
            cliquesize[l] = get!(cliquesize, l, 0) +1
        end
        for (i, (va, size_psd, size_lin)) in enumerate(zip(groups.var_cliques, cliquesizes_psd, cliquesizes_lin))
            print(io, "\n> Clique #", i, ": ", join(va, ", "))
            isempty(size_psd) || print(io, "\n  PSD block sizes:\n    ", sort!(collect(size_psd), rev=true))
            isempty(size_lin) || print(io, "\n  Free block sizes:\n    ", sort!(collect(size_lin), rev=true))
        end
    else
        print(io, "\nVariable cliques:")
        for va in groups.var_cliques
            print(io, "\n  ", join(va, ", "))
        end
        bs = StatsBase.countmap(length.(groups.obj))
        for constr in groups.nonnegs
            mergewith!(+, bs, StatsBase.countmap(length.(constr)))
        end
        for (constr, constr_mat) in zip(groups.psds, poly_problem(x).constr_psd)
            mergewith!(+, bs, StatsBase.countmap(length.(constr) .* size(constr_mat, 1)))
        end
        print(io, "\nPSD block sizes:\n  ", sort!(collect(bs), rev=true))
        if !isempty(groups.zeros)
            empty!(bs)
            for constr in groups.zeros
                mergewith!(+, bs, StatsBase.countmap(length.(constr)))
            end
            print(io, "\nFree block sizes:\n  ", sort!(collect(bs), rev=true))
        end
    end
    return
end

Base.show(io::IO, m::MIME"text/plain", x::AbstractRelaxation) = _show(io, m, x)

# make working with the relaxation as simple as working with the problem itself
Base.getproperty(relaxation::R, f::Symbol) where {R<:AbstractRelaxation} =
    if hasfield(R, f)
        getfield(relaxation, f)
    elseif hasfield(R, :parent)
        getproperty(getfield(relaxation, :parent), f)
    else
        getproperty(getfield(relaxation, :problem), f)
    end
Base.propertynames(relaxation::R) where {P<:Problem,R<:AbstractRelaxation{P}} =
    if hasfield(R, :parent)
        (fieldnames(R)..., propertynames(fieldtype(R, :parent))...)
    else
        (fieldnames(R)..., fieldnames(P)...)
    end
MultivariatePolynomials.variables(relaxation::AbstractRelaxation) = variables(relaxation.problem)
MultivariatePolynomials.nvariables(relaxation::AbstractRelaxation) = nvariables(relaxation.problem)
"""
    degree(problem::AbstractRelaxation)

Returns the degree associated with the relaxation of a polynomial optimization problem.

See also [`poly_problem`](@ref).
"""
function MultivariatePolynomials.degree(relaxation::AbstractRelaxation)
    gr = groupings(relaxation)
    subdegree = v -> maximum(maxdegree_complex, v, init=0)
    return max(
        maximum(maxdegree_complex, gr.obj, init=0),
        maximum(subdegree, gr.zeros, init=0),
        maximum(subdegree, gr.nonnegs, init=0),
        maximum(subdegree, gr.psds, init=0)
    )
end
Base.isreal(relaxation::AbstractRelaxation) = isreal(relaxation.problem)

include("./basis/Basis.jl")
include("./sparse/Sparse.jl")

end
