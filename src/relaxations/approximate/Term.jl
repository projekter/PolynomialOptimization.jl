@enum TermMode begin
    TERM_MODE_DENSE = 0
    TERM_MODE_BLOCK = 1
    TERM_MODE_CLIQUES = 2
    TERM_MODE_CHORDAL_CLIQUES = 3
    TERM_MODE_NONE = 4 # only for subselection
end

@doc """
    @enum TermMode TERM_MODE_DENSE TERM_MODE_BLOCK TERM_MODE_CLIQUES
        TERM_MODE_CHORDAL_CLIQUES TERM_MODE_NONE

Specifies which kind of completion procedure is used for the iteration of term sparsity pattern.
Valid values are `TERM_MODE_DENSE` ([`Dense`](@ref)), `TERM_MODE_BLOCK` ([`SparsityTermBlock`](@ref)),
`TERM_MODE_CHORDAL_CLIQUES` ([`SparsityTermChordal`](@ref)), and `TERM_MODE_CLIQUES` ([`SparsityTermChordal`](@ref) with
`chordal_completion = false`).
`TERM_MODE_NONE` can be used during iteration to disable the iteration of individual constraints.
""" TermMode

const VarcliqueMethods = Union{Missing,<:AbstractVector{<:Union{Missing,TermMode}}}

"""
    SparsityTerm

Common base class that term sparsity methods use or wrap. The [`SparsityTermBlock`](@ref) and [`SparsityTermChordal`](@ref)
constructors are shorthands that create `SparsityTerm` objects with the `method` parameter appropriately set.
[`SparsityCorrelativeTerm`](@ref) is a very thin wrapper around `SparsityTerm`.
"""
mutable struct SparsityTerm{
    I<:Integer,
    P<:Problem{<:(IntPolynomial{<:Any,Nr,Nc,<:IntMonomialVector{Nr,Nc,I}} where {Nr,Nc})},
    PG<:RelaxationGroupings,
    U<:AbstractSet{I},
    S<:IntMonomialVector,
    G<:RelaxationGroupings
} <: AbstractRelaxationSparse{P}
    const problem::P
    const parent # no specialization, therefore keep a typed problem field
    const parentgroupings::PG
    support_union::Set{I}
    const localizing_supports::Vector{S}
    const graphs::Vector{Graphs.SimpleGraph{Int}}
    groupings::G
    const method::AbstractVector{TermMode} # no specialization
    const varclique_method::VarcliqueMethods

    function SparsityTerm(relaxation::AbstractRelaxation{P},
        support_union::AbstractSet{I}; method::Union{TermMode,<:AbstractVector{TermMode}},
        varclique_method::VarcliqueMethods=missing, verbose::Bool=false) where
        {Nr,Nc,I<:Integer,MV<:IntMonomialVector{Nr,Nc,I},P<:Problem{<:IntPolynomial{<:Any,Nr,Nc,MV}}}
        problem = poly_problem(relaxation)
        parent = groupings(relaxation)
        tot = 1 + length(problem.constr_zero) + length(problem.constr_nonneg) + length(problem.constr_psd)
        method isa AbstractVector && length(method) != tot &&
            throw(ArgumentError("The number of term modes specified is incompatible with the number of constraints"))
        (method === TERM_MODE_NONE || (method isa AbstractVector && TERM_MODE_NONE ∈ method)) &&
            throw(ArgumentError("TERM_MODE_NONE is only allowed during iteration"))
        !ismissing(varclique_method) && length(varclique_method) != length(parent.var_cliques) &&
            throw(ArgumentError("The number of clique term modes specified is incompatible with the number of variable cliques"))
        @verbose_info("Generating localizing supports")
        localizing_supports = Vector{MV}(undef, tot)
        @inbounds localizing_supports[1] = monomials(one(problem.objective))
        i = 2
        for constrs in (problem.constr_zero, problem.constr_nonneg, problem.constr_psd)
            for constr in constrs
                @inbounds localizing_supports[i] = monomials(constr)
                i += 1
            end
        end
        methods = method isa TermMode ? ConstantVector(method, tot) : method
        @verbose_info("Converting supports into graphs")
        graphtime = @elapsed begin
            graphs = Vector{Graphs.SimpleGraph{Int}}(undef, length(parent.obj) + sum(length, parent.zeros, init=0) +
                sum(length, parent.nonnegs, init=0) + sum(length, parent.psds, init=0))
            _supports_to_graphs!(graphs, support_union, localizing_supports, parent, methods, varclique_method)
        end
        @verbose_info("Obtained graphs in ", graphtime, " seconds. Generating groupings.")
        grouptime = @elapsed begin
            newgroupings = _extend_graphs!(parent, parent, graphs, methods, varclique_method)
        end
        @verbose_info("Generated new groupings in ", grouptime, " seconds; embedding with old.")
        intersecttime = @elapsed begin
            gr = embed(newgroupings, parent, relaxation isa AbstractRelaxationBasis)
        end
        @verbose_info("Obtained embedding in ", intersecttime, " seconds")

        new{I,P,typeof(parent),typeof(support_union),eltype(localizing_supports),typeof(gr)}(
            problem, relaxation, parent, support_union, localizing_supports, graphs, gr, methods, varclique_method
        )
    end
end

# Only a helper so that we can cheaply obtain the monomial union without having to construct the square basis.
struct _SquareBasis{M<:IntMonomial,B<:AbstractVector{M}}
    parent::B
end

Base.IteratorSize(::Type{<:_SquareBasis{M,B}}) where {M<:IntMonomial,B<:AbstractVector{M}} = Base.IteratorSize(B)
Base.IteratorEltype(::Type{<:_SquareBasis{M}}) where {M} = Base.HasEltype()
Base.eltype(::Type{<:_SquareBasis{M}}) where {M} = IntPolynomials._get_I(M)
Base.length(b::_SquareBasis{<:IntMonomial{<:Any,0}}) = length(b.parent)
Base.length(b::_SquareBasis) = count(isreal, b.parent) # assuming that no conjugates occur, a monomial with complex variables
                                                       # cannot be real-valued
function Base.iterate(b::_SquareBasis{M}, state=missing) where {I<:Integer,M<:IntMonomial{<:Any,<:Any,I}}
    result = ismissing(state) ? iterate(b.parent) : iterate(b.parent, state)
    local parent_mon, parent_state
    while true
        isnothing(result) && return nothing
        parent_mon, parent_state = result
        real_mon = true
        for (v, _) in parent_mon # faster than isreal(parent_mon), as we can assume that no conjugates are present
            if !isreal(v)
                real_mon = false
                break
            end
        end
        real_mon && break
        result = iterate(b.parent, parent_state)
    end
    return monomial_index(parent_mon, parent_mon), parent_state
end

"""
    SparsityTerm(relaxation::AbstractRelaxation; method, varclique_method=missing, verbose=false)

Low-level constructor for `SparsityTerm` objects that provides more nuanced control over which methods are used. `method` must
either be a valid [`TermMode`](@ref) (`TERM_MODE_NONE` is forbidden) or a vector or term modes. In the latter case, the vector
must have exactly as many elements are there are constraints +1. The first element corresponds to the term mode applied to the
objective, the second to the first zero constraints, ..., followed by the nonnegative and psd constraints.
If variable cliques are present, different methods can be assigned to them. Note that a variable clique can cover the objective
and all constraints; in the case of conflicting assignments, the clique assignment takes precedence (but a clique mode may also
be `missing` individually, in which case the default is taken).

See also [`SparsityTermBlock`](@ref), [`SparsityTermChordal`](@ref), [`SparsityCorrelativeTerm`](@ref).
"""
function SparsityTerm(relaxation::AbstractRelaxation{P}; method::Union{TermMode,<:AbstractVector{TermMode}},
    varclique_method::VarcliqueMethods=missing, verbose::Bool=false) where
    {Nr,Nc,I<:Integer,P<:Problem{<:IntPolynomial{<:Any,Nr,Nc,<:IntMonomialVector{Nr,Nc,I}}}}
    # Let 𝒜 = supp(obj) ∪ supp(constrs)
    # support_union corresponds to 𝒮. For real-valued problems, the initialization is
    # 𝒮 = 𝒜 ∪ 2grouping(obj), for complex-valued problems, only 𝒜. There is no good explanation given for why we even use
    # 2grouping(obj) - it seems to be good for convergence, it might be related to the Hankel structure... So what to do with
    # mixed problems? Two immediate possibilities:
    # - form the union with all monomials that are squared in the real variables, provided they have at least one real variable
    #   present
    # - same, but only if the monomials don't have a complex variable
    # Here, we choose the second approach (just because we have to choose).
    @verbose_info("Generating initial support union")
    prob = poly_problem(relaxation)
    parent = groupings(relaxation)
    supptime = @elapsed begin
        # When a Problem is constructed using poly_problem, all polynomials are converted from MP interface; therefore, they
        # will all internally use ExponentsAll as the exponents backend. We can therefore directly access the index.
        # However, just to be on the very safe side (maybe someone somehow constructs a Problem directly, or we want to change
        # the implementation later?), we don't rely on this and use monomial_index instead. This will convert the index to
        # ExponentsAll, but has a short path for recognizing when it already is given in this way, so it can be optimized to
        # just grab the field value if everything is as it should be.
        support_union = Set{I}(Iterators.map(monomial_index, monomials(prob.objective)))
        union!(support_union, Iterators.map(monomial_index, monomials(prob.prefactor)))
        # ^ this is for the minimal value that is subtracted from the objective
        @unroll for constrs in (prob.constr_zero, prob.constr_nonneg, prob.constr_psd)
            for constr in constrs
                union!(support_union, Iterators.map(monomial_index, monomials(constr)))
            end
        end
        if !iszero(Nr)
            for grouping in parent.obj
                union!(support_union, _SquareBasis(grouping))
            end
        end
    end
    @verbose_info("Generated support union in ", supptime, " seconds")
    return SparsityTerm(relaxation, support_union; method, varclique_method, verbose)
end

function _supports_to_graphs!(graphs::Vector{Graphs.SimpleGraph{Int}}, support_union::AbstractSet{I},
    localizing_supports::Vector{MV}, parent::RelaxationGroupings, methods::AbstractVector{TermMode},
    varclique_methods::Union{Missing,<:AbstractVector{Union{TermMode,Missing}}}) where {I<:Integer,MV<:IntMonomialVector}
    # enforce specialization on grouping, which is an abstract type
    c = Channel{Tuple{Any,MV,Int,Int}}(0, spawn=true) do ch
        ipoly = 1
        igroup = 1
        @unroll for constrs in ((parent.obj,), parent.zeros, parent.nonnegs, parent.psds)
            for constr_groupings in constrs
                if methods[ipoly] != TERM_MODE_NONE || !ismissing(varclique_methods)
                    localizing_support = localizing_supports[ipoly]
                    for grouping in constr_groupings
                        put!(ch, (grouping, localizing_support, ipoly, igroup))
                        igroup += 1
                    end
                else
                    igroup += length(constr_groupings)
                end
                ipoly += 1
            end
        end
    end
    Threads.threading_run(_ -> begin
        for (grouping, localizing_support, ipoly, igroup) in c
            if methods[ipoly] == TERM_MODE_NONE
                vcm = varclique_methods[_findclique(grouping, parent.var_cliques)]
                if ismissing(vcm) || vcm == TERM_MODE_NONE
                    return
                end
            end
            graphs[igroup] = graph = Graphs.SimpleGraph(length(grouping))
            for (exp2, g₂) in enumerate(grouping)
                r = isreal(g₂)
                for (exp1, g₁) in enumerate(Iterators.take(grouping, r ? exp2 -1 : typemax(Int)))
                    Graphs.has_edge(graph, exp1, exp2) && continue
                    for supp_constr in localizing_support
                        if monomial_index(g₁, supp_constr, IntConjMonomial(g₂)) ∈ support_union
                            Graphs.add_edge!(graph, Graphs.Edge(exp1, exp2))
                            break
                        end
                    end
                end
            end
            nothing
        end
    end, false)
    return
end
_supports_to_graphs!(relaxation::SparsityTerm, methods, varclique_methods) = _supports_to_graphs!(relaxation.graphs,
    relaxation.support_union, relaxation.localizing_supports, relaxation.parentgroupings, methods, varclique_methods)

function _extend_graphs!(::Val{TERM_MODE_BLOCK}, g::Graphs.SimpleGraph)
    ne = 0
    # connected components are not overlapping
    connections = Graphs.connected_components(g)
    @inbounds for cc in connections
        @assert(issorted(cc))
        for c in cc
            if length(g.fadjlist[c]) != length(cc)
                g.fadjlist[c] = copy(cc) # this will include self-loops, but this is easier and doesn't bother us
            end
        end
        # the edges are bidirectional, but count only once
        ne += (length(cc) * (length(cc) +1)) >> 1
    end
    g.ne = ne
    return connections
end

function _extend_graphs!(::Val{TERM_MODE_CLIQUES}, g::Graphs.SimpleGraph{Int})
    # cliques can be overlapping
    cliques = sort!(sort!.(Graphs.maximal_cliques(g)))
    edge_sets = [Set{Int}() for _ in 1:Graphs.nv(g)]
    @inbounds for mc in cliques
        for c in mc
            union!(edge_sets[c], mc)
        end
    end
    ne = 0
    @inbounds for (i, edges) in enumerate(edge_sets)
        g.fadjlist[i] = sort!(collect(edges))
        ne += (length(edges) * (length(edges) +1)) >> 1
    end
    g.ne = ne
    return cliques
end

_extend_graphs!(::Val{TERM_MODE_CHORDAL_CLIQUES}, g::Graphs.SimpleGraph) = chordal_cliques!(g)

function _extend_graphs!(::Val{TERM_MODE_DENSE}, g::Graphs.SimpleGraph{T}) where {T<:Integer}
    # Make the graph into a complete one. As we know that every modification will always be towards more compelte graphs, we
    # can share the memory for every entry in the adjacency list.
    edges = collect(one(T):T(Graphs.ne(g)))
    fill!(g.fadjlist, edges)
    return [edges]
end

@eval function _extend_graphs!(previous::RelaxationGroupings{Nr,Nc,I}, parent::RelaxationGroupings{Nr,Nc,I}, g::Vector{G},
    methods::AbstractVector{TermMode}, varclique_methods::VarcliqueMethods) where {Nr,Nc,I<:Integer,G<:Graphs.SimpleGraph}
    # better enable bounds checking in this method
    ipoly = 1
    igroup = 1
    $((let
        body = quote
            $(name === :obj ? :newobj : :($(Symbol(:new, name))[i])) = if methods[ipoly] != TERM_MODE_NONE
                newgroup = FastVec{IntMonomialVector{Nr,Nc,I}}()
                for grouping in constr_groupings
                    # while grouping is generic, we don't gain a lot on inference by putting the loop body in a function
                    # barrier
                    if !ismissing(varclique_methods)
                        vcm = varclique_methods[_findclique(grouping, parent.var_cliques)]
                        method = ismissing(vcm) ? methods[ipoly] : vcm
                    else
                        method = methods[ipoly]
                    end
                    cliques = _extend_graphs!(Val(method), g[igroup])
                    for clique in cliques
                        ($(name === :obj) && isone(length(clique)) && isconstant(grouping[first(clique)])) ||
                            push!(newgroup, @view(grouping[clique]))
                    end
                    igroup += 1
                end
                finish!(newgroup)
            elseif ismissing(varclique_methods)
                igroup += length(constr_groupings)
                $(name === :obj ? :(previous.obj) : :(previous.$name[i]))
            else
                newgroup = FastVec{IntMonomialVector{Nr,Nc,I}}()
                missingclique = false
                nonmissingclique = false
                for grouping in constr_groupings
                    vcm = varclique_methods[_findclique(grouping, parent.var_cliques)]
                    missingclique |= ismissing(vcm)
                    nonmissingclique |= !ismissing(vcm)
                    if !ismissing(vcm)
                        cliques = _extend_graphs!(Val(vcm), g[igroup])
                        for clique in cliques
                            ($(name === :obj) && isone(length(clique)) && isconstant(grouping[first(clique)])) ||
                                push!(newgroup, @view(grouping[clique]))
                        end
                    end
                    igroup += 1
                end
                if missingclique
                    # This is an implementation limitation: as we do not know which method was applied in the previous step, we
                    # cannot reconstruct the cliques from the graph alone
                    nonmissingclique &&
                        throw(ArgumentError("Partially missing clique methods cannot be mixed with no known objective method"))
                    $(name === :obj ? :(previous.obj) : :(previous.$name[i]))
                else
                    finish!(newgroup)
                end
            end
            ipoly += 1
        end
        if name === :obj
            Expr(:block,
                :(constr_groupings = parent.obj),
                body
            )
        else
            Expr(:block,
                Expr(:(=), Symbol(:new, name),:(Vector{Vector{IntMonomialVector{Nr,Nc,I}}}(undef, length(parent.$name)))),
                :(for (i, constr_groupings) in enumerate(parent.$name); $body end)
            )
        end
    end for name in (:obj, :zeros, :nonnegs, :psds))...)

    return RelaxationGroupings(newobj, newzeros, newnonnegs, newpsds, previous.var_cliques)
end

function _extend_graphs!(relaxation::SparsityTerm, methods::AbstractVector{TermMode}, varclique_methods::VarcliqueMethods)
    newgroupings = embed(_extend_graphs!(groupings(relaxation), relaxation.parentgroupings, relaxation.graphs, methods,
        varclique_methods), relaxation.parentgroupings, relaxation.parent isa AbstractRelaxationBasis)
    if newgroupings != relaxation.groupings
        relaxation.groupings = newgroupings
        return relaxation
    else
        return nothing
    end
end

function _iterate_supports(parent::RelaxationGroupings, localizing_supports::Vector{MV}, g::Vector{G}) where
    {I<:Integer,Nc,MV<:IntMonomialVector{<:Any,Nc,I},G<:Graphs.SimpleGraph}
    support_union = Set{I}()
    grouping_loop = @capture((grouping, localizing_support, graph) -> begin
        for e in Graphs.edges(graph)
            for α in localizing_support
                push!($support_union, monomial_index(grouping[Graphs.src(e)], α,
                                                     IntConjMonomial(grouping[Graphs.dst(e)])))
                if !iszero(Nc) && Graphs.src(e) != Graphs.dst(e) # or noncommutative
                    push!(support_union, monomial_index(grouping[Graphs.dst(e)], α,
                          IntConjMonomial(grouping[Graphs.src(e)])))
                end
            end
        end
        # TODO: check the complex case
        nothing
    end)
    ipoly = 1
    igroup = 1
    @inbounds @unroll for constrs in ((parent.obj,), parent.zeros, parent.nonnegs, parent.psds)
        for constr_groupings in constrs
            localizing_support = localizing_supports[ipoly]
            for grouping in constr_groupings
                graph = g[igroup]
                sizehint!(support_union, length(support_union) +
                    (iszero(Nc) ? 1 : 2) * Graphs.ne(graph) * length(localizing_support))
                grouping_loop(grouping, localizing_support, graph) # here a function barrier is good
                igroup += 1
            end
            ipoly += 1
        end
    end
    return support_union
end
_iterate_supports!(relaxation::SparsityTerm) = relaxation.support_union =
    _iterate_supports(relaxation.parentgroupings, relaxation.localizing_supports, relaxation.graphs)

"""
    SparsityTermBlock(relaxation::AbstractProblem; verbose::Bool=false)

Analyze the term sparsity of the problem.
[Term sparsity](https://doi.org/10.1137/19M1307871) is a recent iterative sparsity analysis that groups terms with shared
supports. Its last iteration will give the same optimal value as the original problem, although it may still be of a smaller
size. Often, even the uniterated analysis already gives the same bound as the dense problem.
The terms are grouped based on connected components of a graph; this can be improved by using the smallest chordal extension
(see [`SparsityTermChordal`](@ref)), which will lead to even smaller problem sizes, but typically also worse bounds.

If correlative and term sparsity are to be used together, use [`SparsityCorrelativeTerm`](@ref) or nest the sparsity objects.
"""
SparsityTermBlock(args...; kwargs...) = SparsityTerm(args...; method=TERM_MODE_BLOCK, kwargs...)

"""
    SparsityTermChordal(relaxation::AbstractProblem; chordal_completion=true, verbose=false)

Analyze the term sparsity of the problem using chordal cliques.
[Chordal term sparsity](https://doi.org/10.1137/20M1323564) is a recent iterative sparsity analysis that groups terms
with shared supports. Even in its last iteration, it may give strictly smaller values than the dense problem.
The basis elements are grouped in terms of chordal cliques of the term sparsity graph. This uses maximal cliques; as obtaining
maximal cliques of an arbitrary graph is not efficient, the graph is extended to a chordal graph if `chordal_completion` is
`true` using a heuristic. Disabling the chordal completion can lead to smaller problem sizes.

If correlative and term sparsity are to be used together, use [`SparsityCorrelativeTerm`](@ref) or nest the sparsity objects.
"""
SparsityTermChordal(args...; chordal_completion::Bool=true, kwargs...) =
    SparsityTerm(args...; method=chordal_completion ? TERM_MODE_CHORDAL_CLIQUES : TERM_MODE_CLIQUES, kwargs...)

const ModeSpecifiers = Union{Bool,TermMode,<:AbstractSet{<:Integer},<:AbstractVector{TermMode}}

function _jointmethods(problem::Problem, objective::Union{Bool,TermMode}, zero::ModeSpecifiers, nonneg::ModeSpecifiers,
    psd::ModeSpecifiers, basemethod::Union{TermMode,<:AbstractVector{TermMode}})
    len = 1 + length(problem.constr_zero) + length(problem.constr_nonneg) + length(problem.constr_psd)
    basemethod isa AbstractVector && @assert(length(basemethod) == len)
    if objective isa TermMode || zero isa AbstractVector || zero isa AbstractSet || nonneg isa AbstractVector ||
        nonneg isa AbstractSet || psd isa AbstractVector || psd isa AbstractSet
        if basemethod isa AbstractVector
            methods = collect(basemethod) # not copy - if it is a ConstantVector, we need to make it a Vector
        else
            methods = fill(basemethod, len)
        end
    elseif basemethod isa AbstractVector
        methods = copy(basemethod)
    else
        methods = ConstantVector(basemethod, len)
    end
    @inbounds if objective isa TermMode
        methods[1] = objective
    elseif objective === false
        methods[1] = TERM_MODE_NONE
    elseif basemethod === TERM_MODE_NONE
        throw(ArgumentError("Default mode TERM_MODE_NONE cannot be combined with `true` iteration value"))
    end
    if zero !== true
        @inbounds for i in 2:length(problem.constr_zero)+1
            if zero isa TermMode
                methods[i] = zero
            elseif zero isa Bool || !(i -1 ∈ zero)
                methods[i] = TERM_MODE_NONE
            end
        end
    elseif basemethod === TERM_MODE_NONE
        throw(ArgumentError("Default mode TERM_MODE_NONE cannot be combined with `true` iteration value"))
    end
    δ = 1 + length(problem.constr_zero)
    if nonneg !== true
        @inbounds for i in δ+1:length(problem.constr_nonneg)+δ
            if nonneg isa TermMode
                methods[i] = nonneg
            elseif nonneg isa Bool || !(i - δ ∈ nonneg)
                methods[i] = TERM_MODE_NONE
            end
        end
    elseif basemethod === TERM_MODE_NONE
        throw(ArgumentError("Default mode TERM_MODE_NONE cannot be combined with `true` iteration value"))
    end
    δ += length(problem.constr_nonneg)
    if psd !== true
        @inbounds for i in δ+1:length(problem.constr_psd)+δ
            if psd isa TermMode
                methods[i] = nonneg
            elseif psd isa Bool || !(i - δ ∈ psd)
                methods[i] = TERM_MODE_NONE
            end
        end
    elseif basemethod === TERM_MODE_NONE
        throw(ArgumentError("Default mode TERM_MODE_NONE cannot be combined with `true` iteration value"))
    end
    return methods
end

"""
    iterate!(relaxation::Union{SparsityTerm,SparsityCorrelativeTerm}; [method,]
        objective=true, zero=true, nonneg=true, psd=true, varclique_methods=missing)

[`SparsityTerm`](@ref) implementations allow to customize the iteration procedure by the keyword arguments. The arguments
`objective`, `zero`, `nonneg`, and `psd` can be boolean values (`false` means that these elements will not contribute to the
iteration, `true` that they will). `zero`, `nonneg`, and `psd` can also be `AbstractSet`s of integers, indicating that only the
constraints with the indices specified in the set will contribute to the iteration. This all implies that the method used for
their iteration will be given by `method`. Custom methods can be assigned if the parameters are set to a [`TermMode`](@ref) or
a vector of `TermMode`s.

The parameter `method` therefore determines the default that is assigned to the elements, and if not specified, it will be
determined by the default method with which `relaxation` was constructed. Instead of a single `TermMode`, `method` may also be
a vector successively assigning modes to the objective, the first zero constraints, ..., the nonnegative constraints, the psd
constraints (eliminating the need for the other keywords).

The parameter `varclique_methods` instead allows to assign custom methods to individual variable cliques. Note that a variable
clique can cover the objective and all constraints; in the case of conflicting assignments, the clique assignment takes
precedence (but a clique mode may also be `missing` individually, in which case the default is taken).
"""
function iterate!(relaxation::SparsityTerm; method::Union{TermMode,<:AbstractVector{TermMode}}=relaxation.method,
    objective::Union{Bool,TermMode}=true, zero::ModeSpecifiers=true, nonneg::ModeSpecifiers=true, psd::ModeSpecifiers=true,
    varclique_methods::VarcliqueMethods=missing)
    problem = poly_problem(relaxation)
    _iterate_supports!(relaxation)
    # Even if new_supports is the same as ts.support_union, we might run into the unlikely scenario that the user changed the
    # type of term sparsity, so that now extending with a different type might lead to a different result despite having the
    # same support.
    methods = _jointmethods(problem, objective, zero, nonneg, psd, method)
    _supports_to_graphs!(relaxation, methods, varclique_methods)
    return _extend_graphs!(relaxation, methods, varclique_methods)
end

default_solution_method(::SparsityTerm) = :heuristic