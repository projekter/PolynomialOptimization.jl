struct SparsityCorrelative{P<:Problem,G<:RelaxationGroupings} <: AbstractRelaxationSparse{P}
    parent::AbstractRelaxation{P} # no specialization
    groupings::G

    @doc """
        SparsityCorrelative(relaxation::AbstractRelaxation; [high_order_zero,]
            [high_order_nonneg,] [high_order_psd,] [low_order_zero,] [low_order_nonneg,]
            [low_order_psd,] chordal_completion=CliqueTrees.MF(),
            [split_psd,] [join_psd,] verbose::Bool=false)

    Analyze the correlative sparsity of a problem.
    Correlative sparsity is a variable-based sparsity analysis. It was first defined by
    [Waki et al.](https://doi.org/10.1137/050623802) in 2006 and extended by
    [Josz and Molzahn](https://doi.org/10.1137/15M1034386) in 2018.
    Variables are grouped into cliques based on the terms in the objective in which they appear together.
    Additional grouping is induced by variables that occur anywhere in a constraint of high order; or by variables that occur
    in a term in a constraint of low order.
    The parameters `high_order_...` allow to specify which constraints - identified by their indices - are of high order. If
    the parameter is omitted, all such constraints are of high order. Conversely, `low_order_...` can be used to specify that
    all _but_ the listed constraints are of high order. Both parameters cannot be used simultaneously for the same set of
    constraints.
    Note that the order of the constraints is also influenced by the parent relaxation. If a correlative sparsity relaxation is
    applied to another relaxation that already limited the prefactor of a constraint to be of degree zero, it must necessarily
    be of low order.
    A PSD constraint is by default handled as one monolithic constraint. However, the groupings can also be constructed
    separately per index; then, only variables that occur in the kᵗʰ row or kᵗʰ column are put into a clique. This corresponds
    to setting the `split_psd` parameters to the indices of the constraints that should be separated in this way; the opposite
    effect (default) is achieved by `join_psd`.
    All values may also be `true` (to indicate all constraints) or `false` (to indicate none)

    By default, the correlative sparsity graph is completed to a chordal graph before the cliques are determined, which
    guarantees that the maximal cliques can be determined quickly; however, this may degrade the sparsity and it may be
    favorable not to carry out the completion. The argument `chordal_completion` may be set to `false` for this; any of the
    allowed elimination algorithms from `CliqueTrees.jl` are other possible arguments. `true` (deprecated) defaults to `MF()`.

    If correlative and term sparsity are to be used together, use [`SparsityCorrelativeTerm`](@ref) instead of
    nesting the sparsity objects.
    """
    function SparsityCorrelative(relaxation::AbstractRelaxation{P}; high_order_zero=missing,
        high_order_nonneg=missing, high_order_psd=missing, low_order_zero=missing, low_order_nonneg=missing,
        low_order_psd=missing, chordal_completion::Union{Bool,<:CliqueTrees.EliminationAlgorithm}=CliqueTrees.MF(),
        split_psd=missing, join_psd=missing, verbose::Bool=false) where # sync with SparsityCorrelativeTerm
        {Nr,Nc,I<:Integer,P<:Problem{<:IntPolynomial{<:Any,Nr,Nc,<:IntMonomialVector{Nr,Nc,I}}}}
        ((!ismissing(high_order_zero) && !ismissing(low_order_zero)) ||
            (!ismissing(high_order_nonneg) && !ismissing(low_order_nonneg)) ||
            (!ismissing(high_order_psd) && !ismissing(low_order_psd))) &&
            error("high_order_... and low_order_... specifications cannot be used at the same time")
        !ismissing(split_psd) && !ismissing(join_psd) &&
            error("split_psd and join_psd specifications cannot be used at the same time")
        problem = poly_problem(relaxation)
        if high_order_zero isa Bool || low_order_zero isa Bool
            high_order_zero = 1:(high_order_zero === true || low_order_zero === false ? length(problem.constr_zero) : 0)
            low_order_zero = missing
        end
        if high_order_nonneg isa Bool || low_order_nonneg isa Bool
            high_order_nonneg = 1:(high_order_nonneg === true || low_order_nonneg === false ?
                                   length(problem.constr_nonneg) : 0)
            low_order_nonneg = missing
        end
        if high_order_psd isa Bool || low_order_psd isa Bool
            high_order_psd = 1:(high_order_psd === true || low_order_psd === false ? length(problem.constr_psd) : 0)
            low_order_psd = missing
        end
        if split_psd isa Bool || join_psd isa Bool
            split_psd = 1:(split_psd === true || join_psd === false ? length(problem.constr_psd) : 0)
            join_psd = missing
        end
        ((ismissing(high_order_zero) || high_order_zero ⊆ 1:length(problem.constr_zero)) &&
            (ismissing(low_order_zero) || low_order_zero ⊆ 1:length(problem.constr_zero)) &&
            (ismissing(high_order_nonneg) || high_order_nonneg ⊆ 1:length(problem.constr_nonneg)) &&
            (ismissing(low_order_nonneg) || low_order_nonneg ⊆ 1:length(problem.constr_nonneg)) &&
            (ismissing(high_order_psd) || high_order_psd ⊆ 1:length(problem.constr_psd)) &&
            (ismissing(low_order_psd) || low_order_psd ⊆ 1:length(problem.constr_psd)) &&
            (ismissing(split_psd) || split_psd ⊆ 1:length(problem.constr_psd)) &&
            (ismissing(join_psd) || join_psd ⊆ 1:length(problem.constr_psd))) ||
            error("Unknown constraint index specified")

        parent = groupings(relaxation)
        parentmaxobjdeg = maximum(maxdegree, parent.obj, init=0)::Int
        parentmaxzerodeg = Int[maximum(maxdegree, z, init=0) for z in parent.zeros]
        parentmaxnonnegdeg = Int[maximum(maxdegree, n, init=0) for n in parent.nonnegs]
        parentmaxpsddeg = Int[maximum(pᵢ -> maximum(maxdegree, pᵢ, init=0), p, init=0) for p in parent.psds]

        @verbose_info("Constructing correlative sparsity graph")
        g = Graphs.SimpleGraph(Nr + Nc)
        # objective: check which pairs of exponents appear together in a term
        # effective objective due to polynomial optimization: it is actually checked whether
        # prefactor * (obj - lowerbound) is SOS. The prefactor * obj multiplication was already done during problem
        # construction, but prefactor * lowerbound is still missing.
        for poly in (problem.objective, problem.prefactor)
            for term in poly
                mon = monomial(term)
                for (i, (var1, _)) in Iterators.drop(enumerate(mon), 1)
                    var1i = variable_index(var1)
                    for (var2, _) in Iterators.take(mon, i -1)
                        Graphs.add_edge!(g, Graphs.Edge(var1i, variable_index(var2)))
                    end
                end
            end
        end
        # This parent grouping will put some constraints on the form of the prefactor(s). With correlative sparsity, we
        # additionally enforce that no variable occurs in a grouping that is not already present in the constraint.
        _correlative_addedges!(g, problem.constr_zero, parent.zeros, high_order_zero, low_order_zero, parentmaxzerodeg)
        _correlative_addedges!(g, problem.constr_nonneg, parent.nonnegs, high_order_nonneg, low_order_nonneg,
            parentmaxnonnegdeg)
        _correlative_addedges!(g, problem.constr_psd, parent.psds, high_order_psd, low_order_psd, parentmaxpsddeg, split_psd,
            join_psd)

        @verbose_info("Determining cliques")
        if chordal_completion === true
            chordal_completion = CliqueTrees.MF()
        end
        gentime = @elapsed(cliques = sort!(chordal_completion !== false ? chordal_cliques(g, alg=chordal_completion) :
                                                                          Graphs.maximal_cliques(g), by=length))
        sort!.(cliques)
        if Nc > 1
            # The graph needs continuous indices; but the complex variables are not indexed contiuously. Normal and conjugated
            # variables alternate. The first complex variable doesn't need to be touched.
            for clique in cliques
                @inbounds for i in searchsortedfirst(clique, Nr +2):length(clique)
                    clique[i] = ((clique[i] - Nr) << 1) + Nr -1
                end
            end
        end
        @verbose_info("Obtained ", length(cliques), " clique", length(cliques) > 1 ? "s" : "", " in ", gentime,
            " seconds. Generating groupings.")
        # The correlative iterator could potentially be made even smaller by determining all the multideg boundaries, but would
        # this be worth the effort?
        minmultideg = IntPolynomials.ConstantVector(0, Nr + 2Nc)
        newobj = Vector{IntMonomialVector{Nr,Nc,I}}(undef, length(cliques))
        newzero = [IntMonomialVector{Nr,Nc,I}[] for _ in 1:length(problem.constr_zero)]
        newnonneg = [IntMonomialVector{Nr,Nc,I}[] for _ in 1:length(problem.constr_nonneg)]
        newpsd = [AbstractVector{IntMonomialVector{Nr,Nc,I}}[] for _ in 1:length(problem.constr_psd)]
        @inbounds for (i, clique) in enumerate(cliques)
            maxmultideg = zeros(Int, Nr + 2Nc)
            fill!(@view(maxmultideg[clique]), parentmaxobjdeg)
            newobj[i] = IntMonomialVector{Nr,Nc}(ExponentsMultideg{Nr+2Nc,I}(0, parentmaxobjdeg, minmultideg, maxmultideg))
        end
        _correlative_gengroupings!(newzero, cliques, problem.constr_zero, parentmaxzerodeg)
        _correlative_gengroupings!(newnonneg, cliques, problem.constr_nonneg, parentmaxnonnegdeg)
        _correlative_gengroupings!(newpsd, cliques, problem.constr_psd, parentmaxpsddeg, split_psd, join_psd)
        @verbose_info("Generated new groupings; embedding in old.")
        gentime = @elapsed(gr = embed!(
            RelaxationGroupings(newobj, newzero, newnonneg, newpsd, map.(IntVariable{Nr,Nc}, cliques)),
            parent,
            relaxation isa AbstractRelaxationBasis
        ))
        @verbose_info("Obtained embedding in ", gentime, " seconds")

        return new{P,typeof(gr)}(relaxation, gr)
    end
end

function _correlative_addedges!(g::Graphs.SimpleGraph, constrs,
    parentgroupings::Union{Vector{Vector{IntMonomialVector{Nr,Nc,I}}},
                           Vector{Vector{AbstractVector{IntMonomialVector{Nr,Nc,I}}}}} where {Nr,Nc,I<:Integer}, h, l, md,
    splits=nothing, joins=nothing)
    @assert isnothing(splits) == isnothing(joins)
    @assert isnothing(splits) == (parentgroupings isa (Vector{Vector{IntMonomialVector{Nr,Nc,I}}} where {Nr,Nc,I<:Integer}))
    if !ismissing(h) && !(h isa AbstractSet)
        h = Set(h)
    end
    if !ismissing(l) && !(l isa AbstractSet)
        l = Set(l)
    end
    for (i, (constr, groupings)) in enumerate(zip(constrs, parentgroupings))
        split = false
        if !isnothing(splits)
            if !ismissing(splits)
                split = i ∈ splits
            elseif !ismissing(joins)
                split = !(i ∈ joins)
            end
        end

        low_deg = false
        if !ismissing(h)
            low_deg = !(i ∈ h)
        elseif !ismissing(l)
            low_deg = i ∈ l
        end
        if !low_deg && isone(length(groupings))
            low_deg = isnothing(splits) ? isconstant(last(first(groupings))) :
                                          all(∘(isconstant, last), first(groupings))
        end
        if low_deg
            # we must make sure that the grouping only contains the constant, else mixing will occur
            md[i] = 0
            # then we can go back to considering only variables in terms instead of all variables
            for mon in monomials(constr)
                for (i, (var1, _)) in Iterators.drop(enumerate(mon), 1)
                    var1i = variable_index(var1)
                    for (var2, _) in Iterators.take(mon, i -1)
                        Graphs.add_edge!(g, Graphs.Edge(var1i, variable_index(var2)))
                    end
                end
            end
        elseif split
            for constr_col in eachcol(constr) # constr is Hermitian, so no need to also consider eachrow
                vars = effective_variables(constr_col, rettype=Set, by=ordinary_variable)
                @inbounds for (i, var1) in Iterators.drop(enumerate(vars), 1)
                    var1i = variable_index(var1)
                    for var2 in Iterators.take(vars, i -1)
                        Graphs.add_edge!(g, Graphs.Edge(var1i, variable_index(var2)))
                    end
                end
            end
        else
            vars = effective_variables(constr, rettype=Set, by=ordinary_variable)
            @inbounds for (i, var1) in Iterators.drop(enumerate(vars), 1)
                var1i = variable_index(var1)
                for var2 in Iterators.take(vars, i -1)
                    Graphs.add_edge!(g, Graphs.Edge(var1i, variable_index(var2)))
                end
            end
        end
    end
    return g
end

function _correlative_gengroupings!(news::Vector{Vector{IntMonomialVector{Nr,Nc,I}}}, cliques, constrs, parentdeg) where {Nr,Nc,I<:Integer}
    minmultideg = IntPolynomials.ConstantVector(0, Nr + 2Nc)
    for (constr, maxdeg, newel) in zip(constrs, parentdeg, news)
        constrvars = effective_variables(constr, rettype=Vector,
            by=∘(Base.Fix2(getproperty, :index), ordinary_variable))::Vector{IntPolynomials.smallest_unsigned(Nr+2Nc)}
        @assert(issorted(constrvars))
        # There will be at least one clique that contains all the variables in the constraint; take the smallest one...
        cliqueᵢ = findfirst(Base.Fix1(issubset_sorted, constrvars), cliques)
        @inbounds if !isnothing(cliqueᵢ)
            clique = cliques[cliqueᵢ]
            maxmultideg = zeros(Int, Nr + 2Nc)
            fill!(@view(maxmultideg[clique]), maxdeg)
            push!(newel,
                IntMonomialVector{Nr,Nc}(
                    ExponentsMultideg{Nr+2Nc,I}(0, min(maxdeg, length(clique) * maxdeg),
                    minmultideg, maxmultideg)
                )
            )
        else
            # ...unless the prefactor was of low order, as then we didn't necessarily introduce all the couplings.
            # But low degree = prefactor is 1
            iszero(maxdeg) || error("Something is wrong with graph theory")
            push!(newel, IntMonomialVector{Nr,Nc}(ExponentsDegree{Nr+2Nc,I}(0, 0)))
        end
    end
end

function _correlative_gengroupings!(news::Vector{Vector{AbstractVector{IntMonomialVector{Nr,Nc,I}}}}, cliques, constrs,
    parentdeg, splits, joins) where {Nr,Nc,I<:Integer}
    minmultideg = IntPolynomials.ConstantVector(0, Nr + 2Nc)
    for (i, (constr, maxdeg, newel)) in enumerate(zip(constrs, parentdeg, news))
        split = false
        if !isnothing(splits)
            if !ismissing(splits)
                split = i ∈ splits
            elseif !ismissing(joins)
                split = !(i ∈ joins)
            end
        end

        if split
            newsub = FastVec{IntMonomialVector{Nr,Nc,I}}(buffer=size(constr, 2))
            for constr_col in eachcol(constr)
                constrvars = effective_variables(constr_col, rettype=Vector,
                    by=∘(Base.Fix2(getproperty, :index), ordinary_variable))::Vector{IntPolynomials.smallest_unsigned(Nr+2Nc)}
                @assert(issorted(constrvars))
                # There will be at least one clique that contains all the variables in the constraint; take the smallest one...
                cliqueᵢ = findfirst(Base.Fix1(issubset_sorted, constrvars), cliques)
                @inbounds if !isnothing(cliqueᵢ)
                    clique = cliques[cliqueᵢ]
                    maxmultideg = zeros(Int, Nr + 2Nc)
                    fill!(@view(maxmultideg[clique]), maxdeg)
                    unsafe_push!(newsub,
                        IntMonomialVector{Nr,Nc}(
                            ExponentsMultideg{Nr+2Nc,I}(0, min(maxdeg, length(clique) * maxdeg),
                            minmultideg, maxmultideg)
                        )
                    )
                else
                    # ...unless the prefactor was of low order, as then we didn't necessarily introduce all the couplings.
                    # But low degree = prefactor is 1
                    iszero(maxdeg) || error("Something is wrong with graph theory")
                    push!(newsub, IntMonomialVector{Nr,Nc}(ExponentsDegree{Nr+2Nc,I}(0, 0)))
                end
            end
            push!(newel, allequal(newsub) ? ConstantVector{IntMonomialVector{Nr,Nc,I}}(first(newsub), length(newsub)) :
                                            finish!(newsub))
        else
            constrvars = effective_variables(constr, rettype=Vector,
                by=∘(Base.Fix2(getproperty, :index), ordinary_variable))::Vector{IntPolynomials.smallest_unsigned(Nr+2Nc)}
            @assert(issorted(constrvars))
            # There will be at least one clique that contains all the variables in the constraint; take the smallest one...
            cliqueᵢ = findfirst(Base.Fix1(issubset_sorted, constrvars), cliques)
            @inbounds if !isnothing(cliqueᵢ)
                clique = cliques[cliqueᵢ]
                maxmultideg = zeros(Int, Nr + 2Nc)
                fill!(@view(maxmultideg[clique]), maxdeg)
                push!(newel,
                    ConstantVector{IntMonomialVector{Nr,Nc,I}}(IntMonomialVector{Nr,Nc}(
                        ExponentsMultideg{Nr+2Nc,I}(0, min(maxdeg, length(clique) * maxdeg),
                        minmultideg, maxmultideg)
                    ), size(constr, 2))
                )
            else
                # ...unless the prefactor was of low order, as then we didn't necessarily introduce all the couplings.
                # But low degree = prefactor is 1
                iszero(maxdeg) || error("Something is wrong with graph theory")
                push!(newel, ConstantVector{IntMonomialVector{Nr,Nc,I}}(
                    IntMonomialVector{Nr,Nc}(ExponentsDegree{Nr+2Nc,I}(0, 0)),
                    size(constr, 2)
                ))
            end
        end
    end
end

default_solution_method(::SparsityCorrelative) = :mvhankel