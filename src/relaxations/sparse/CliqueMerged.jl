function merge_cliques!(_cliques::AbstractVector{<:AbstractSet{T}}) where {T}
    # directly drop cliques of length 1 and 2. They are so efficient (linear vs. quadratic constraints) that we don't even
    # consider them in the merge process
    smallcliques = @view _cliques[length.(_cliques).≤2]
    cliques = @view _cliques[length.(_cliques).≥3]
    # first form the clique graph; this time, we work with the adjacency matrix
    n = length(cliques)
    n ≤ 1 && return [smallcliques; cliques]
    @inbounds adjmOwn = collect(@capture(i > j ? length($cliques[i])^3 + length(cliques[j])^3 -
                                                 length(cliques[i] ∪ cliques[j])^3 : 0
                                         for i in 1:n, j in 1:n))::Matrix{Int} # capture n?
    idxOwn = fill(true, n)
    GC.@preserve adjmOwn idxOwn begin
        adjm = unsafe_wrap(Array, pointer(adjmOwn), (n, n), own=false)
        idx = unsafe_wrap(Array, pointer(idxOwn), n, own=false)
        deleted = 0
        @inbounds while true
            # select two permissible cliques with the highest weight
            w, maxidx = findmax(adjm)
            i = maxidx[1]
            j = maxidx[2]
            # while clique graph contains positive weights
            w ≤ 0 && break
            # merge cliques
            union!(cliques[i], cliques[j])
            idx[j] = false
            # In every iteration, we have to search through n^2 elements, we have to reset n-1 elements and also recalculate
            # n-1 elements. The calculation a^3+b^3-c^3 itself consists of 6 multiplications, one addition and one subtraction
            # (and a couple of loads). Multiplication cost is about 3*addition cost, but with out-of-order execution this may
            # be reduced again, but assuming 20 cycles per calcuation (which is also conditional) is fair. Resetting to zero
            # can be done using ymm registers for four items in one step, but vmovups also has a latency of ~4 (varies greatly
            # between the architectures). But it is unrolled with four instructions following each other, which may partially
            # compensate. And we have the counter add, check and jump, so about 4.25(n-1) clock cycles for the reset would be
            # fair. In comparison, the @inbounds findmax requires about 8 clock cycles per element.
            # So the total cost is 8n^2 [search] + 4.25(n-1) [reset] + 20(n-1) [recalculate] for every iteration, until we are
            # done. If we instead restart the evaluation, we have to recalculate everything: 20n^2, but then the n is smaller.
            # So when is n * (8n^2 + 24.25(n +1)) ≥ deleted * (8n^2 + 24.25(n +1)) + 20(n - deleted)^2 +
            #                                       (n - deleted)*(8(n - deleted)^2 + 24.25(n - deleted +1))
            # As soon as n > deleted ≥ 2, we find this to be fulfilled. However, these theoretical considerations do not seem
            # to be particularly successful.
            deleted += 1
            if n > 100 && deleted == 50
                cliques = @view cliques[idx]
                n -= deleted
                # we already have enough space allocated at adjm; we will now simply overwrite it.
                adjm = unsafe_wrap(Array, pointer(adjmOwn), (n, n), own=false)
                adjm .= collect(@capture(i > j ? length($cliques[i])^3 + length(cliques[j])^3 -
                                                 length(cliques[i] ∪ cliques[j])^3 : 0
                                         for i in 1:n, j in 1:n))
                # same for idx
                idx = unsafe_wrap(Array, pointer(idxOwn), n, own=false)
                idx .= true
                deleted = 0
            else
                # update clique graph
                adjm[j, 1:j-1] .= 0
                adjm[j+1:n, j] .= 0
                # recompute weights for updated clique graph
                newLenCube = length(cliques[i])^3
                for k in 1:i-1
                    idx[k] && (adjm[i, k] = length(cliques[k])^3 + newLenCube - length(cliques[k] ∪ cliques[i])^3)
                end
                for k in i+1:n
                    idx[k] && (adjm[k, i] = newLenCube + length(cliques[k])^3 - length(cliques[k] ∪ cliques[i])^3)
                end
            end
        end
    end
    return [smallcliques; cliques[idx]]
end

function _convert_cliques(cl::Vector{<:SimpleMonomialVector{Nr,Nc,I}}) where {Nr,Nc,I<:Integer}
    # best situation (which should always be the case: all vectors have the same set of exponents)
    a, rest = Iterators.peel(cl)
    required_e = a.e
    for item in rest
        if required_e != item.e
            use_e = ExponentsAll{Nr+2Nc,I}()
            return let E=ExponentsAll{Nr+2Nc,I}()
                E, [v.e == use_e ? Set(v.indices) : Set(monomial_index.((E,), v)) for v in cl]
            end
        end
    end
    return required_e, Set.(getproperty.(cl, :indices))
end
function merge_cliques(cl::Vector{<:SimpleMonomialVector{Nr,Nc,I}}) where {Nr,Nc,I<:Integer}
    converted = _convert_cliques(cl)
    merged = merge_cliques!(converted[2])
    return SimpleMonomialVector{Nr,Nc,I}[SimpleMonomialVector{Nr,Nc}(converted[1], collect(mergedᵢ)) for mergedᵢ in merged]
end

function merge_cliques(groupings::Relaxation.RelaxationGroupings{Nr,Nc}) where {Nr,Nc}
    outtype_part = Vector{<:SimpleMonomialVector{Nr,Nc}}
    newobj = merge_cliques(groupings.obj)
    # it does not make sense to merge the zero constraints, as every unique product corresponds to a separate variable anyway.
    newnonnegs = similar(groupings.nonnegs)
    newnonnegs .= merge_cliques.(groupings.nonnegs)
    newpsds = similar(groupings.psds)
    newpsds .= merge_cliques.(groupings.psds)
    return Relaxation.RelaxationGroupings(newobj, groupings.zeros, newnonnegs, newpsds, groupings.var_cliques)
end

struct CliqueMerged{P<:Problem,R<:AbstractRelaxation{P},G<:RelaxationGroupings} <: AbstractRelaxationSparse{P}
    problem::P
    parent::R
    groupings::G

    @doc """
        CliqueMerged(relaxation::AbstractRelaxation)

    Performs clique merging on the parent `relaxation`. Clique merging may allow to reduce the solver time by merging together
    smaller blocks of variables with huge overlap into a single larger one; however, it comes at a significant cost itself.
    Basically, clique merging undoes some of the sparsity analysis performed before when it might be too excessive. This is
    is an equivalent reformulation that is as exact as `relaxation` itself.
    """
    function CliqueMerged(relaxation::AbstractRelaxation{P}) where {P}
        relaxation isa CliqueMerged && throw(MethodError(CliqueMerged, (relaxation,)))
        g = merge_cliques(groupings(relaxation))
        new{P,typeof(relaxation),typeof(g)}(poly_problem(relaxation), relaxation, g)
    end
end

function Base.show(io::IO, m::MIME"text/plain", relaxation::CliqueMerged)
    print(io, "Clique-merged version of ")
    return show(io, m, relaxation.parent)
end

MultivariatePolynomials.degree(c::CliqueMerged) = degree(c.parent)

basis(c::(CliqueMerged{P,<:AbstractRelaxationBasis{P}} where {P<:Problem})) = basis(c.parent)

default_solution_method(c::CliqueMerged) = default_solution_method(c.parent)