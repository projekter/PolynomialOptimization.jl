export monomial_count, trisize, count_uniques, SparseMatrixCOO

"""
    monomial_count(n, d)

Short helper function that allows to determine the number of monomials in `n` variables up to degree `d`,
``\\binom{n + d}{n}``. Uses the underlying cache machinery provided by [`ExponentsDegree`](@ref).
"""
monomial_count(n, d) = length(ExponentsDegree{n,UInt}(0:d))

macro twice(symb::Symbol, condition, body)
    esc(quote
        let $symb=false
            $body
        end
        $condition && let $symb=true
            $body
        end
    end)
end

# For complex-valued problems, there are m and conj(m); all solvers require real-valued inputs.
# Therefore, we first define the canonicalized monomial m̃ as the one of m, conj(m) that has the smaller index.
iscanonical(::IntMonomial{<:Any,0}) = true
Base.@assume_effects :consistent iscanonical(m::IntMonomialOrConj) = m ≤ IntConjMonomial(m)
canonicalize(m::IntMonomialOrConj) = iscanonical(m) ? m : IntConjMonomial(m)

@inline function getreim(state, args::IntMonomialOrConj{Nr,Nc}...) where {Nr,Nc}
    idx₁ = mindex(state, args...)
    idx₂ = mindex(state, IntConjMonomial.(reverse(args))...) # reverse is rather unnecessary (commuting)
    if idx₁ ≤ idx₂
        return idx₁, idx₂, true # make sure to test for equality whenever the imaginary part is involved
    else
        return idx₂, idx₁, false
    end
end

@inline @generated function findin(haystacks::NTuple{N,Any}, needles::NTuple{N,Any}, startfrom::Integer=firstindex(haystacks[1])) where {N}
    quote
        # should we assert that all lengths be equal?
        @inbounds for i in startfrom:length(haystacks[1])
            # Assume haystacks and needles are short and the comparison itself cheap, but jumping after the comparison not. So
            # to minimize the number of jumps, we compare everything always.
            if +($((:(haystacks[$i][i] == needles[$i]) for i in 1:N)...),) == $N
                return i
            end
        end
        return nothing
    end
end

macro pushorupdate!(part::Symbol, icoeff...)
    iseven(length(icoeff)) || throw(MethodError(var"@pushorupdate!", (__source__, __module__, part, icoeff...)))
    n = length(icoeff) ÷ 2
    esc(quote
        let $((Expr(:(=), Symbol(:coeff, i), icoeff[2i]) for i in 1:n)...),
            partidx=@isdefined(rows) ? findin((indices, rows), ($part, $(icoeff[1])), length(indices) - curlen +1) :
                                       findin((indices,), ($part,), length(indices) - curlen +1)
            if isnothing(partidx)
                @isdefined(rows) && unsafe_push!(rows, $((icoeff[2i-1] for i in 1:n)...))
                unsafe_push!(indices, $((part for _ in 1:n)...))
                unsafe_push!(values, $((Symbol(:coeff, i) for i in 1:n)...))
                curlen += $n
            else
                $((:(values[partidx+$(i-1)] += $(Symbol(:coeff, i))) for i in 1:n)...)
            end
        end
    end)
end

@inline function pushorupdate!(idxvec::FastVec, index, valuevec::FastVec, value)
    i = findfirst(isequal(index), idxvec)
    if isnothing(i)
        unsafe_push!(idxvec, index)
        unsafe_push!(valuevec, value)
    else
        @inbounds valuevec[i] += value
    end
end

struct ScalarArray{X,N} <: AbstractArray{X,N} # we could make it mutable and fully implement the AbstractArray interface, but
                                              # we don't need it
    x::X
end

ScalarArray{<:Any,N}(x::X) where {X,N} = ScalarArray{X,N}(x)

Base.size(::ScalarArray{<:Any,N}) where {N} = ntuple(_ -> 1, Val(N))
Base.length(::ScalarArray) = 1
@inline function Base.getindex(m::ScalarArray, args...)
    @boundscheck checkbounds(m, args...)
    return m.x
end
Base.iterate(m::ScalarArray) = m.x, nothing
Base.iterate(::ScalarArray, ::Nothing) = nothing

const ScalarVector{X} = ScalarArray{X,1}
const ScalarMatrix{X} = ScalarArray{X,2}

"""
    count_uniques(vec::AbstractVector[, callback])
    count_uniques(vec₁::AbstractVector, vec₂::AbstractVector[, callback])

Return the number of unique elements in the vector(s), which must be sorted but may possibly contain duplicates. The callback
is invoked once for every unique entry. Its first parameter is the index of the element in the unique total vector, its second
(and third) is the last index/indices correponding to the element in the input vector(s). In the second form which allows to
check for two vectors jointly, one of the callback parameters can be `missing` if the element is present only in one of the two
vectors.
"""
function count_uniques(vec::AbstractVector, callback::F=(_, _) -> nothing) where {F}
    @assert(issorted(vec))
    i = 1
    remaining = length(vec)
    index = 1
    @inbounds while !iszero(remaining)
        cur = vec[i]
        # skip over duplicates
        while remaining > 1 && vec[i+1] == cur
            i += 1; remaining -= 1
        end
        @inline callback(index, i)
        i += 1; remaining -= 1
        index += 1
    end
    return index -1 # return the count
end

function count_uniques(vec₁::AbstractVector{I}, vec₂::AbstractVector{I}, callback::F=(_, _, _) -> nothing) where {I,F}
    @assert(issorted(vec₁) && issorted(vec₂))
    # vec₁ and vec₂ are sorted vectors with possible duplicates. Iterate through all of them, count the unique ones, call the
    # callback for every unique entry with the respective last indices that correspond to this element.
    i₁ = 1
    i₂ = 1
    remaining₁ = length(vec₁)
    remaining₂ = length(vec₂)
    index = 1
    @inbounds while !iszero(remaining₁) && !iszero(remaining₂)
        cur₁ = vec₁[i₁]
        cur₂ = vec₂[i₂]
        # skip over duplicates
        while remaining₁ > 1 && vec₁[i₁+1] == cur₁
            i₁ += 1; remaining₁ -= 1
        end
        while remaining₂ > 1 && vec₂[i₂+1] == cur₂
            i₂ += 1; remaining₂ -= 1
        end
        # and work with the smaller one until it is no longer the smaller one
        if cur₁ == cur₂
            @inline callback(index, i₁, i₂)
            i₁ += 1; remaining₁ -= 1
            i₂ += 1; remaining₂ -= 1
            index += 1
        elseif cur₁ < cur₂
            @inline callback(index, i₁, missing)
            i₁ += 1; remaining₁ -= 1
            index += 1
            while !iszero(remaining₁)
                cur₁ = vec₁[i₁]
                cur₁ ≥ cur₂ && break
                # skip over duplicates
                while remaining₁ > 1 && vec₁[i₁+1] == cur₁
                    i₁ += 1; remaining₁ -= 1
                end
                @inline callback(index, i₁, missing)
                i₁ += 1; remaining₁ -= 1
                index += 1
            end
        else
            @inline callback(index, missing, i₂)
            i₂ += 1; remaining₂ -= 1
            index += 1
            while remaining₂ > 1
                cur₂ = vec₂[i₂]
                cur₂ ≥ cur₁ && break
                # skip over duplicates
                while remaining₂ > 1 && vec₂[i₂+1] == cur₂
                    i₂ += 1; remaining₂ -= 1
                end
                @inline callback(index, missing, i₂)
                i₂ += 1; remaining₂ -= 1
                index += 1
            end
        end
    end
    # tail checks. At most one of the two still has elements.
    @inbounds while !iszero(remaining₁)
        cur₁ = vec₁[i₁]
        # skip over duplicates
        while remaining₁ > 1 && vec₁[i₁+1] == cur₁
            i₁ += 1; remaining₁ -= 1
        end
        @inline callback(index, i₁, missing)
        i₁ += 1; remaining₁ -= 1
        index += 1
    end
    @inbounds while !iszero(remaining₂)
        cur₂ = vec₂[i₂]
        # skip over duplicates
        while remaining₂ > 1 && vec₂[i₂+1] == cur₂
            i₂ += 1; remaining₂ -= 1
        end
        @inline callback(index, missing, i₂)
        i₂ += 1; remaining₂ -= 1
        index += 1
    end
    return index -1 # return the count
end

"""
    SparseMatrixCOO{I<:Integer,K<:Integer,V<:Real,Offset}

Representation of a sparse matrix in COO form. Fields are `rowinds::FastVec{I}`, `moninds::FastVec{K}` (where `K` is of the
type returned by `monomial_index`), and `nzvals::FastVec{V}`. The first row/column for the solver has index `Offset` (of type
`I`).

If this matrix is used in the context of a vectorized COO matrix, the fields are `rows`, `cols`, and `vals` with `K === I`.
"""
struct SparseMatrixCOO{I<:Integer,K<:Integer,V<:Real,Offset}
    rowinds::FastVec{I}
    moninds::FastVec{K}
    nzvals::FastVec{V}

    function SparseMatrixCOO{I,K,V,Offset}() where {I<:Integer,K<:Integer,V<:Real,Offset}
        Offset isa I || throw(MethodError(SparseMatrixCOO{I,K,V,Offset}, ()))
        new{I,K,V,Offset}(FastVec{I}(), FastVec{K}(), FastVec{V}())
    end

    SparseMatrixCOO(rows::FastVec{I}, cols::FastVec{I}, vals::FastVec{V}, offset::I) where {I,V<:Real} =
        new{I,I,V,offset}(rows, cols, vals)
end

function Base.getproperty(smc::SparseMatrixCOO, f::Symbol)
    (f === :rows || f === :rowinds) && return getfield(smc, :rowinds)
    (f === :cols || f === :moninds) && return getfield(smc, :moninds)
    (f === :vals || f === :nzvals) && return getfield(smc, :nzvals)
    return getfield(smc, f)
end

Base.length(smc::SparseMatrixCOO) = length(smc.rowinds)
@inline function Base.size(smc::SparseMatrixCOO{<:Integer,<:Integer,<:Real,Offset}, dim) where {Offset}
    dim == 1 || error("Not implemented")
    @inbounds return isempty(smc.rowinds) ? 0 : Int(smc.rowinds[end]) + (1 - Int(Offset))
end

"""
    unique_outer_groupings(grouping::AbstractVector{<:IntMonomial}[; result::Set]
        [, exponents::AbstractExponents])

Computes a set of all unique groupings that are produced via `groupings * groupings'`. For purely complex-valued groupings,
this is the full list; as soon as we have a real variable present, it is smaller.
Returns a set of `FastKey`-wrapped monomial indices, a boolean indicating whether the grouping is effectively real-valued, and
the number of (real-valued) equalities that arise from this grouping (which should be multiplied by 2 if the groupings are to
be multiplied by a complex-valued constraint).
If present, the data will be appended to the `result` parameter.
If present, the output index will be in the exponent set `e`.
"""
function unique_outer_groupings(grouping::AbstractVector{M}; e::Union{Nothing,AbstractExponents}=nothing,
    result::Set{FastKey{I}}=Set{FastKey{isnothing(e) ? IntPolynomials._get_I(M) : indextype(e)}}()) where
    {Nr,Nc,I<:Integer,M<:IntMonomial{Nr,Nc}}
    isempty(grouping) && return result, true, 0
    if isnothing(e)
        e = @inbounds grouping[begin].e
    end
    # To avoid rehashings, get an overestimator of the total grouping size first.
    # TODO (maybe): In the first loop to populate unique_groupings, we determine whether the grouping is real-valued. So we
    # could instead populate two sets, saving isreal and a lot of conditionals in the second loop.
    result = sizehint_!(result, length(result) + (iszero(Nc) ? trisize(length(grouping)) : length(grouping)^2), shrink=false)
    real_grouping = true
    totalsize = 0
    for (i, g₁) in enumerate(grouping)
        if !iszero(Nc)
            g₁real = !iszero(Nr) && isreal(g₁)
            # Consider the g₂ = ḡ₁ case separately in the complex case. Explanations below.
            let g₂=g₁
                prodidx = FastKey(monomial_index(e, g₁, IntConjMonomial(g₂)))
                indexug, sh = Base.ht_keyindex2_shorthash!(result.dict, prodidx)
                if indexug ≤ 0
                    @inbounds Base._setindex!(result.dict, nothing, prodidx, -indexug, sh)
                    totalsize += 1
                end
            end
        end
        # In the real case, we can skip the entries behind i as they would lead to duplicates.
        # In the complex case, we can also skip them, as they would lead to exact conjugates, which in the end give rise to the
        # same conditions (but note that i is already handled above).
        for g₂ in Iterators.take(grouping, iszero(Nc) ? i : i -1)
            # We don't use mindex, as this can have unintended side-effects on the solver state (such as creating a
            # representation for this monomial, although we don't even know whether we need it - if constraint does not contain
            # a constant term, this function must not automatically add all the squared groupings as monomials, even if they
            # will probably appear at some place).
            prodidx = FastKey(monomial_index(e, g₁, IntConjMonomial(g₂)))
            # We need to add the product to the set if it does not exists; we also need to count the number of conditions that
            # we get out of it.
            indexug, sh = Base.ht_keyindex2_shorthash!(result.dict, prodidx)
            if indexug ≤ 0
                # It does not exist.
                @inbounds Base._setindex!(result.dict, nothing, prodidx, -indexug, sh)
                # Assume we have a grouping g = (gᵣ + im*gᵢ) and a polynomial p = pᵣ + im*pᵢ, where the individual parts are
                # real-valued. Then, add_equality! means that g*p = 0 and ḡ*p = 0. Of course we can also conjugate everything.
                # We must split each constraint into its real and imaginary parts:
                # (I)   Re(g*p) = gᵣ*pᵣ - gᵢ*pᵢ
                # (II)  Im(g*p) = gᵣ*pᵢ + gᵢ*pᵣ
                # (III) Re(ḡ*p) = gᵣ*pᵣ + gᵢ*pᵢ
                # (IV)  Im(ḡ*p) = gᵣ*pᵢ - gᵢ*pᵣ
                # To analyze this (which would be easier if we added and subtracted the equalities, but in the
                # IntPolynomials setup, the given form is most easy to handle), let's consider linear dependencies.
                # - If the constraint is real-valued, (III) is equal to (I) and (IV) is -(II), so we only take (I) and (II).
                # - If the grouping is real-valued, (III) is equal to (I) and (IV) is equal to (II); we only take (I) and (II).
                # - If both are real-valued, (III) is equal to (I) while (II) and (IV) are zero, so we only take (I).
                # - If both are complex-valued, all constraints are linearly independent.
                # Rearranging this, we always take (I); if at least one is complex-valued, we also take (II); if both are, we
                # take all. Note that we don't have to consider the conjugates of the groupings separately, as they only yield
                # a global sign in the zero-equality.
                # For this loop, this means that we will only check whether g₁*ḡ₂ belongs to a real-valued monomial, in which
                # case we add 1; or to a complex-valued monomial, in which case we add 2. After the loop, we multiply by 2 if
                # the constraint was also complex-valued.
                if iszero(Nc) || (!iszero(Nr) && g₁real && isreal(g₂)) # note that g₁ ≠ ḡ₂
                    totalsize += 1
                else
                    totalsize += 2
                    real_grouping = false
                end
            end
        end
    end
    return result, real_grouping, totalsize
end