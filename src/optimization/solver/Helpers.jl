export monomial_count, trisize, count_uniques, SparseMatrixCOO

"""
    monomial_count(n, d)

Short helper function that allows to determine the number of monomials in `n` variables up to degree `d`,
``\\binom{n + d}{n}``. Uses the underlying cache machinery provided by [`ExponentsDegree`](@ref).
"""
monomial_count(n, d) = length(ExponentsDegree{n,UInt}(0:d))

"""
    trisize(n)

Returns the number of items in the triangle of a matrix of side dimension `n`, ``\\frac{n(n +1)}{2}``.
"""
trisize(n) = (n * (n +1)) >> 1

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
iscanonical(::SimpleMonomial{<:Any,0}) = true
Base.@assume_effects :consistent iscanonical(m::SimpleMonomialOrConj) = m ≤ SimpleConjMonomial(m)
canonicalize(m::SimpleMonomialOrConj) = iscanonical(m) ? m : SimpleConjMonomial(m)

@inline function getreim(state, args::SimpleMonomialOrConj{Nr,Nc}...) where {Nr,Nc}
    idx₁ = mindex(state, args...)
    idx₂ = mindex(state, SimpleConjMonomial.(reverse(args))...) # reverse is rather unnecessary (commuting)
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

struct ScalarMatrix{X} <: AbstractMatrix{X} # we could make it mutable and fully implement the AbstractArray interface, but we
                                            # don't need it
    x::X
end

Base.size(::ScalarMatrix) = (1, 1)
Base.length(::ScalarMatrix) = 1
@inline function Base.getindex(m::ScalarMatrix, args...)
    @boundscheck checkbounds(m, args...)
    return m.x
end
Base.iterate(m::ScalarMatrix) = m.x, nothing
Base.iterate(::ScalarMatrix, ::Nothing) = nothing

collect_grouping(g::AbstractVector{M} where M<:SimpleMonomial) = g
collect_grouping(g) = collect(g)

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