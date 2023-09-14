# Here, we provide the helper function sort_along!, which allows to sort! a vector (the first argument) while at the same time
# sorting various other vectors using the exact same permutation. Hence,
#    sort_along!(v, a, b, c)
# is equivalent to
#    p = sortperm(v)
#    permute!.((v, a, b, c), (p,))
# but does not require any new intermediate allocations.
# Note that the sorting functions are exactly identical to the ones used in Julia 1.8 (there was a big overhaul in 1.9 which
# made everything extremely complicated; we don't reproduce the new behavior).

function sort_along!(v::AbstractVector, along::AbstractVector...; lo::Integer=1, hi::Integer=length(v),
    o::Base.Ordering=Base.Order.Forward)
    lv = length(v)
    all(x -> length(x) == lv, along) || error("sort_along! requires vectors of the same length")
    lo ≥ 1 || error("sort_along! parameter lo must be strictly positive")
    hi ≤ lv || error("sort_along! parameter hi must not exceed length of the vectors")
    return sort_along!(v, lo, hi, o, along...)
end

function sort_along!(v::AbstractVector, lo::Integer, hi::Integer, o::Base.Ordering, along::Vararg{<:AbstractVector,N}) where {N}
    while lo < hi
        if hi - lo <= Base.SMALL_THRESHOLD
            # function sort!(v::AbstractVector, lo::Integer, hi::Integer, ::InsertionSortAlg, o::Ordering)
            lo_plus_1 = (lo +1)::Integer
            for i in lo_plus_1:hi
                j = i
                x = v[i]
                if @generated
                    Expr(:block, [:($(Symbol("x", i)) = along[$i][i]) for i in 1:N]...)
                else
                    x_along = [alongᵢ[i] for alongᵢ in along]
                end
                while j > lo
                    y = v[j-1]
                    Base.lt(o, x, y)::Bool || break
                    v[j] = y
                    if @generated
                        Expr(:block, [:(along[$i][j] = along[$i][j-1]) for i in 1:N]...)
                    else
                        for alongᵢ in along
                            alongᵢ[j] = alongᵢ[j-1]
                        end
                    end
                    j -= 1
                end
                v[j] = x
                if @generated
                    Expr(:block, [:(along[$i][j] = $(Symbol("x", i))) for i in 1:N]...)
                else
                    for (alongᵢ, x_alongᵢ) in zip(along, x_along)
                        alongᵢ[j] = x_alongᵢ
                    end
                end
            end
            return v, along...
            # end
        end
        # function partition!(v::AbstractVector, lo::Integer, hi::Integer, o::ordering)
            # function selectpivot!(v::AbstractVector, lo::Integer, hi::Integer, o::Ordering)
            mi = Base.midpoint(lo, hi)
            if Base.lt(o, v[lo], v[mi])
                v[mi], v[lo] = v[lo], v[mi]
                if @generated
                    Expr(:block, [:((along[$i][mi], along[$i][lo]) = (along[$i][lo], along[$i][mi])) for i in 1:N]...)
                else
                    for alongᵢ in along
                        alongᵢ[mi], alongᵢ[lo] = alongᵢ[lo], alongᵢ[mi]
                    end
                end
            end
            if Base.lt(o, v[hi], v[lo])
                if Base.lt(o, v[hi], v[mi])
                    v[hi], v[lo], v[mi] = v[lo], v[mi], v[hi]
                    if @generated
                        Expr(:block, [:((along[$i][hi], along[$i][lo], along[$i][mi]) =
                                        (along[$i][lo], along[$i][mi], along[$i][hi])) for i in 1:N]...)
                    else
                        for alongᵢ in along
                            alongᵢ[hi], alongᵢ[lo], alongᵢ[mi] = alongᵢ[lo], alongᵢ[mi], alongᵢ[hi]
                        end
                    end
                else
                    v[hi], v[lo] = v[lo], v[hi]
                    if @generated
                        Expr(:block, [:((along[$i][hi], along[$i][lo]) = (along[$i][lo], along[$i][hi])) for i in 1:N]...)
                    else
                        for alongᵢ in along
                            alongᵢ[hi], alongᵢ[lo] = alongᵢ[lo], alongᵢ[hi]
                        end
                    end
                end
            end
            pivot = v[lo]
            if @generated
                Expr(:block, [:($(Symbol("pivot", i)) = along[$i][lo]) for i in 1:N]...)
            else
                pivot_along = [alongᵢ[lo] for alongᵢ in along]
            end
            # end
            i, j = lo, hi
            while true
                i += 1; j -= 1
                while Base.lt(o, v[i], pivot); i += 1; end;
                while Base.lt(o, pivot, v[j]); j -= 1; end;
                i >= j && break
                v[i], v[j] = v[j], v[i]
                if @generated
                    Expr(:block, [:((along[$i][i], along[$i][j]) = (along[$i][j], along[$i][i])) for i in 1:N]...)
                else
                    for alongᵢ in along
                        alongᵢ[i], alongᵢ[j] = alongᵢ[j], alongᵢ[i]
                    end
                end
            end
            v[j], v[lo] = pivot, v[j]
            if @generated
                Expr(:block, [:((along[$i][j], along[$i][lo]) = ($(Symbol("pivot", i)), along[$i][j])) for i in 1:N]...)
            else
                for (alongᵢ, pivot_alongᵢ) in zip(along, pivot_along)
                    alongᵢ[j], alongᵢ[lo] = pivot_alongᵢ, alongᵢ[j]
                end
            end
        # end
        if j - lo < hi - j
            lo < (j -1) && sort_along!(v, lo, j -1, o, along...)
            lo = j +1
        else
            j +1 < hi && sort_along!(v, j +1, hi, o, along...)
            hi = j -1
        end
    end
    return v, along...
end