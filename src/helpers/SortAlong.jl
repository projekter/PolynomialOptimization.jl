"""
    sort_along!(v::AbstractVector, along::AbstractVector...; lo=1, hi=length(v),
        o=Base.Order.Forward, relevant=1)

This helper function sorts the vector `v` as `sort!` would do, but at the same time puts all vectors in `along` in the same
order as `v`. Therefore, `sort_along!(v, a, b, c)` is equivalent to
```julia
p = sortperm(v)
permute!.((v, a, b, c), (p,))
```
but does not require any new intermediate allocations.
Note that by increasing `relevant` (to at most `1 + length(along)`), additional vectors in `along` can be taken into account
when breaking ties in `v`.

!!! info "Implementation"
    Note that the sorting functions are exactly identical to the ones used in Julia 1.8 (there was a big overhaul in 1.9 which
    made everything extremely complicated; we don't reproduce the new behavior). This means that either InsertionSort or
    QuickSort are used.
"""
function sort_along!(v::AbstractVector, along::AbstractVector...; lo::Integer=1, hi::Integer=length(v),
    o::Base.Ordering=Base.Order.Forward, relevant::Integer=1)
    lv = length(v)
    all(x -> length(x) == lv, along) || error("sort_along! requires vectors of the same length")
    lo ≥ 1 || error("sort_along! parameter lo must be strictly positive")
    hi ≤ lv || error("sort_along! parameter hi must not exceed length of the vectors")
    1 ≤ relevant ≤ 1 + length(along) || error("sort_along! parameter relevant must be a valid number of vectors")
    return sort_along!((v, along...), lo, hi, o, Val(relevant))
end

function sort_along_expr(v1, i1, v2, i2, relevant)
    check = Expr(:||, :(Base.lt(o, $(isnothing(i1) ? :($v1[1]) : :($v1[1][$i1])),
                                   $(isnothing(i2) ? :($v2[1]) : :($v2[1][$i2])))::Bool))
    curlog = check
    for i in 2:relevant
        newor = Expr(:||, :(Base.lt(o, $(isnothing(i1) ? :($v1[$i]) : :($v1[$i][$i1])),
                                       $(isnothing(i2) ? :($v2[$i]) : :($v2[$i][$i2])))::Bool))
        newand = :($(isnothing(i1) ? :($v1[$(i-1)]) : :($v1[$(i-1)][$i1])) ==
                   $(isnothing(i2) ? :($v2[$(i-1)]) : :($v2[$(i-1)][$i2])) && $newor)
        push!(curlog.args, newand)
        curlog = newor
    end
    return check
end

@generated function sort_along!(v::NTuple{N,AbstractVector}, lo::Integer, hi::Integer, o::Base.Ordering, ::Val{relevant}) where {N,relevant}
    quote
        @inbounds while lo < hi
            if hi - lo <= Base.SMALL_THRESHOLD
                # function sort!(v::AbstractVector, lo::Integer, hi::Integer, ::InsertionSortAlg, o::Ordering)
                lo_plus_1 = (lo +1)::Integer
                for i in lo_plus_1:hi
                    j = i
                    x = ($((:(v[$i][i]) for i in 1:N)...),)
                    while j > lo
                        y = ($((:(v[$i][j-1]) for i in 1:relevant)...),)
                        $(sort_along_expr(:x, nothing, :y, nothing, relevant)) || break
                        $((:(v[$i][j] = y[$i]) for i in 1:relevant)...)
                        $((:(v[$i][j] = v[$i][j-1]) for i in relevant+1:N)...)
                        j -= 1
                    end
                    $((:(v[$i][j] = x[$i]) for i in 1:N)...)
                end
                return v
                # end
            end
            # function partition!(v::AbstractVector, lo::Integer, hi::Integer, o::ordering)
                # function selectpivot!(v::AbstractVector, lo::Integer, hi::Integer, o::Ordering)
                mi = Base.midpoint(lo, hi)
                if $(sort_along_expr(:v, :lo, :v, :mi, relevant))
                    $((:((v[$i][mi], v[$i][lo]) = (v[$i][lo], v[$i][mi])) for i in 1:N)...)
                end
                if $(sort_along_expr(:v, :hi, :v, :lo, relevant))
                    if $(sort_along_expr(:v, :hi, :v, :mi, relevant))
                        $((:((v[$i][hi], v[$i][lo], v[$i][mi]) = (v[$i][lo], v[$i][mi], v[$i][hi])) for i in 1:N)...)
                    else
                        $((:((v[$i][hi], v[$i][lo]) = (v[$i][lo], v[$i][hi])) for i in 1:N)...)
                    end
                end
                pivot = ($((:(v[$i][lo]) for i in 1:N)...),)
                # end
                i, j = lo, hi
                while true
                    i += 1; j -= 1
                    while $(sort_along_expr(:v, :i, :pivot, nothing, relevant)); i += 1; end;
                    while $(sort_along_expr(:pivot, nothing, :v, :j, relevant)); j -= 1; end;
                    i >= j && break
                    $((:((v[$i][i], v[$i][j]) = (v[$i][j], v[$i][i])) for i in 1:N)...)
                end
                $((:((v[$i][j], v[$i][lo]) = (pivot[$i], v[$i][j])) for i in 1:N)...)
            # end
            if j - lo < hi - j
                lo < (j -1) && sort_along!(v, lo, j -1, o, Val(relevant))
                lo = j +1
            else
                j +1 < hi && sort_along!(v, j +1, hi, o, Val(relevant))
                hi = j -1
            end
        end
        return v
    end
end