include("./Macros.jl")
include("./Mutation.jl")
include("./FastVector.jl")
using .FastVector
include("./StackVector.jl")
include("./SortAlong.jl")
include("./FastKey.jl")
include("./Allocations.jl")

# we must load this before MatrixPolynomials (so that effective_variables_in is known), but after FastVector, which is required
# by IntPolynomials.
include("../poly/IntPolynomials.jl")
using .IntPolynomials, .IntPolynomials.MultivariateExponents
using .IntPolynomials: IntMonomialOrConj
import .IntPolynomials.MultivariateExponents: iterate! # be careful to avoid duplication of methods; let's reuse this one

include("./MatrixPolynomials.jl")

"""
    issubset_sorted(a, b)

Equivalent to `issubset(a, b)`, but assumes that both `a` and `b` are sorted vectors.
"""
function issubset_sorted(a, b)
    i = firstindex(b)
    @inbounds for x in a
        i += searchsortedlast(@view(b[i:end]), x) -1
        if i < firstindex(b) || b[i] != x
            return false
        end
    end
    return true
end

if VERSION < v"1.11"
    # Julia 1.10 can be very slow since the union will always try to shrink. To avoid piracy, copy the shrink-customizable
    # definitions from 1.11 under new names.
    function sizehint_!(d::Dict{T}, newsz; shrink::Bool=true) where T
        oldsz = length(d.slots)
        # limit new element count to max_values of the key type
        newsz = min(max(newsz, length(d)), Base.max_values(T)::Int)
        # need at least 1.5n space to hold n elements
        newsz = Base._tablesz(cld(3 * newsz, 2))
        return (shrink ? newsz == oldsz : newsz <= oldsz) ? d : Base.rehash!(d, newsz)
    end
    sizehint_!(s::Set, newsz; shrink::Bool=true) = (sizehint_!(s.dict, newsz; shrink); s)
    sizehint_!(x, y; shrink::Bool=true) = sizehint!(x, y) # not implemented, but provide for upwards-compatibility

    function union_!(s::AbstractSet{T}, itr) where T
        Base.haslength(itr) && sizehint_!(s, length(s) + Int(length(itr))::Int; shrink = false)
        for x in itr
            push!(s, x)
            length(s) == Base.max_values(T) && break
        end
        return s
    end
else
    sizehint_! = sizehint!
    union_! = union!
end