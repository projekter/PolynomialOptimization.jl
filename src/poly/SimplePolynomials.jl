module SimplePolynomials

using MultivariatePolynomials, ..FastVector
using ..PolynomialOptimization: @assert, @inbounds, sort_along!

include("./Utils.jl")
include("./exponents/MultivariateExponents.jl")
using .MultivariateExponents
using .MultivariateExponents: ExponentIndices, Unsafe, index_counts
import .MultivariateExponents: veciter
include("./Variable.jl")
include("./Monomial.jl")
include("./MonomialVector.jl")
include("./Polynomial.jl")
include("./TypeHandling.jl")

end