module MultivariateExponents

using ...PolynomialOptimization: @assert, @inbounds

export unsafe

struct Unsafe end
const unsafe = Unsafe()

include("./AbstractExponents.jl")
include("./ExponentsIndices.jl")
include("./AbstractExponentsUnbounded.jl")
include("./AbstractExponentsDegreeBounded.jl")


end