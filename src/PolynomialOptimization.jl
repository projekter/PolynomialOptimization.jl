__precompile__()
module PolynomialOptimization

using MultivariatePolynomials
using SparseArrays
using LinearAlgebra
using Preferences
using Printf
import MutableArithmetics
import StatsBase

export Newton, Relaxation

const haveMPI = Ref{Bool}(false)
const debugging = @load_preference("debugging", false) # only for testing

include("./helpers/Helpers.jl")

include("./Problem.jl")
include("./relaxations/Relaxation.jl")
using .Relaxation
include("./optimization/Optimization.jl")
include("./newton/Newton.jl")
import .Newton
include("./solutions/SolutionExtraction.jl")
include("./Tightening.jl")
include("./solvers/Solvers.jl")

end