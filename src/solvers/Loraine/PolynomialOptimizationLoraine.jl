module PolynomialOptimizationLoraine

using MultivariatePolynomials, LinearAlgebra, SparseArrays, PolynomialOptimization.Solver,
    PolynomialOptimization.SimplePolynomials, PolynomialOptimization.SimplePolynomials.MultivariateExponents
import ..Solvers.Loraine
using PolynomialOptimization: @assert, @inbounds

include("./LoraineMoment.jl")

__init__() = push!(solver_methods, :Loraine, :LoraineMoment)

@solver_alias Loraine LoraineMoment

end