module PolynomialOptimizationSketchyCGAL

using MultivariatePolynomials, LinearAlgebra, SparseArrays, PolynomialOptimization.Solver,
    PolynomialOptimization.SimplePolynomials, PolynomialOptimization.SimplePolynomials.MultivariateExponents
import ..Solvers.SketchyCGAL
using PolynomialOptimization: @assert, @inbounds

include("./SketchyCGALMoment.jl")

__init__() = push!(solver_methods, :SketchyCGAL, :SketchyCGALMoment)

@solver_alias SketchyCGAL SketchyCGALMoment

end