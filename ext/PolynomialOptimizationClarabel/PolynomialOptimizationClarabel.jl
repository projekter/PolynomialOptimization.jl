module PolynomialOptimizationClarabel

using Clarabel, MultivariatePolynomials, LinearAlgebra, SparseArrays, PolynomialOptimization.Solver
using PolynomialOptimization: @assert, @inbounds

include("./ClarabelMoment.jl")

__init__() = push!(solver_methods, :Clarabel, :ClarabelMoment)

@solver_alias Clarabel ClarabelMoment

end