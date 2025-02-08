module PolynomialOptimizationStride

using ..Solvers.Stride, MultivariatePolynomials, LinearAlgebra, SparseArrays, StandardPacked, ..Solver
using ..Solvers: EfficientCholmod
using ..PolynomialOptimization: @assert, @inbounds

include("./StrideMoment.jl")

__init__() = push!(solver_methods, :Stride, :StrideMoment)

@solver_alias Stride StrideMoment

end