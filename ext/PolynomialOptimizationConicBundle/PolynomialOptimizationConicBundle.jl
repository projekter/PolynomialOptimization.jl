module PolynomialOptimizationConicBundle

using ConicBundle, MultivariatePolynomials, LinearAlgebra, OffsetArrays, SparseArrays, StandardPacked,
    PolynomialOptimization.Solver
using PolynomialOptimization: @assert, @inbounds

include("./ConicBundleMoment.jl")

__init__() = push!(solver_methods, :ConicBundle, :ConicBundleMoment)

@solver_alias ConicBundle ConicBundleMoment

end