module PolynomialOptimizationProxSDP

using ProxSDP, MultivariatePolynomials, LinearAlgebra, SparseArrays, StandardPacked, PolynomialOptimization.Solver
using PolynomialOptimization: @assert, @inbounds

include("./ProxSDPMoment.jl")

__init__() = push!(solver_methods, :ProxSDP, :ProxSDPMoment)

@solver_alias ProxSDP ProxSDPMoment

end