# while this is not a weak dependency as the solver is shipped with PO, we follow the naming convention
module PolynomialOptimizationSpecBM

using ..Solvers.SpecBM, MultivariatePolynomials, LinearAlgebra, SparseArrays, ..Solver
using ..PolynomialOptimization: @assert, @inbounds

include("./SpecBMSOS.jl")

__init__() = push!(solver_methods, :SpecBM, :SpecBMSOS)

@solver_alias SpecBM SpecBMSOS

end