# while this is not a weak dependecy as the solver interface is shipped with PO, we follow the naming convention
module PolynomialOptimizationLoRADS

using ..Solvers.LoRADS, MultivariatePolynomials, LinearAlgebra, SparseArrays, ..Solver
using ..PolynomialOptimization: @assert, @inbounds

include("./LoRADSMoment.jl")

__init__() = !isempty(LoRADS.solverlib) && push!(solver_methods, :LoRADS, :LoRADSMoment)

@solver_alias LoRADS LoRADSMoment

end