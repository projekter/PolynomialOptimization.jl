module PolynomialOptimizationSCS

using SCS, MultivariatePolynomials, LinearAlgebra, SparseArrays, PolynomialOptimization.Solver
using PolynomialOptimization: @assert, @inbounds
using SCS: ScsCone, ScsData, ScsSettings, ScsSolution, ScsInfo, ScsMatrix, scsint_t, scs_init, scs_solve, scs_finish,
    LinearSolver

include("./SCSMoment.jl")

__init__() = push!(solver_methods, :SCS, :SCSMoment)

@solver_alias SCS SCSMoment

end