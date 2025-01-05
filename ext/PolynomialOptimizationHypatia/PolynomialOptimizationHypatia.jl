module PolynomialOptimizationHypatia

using Hypatia, MultivariatePolynomials, LinearAlgebra, SparseArrays, PolynomialOptimization.Solver,
    PolynomialOptimization.Solvers.SpecBM
using PolynomialOptimization: @assert, @inbounds
using Hypatia: Cones, Models, Solvers

include("./HypatiaMoment.jl")
include("./SpecBM.jl")

function __init__()
    push!(solver_methods, :Hypatia, :HypatiaMoment)
    push!(SpecBM.specbm_methods, :Hypatia)
end

@solver_alias Hypatia HypatiaMoment

end