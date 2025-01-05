module PolynomialOptimizationMosek

using Mosek, MultivariatePolynomials, LinearAlgebra, SparseArrays, PolynomialOptimization.Solver,
    PolynomialOptimization.Newton, PolynomialOptimization.Solvers.SpecBM, StandardPacked
using PolynomialOptimization: @assert, @inbounds, @allocdiff
using PolynomialOptimization.SimplePolynomials: veciter
using Mosek: msk_global_env, Env, deletetask
using StandardPacked: SPMatrix
import PolynomialOptimization

printstream(msg::String) = (print(msg); flush(stdout))

include("./MosekSOS.jl")
# Do we have Mosek version at least 10?
if isdefined(Mosek, :appendafes)
    include("./MosekMoment.jl")
    if VersionNumber(Mosek.getversion()) ≥ v"10.1.11"
        isdefined(Mosek, :appendafes) && include("./SpecBM.jl")
        VersionNumber(Mosek.getversion()) < v"10.1.13" &&
            @warn("Consider upgrading your version of Mosek to avoid rare crashes.")
    else
        @warn("The SpecBM method Mosek is not available: upgrade your Mosek distribution to at least version 10.1.11.")
        # Don't issue this warning if <10 is installed - most likely, the license doesn't allow for 10. But if it does, then
        # we'll need the user to fix the bugs.
    end
end
include("./Newton.jl")
include("./Tightening.jl")

function __init__()
    if isdefined(Mosek, :appendafes)
        push!(solver_methods, :MosekMoment)
        if VersionNumber(Mosek.getversion()) < v"10.1.13"
            push!(SpecBM.specbm_methods, :Mosek)
        else
            pushfirst!(SpecBM.specbm_methods, :Mosek)
        end
    end
    pushfirst!(solver_methods, :Mosek, :MosekSOS)
    pushfirst!(Newton.newton_methods, :Mosek)
    pushfirst!(PolynomialOptimization.tightening_methods, :Mosek)
end

@solver_alias Mosek MosekSOS

end