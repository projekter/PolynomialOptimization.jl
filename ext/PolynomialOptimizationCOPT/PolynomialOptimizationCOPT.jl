module PolynomialOptimizationCOPT

using COPT, MultivariatePolynomials, SparseArrays, PolynomialOptimization.Solver, PolynomialOptimization.Newton
using PolynomialOptimization: @assert, @inbounds, @allocdiff
using PolynomialOptimization.SimplePolynomials: veciter
using COPT: _check_ret, Env#, libcopt
import PolynomialOptimization

global copt_env::Env

mutable struct COPTProb
    ptr::Ptr{copt_prob}
    finalize_called::Bool

    function COPTProb(env::Env)
        p_ptr = Ref{Ptr{copt_prob}}(C_NULL)
        _check_ret(env, COPT_CreateProb(env.ptr, p_ptr))
        problem = new(p_ptr[], false)
        finalizer(problem) do p
            p.finalize_called = true
            if p.ptr != C_NULL
                COPT_DeleteProb(Ref(p.ptr))
                p.ptr = C_NULL
            end
        end
        return problem
    end
end

Base.unsafe_convert(::Type{Ptr{copt_prob}}, problem::COPTProb) = problem.ptr
Base.unsafe_convert(::Type{Ptr{Ptr{copt_prob}}}, problem::COPTProb) = Ptr{Ptr{copt_prob}}(pointer_from_objref(problem))

include("./COPTMoment.jl")
include("./Newton.jl")
include("./Tightening.jl")

function __init__()
    global copt_env = Env()
    pushfirst!(solver_methods, :COPT, :COPTMoment)
    pushfirst!(Newton.newton_methods, :COPT)
    pushfirst!(PolynomialOptimization.tightening_methods, :COPT)
end

@solver_alias COPT COPTMoment

end