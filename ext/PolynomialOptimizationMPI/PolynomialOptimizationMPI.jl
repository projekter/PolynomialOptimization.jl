module PolynomialOptimizationMPI

using MultivariatePolynomials, PolynomialOptimization.SimplePolynomials, PolynomialOptimization.FastVector,
    PolynomialOptimization.SimplePolynomials.MultivariateExponents, Printf
import MPI, Random
using PolynomialOptimization: @assert, @inbounds, @verbose_info, @capture, haveMPI
using PolynomialOptimization.Relaxation: RelaxationGroupings
import PolynomialOptimization.Newton: execute_taskfun, execute, verbose_worker, halfpolytope
using PolynomialOptimization.Newton: preproc, analyze, prepare, alloc_global, alloc_local, clonetask, work, step_callback

__init__() = haveMPI[] = true

abstract type MPIRank end

struct RootRank <: MPIRank end
struct OtherRank <: MPIRank
    rank::Int
end

const root = 0

Base.convert(::Type{<:Integer}, ::RootRank) = root
Base.convert(::Type{<:Integer}, rank::OtherRank) = rank.rank
Base.:(==)(rank::MPIRank, with::Integer) = convert(typeof(with), rank) == with
Base.:(==)(with::Integer, rank::MPIRank) = with == convert(typeof(with), rank)
isroot(::MPIRank) = false
isroot(::RootRank) = true
MPIRank(rank::Integer) = iszero(rank) ? RootRank() : OtherRank(rank)

include("./Newton.jl")
include("./Worker.jl")

end