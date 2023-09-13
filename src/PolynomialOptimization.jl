__precompile__()
module PolynomialOptimization

using MultivariatePolynomials
using SparseArrays
using LinearAlgebra
using Printf
import SemialgebraicSets
import Graphs
import Mosek
import COSMO
import Hypatia
import ConicBundle
import Combinatorics
import DynamicPolynomials
import MutableArithmetics
import StatsBase

const sqrt2 = sqrt(2.0)

macro verbose_info(str...)
    quote
        if $(esc(:verbose))
            println($(esc.(str)...))
            flush(stdout)
        end
    end
end

include("./helpers/FastVector.jl")
include("./helpers/ComplexPolynomials.jl")
include("./helpers/MatrixPolynomials.jl")

include("./sparsity/Chordal.jl")
include("./Newton.jl")
include("./Problem.jl")
include("./Tightening.jl")
include("./sparsity/SparseAnalysis.jl")
include("./SolutionExtraction.jl")
include("./SolutionExtractionHeuristic.jl")
# Do we have Mosek version at least 10?
isdefined(Mosek, :appendafes) && include("./solvers/MosekMoment.jl")
include("./solvers/MosekSOS.jl")
include("./solvers/COSMOMoment.jl")
include("./solvers/HypatiaMoment.jl")
include("./solvers/ConicBundleSOS.jl")
include("./sparsity/SparsityNone.jl")
include("./sparsity/SparsityCorrelative.jl")
include("./sparsity/SparsityTerm.jl")
include("./sparsity/SparsityTermBlock.jl")
include("./sparsity/SparsityTermCliques.jl")
include("./sparsity/SparsityCorrelativeTerm.jl")

end