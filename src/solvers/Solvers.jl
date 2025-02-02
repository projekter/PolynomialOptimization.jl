module Solvers

using LinearAlgebra
using ..PolynomialOptimization: @assert, @inbounds

include("helpers/Helpers.jl")

include("SpecBM/SpecBM.jl")
include("Lancelot/Lancelot.jl")
include("Loraine/Loraine.jl")
include("LoRADS/LoRADS.jl")
#include("SketchyCGAL/SketchyCGAL.jl") << this is the single-matrix variant that directly follows the paper
include("SketchyCGAL/SketchyCGALBlock.jl")

end

include("SpecBM/PolynomialOptimizationSpecBM.jl")
include("Loraine/PolynomialOptimizationLoraine.jl")
include("LoRADS/PolynomialOptimizationLoRADS.jl")