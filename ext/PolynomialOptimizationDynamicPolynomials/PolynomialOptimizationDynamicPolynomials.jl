module PolynomialOptimizationDynamicPolynomials

import DynamicPolynomials

function __init__()
    # this implementation is buggy; we define one that should supersede it in the ArrayPolynomials.jl helper
    # https://github.com/JuliaAlgebra/DynamicPolynomials.jl/issues/171
    try
        m = which(DynamicPolynomials.variables, (AbstractArray{<:DynamicPolynomials.PolyType},))
        m.module === DynamicPolynomials && Base.delete_method(m)
    catch e
        e isa MethodError || rethrow(e)
    end
    return
end

end