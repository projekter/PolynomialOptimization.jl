using Test, Random
using PolynomialOptimization, PolynomialOptimization.SimplePolynomials, PolynomialOptimization.Relaxation
# ^ we don't require SimplePolynomials in the namespace, but for printing we want the guarantee that it is there
using MultivariatePolynomials
import DynamicPolynomials
import Clarabel, COPT, #= COSMO, =# Hypatia, Mosek, SCS
import StatsBase

if !@isdefined(optimize)
    optimize = true

    const all_solvers = filter(s -> occursin("Moment", string(s)) || occursin("SOS", string(s)),
        PolynomialOptimization.Solver.solver_methods) # remove all shorthands
    # Mosek requires a license to work at all. COPT will work without a license for small problems.
    try
        Mosek.maketask() do t
            Mosek.optimize(t)
        end
    catch e
        if e isa Mosek.MosekError && Mosek.MSK_RES_ERR_LICENSE.value ≤ e.rcode ≤ Mosek.MSK_RES_ERR_LICENSE_NO_SERVER_LINE.value
            filter!(s -> !occursin("Mosek", string(s)), all_solvers)
        else
            rethrow(e)
        end
    end

    # allow COPT to fail: if no license is present, the larger tests won't work
    function Test.do_test(result::Test.Threw, orig_expr)
        if result.exception isa ErrorException && result.exception.msg == "API call failed: invalid license"
            testres = Test.Broken(:test, orig_expr)
            Test.record(Test.get_testset(), testres)
        else
            @invoke Test.do_test(result::Test.ExecutionResult, orig_expr::Any)
        end
    end

    function strRep(x)
        io = IOBuffer()
        show(IOContext(io, :limit => false, :displaysize => (10, 10)), "text/plain", x)
        return String(take!(io))
    end
end

solvers = copy(all_solvers)