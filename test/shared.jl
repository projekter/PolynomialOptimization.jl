using Test, Random, MKL
using PolynomialOptimization, PolynomialOptimization.IntPolynomials, PolynomialOptimization.Relaxation
# ^ we don't require IntPolynomials in the namespace, but for printing we want the guarantee that it is there
using MultivariatePolynomials
import DynamicPolynomials
import Clarabel, COPT, Hypatia, GALAHAD, Mosek, ProxSDP, SCS
import StatsBase

if !@isdefined(strRep)
    function strRep(x)
        io = IOBuffer()
        show(IOContext(io, :limit => false, :displaysize => (10, 10)), "text/plain", x)
        return String(take!(io))
    end

    # Mosek requires a license to work at all. COPT will work without a license for small problems.
    const have_mosek = try
        Mosek.maketask() do t
            Mosek.optimize(t)
        end
        true
    catch e
        if e isa Mosek.MosekError && Mosek.MSK_RES_ERR_LICENSE.value ≤ e.rcode ≤ Mosek.MSK_RES_ERR_LICENSE_NO_SERVER_LINE.value
            false
        else
            rethrow(e)
        end
    end
    if !have_mosek
        findmosek = isequal(:Mosek)
        deleteat!(Newton.newton_methods, findfirst(findmosek, Newton.newton_methods))
        deleteat!(PolynomialOptimization.tightening_methods, findfirst(findmosek, PolynomialOptimization.tightening_methods))
    end

    # LoRADS requires compiled binaries.
    const have_lorads = isfile(PolynomialOptimization.Solvers.LoRADS.solverlib)

    function skipsolver(solver::Symbol)
        have_lorads || solver !== :LoRADSMoment || return true
        have_mosek || !startswith(string(solver), "Mosek") || return true
        return false
    end
end