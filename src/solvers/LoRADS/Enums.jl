@enum Retcode RETCODE_OK RETCODE_FAILED RETCODE_MEMORY RETCODE_EXIT RETCODE_RANK RETCODE_SPLIT

@enum Status begin
    ASDP_UNKNOWN
    ASDP_DUAL_FEASIBLE
    ASDP_DUAL_OPTIMAL
    ASDP_PRIMAL_DUAL_OPTIMAL
    ASDP_PRIMAL_OPTIMAL
    ASDP_MAXITER
    ASDP_SUSPECT_INFEAS_OR_UNBOUNDED
    ASDP_INFEAS_OR_UNBOUNDED
    ASDP_TIMELIMIT
    ASDP_USER_INTERRUPT
    ASDP_INTERNAL_ERROR
    ASDP_NUMERICAL
end

@enum ConstrType ASDP_DENSE ASDP_SPARSE

@enum Method ADMMMethod BMMethod

"""
    Strategy

A strategy to employ for the adaptive change in rho. Possible values are:
- `STRATEGY_DEFAULT`
- `STRATEGY_MIN_BISECTION`: divide ``\\rho`` by two if the constraint violations are sufficiently small. This is the official
  termination scheme, which is currently disabled.
- `STRATEGY_MAX_CUT`
- `STRATEGY_SNL`
- `STRATEGY_OPF`
- `STRATEGY_GENERAL`
- `STRATEGY_GAMMA_DYNAMIC`: multiply ``\\rho`` by ``\\gamma`` if the constraint violations exceed ``tau`` times the gradient
  norm. This is applied in the ADMM optimizer if it is given as a strategy for [`solve`](@ref).

Apart from the noted behavior, all strategies behave the same.
"""
@enum Strategy begin
    STRATEGY_DEFAULT
    STRATEGY_MIN_BISECTION
    STRATEGY_MAX_CUT
    STRATEGY_SNL
    STRATEGY_OPF
    STRATEGY_GENERAL
    STRATEGY_GAMMA_DYNAMIC
end

"""
    ConeType

The type of the current cone. Only `ASDP_CONETYPE_DENSE_SDP` and `ASDP_CONETYPE_SPARSE_SDP` are implemented.
"""
@enum ConeType begin
    ASDP_CONETYPE_UNKNOWN
    ASDP_CONETYPE_LP
    ASDP_CONETYPE_BOUND
    ASDP_CONETYPE_SCALAR_BOUND
    ASDP_CONETYPE_DENSE_SDP
    ASDP_CONETYPE_SPARSE_SDP
    ASDP_CONETYPE_SOCP
end

@enum SDPCoeffType begin
    SDP_COEFF_ZERO
    SDP_COEFF_SPARSE
    SDP_COEFF_DENSE
    SDP_COEFF_SPR1
    SDP_COEFF_DSR1
end