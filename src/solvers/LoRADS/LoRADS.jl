# This is an interface to some of the functions of the experimental LoRADS solver, https://github.com/COPT-Public/LoRADS,
# tightly integrated with the PolynomialOptimization framework
module LoRADS

export set_solverlib, init_solver, set_dual_objective, conedata_to_userdata, set_cone, set_lp_cone, init_cone_data, preprocess,
    load_sdpa, solver, get_X, get_Xlin

using LinearAlgebra: chkstride1, Transpose, Symmetric
using SparseArrays, Preferences
using LinearAlgebra.BLAS: axpy!, syrk!

const solverlib = @load_preference("lorads-solver", "")

function __init__()
    !isempty(solverlib) && let dl=Libc.dlopen(solverlib, throw_error=false)
        if isnothing(dl)
            @warn("The LoRADS library is configured to $solverlib, but it could not be opened. Call \
                   `PolynomialOptimization.Solvers.LoRADS.set_solverlib` to change the configuration; set it to an empty \
                   value to disable the solver.")
        else
            isnothing(Libc.dlsym(dl, :ASDPSetLpCone, throw_error=false)) &&
                @warn("The unpatched version of the LoRADS library is used. Expect segfaults.")
            Libc.dlclose(dl)
        end
    end
    return
end

"""
    set_solverlib(path)

Changes the path to the LoRADS library, which takes effect after restarting Julia.
"""
function set_solverlib(path::String)
    if isfile(path)
        @set_preferences!("lorads-solver" => path)
        @info("New path to solver library set; restart you Julia session for this change to take effect!")
        return
    else
        throw(ArgumentError("Library not found: $path"))
    end
end

include("./Enums.jl")
include("./DataTypes.jl")

check(r::Retcode) = r === RETCODE_OK || error("ASDP error: $r")

"""
    init_solver(solver, nConstrRows, coneDims::AbstractVector{Cint}, nLpCols;
        rho=1.5, rhoMax=5000., maxIter=10_000, strategy=STRATEGY_DEFAULT)

Initializes a fresh `ASDP` object with `nConstrRows` constraints, positive semidefinite variables of side dimension `coneDims`
(a vector of `Cint`), `nLpCols` scalar nonnegative variables, penalty parameter `rho` (zero will set it to ``\frac{1}{n}``,
where `n` is the sum of all side dimensions plus the number of linear variables) which can be dynamically increased up to at
most `rhoMax`, total number of iterations `maxIter`, and strategy to adapt `rho`.
"""
function init_solver(solver::ASDP, nConstrRows::Integer, coneDims::AbstractVector{Cint}, nLpCols::Integer; rho::Real=1.5,
    rhoMax::Real=5000., maxIter::Integer=10_000, strategy::Strategy=STRATEGY_DEFAULT)
    getfield(solver, :init_called) && error("Double initialization")
    chkstride1(coneDims)
    setfield!(solver, :init_called, true)
    check(@ccall solverlib.ASDPInitSolver(solver::Ptr{Cvoid}, nConstrRows::Cint, length(coneDims)::Cint, coneDims::Ptr{Cint},
        nLpCols::Cint, rho::Cdouble, rhoMax::Cdouble, maxIter::Cint, strategy::Cint)::Retcode)
    setfield!(solver, :init_success, true)
    return solver
end

function cleanup(solver::ASDP)
    getfield(solver, :init_success) || return solver
    @ccall solverlib.ASDPDestroyADMMVars(solver::Ptr{Cvoid})::Cvoid
    @ccall solverlib.ASDPDestroyBMVars(solver::Ptr{Cvoid})::Cvoid
    @ccall solverlib.freeDetectSparsitySDPCoeff(solver::Ptr{Cvoid})::Cvoid
    @ccall solverlib.ASDPDestroyRankElem(solver::Ptr{Cvoid})::Cvoid
    @ccall solverlib.destroyPreprocess(solver::Ptr{Cvoid})::Cvoid
    @ccall solverlib.ASDPDestroyConeData(solver::Ptr{Cvoid})::Cvoid
    setfield!(solver, :init_success, false)
    return solver
end

function destroy_solver(solver::ASDP)
    getfield(solver, :init_success) && cleanup(x)
    @ccall solverlib.ASDPDestroySolver(solver::Ptr{Cvoid})::Cvoid
    setfield!(solver, :init_called, false)
    return solver
end

"""
    set_dual_objective(solver, dObj::AbstractVector{Cdouble})

Sets the dual objective, i.e., the right-hand side of the constraints, in an initialized `ASDP` object.
"""
function set_dual_objective(solver::ASDP, dObj::AbstractVector{Cdouble})
    chkstride1(dObj)
    @ccall solverlib.ASDPSetDualObjective(solver::Ptr{Cvoid}, dObj::Ptr{Cdouble})::Cvoid
    return solver
end

"""
    conedata_to_userdata(cone::ConeType, nConstrRows, dim, coneMatBeg::AbstractVector{Cint},
        coneMatIdx::AbstractVector{Cint}, coneMatElem::AbstractVector{Cdouble})

Allocates a new user data objects and sets its conic data. This consists of a cone type (only `ASDP_CONETYPE_DENSE_SDP` and
`ASDP_CONETYPE_SPARSE_SDP` are supported), the number of rows (which is the same as the number of constraints in the solver)
and the side dimension of the semidefinite variable, followed by the constraint matrices in zero-indexed CSR format. Every row
corresponds to the vectorized lower triangle of the column of a constraint matrix. The zeroth row is the coefficient matrix for
the objective.
Therefore, `nConstrRows +1 = length(coneMatBeg) -1` should hold (`+1` for the objective; `-1` for CSR).

The returned userdata pointer should be assigned to a solver, which will take care of freeing the allocated data. Note that the
vectors passed to this function must be preserved until the [`preprocess`](@ref) function was called, after which they can be
freed.

See also [`set_cone`](@ref), [`init_cone_data`](@ref).
"""
function conedata_to_userdata(cone::ConeType, nConstrRows::Integer, dim::Integer, coneMatBeg::AbstractVector{Cint},
    coneMatIdx::AbstractVector{Cint}, coneMatElem::AbstractVector{Cdouble})
    chkstride1(coneMatBeg)
    chkstride1(coneMatIdx)
    chkstride1(coneMatElem)
    nConstrRows == length(coneMatBeg) -2 ||
        throw(ArgumentError("The number of constraint rows is not compatible with the given matrix."))
    (length(coneMatIdx) == length(coneMatElem) && iszero(coneMatBeg[begin]) && coneMatBeg[end] == length(coneMatIdx) &&
        issorted(coneMatBeg)) || throw(ArgumentError("Error in the CSR format"))
    result = Ref{Ptr{Cvoid}}()
    check(@ccall solverlib.AUserDataCreate(result::Ref{Ptr{Cvoid}})::Retcode)
    @ccall solverlib.AUserDataSetConeData(result[]::Ptr{Cvoid}, cone::Cint, nConstrRows::Cint, dim::Cint,
        coneMatBeg::Ptr{Cint}, coneMatIdx::Ptr{Cint}, coneMatElem::Ptr{Cdouble})::Cvoid
    return result[]
end

"""
    set_cone(solver, iCone, userCone)

Sets the `iCone`th cone to the data previously defined using [`conedata_to_userdata`](@ref).

See also [`init_cone_data`](@ref).
"""
function set_cone(solver::ASDP, iCone::Integer, userCone::Ptr{Cvoid})
    check(@ccall solverlib.ASDPSetCone(solver::Ptr{Cvoid}, iCone::Cint, userCone::Ptr{Cvoid})::Retcode)
    return solver
end

"""
    set_lp_cone(solver, nConstrRows, nLpCols, lpMatBeg::AbstractVector{Cint},
        lpMatIdx::AbstractVector{Cint}, lpMatElem::AbstractVector{Cdouble})

Set the data of the constraint matrix for the linear variables according to the CSR data specified in the parameters.

!!! warning
    This function is not exported on the original code release and can therefore not be used. However, only the patched version
    should be used, as it fixes heap corruption errors that can arise during the optimization.
"""
function set_lp_cone(solver::ASDP, nConstrRows::Integer, nLpCols::Integer, lpMatBeg::AbstractVector{Cint},
    lpMatIdx::AbstractVector{Cint}, lpMatElem::AbstractVector{Cdouble})
    chkstride1(lpMatBeg)
    chkstride1(lpMatIdx)
    chkstride1(lpMatElem)
    nConstrRows == length(lpMatBeg) -2 ||
        throw(ArgumentError("The number of constraint rows is not compatible with the given matrix."))
    (length(lpMatIdx) == length(lpMatElem) && iszero(lpMatBeg[begin]) && lpMatBeg[end] == length(lpMatIdx) &&
        issorted(lpMatBeg)) || throw(ArgumentError("Error in the CSR format"))
    check(@ccall solverlib.ASDPSetLpCone(solver.lpCone::Ptr{Cvoid}, nConstrRows::Cint, nLpCols::Cint, lpMatBeg::Ptr{Cint},
        lpMatIdx::Ptr{Cint}, lpMatElem::Ptr{Cdouble})::Retcode)
    return solver
end

@doc raw"""
    init_cone_data(solver, coneMat, coneDims, lpMat)

Initializes the solver for a problem in the form
```math
   \min \vec a_0 \cdot \vec x + \sum_j \langle\operatorname{mat}(G_{j, 0}), Z_j\rangle \\
   \text{such that} \\
   x_i \geq 0 \ \forall i \\
   Z_j \succeq 0 \ \forall j \\
   \vec a_k \cdot \vec x - \sum_j \langle\operatorname{mat}(G_{j, k}, Z_j)\rangle = c_k \ \forall k
```
with the following representation in the variables:
- `1 ≤ j ≤ length(coneDims) = length(coneMat)`
- `coneMat` is a vector of matrices, `lpMat` is a matrix. They should be in CSR storage, where the row index (starting at 0 for
  the objective, then `k` for the `k`th constraint) is the constraint. Since CSR is not natively supported by Julia, the
  transpose of a `SparseMatrixCSC{Cdouble,Cint}` is expected.
- `mat` makes the unscaled lower triangle into a full matrix

This is a convenience function that does the job of [`conedata_to_userdata`](@ref), [`set_cone`](@ref), and
[`preprocess`](@ref) in one step. However, note that it is more efficient to call these functions individually.
"""
function init_cone_data(solver::ASDP, coneMat::AbstractVector{Transpose{Cdouble,SparseMatrixCSC{Cdouble,Cint}}},
    coneDims::AbstractVector{Cint}, lpMat::Union{Nothing,Transpose{Cdouble,SparseMatrixCSC{Cdouble,Cint}}})
    chkstride1(coneDims)
    if isempty(coneMat)
        isnothing(lpMat) && throw(ArgumentError("No data"))
        nConstrs = size(lpMat, 1)
    else
        nConstrs = size(first(coneMat), 1)
    end
    length(coneDims) == length(coneMat) || throw(ArgumentError("Lengths of cones are inconsistent"))
    coneMatElem = Vector{Ptr{Cdouble}}(undef, length(coneMat))
    coneMatBeg = similar(coneMatElem, Ptr{Cint})
    coneMatIdx = similar(coneMatElem, Ptr{Cint})
    @inbounds for (i, cm) in enumerate(coneMat)
        size(cm, 1) == nConstrs || throw(ArgumentError("Number of constraints inconsistent"))
        # cm is the transpose of CSC, i.e., CSR. Each row in cm is a constraint, each column an entry in the lower triangle
        coneDims[i] * (coneDims[i] +1) ÷ 2 == size(cm, 2) ||
            throw(ArgumentError(lazy"Length of cone $i is wrong (is $(coneDims[i])), should be $((isqrt(1 + 8size(cm, 2)) -1) ÷ 2)"))
        colptr = SparseArrays.getcolptr(parent(cm))
        colptr .-= one(Cint)
        rowval = SparseArrays.getrowval(parent(cm))
        rowval .-= one(Cint)
        coneMatBeg[i] = pointer(colptr)
        coneMatIdx[i] = pointer(rowval)
        coneMatElem[i] = pointer(SparseArrays.getnzval(parent(cm)))
    end
    if !isnothing(lpMat)
        size(lpMat, 1) == nConstrs || throw(ArgumentError("Number of constraints inconsistent"))
        nLpCols = size(lpMat, 2)
        LpMatBeg = SparseArray.getcolptr(lpMat)
        LpMatIdx = SparseArray.getrowval(lpMat)
        LpMatElem = SparseArray.getnzval(lpMat)
        LpMatBeg .-= one(Cint)
        LpMatIdx .-= one(Cint)
    else
        nLpCols = 0
        LpMatBeg = C_NULL
        LpMatIdx = C_NULL
        LpMatElem = C_NULL
    end
    sdpDatas = similar(coneMatElem, Ptr{Cvoid})
    check(@ccall solverlib.ASDPInitConeData(
        solver::Ptr{Cvoid},
        sdpDatas::Ptr{Ptr{Cvoid}},
        coneMatElem::Ptr{Ptr{Cdouble}},
        coneMatBeg::Ptr{Ptr{Cint}},
        coneMatIdx::Ptr{Ptr{Cint}},
        coneDims::Ptr{Cint},
        nConstrs::Cint,
        length(coneMat)::Cint,
        nLpCols::Cint,
        LpMatBeg::Ptr{Cint},
        LpMatIdx::Ptr{Cint},
        LpMatElem::Ptr{Cdouble}
    )::Retcode)
    # sdpDatas is filled. Every cone in solver has a field usrData that contains the same pointer. And when the cones are
    # destroyed when finalizing solver, this userdata is also freed (just the struct that holds the references, not the arrays
    # themselves). So we can just forget about our sdpDatas.
    preprocess(solver, coneDims)
    # Now all the data was converted to internal data. We don't need any of them any more.
    @inbounds for (i, cm) in enumerate(coneMat)
        SparseArrays.getcolptr(parent(cm)) .+= one(Cint)
        SparseArrays.getrowval(parent(cm)) .+= one(Cint)
    end
    if !isnothing(lpMat)
        LpMatBeg .+= one(Cint)
        LpMatIdx .+= one(Cint)
    end
    return solver
end

"""
    preprocess(solver, coneDims::AbstractVector{Cint})

Invokes the preprocessor. This should be called after all cones were set up, after which their original data may be reused or
destroyed.

See also [`conedata_to_userdata`](@ref), [`set_cone`](@ref).
"""
function preprocess(solver::ASDP, coneDims::AbstractVector{Cint})
    chkstride1(coneDims)
    check(@ccall solverlib.ASDPPreprocess(solver::Ptr{Cvoid}, coneDims::Ptr{Cint})::Retcode)
    return solver
end

function determine_rank(solver::ASDP, coneDims::AbstractVector{Cint}, timesRank::Real, rankSpecify::Integer)
    chkstride1(coneDims)
    check(@ccall solverlib.ASDPDetermineRank(solver::Ptr{Cvoid}, coneDims::Ptr{Cint}, timesRank::Cdouble,
        rankSpecify::Cint)::Retcode)
    return solver
end

function detect_max_cut_prob(solver::ASDP, coneDims::AbstractVector{Cint})
    chkstride1(coneDims)
    maxCut = Ref{Cint}()
    @ccall solverlib.detectMaxCutProb(solver::Ptr{Cvoid}, coneDims::Ptr{Cint}, maxCut::Ref{Cint})::Cvoid
    return maxCut[]
end

function detect_sparsity_sdp_coeff(solver::ASDP)
    @ccall solverlib.detectSparsitySDPCoeff(solver::Ptr{Cvoid})::Cvoid
    return solver
end

function init_bm_vars(solver::ASDP, coneDims::AbstractVector{Cint}, nLpCols::Integer)
    chkstride1(coneDims)
    check(@ccall solverlib.ASDPInitBMVars(solver::Ptr{Cvoid}, coneDims::Ptr{Cint}, length(coneDims)::Cint,
        nLpCols::Cint)::Retcode)
    return solver
end

function init_admm_vars(solver::ASDP, coneDims::AbstractVector{Cint}, nLpCols::Integer)
    chkstride1(coneDims)
    check(@ccall solverlib.ASDPInitADMMVars(solver::Ptr{Cvoid}, coneDims::Ptr{Cint}, length(coneDims)::Cint,
        nLpCols::Cint)::Retcode)
    return solver
end

function scale(solver::ASDP)
    @ccall solverlib.ASDP_SCALE(solver::Ptr{Cvoid})::Cvoid
    return solver
end

bm_optimize(solver::ASDP, endBMTol::Real, endBMTol_pd::Real, endTauTol::Real, endBMALSub::Real, ori_start::Real,
    is_rank_max::Bool, pre_mainiter::Ref{Cint}, pre_miniter::Ref{Cint}, timeLimit::Real) =
    @ccall solverlib.ASDP_BMOptimize(solver::Ptr{Cvoid}, endBMTol::Cdouble, endBMTol_pd::Cdouble, endTauTol::Cdouble,
        endBMALSub::Cdouble, ori_start::Cdouble, is_rank_max::Cint, pre_mainiter::Ref{Cint}, pre_miniter::Ref{Cint},
        timeLimit::Cdouble)::Retcode

check_all_rank_max(solver::ASDP, aug_factor::Real) =
    !iszero(@ccall solverlib.CheckAllRankMax(solver::Ptr{Cvoid}, aug_factor::Cdouble)::Cint)

function aug_rank(solver::ASDP, coneDims::AbstractVector{Cint}, aug_factor::Real)
    chkstride1(coneDims)
    return !iszero(@ccall solverlib.AUG_RANK(solver::Ptr{Cvoid}, coneDims::Ptr{Cint}, length(coneDims)::Cint,
        aug_factor::Cdouble)::Cint)
end

function bm_to_admm(solver::ASDP, heuristic::Real)
    @ccall solverlib.ASDP_BMtoADMM(solver::Ptr{Cvoid}, heuristic::Cdouble)::Cvoid
    return solver
end

optimize(solver::ASDP, rhoFreq::Integer, rhoFactor::Real, rhoStrategy::Strategy, tau::Real, gamma::Real, rhoMin::Real,
    orig_start::Real, timeLimit::Real) =
    @ccall solverlib.ASDPOptimize(solver::Ptr{Cvoid}, rhoFreq::Cint, rhoFactor::Cdouble, rhoStrategy::Cint, tau::Cdouble,
        gamma::Cdouble, rhoMin::Cdouble, orig_start::Cdouble, timeLimit::Cdouble)::Retcode

function end_program(solver::ASDP)
    @ccall solverlib.ASDPEndProgram(solver::Ptr{Cvoid})::Cvoid
    return solver
end

get_timestamp() = @ccall solverlib.AUtilGetTimeStamp()::Cdouble

"""
    load_sdpa(fn; maxIter=10_000, rho=0., rhoMax=5000., strategy=STRATEGY_DEFAULT)

Loads a problem from a file `fn` in SDPA format and returns a preprocessed `ASDP` instance, a vector containing the cone
dimensions, and the number of nonnegative scalar variables.

!!! warning
    This function will produce memory leaks. The data that is allocated by the LoRADS library cannot be freed by Julia, as no
    corresponding functions are exported. Only use it for quick tests, not in production.
"""
function load_sdpa(fn::AbstractString; maxIter::Integer=10_000, rho::Real=0., rhoMax::Real=5000.,
    strategy::Strategy=STRATEGY_DEFAULT)
    isempty(solverlib) &&
        error("LoRADS is not configured. Call `PolynomialOptimization.Solvers.LoRADS.set_solverlib` to specify the location \
               of the solver library")
    fn_str = Vector{UInt8}(fn)
    fn_str[end] == 0x00 || push!(fn_str, 0x00)
    pnConstrs = Ref{Cint}()
    pnBlks = Ref{Cint}()
    pblkDims = Ref{Ptr{Cint}}()
    prowRHS = Ref{Ptr{Cdouble}}()
    pconeMatBeg = Ref{Ptr{Ptr{Cint}}}()
    pconeMatIdx = Ref{Ptr{Ptr{Cint}}}()
    pconeMatElem = Ref{Ptr{Ptr{Cdouble}}}()
    pnLpCols = Ref{Cint}()
    pLpMatBeg = Ref{Ptr{Cint}}()
    pLpMatIdx = Ref{Ptr{Cint}}()
    pLpMatElem = Ref{Ptr{Cdouble}}()

    check(@ccall solverlib.AReadSDPA(fn_str::Ptr{Cvoid}, pnConstrs::Ptr{Cvoid}, pnBlks::Ptr{Cvoid}, pblkDims::Ptr{Cvoid},
        prowRHS::Ptr{Cvoid}, pconeMatBeg::Ptr{Cvoid}, pconeMatIdx::Ptr{Cvoid}, pconeMatElem::Ptr{Cvoid},
        Ref{Cint}()::Ptr{Cvoid}, pnLpCols::Ptr{Cvoid}, pLpMatBeg::Ptr{Cvoid}, pLpMatIdx::Ptr{Cvoid},
        pLpMatElem::Ptr{Cvoid}, Ref{Cint}()::Ptr{Cvoid})::Retcode)

    coneDims = unsafe_wrap(Vector{Cint}, pblkDims[], pnBlks[], own=false)
    rowRHS = unsafe_wrap(Vector{Cdouble}, prowRHS[], pnConstrs[], own=false)

    solver = ASDP()
    init_solver(solver, pnConstrs[], coneDims, pnLpCols[]; rho, rhoMax, maxIter, strategy)
    set_dual_objective(solver, rowRHS)
    sdpDatas = Vector{Ptr{Cvoid}}(undef, pnBlks[])
    check(@ccall solverlib.ASDPInitConeData(
        solver::Ptr{Cvoid},
        sdpDatas::Ptr{Ptr{Cvoid}},
        pconeMatElem[]::Ptr{Ptr{Cdouble}},
        pconeMatBeg[]::Ptr{Ptr{Cint}},
        pconeMatIdx[]::Ptr{Ptr{Cint}},
        coneDims::Ptr{Cint},
        pnConstrs[]::Cint,
        pnBlks[]::Cint,
        pnLpCols[]::Cint,
        pLpMatBeg[]::Ptr{Cint},
        pLpMatIdx[]::Ptr{Cint},
        pLpMatElem[]::Ptr{Cdouble}
    )::Retcode)
    preprocess(solver, coneDims)
    return solver, coneDims, pnLpCols[]
end

"""
    solve(solver, coneDims, nLpCols; timesLogRank=2., phase1Tol=1e-3, time_limit=10_000.,
        rhoFreq=5, rhoFactor=1.2, admmStrategy=STRATEGY_MIN_BISECTION, tau=0., gamma=0.,
        rankFactor=1.5)

Solves a preprocessed `ASDP` instance.
"""
function solve(solver::ASDP, coneDims::AbstractVector{Cint}, nLpCols::Integer; timesLogRank::Real=2., phase1Tol::Real=1e-3,
    time_limit::Real=10_000., rhoFreq::Integer=5, rhoFactor::Real=1.2, admmStrategy::Strategy=STRATEGY_MIN_BISECTION,
    tau::Real=0., gamma::Real=0., rankFactor::Real=1.5)
    isempty(solverlib) &&
        error("LoRADS is not configured. Call `PolynomialOptimization.Solvers.LoRADS.set_solverlib` to specify the location of the solver library")
    chkstride1(coneDims)
    determine_rank(solver, coneDims, timesLogRank, 0)
    if detect_max_cut_prob(solver, coneDims) != -1
        println("**Detected MaxCut problem: set phase1Tol -> 1e-2 and heuristicFactor -> 10")
        phase1Tol = 1e-2
        heuristic = 10.
    else
        heuristic = 1.
    end
    detect_sparsity_sdp_coeff(solver)
    init_bm_vars(solver, coneDims, nLpCols)
    init_admm_vars(solver, coneDims, nLpCols)
    scale(solver)
    println("**First using BM method as warm start")
    pre_mainiter = Ref(zero(Cint))
    pre_miniter = Ref(zero(Cint))
    time = get_timestamp()
    is_rank_max = 1e-10 > timesLogRank
    local bm_ret
    while true
        bm_ret = bm_optimize(solver, phase1Tol, -.001, 1e-16, 1e-10, time, is_rank_max, pre_mainiter, pre_miniter,
            time_limit)
        bm_ret == RETCODE_RANK || break
        if !check_all_rank_max(solver, rankFactor)
            if (is_rank_max = aug_rank(solver, coneDims, rankFactor))
                println("Restarting BM with maximum rank")
            else
                timesLogRank *= rankFactor
                println("Restarting BM with updated rank (now ", timesLogRank, " logm)")
            end
        end
    end
    if bm_ret != RETCODE_EXIT
        bm_to_admm(solver, heuristic)
        optimize(solver, rhoFreq, rhoFactor, admmStrategy, tau, gamma, 0., time, time_limit)
    end
    end_program(solver)

    return solver
end

"""
    get_X(solver, i, admm::Bool=true)

Returns the `i`th PSD solution matrix ``X_i``. The result will be a freshly allocated symmetric view of a dense matrix.
Indicate with `admm` whether the ADMM split optimization has run or if only the BM optimization finished.

!!! warning
    This method may only be called once per `i`. All further calls with the same `i` will give wrong output, as the internal
    solver data is modified.
"""
function get_X(solver::ASDP, i::Integer; admm::Bool=true)
    # The solver doesn't contain the primal matrix in full form, but in a low-rank factorization.
    V = unsafe_load(unsafe_load(solver.V, i))
    if admm
        # Since the ADMM approach uses a U Vᵀ, V = U rewrite with additional constraints, even then it is not just U Uᵀ. So the
        # recommended reconstruction is Û = (U + V)/2, X = Û Ûᵀ.
        U = unsafe_load(unsafe_load(solver.U, i))
        @assert U.nRows == V.nRows && U.rank == V.rank
        axpy!(V.nRows * V.rank, one(Cdouble), U.matElem, 1, V.matElem, 1)
    end
    result = Matrix{Cdouble}(undef, V.nRows, V.nRows)
    Vwrap = unsafe_wrap(Array, V.matElem, (V.nRows, V.rank))
    syrk!('L', 'N', admm ? .25 : 1., Vwrap, false, result)
    return Symmetric(result, :L)
end

"""
    get_Xlin(solver, admm::Bool=true)

Returns the linear solution vector ``x``. The result will be a vector backed by internal solver data and will be invalidated if
the solver is destroyed. Copy it if desired.
Indicate with `admm` whether the ADMM split optimization has run or if only the BM optimization finished.

!!! warning
    This method may only be called once. All further calls will give wrong output, as the internal solver data is modified.
"""
function get_Xlin(solver::ASDP; admm::Bool=true)
    v = unsafe_load(solver.vLp)
    vWrap = unsafe_wrap(Array, v.matElem, v.nLPCols)
    if admm
        u = unsafe_load(solver.uLp)
        @assert u.nLPCols == v.nLPCols
        uWrap = unsafe_wrap(Array, u.matElem, u.nLPCols)
        vWrap .= .25 .* (uWrap .+ vWrap) .^ 2
    else
        vWrap .^= 2
    end
    return vWrap
end

"""
    get_S(solver, i)

Returns the `i`th slack variable for the PSD solution matrix ``S_i``. The result will be a freshly allocated symmetric view of
a dense matrix.
"""
function get_S(solver::ASDP, i::Integer)
    cone = unsafe_load(unsafe_load(solver.ACones, i))
    dim = unsafe_load(Ptr{Cint}(cone.sdp_slack_var)) # first field is nSDPCol
    S = Matrix{Cdouble}(undef, dim, dim)
    @ccall solverlib.reconstructSymmetricMatrix(cone.sdp_slack_var::Ptr{Cvoid}, S::Ptr{Cdouble}, dim::Cint)::Cvoid
    return Symmetric(S, :L)
end

"""
    get_Slin(solver[, i::AbstractUnitRange])

Returns the slack variables to the nonnegative variables of index `i` (by default, all). The result will be a freshly allocated
vector.
"""
function get_Slin(solver::ASDP, i::AbstractUnitRange=1:solver.nLpCols)
    cone = unsafe_load(solver.lpCone)
    # There's no dual data available, we have to construct it ourselves
    s = Vector{Cdouble}(undef, length(i))
    unsafe_copyto!(pointer(s), unsafe_load(cone.coneData).objMatElem + (first(i) -1) * sizeof(Cdouble), length(s))
    # ^ so that we don't have to call objCoeffSum in a loop
    negLambd = solver.negLambd
    for idx in i
        ccall(cone.lpDataWSum, Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Cint), cone.coneData, negLambd, s, idx -1)
    end
    return s
end

end