struct ASDPRkMatDense
    rank::Cint # rank of the matrix, i.e., column number of the matrix
    nRows::Cint
    matElem::Ptr{Cdouble}
end

struct ASDPRkMatLp
    nLPCols::Cint
    matElem::Ptr{Cdouble}
end

struct ConstrValSparse
    nnz::Cint
    nnzIdx::Ptr{Cint}
    val::Ptr{Cdouble}
end

struct ConstrValDense
    nnz::Cint
    val::Ptr{Cdouble}
end

struct ConstrValStruct
    constrVal::Ptr{Cvoid}
    type::ConstrType
end

struct ASDPCone
    cone::ConeType
    sdp_coeff_w_sum::Ptr{Cvoid}
    sdp_obj_sum::Ptr{Cvoid}
    sdp_slack_var::Ptr{Cvoid}
    UVt_w_sum::Ptr{Cvoid}
    UVt_obj_sum::Ptr{Cvoid}

    usrData::Ptr{Cvoid}
    coneData::Ptr{Cvoid}

    sdp_coeff_w_sum_sp_ratio::Cdouble
    nConstr::Cint
end

struct ASDPConeLp
    nRow::Cint
    nCol::Cint
    objMatElem::Ptr{Cdouble}
    rowMatBeg::Ptr{Cint}
    rowMatIdx::Ptr{Cint}
    rowMatElem::Ptr{Cdouble}
    lpCol::Ptr{Ptr{Cvoid}}
    nrm2Square::Ptr{Cdouble}
end

struct ASDPLpCone
    coneData::Ptr{ASDPConeLp}
    nCol::Cint
    coneCreate::Ptr{Cvoid}
    coneProcData::Ptr{Cvoid}
    conePresolveData::Ptr{Cvoid}
    coneDestroyData::Ptr{Cvoid}
    coneView::Ptr{Cvoid}
    coneAUV::Ptr{Cvoid}
    coneAUV2::Ptr{Cvoid}
    objAUV::Ptr{Cvoid}
    coneObjNrm1::Ptr{Cvoid}
    coneObjNrm2Square::Ptr{Cvoid}
    coneObjNrmInf::Ptr{Cvoid}
    lpDataWSum::Ptr{Cvoid}
    objCoeffSum::Ptr{Cvoid}
    scaleData::Ptr{Cvoid}
end

struct _ASDPSolverInternal
    # User data
    nRows::Cint # constraint number
    rowRHS::Ptr{Cdouble} # b of Ax = b

    # Cones
    nCones::Cint # sdp cones block number
    ACones::Ptr{Ptr{ASDPCone}}
    CGLinsys::Ptr{Cvoid}
    LanczosSys::Ptr{Cvoid}
    LanczosStart::Ptr{Ptr{Cdouble}}

    # Auxiliary variable SDPCone
    constrVal::Ptr{Ptr{ConstrValStruct}} # constraint violation [iCone][iConstr]
    constrValSum::Ptr{Cdouble}           # constraint violation [iConstr]
    ARDSum::Ptr{Cdouble}                 # q1 in line search of BM
    ADDSum::Ptr{Cdouble}                 # q2 in line search of BM
    bLinSys::Ptr{Ptr{Cdouble}}           # for solving linear system, bInLinSys[iCone]
    rankElem::Ptr{Cint}                  # all rank
    M1temp::Ptr{Cdouble}                 # M1 for solving linear system
    bestDualVar::Ptr{Cdouble}
    M2temp::Ptr{Ptr{ASDPRkMatDense}}

    # Variables SDPCone
    U::Ptr{Ptr{ASDPRkMatDense}}    # admm variable, and lbfgs descent direction D
    V::Ptr{Ptr{ASDPRkMatDense}}    # admm variable only
    R::Ptr{Ptr{ASDPRkMatDense}}    # average variable for storage and BM variable
    Grad::Ptr{Ptr{ASDPRkMatDense}} # grad of R

    # Auxiliary variable LPCone (SDPCone but rank is 1, dim is 1)
    rLp::Ptr{ASDPRkMatLp}
    uLp::Ptr{ASDPRkMatLp}
    vLp::Ptr{ASDPRkMatLp}
    vlagLp::Ptr{ASDPRkMatLp}
    lagLp::Ptr{ASDPRkMatLp}
    gradLp::Ptr{ASDPRkMatLp}
    lpCone::Ptr{ASDPLpCone}
    nLpCols::Cint
    constrValLP::Ptr{Ptr{ConstrValStruct}}

    # BM lbfgs and ADMM Variables
    dualVar::Ptr{Cdouble}
    Vlag::Ptr{Ptr{ASDPRkMatDense}} # ADMM dual variable and lbfgs gradient
    lag::Ptr{Ptr{ASDPRkMatDense}}

    # BM lbfgs Variables
    hisRecT::Cint          # all node number is hisRecT + 1
    lbfgsHis::Ptr{Cvoid}   # record difference of primal variable history and gradient history, lbfgsHis[iCone]
    maxBMInIter::Cint      # inner iteration
    maxBMOutIter::Cint     # outer interation, increase penalty

    # ADMM Parameters
    maxIter::Cint

    # Monitor
    whichMethod::Method
    nIterCount::Cint
    cgTime::Cdouble
    cgIter::Cint
    checkSolTimes::Cint
    traceSum::Cdouble

    # Convergence criterion
    pObjVal::Cdouble
    dObjVal::Cdouble
    pInfeas::Cdouble
    dInfeas::Cdouble
    cObjNrm1::Cdouble
    cObjNrm2::Cdouble
    cObjNrmInf::Cdouble
    bRHSNrm1::Cdouble
    bRHSNrmInf::Cdouble
    bRHSNrm2::Cdouble
    dimacError::Ptr{Cdouble}
    constrVio::Ptr{Cdouble}
    UsubV::Ptr{Ptr{Cdouble}}
    uSubvLp::Ptr{Cdouble}
    negLambd::Ptr{Cdouble}

    AStatus::Status

    # Starting time
    dTimeBegin::Cdouble

    # Parameters
    rho::Cdouble # BM and ADMM
    rhoMax::Cdouble
    strategy::Cint

    # scale
    cScaleFactor::Cdouble
    bScaleFactor::Cdouble

    # check exit bm
    rank_max::Ptr{Cint}
    sparsitySDPCoeff::Ptr{Cdouble}
    overallSparse::Cdouble
    nnzSDPCoeffSum::Cint
    SDPCoeffSum::Cint
end

mutable struct ASDP
    ptr::Ptr{Cvoid}
    init_called::Bool
    init_success::Bool

    @doc """
        ASDP()

    Creates a new empty LoRADS solver object. This has to be initialized by called [`init_solver`](@ref).
    All internal properties of the solver can be retrieved by accessing the properties of this object. Note that this assumes
    the solver interface of version 1.0.0. Any other version might be broken.
    """
    function ASDP()
        result = new(C_NULL, false, false)
        check(@ccall solverlib.ASDPSolverCreate(result::Ptr{Ptr{Cvoid}})::Retcode)
        finalizer(result) do x
            getfield(x, :init_success) && cleanup(x)
            getfield(x, :init_called) && destroy_solver(x)
            @ccall solverlib.ASDPDestroy(x::Ptr{Ptr{Cvoid}})::Cvoid
        end
        result
    end
end

Base.unsafe_convert(::Type{Ptr{Cvoid}}, x::ASDP) = getfield(x, :ptr)
Base.unsafe_convert(::Type{Ptr{Ptr{Cvoid}}}, x::ASDP) = Ptr{Ptr{Cvoid}}(pointer_from_objref(x))

Base.getproperty(solver::ASDP, f::Symbol) =
    unsafe_load(
        Ptr{fieldtype(_ASDPSolverInternal, f)}(getfield(solver, :ptr)) +
        fieldoffset(_ASDPSolverInternal, Base.fieldindex(_ASDPSolverInternal, f)),
    )
Base.propertynames(solver::ASDP) = fieldnames(_ASDPSolverInternal)
Base.setproperty!(solver::ASDP, f::Symbol, v) =
    unsafe_store!(
        Ptr{fieldtype(_ASDPSolverInternal, f)}(getfield(solver, :ptr)) +
        fieldoffset(_ASDPSolverInternal, Base.fieldindex(_ASDPSolverInternal, f)),
        v
    )