pushnm!(a::AbstractVector, ::Missing) = missing
pushnm!(a::AbstractVector, x) = (push!(a, x); x)

mutable struct LANCELOT_problem_type
    n::Cint
    ng::Cint
    nel::Cint
    IELING::GFortranVectorDescriptorRaw{Cint}
    ISTADG::GFortranVectorDescriptorRaw{Cint}
    IELVAR::GFortranVectorDescriptorRaw{Cint}
    ISTAEV::GFortranVectorDescriptorRaw{Cint}
    INTVAR::GFortranVectorDescriptorRaw{Cint}
    ISTADH::GFortranVectorDescriptorRaw{Cint}
    ICNA::GFortranVectorDescriptorRaw{Cint}
    ISTADA::GFortranVectorDescriptorRaw{Cint}
    KNDOFG::GFortranVectorDescriptorRaw{Cint}
    ITYPEE::GFortranVectorDescriptorRaw{Cint}
    ISTEPA::GFortranVectorDescriptorRaw{Cint}
    ITYPEG::GFortranVectorDescriptorRaw{Cint}
    ISTGPA::GFortranVectorDescriptorRaw{Cint}
    A::GFortranVectorDescriptorRaw{Cdouble}
    B::GFortranVectorDescriptorRaw{Cdouble}
    BL::GFortranVectorDescriptorRaw{Cdouble}
    BU::GFortranVectorDescriptorRaw{Cdouble}
    X::GFortranVectorDescriptorRaw{Cdouble}
    C::GFortranVectorDescriptorRaw{Cdouble}
    Y::GFortranVectorDescriptorRaw{Cdouble}
    GSCALE::GFortranVectorDescriptorRaw{Cdouble}
    ESCALE::GFortranVectorDescriptorRaw{Cdouble}
    VSCALE::GFortranVectorDescriptorRaw{Cdouble}
    EPVALU::GFortranVectorDescriptorRaw{Cdouble}
    GPVALU::GFortranVectorDescriptorRaw{Cdouble}
    GXEQX::GFortranVectorDescriptorRaw{FortranBool}
    INTREP::GFortranVectorDescriptorRaw{FortranBool}
    VNAMES::GFortranVectorDescriptorRaw{FortranString{10}}
    GNAMES::GFortranVectorDescriptorRaw{FortranString{10}}
    _gc::Vector{AbstractVector}

    @doc raw"""
    LANCELOT_problem_type(; IELING, ISTADG, IELVAR, ISTAEV, INTVAR, ISTADH, ICNA, ISTADA
    [, KNDOFG], ITYPEE[, ISTEPA][, ITYPEG][, ISTGPA], A, B, BL, BU[, X][, C][, Y],
    GSCALE, ESCALE, VSCALE,[ EPVALU,][ GPVALU,] GXEQX, INTREP[, VNAMES][, GNAMES])

Creates a new problem for LANCELOT.
# Arguments
- `IELING::AbstractVector{Cint}` is a vector of dimension `ISTADG[end] -1` that holds the indices of the nonlinear elements
  ``\mathcal E_i`` used by each group. The indices for group ``i`` must immediately precede those for group ``i + 1`` and each
  group's indices must appear in a contiguous list.
- `ISTADG::AbstractVector{Cint}` is a vector of dimension ``n_g+1`` whose ``i``-th value gives the position in `IELING` of the
  first nonlinear element in group function ``i``. In addition, `ISTADG[end]` should be equal to the position in `IELING` of
  the last nonlinear element in group ``n_g`` plus one.
- `IELVAR::AbstractVector{Cint}` is a vector of dimension `ISTAEV[end] -1` that holds the indices of the variables in the first
  nonlinear element ``e_1``, followed by those in the second nonlinear element ``e_2``, ....
- `ISTAEV::AbstractVector{Cint}` is a vector of dimension ``n_e +1`` whose ``k``-th value is the position of the first variable
  of the ``k``-th nonlinear element function, in the list `IELVAR`. In addition, `ISTAEV[end]` must be equal to the position of
  the last variable of element ``n_e`` in `IELVAR` plus one.
- `INTVAR::AbstractVector{Cint}` is a vector of dimension ``n_e +1`` whose ``i``-th value must be set to the number of internal
  variables required for the ``i``-th nonlinear element function ``e_i`` on initial entry.
- `ICNA::AbstractVector{Cint}` is a vector of dimension `ISTADA[end] -1` that holds the indices of the nonzero components of
  ``a_1``, the gradient of the first linear element, in any order, followed by those in ``a_2``, etc.
- `ISTADA::AbstractVector{Cint}` is a vector of dimension ``n_g+1`` whose ``i``-th value is the position of the first nonzero
  component of the ``i``-th linear element gradient, ``a_i``, in the list `ICNA`. In addition, `ISTADA[end]` must be equal to
  the position of the last nonzero component of ``a_{n_g}`` in `ICNA` plus one.
- `KNDOFG::AbstractVector{Cint}` is a vector of dimension ``n_g`` that is used to indicate which of the groups are to be
  included in the objective function, which define equality constraints, and which are to be ignored. If KNDOFG is omitted, all
  groups will be included in the objective function, and it will be assumed that there are no general constraints. If `KNDOFG`
  is present, each of the first ``n_g`` entries of `KNDOFG` must be set to `0`, `1` or `2` on initial entry. If `KNDOFG[i]` has
  the value `1`, ``i`` lies in the set ``\mathcal G_O`` and the group will be included in the objective function. If
  `KNDOFG[i]` has the value `2`, ``i`` lies in the complement set ``\mathcal G_C``, and the ``i``-th group defines an equality
  constraint. Finally, if `KNDOFG[i]` has the value `0`, group ``i`` will be ignored; this is useful when, for example, many
  optimizations are required with different subsets of constraints, or when a feasible point is sought without reference to the
  objective.
- `ITYPEE::AbstractVector{Cint}` is a vector dimension ``n_e`` that is used to indicate the types of the nonlinear element
  functions. The ``i``-th component of `ITYPEE` specifies the type of element ``i``.
- `ISTEPA::AbstractVector{Cint}` is a vector of dimension ``n_e +1``, whose ``i``-th component gives the position in `EPVALU`
  of the first parameter for element function ``i``. In addition, `ISTEPA[end]` is the position in `EPVALU` of the last
  parameter for element function ``n_e`` plus one.
- `ITYPEG::AbstractVector{Cint}` is a vector of dimension ``n_g`` that is used to indicate the types of the group functions.
  The ``i``-th component of `ITYPEG` specifies the type of group ``i``.
- `ISTGPA::AbstractVector{Cint}` is a vector of dimension ``n_g +1`` whose ``i``-th component gives the position in `GPVALU` of
  the first parameter for group function ``i``. In addition, `ISTGPA[end]` is the position in `GPVALU` of the last parameter
  for element function ``n_g`` plus one.
- `A::AbstractVector{Cdouble}` is a vector of dimension `ISTADA[end]` that holds the values of the nonzero components of the
  gradients of the linear element functions, ``a_i``, ``i = 1, ..., n_g``. The values must appear in the same order as their
  indices appear in `ICNA`, i.e., the nonzero from element ``i``, whose index is, say, `ICNA[k]` will have value `A[k]`.
- `B::AbstractVector{Cdouble}` is a vector of dimension ``n_g` whose ``i``-th entry holds the value of the constant ``b_i`` for
  each group.
- `BL::AbstractVector{Cdouble}` is a vector of dimension `n` whose ``i``-th entry must be set to the value of the lower bound
  ``l_i`` on the ``i``-th variable. If the ``i``-th variable has no lower bound, `BL[i]` should be set to a large negative
  number.
- `BU::AbstractVector{Cdouble}` is a vector of dimension ``n`` whose ``i``-th entry must be set to the value of the upper bound
  ``u_i`` on the ``i``-th variable. If the ``i``-th variable has no upper bound, `BU[i]` should be set to a large positive
  number.
- `X::AbstractVector{Cdouble}` is a vector of dimension ``n`` that holds the current values of the minimization variables,
  ``\vec x``.
- `C::AbstractVector{Cdouble}` is a vector of dimension ``n_g`` that holds the current estimates of the values of the equality
  constraints for the problem. If `KNDOFG[i] < 2`, `C[i]` will not be set, but if `KNDOFG[i] = 2`, `C[i]` contains the
  constraint value ``c_i(x)`` of the ``i``-th constraint. `C` need not be present if `KNDOFG` is not.
- `Y::AbstractVector{Cdouble}` is a vector of dimension ``n_g`` that holds the current estimates of the Lagrange multipliers,
  ``y``, for the problem. If `KNDOFG[i] < 2`, `Y[i]` will not be set, but if `KNDOFG[i] = 2`, `Y[i]` contains the multiplier
  estimate ``y_i`` for the ``i``-th constraint. `Y` need not be present if `KNDOFG` is not.
- `GSCALE::AbstractVector{Cdouble}` is a vector of dimension ``n_g`` whose ``i``-th entry holds the value of ``i``-th group
  weight ``w^g_i``.
- `ESCALE::AbstractVector{Cdouble}` is a vector of dimension `ISTADG[end] -1` whose entries hold the values of element weights
  ``w^e_{i, j}``. The weights must occur in the same order as the indices of the elements assigned to each group in `IELING`,
  with the weights for the elements in group ``i`` preceding those in group ``i + 1``, ``i = 1, \dots, n_g − 1``.
- `VSCALE::AbstractVector{Cdouble}` is a vector of dimension ``n`` that holds suitable positive scale factors for the problem
  variables ``x``. The ``i``-th variable ``x_i`` will implicitly be divided by `VSCALE[i]` within `LANCELOT_solve`.
  The scale factors should ideally be chosen so that the rescaled variables are of order one at the solution to the
  minimization problem. If the user does not know suitable scalings, each component of `VSCALE` should be set to `1.0`. Good
  variable scalings can result in considerable savings in computing times.
- `EPVALU::AbstractVector{Cdouble}` is a vector of dimension `ISTEPA[end] -1` that holds the values of the element parameters
  ``\vec p^e_i``, ``i = 1, \dots, n_e``. The indices for the parameters for element ``i`` immediately precede those for
  element ``i + 1``, and each element's parameters appear in a contiguous list.
- `GPVALU::AbstractVector{Cdouble}` is a vector of dimension `ISTGPA[end] -1` that holds the values of the group parameters
  ``\vec p^g_i``, ``i = 1, \dots, n_g``. The indices for the parameters for group ``i`` immediately precede those for group
  ``i + 1``, and each group's parameters appear in a contiguous list.
- `GXEQX::AbstractVector{<:Union{Bool,FortranBool}}` is a vector of dimension ``n_g`` whose ``i``-th entry must be set `true`
  if the ``i``-th group function is the trivial function ``g(x) = x`` and `false` otherwise.
- `INTREP::AbstractVector{<:Union{Bool,FortranBool}}` is a vector of dimension ``n_e`` whose ``i``-th entry must be set `true`
  if the ``i``-th nonlinear element function has a useful transformation between elemental and internal variables and `false`
  otherwise.
- `VNAMES::AbstractVector{<:AbstractString}` is a vector of dimension ``n`` whose ``j``-th entry contains the "name" of the
  ``j``-th variable (at most 10 chars).
- `GNAMES::AbstractVector{<:AbstractString}` is a vector of dimension ``n_g`` whose ``i``-th entry contains the "name" of the
  ``i``-th group (at most 10 chars).

# More fields
- `n::Cint` holds the number of optimization variables, ``n``.
- `ng::Cint` holds the number of group functions, ``n_g``.
- `nel::Cint` that holds the number of nonlinear element functions, ``n_e``.
- `ISTADH::AbstractVector{Cint}` is a vector of dimension ``n_e +1`` that is set in [`LANCELOT_solve`](@ref) so that its
  ``i``-th value (``1 \leq i \leq n_e``) gives the position in the array `FUVALS` of the first component of the Hessian matrix
  of the ``i``-th nonlinear element function ``e_i``, with respect to its internal variables. Only the upper triangular part of
  each Hessian matrix is stored and the storage is by columns. The element `ISTADH[end]` gives the position in `FUVALS` of the
  first component of the gradient of the objective function.
    """
    function LANCELOT_problem_type(;
        IELING::AbstractVector{Cint}, ISTADG::AbstractVector{Cint}, IELVAR::AbstractVector{Cint}, ISTAEV::AbstractVector{Cint},
        INTVAR::AbstractVector{Cint}, ICNA::AbstractVector{Cint}, ISTADA::AbstractVector{Cint},
        KNDOFG::Union{<:AbstractVector{Cint},Missing}=missing, ITYPEE::AbstractVector{Cint},
        ISTEPA::Union{<:AbstractVector{Cint},Missing}=missing, ITYPEG::Union{<:AbstractVector{Cint},Missing}=missing,
        ISTGPA::Union{<:AbstractVector{Cint},Missing}=missing, A::AbstractVector{Cdouble},
        B::AbstractVector{Cdouble}, BL::AbstractVector{Cdouble}, BU::AbstractVector{Cdouble}, X::AbstractVector{Cdouble},
        C::Union{<:AbstractVector{Cdouble},Missing}=missing, Y::Union{<:AbstractVector{Cdouble},Missing}=missing,
        GSCALE::AbstractVector{Cdouble}, ESCALE::AbstractVector{Cdouble}, VSCALE::AbstractVector{Cdouble},
        EPVALU::Union{<:AbstractVector{Cdouble},Missing}=missing, GPVALU::Union{<:AbstractVector{Cdouble},Missing}=missing,
        GXEQX::Union{<:AbstractVector{Bool},<:AbstractVector{FortranBool},Missing},
        INTREP::Union{<:AbstractVector{Bool},<:AbstractVector{FortranBool},Missing},
        VNAMES::AbstractVector{<:AbstractString}=["x" * i for i in 1:length(X)],
        GNAMES::AbstractVector{<:AbstractString}=["g" * i for i in 1:length(ISTADG)-1])
        n = length(X)
        ng = length(ISTADG) -1
        nel = length(ISTAEV) -1
        ((ismissing(KNDOFG) || (!ismissing(C) && !ismissing(Y))) &&
            ismissing(ISTEPA) == ismissing(EPVALU) &&
            ismissing(ISTGPA) == ismissing(ITYPEG) == ismissing(GPVALU)) || throw(MethodError(LANCELOT_problem_type, ()))
        (length(IELING) == ISTADG[end] -1 &&
            length(IELVAR) == ISTAEV[end] -1 &&
            length(INTVAR) == nel +1 &&
            length(ISTADA) == ng +1 &&
            length(ICNA) == ISTADA[end] -1 &&
            (ismissing(KNDOFG) || length(KNDOFG) == ng) &&
            length(ITYPEE) == nel &&
            (ismissing(ISTEPA) || length(ISTEPA) == nel +1) &&
            (ismissing(ITYPEG) || length(ITYPEG) == ng) &&
            (ismissing(ISTGPA) || length(ISTGPA) == ng +1) &&
            length(A) == ISTADA[end] -1 &&
            length(B) == ng &&
            length(BL) == n &&
            length(BU) == n &&
            (ismissing(C) || length(C) == ng) &&
            (ismissing(Y) || length(Y) == ng) &&
            length(GSCALE) == ng &&
            length(ESCALE) == ISTADG[end] -1 &&
            length(VSCALE) == n &&
            (ismissing(EPVALU) || length(EPVALU) == ISTEPA[end] -1) &&
            (ismissing(GPVALU) || length(GPVALU) == ISTGPA[end] -1) &&
            length(GXEQX) == ng &&
            length(INTREP) == nel &&
            length(VNAMES) == n &&
            length(GNAMES) == ng) || throw(ArgumentError("The input parameter sizes were not suitable"))
        _gc = AbstractVector[]
        new(
            n, ng, nel,
            GFortranVectorDescriptorRaw{Cint}(pushnm!(_gc, IELING)),
            GFortranVectorDescriptorRaw{Cint}(pushnm!(_gc, ISTADG)),
            GFortranVectorDescriptorRaw{Cint}(pushnm!(_gc, IELVAR)),
            GFortranVectorDescriptorRaw{Cint}(pushnm!(_gc, ISTAEV)),
            GFortranVectorDescriptorRaw{Cint}(pushnm!(_gc, INTVAR)),
            GFortranVectorDescriptorRaw{Cint}(pushnm!(_gc, Vector{Cint}(undef, nel +1))),
            GFortranVectorDescriptorRaw{Cint}(pushnm!(_gc, ICNA)),
            GFortranVectorDescriptorRaw{Cint}(pushnm!(_gc, ISTADA)),
            GFortranVectorDescriptorRaw{Cint}(pushnm!(_gc, KNDOFG)),
            GFortranVectorDescriptorRaw{Cint}(pushnm!(_gc, ITYPEE)),
            GFortranVectorDescriptorRaw{Cint}(pushnm!(_gc, ISTEPA)),
            GFortranVectorDescriptorRaw{Cint}(pushnm!(_gc, ITYPEG)),
            GFortranVectorDescriptorRaw{Cint}(pushnm!(_gc, ISTGPA)),
            GFortranVectorDescriptorRaw{Cdouble}(pushnm!(_gc, A)),
            GFortranVectorDescriptorRaw{Cdouble}(pushnm!(_gc, B)),
            GFortranVectorDescriptorRaw{Cdouble}(pushnm!(_gc, BL)),
            GFortranVectorDescriptorRaw{Cdouble}(pushnm!(_gc, BU)),
            GFortranVectorDescriptorRaw{Cdouble}(pushnm!(_gc, X)),
            GFortranVectorDescriptorRaw{Cdouble}(pushnm!(_gc, C)),
            GFortranVectorDescriptorRaw{Cdouble}(pushnm!(_gc, Y)),
            GFortranVectorDescriptorRaw{Cdouble}(pushnm!(_gc, GSCALE)),
            GFortranVectorDescriptorRaw{Cdouble}(pushnm!(_gc, ESCALE)),
            GFortranVectorDescriptorRaw{Cdouble}(pushnm!(_gc, VSCALE)),
            GFortranVectorDescriptorRaw{Cdouble}(pushnm!(_gc, EPVALU)),
            GFortranVectorDescriptorRaw{Cdouble}(pushnm!(_gc, GPVALU)),
            GFortranVectorDescriptorRaw{FortranBool}(pushnm!(_gc, eltype(GXEQX) === FortranBool ? GXEQX :
                                                                                convert(AbstractVector{FortranBool}, GXEQX))),
            GFortranVectorDescriptorRaw{FortranBool}(pushnm!(_gc, eltype(INTREP) === FortranBool ? INTREP :
                                                                                convert(AbstractVector{FortranBool}, INTREP))),
            GFortranVectorDescriptorRaw{FortranString{10}}(pushnm!(_gc, FortranString{10}.(VNAMES))),
            GFortranVectorDescriptorRaw{FortranString{10}}(pushnm!(_gc, FortranString{10}.(GNAMES))),
            _gc
        )
    end
end

struct EXTEND_save_type
    lirnh::Cint
    ljcnh::Cint
    llink_min::Cint
    lirnh_min::Cint
    ljcnh_min::Cint
    lh_min::Cint
    lh::Cint
    litran_min::Cint
    lwtran_min::Cint
    lwtran::Cint
    litran::Cint
    l_link_e_u_v::Cint
    llink::Cint
    lrowst::Cint
    lpos::Cint
    lused::Cint
    lfilled::Cint

    EXTEND_save_type() = new()
end

struct CAUCHY_save_type
    iterca::Cint
    iter::Cint
    itmax::Cint
    nfreed::Cint
    nbreak::Cint
    nzero::Cint
    tk::Cdouble
    gxt::Cdouble
    hxt::Cdouble
    epstl2::Cdouble
    tpttp::Cdouble
    tcauch::Cdouble
    tbreak::Cdouble
    deltat::Cdouble
    epsqrt::Cdouble
    gxtold::Cdouble
    g0tp::Cdouble
    t::Cdouble
    tamax::Cdouble
    ptp::Cdouble
    gtp::Cdouble
    flxt::Cdouble
    tnew::Cdouble
    prnter::FortranBool
    pronel::FortranBool
    recomp::FortranBool

    CAUCHY_save_type() = new()
end

struct CG_save_type
    iter::Cint
    itsle::Cint
    alpha::Cdouble
    oldgns::Cdouble
    onepep::Cdouble
    prnter::FortranBool
    pronel::FortranBool

    CG_save_type() = new()
end

struct ASMBL_save_type
    ptr_status::FortranBool
    ICNTL::NTuple{30,Cint}
    INFO::NTuple{20,Cint}
    CNTL::NTuple{5,Cdouble}

    ASMBL_save_type() = new()
end

struct PRECN_save_type
    liw::Cint
    lw::Cint
    nsemiw::Cint
    nupdat::Cint
    liccgg::Cint
    nextra::Cint
    nz01::Cint
    iaj::Cint
    tfactr::Cfloat
    t1stsl::Cfloat
    tupdat::Cfloat
    tsolve::Cfloat
    ICNTL_iccg::NTuple{5,Cint}
    KEEP_iccg::NTuple{12,Cint}
    INFO_iccg::NTuple{10,Cint}
    CNTL_iccg::NTuple{3,Cdouble}

    PRECN_save_type() = new()
end

struct OTHERS_fdgrad_save_type
    backwd::FortranBool

    OTHERS_fdgrad_save_type() = new()
end

@kwdef struct SCU_matrix_type
    n::Cint = 0 # uninitialized
    m::Cint = 0 # uninitialized
    m_max::Cint = 0 # uninitialized
    class::Cint = 0
    BD_row::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    BD_col_start::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    CD_col::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    CD_row_start::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    BD_val::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    CD_val::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
end

@kwdef struct SCU_data_type
    m::Cint = 0 # uninitialized
    m_max::Cint = 0 # uninitialized
    jumpto::Cint = 0 # uninitialized
    jcol::Cint = 0 # uninitialized
    newdia::Cint = 0 # uninitialized
    sign_determinant::Cint = 0 # uninitialized
    class::Cint = 3
    got_factors::FortranBool = false # uninitialized
    dianew::Cdouble = 0. # uninitialized
    R::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    W::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    Q::GFortranMatrixDescriptorRaw{Cdouble} = GFortranMatrixDescriptorRaw{Cdouble}()
end

@kwdef struct SMT_type
    m::Cint = 0 # uninitialized
    n::Cint = 0 # uninitialized
    ne::Cint = 0 # uninitialized
    id::GFortranVectorDescriptorRaw{FortranChar} = GFortranVectorDescriptorRaw{FortranChar}()
    type::GFortranVectorDescriptorRaw{FortranChar} = GFortranVectorDescriptorRaw{FortranChar}()
    row::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    col::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    ptr::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    val::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
end

@kwdef struct SILS_factors
    keep::GFortranMatrixDescriptorRaw{Cint} = GFortranMatrixDescriptorRaw{Cint}()
    # keep::GFortranVectorDescriptorRaw{Cint}
    iw::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    iw1::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    iw2::GFortranMatrixDescriptorRaw{Cint} = GFortranMatrixDescriptorRaw{Cint}()
    val::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    w::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    r::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    n::Cint = -1
    nrltot::Cint = -1
    nirtot::Cint = -1
    nrlnec::Cint = -1
    nirnec::Cint = -1
    nsteps::Cint = -1
    maxfrt::Cint = -1
    latop::Cint = -1
    dim_iw1::Cint = -1
    pivoting::Cint = -1
    ops::Cdouble = -1
end

struct LANCELOT_save_type
    full_solution::FortranBool
    igetfd::Cint
    unsucc::FortranBool
    nobjgr::Cint
    m::Cint
    icrit::Cint
    ncrit::Cint
    p_type::Cint
    ocnorm::Cdouble
    cnorm_major::Cdouble
    etak::Cdouble
    eta0::Cdouble
    omegak::Cdouble
    omega0::Cdouble
    tau::Cdouble
    tau_steering::Cdouble
    gamma1::Cdouble
    alphae::Cdouble
    betae::Cdouble
    alphak::Cdouble
    alphao::Cdouble
    betao::Cdouble
    omega_min::Cdouble
    eta_min::Cdouble
    epstol::Cdouble
    epsgrd::Cdouble
    cnorm::Cdouble
    STATE::NTuple{5,FortranString{5}}
    itzero::FortranBool
    reeval::FortranBool
    ifactr::Cint
    ldx::Cint
    lfxi::Cint
    lgxi::Cint
    lhxi::Cint
    lggfx::Cint
    nvar2::Cint
    nfreef::Cint
    nfree::Cint
    nnonnz::Cint
    nadd::Cint
    icfact::Cint
    jumpto::Cint
    nbprod::Cint
    infor::Cint
    number::Cint
    nfixed::Cint
    ibqpst::Cint
    nmhist::Cint
    maxsel::Cint
    ntotin::Cint
    nfreec::Cint
    lnguvl::Cint
    lnhuvl::Cint
    ntype::Cint
    nsets::Cint
    nvargp::Cint
    l_suc::Cint
    msweep::Cint
    nbnd::Cint
    mortor_its::Cint
    ntotel::Cint
    inform_status::Cint
    nvrels::Cint
    nnza::Cint
    error::Cint
    out::Cint
    print_level::Cint
    start_print::Cint
    stop_print::Cint
    print_gap::Cint
    n_steering::Cint
    n_steering_this_iteration::Cint
    first_derivatives::Cint
    second_derivatives::Cint
    epstlp::Cdouble
    gmodel::Cdouble
    vscmax::Cdouble
    rad::Cdouble
    maximum_radius::Cdouble
    epsrcg::Cdouble
    fnew::Cdouble
    radmin::Cdouble
    cgstop::Cdouble
    diamin::Cdouble
    diamax::Cdouble
    ared::Cdouble
    prered::Cdouble
    rho::Cdouble
    fmodel::Cdouble
    curv::Cdouble
    dxsqr::Cdouble
    fcp::Cdouble
    f0::Cdouble
    stepmx::Cdouble
    smallh::Cdouble
    resmin::Cdouble
    qgnorm::Cdouble
    oldrad::Cdouble
    epscns::Cdouble
    radtol::Cdouble
    fill::Cdouble
    step::Cdouble
    teneps::Cdouble
    stpmin::Cdouble
    epstln::Cdouble
    f_min::Cdouble
    f_r::Cdouble
    f_c::Cdouble
    sigma_r::Cdouble
    sigma_c::Cdouble
    findmx::Cdouble
    f_min_lag::Cdouble
    f_r_lag::Cdouble
    f_c_lag::Cdouble
    f_min_viol::Cdouble
    f_r_viol::Cdouble
    f_c_viol::Cdouble
    violation::Cdouble
    delta_qv::Cdouble
    delta_qv_steering::Cdouble
    alllin::FortranBool
    altriv::FortranBool
    next::FortranBool
    second::FortranBool
    print_header::FortranBool
    modchl::FortranBool
    iprcnd::FortranBool
    munks::FortranBool
    seprec::FortranBool
    densep::FortranBool
    calcdi::FortranBool
    xactcp::FortranBool
    reusec::FortranBool
    gmpspr::FortranBool
    slvbqp::FortranBool
    refact::FortranBool
    fdgrad::FortranBool
    centrl::FortranBool
    dprcnd::FortranBool
    strctr::FortranBool
    use_band::FortranBool
    icfs::FortranBool
    mortor::FortranBool
    firsup::FortranBool
    twonrm::FortranBool
    direct::FortranBool
    myprec::FortranBool
    prcond::FortranBool
    firstc::FortranBool
    nobnds::FortranBool
    getders::FortranBool
    save_c::FortranBool
    printt::FortranBool
    printi::FortranBool
    printm::FortranBool
    printw::FortranBool
    printd::FortranBool
    printe::FortranBool
    set_printe::FortranBool
    set_printt::FortranBool
    set_printi::FortranBool
    set_printm::FortranBool
    set_printw::FortranBool
    set_printd::FortranBool
    skipg::FortranBool
    steering::FortranBool
    new_major::FortranBool
    cgend::FortranString{6}
    lisend::FortranString{6}
    cgend1::FortranChar
    lisend1::FortranChar
    t::Cfloat
    time::Cfloat
    tmv::Cfloat
    tca::Cfloat
    tls::Cfloat
    tup::Cfloat
    ISYS::NTuple{5,Cint}
    CGENDS::NTuple{6,FortranString{6}}
    LSENDS::NTuple{5,FortranString{6}}
    CGENDS1::NTuple{6,FortranChar}
    LSENDS1::NTuple{5,FortranChar}
    CAUCHY::CAUCHY_save_type
    CG::CG_save_type
    ASMBL::ASMBL_save_type
    PRECN::PRECN_save_type
    OTHERS::OTHERS_fdgrad_save_type
    EXTEND::EXTEND_save_type

    LANCELOT_save_type() = new()
end

@kwdef mutable struct LANCELOT_data_type
    S::LANCELOT_save_type = LANCELOT_save_type()
    ITRANS::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    ROW_start::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    POS_in_H::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    USED::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    FILLED::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    LINK_elem_uses_var::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    WTRANS::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    ISYMMD::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    ISWKSP::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    ISTAJC::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    ISTAGV::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    ISVGRP::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    ISLGRP::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    IGCOLJ::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    IVALJR::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    IUSED::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    ITYPER::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    ISSWTR::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    ISSITR::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    ISET::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    ISVSET::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    INVSET::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    IFREE::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    INDEX::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    IFREEC::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    INNONZ::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    LIST_elements::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    ISYMMH::GFortranMatrixDescriptorRaw{Cint} = GFortranMatrixDescriptorRaw{Cint}()
    FUVALS_temp::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    P::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    X0::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    XCP::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    GX0::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    RADII::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    DELTAX::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    QGRAD::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    GRJAC::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    CDASH::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    C2DASH::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    GV_old::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    BND::GFortranMatrixDescriptorRaw{Cdouble} = GFortranMatrixDescriptorRaw{Cdouble}()
    BND_radius::GFortranMatrixDescriptorRaw{Cdouble} = GFortranMatrixDescriptorRaw{Cdouble}()
    IW_asmbl::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    NZ_comp_w::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    W_ws::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    W_el::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    W_in::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    H_el::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    H_in::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    IKEEP::GFortranMatrixDescriptorRaw{Cint} = GFortranMatrixDescriptorRaw{Cint}()
    IW1::GFortranMatrixDescriptorRaw{Cint} = GFortranMatrixDescriptorRaw{Cint}()
    IW::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    IVUSE::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    H_col_ptr::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    L_col_ptr::GFortranVectorDescriptorRaw{Cint} = GFortranVectorDescriptorRaw{Cint}()
    W::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    RHS::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    RHS2::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    P2::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    G::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    DIAG::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    BREAKP::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    GRAD::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    W1::GFortranMatrixDescriptorRaw{Cdouble} = GFortranMatrixDescriptorRaw{Cdouble}()
    OFFDIA::GFortranMatrixDescriptorRaw{Cdouble} = GFortranMatrixDescriptorRaw{Cdouble}()
    GROUP_SCALING::GFortranVectorDescriptorRaw{Cdouble} = GFortranVectorDescriptorRaw{Cdouble}()
    GXEQX_AUG::GFortranVectorDescriptorRaw{FortranBool} = GFortranVectorDescriptorRaw{FortranBool}()
    SCU_matrix::SCU_matrix_type = SCU_matrix_type()
    SCU_data::SCU_data_type = SCU_data_type()
    matrix::SMT_type = SMT_type()
    SILS_data::SILS_factors = SILS_factors()
end

@kwdef struct SILS_control
    ICNTL::NTuple{30,Cint} = (6, 6, 0, 2139062143, 1, 32639, 32639, 32639, 32639, 14, 9, 8, 8, 9, 10, 32639, 32639, 32639,
    32689, 24, 11, 9, 8, 9, 10, 0, 0, 0, 0, 0)
    lp::Cint = 6
    wp::Cint = 6
    mp::Cint = 6
    sp::Cint = -1
    ldiag::Cint = 0
    la::Cint = 0
    liw::Cint = 0
    maxla::Cint = typemax(Cint)
    maxliw::Cint = typemax(Cint)
    pivoting::Cint = 1
    nemin::Cint = 1
    factorblocking::Cint = 16
    solveblocking::Cint = 16
    thresh::Cint = 50
    ordering::Cint = 3
    scaling::Cint = 0
    CNTL::NTuple{5,Cdouble} = (0.1, 1.0, 0.0, 0.0, 0.0)
    multiplier::Cdouble = 2.0
    reduce::Cdouble = 2.0
    u::Cdouble = 0.1
    static_tolerance::Cdouble = 0.0
    static_level::Cdouble = 0.0
    tolerance::Cdouble = 0.0
    convergence::Cdouble = 0.5
end

@kwdef mutable struct LANCELOT_control_type
    # This is also nasty, but less so. And we have lots and lots of default values. So let's actually define this type instead
    # of extracting offsets in the assembly.
    error::Cint = 6
    out::Cint = 6
    alive_unit::Cint = 60
    alive_file::NTuple{30,Cchar} = ('A', 'L', 'I', 'V', 'E', '.', 'd', ' ', ' ', ' ',
                                    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
                                    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ')
    print_level::Cint = 0
    maxit::Cint = 1000
    start_print::Cint = -1
    stop_print::Cint = -1
    print_gap::Cint = 1
    linear_solver::Cint = 8
    icfact::Cint = 5
    semibandwidth::Cint = 5
    max_sc::Cint = 100
    io_buffer::Cint = 75
    more_toraldo::Cint = 0
    non_monotone::Cint = 1
    first_derivatives::Cint = 0
    second_derivatives::Cint = 0
    stopc::Cdouble = 1e-5
    stopg::Cdouble = 1e-5
    min_aug::Cdouble = floatmax(Cdouble) / 8
    acccg::Cdouble = 0.01
    initial_radius::Cdouble = -1.
    maximum_radius::Cdouble = 1e20
    eta_successful = 0.01
    eta_very_successful::Cdouble = 0.9
    eta_extremely_successful::Cdouble = 0.95
    gamma_smallest::Cdouble = 0.0625
    gamma_decrease::Cdouble = 0.25
    gamma_increase::Cdouble = 2.0
    mu_meaningful_model::Cdouble = 0.01
    mu_meaningful_group::Cdouble = 0.1
    initial_mu::Cdouble = 0.1
    mu_decrease::Cdouble = 0.1
    mu_steering_decrease::Cdouble = 0.7
    mu_tol::Cdouble = 0.1
    firstg::Cdouble = 0.1
    firstc::Cdouble = 0.1
    num_mudec::Cint = typemax(Cint)
    num_mudec_per_iteration::Cint = typemax(Cint)
    kappa_3::Cdouble = 1e-5
    kappa_t::Cdouble = 0.9
    mu_min::Cdouble = 0.
    cpu_time_limit::Cdouble = -1.
    quadratic_problem::FortranBool = false
    steering::FortranBool = false
    two_norm_tr::FortranBool = false
    exact_gcp::FortranBool = true
    gn_model::FortranBool = false
    gn_model_after_cauchy::FortranBool = false
    magical_steps::FortranBool = false
    accurate_bqp::FortranBool = false
    structured_tr::FortranBool = false
    print_max::FortranBool = false
    full_solution::FortranBool = true
    space_critical::FortranBool = false
    deallocate_error_fatal::FortranBool = false
    SILS_cntl::SILS_control = SILS_control()
end

@kwdef struct SCU_inform_type
    status::Cint = 0
    alloc_status::Cint = 0
    intertia::NTuple{3,Cint} = (0, 0, 0)
end

@kwdef struct SILS_ainfo
    flag::Cint = 0
    more::Cint = 0
    nsteps::Cint = 0
    nrltot::Cint = -1
    nirtot::Cint = -1
    nrlnec::Cint = -1
    nirnec::Cint = -1
    nrladu::Cint = -1
    niradu::Cint = -1
    ncmpa::Cint = 0
    oor::Cint = 0
    dup::Cint = 0
    maxfrt::Cint = -1
    stat::Cint = 0
    faulty::Cint = 0
    opsa::Cdouble = -1.0
    opse::Cdouble = -1.0
end

@kwdef struct SILS_finfo
    flag::Cint = 0
    more::Cint = 0
    maxfrt::Cint = -1
    nebdu::Cint = -1
    nrlbdu::Cint = -1
    nirbdu::Cint = -1
    nrltot::Cint = -1
    nirtot::Cint = -1
    nrlnec::Cint = -1
    nirnec::Cint = -1
    ncmpbr::Cint = -1
    ncmpbi::Cint = -1
    ntwo::Cint = -1
    neig::Cint = -1
    delay::Cint = -1
    signc::Cint = -1
    static::Cint = -1
    modstep::Cint = -1
    rank::Cint = -1
    stat::Cint = -1
    faulty::Cint = -1
    step::Cint = -1
    opsa::Cdouble = -1.0
    opse::Cdouble = -1.0
    opsb::Cdouble = -1.0
    maxchange::Cdouble = -1.0
    smin::Cdouble = 0.0
    smax::Cdouble = 0.0
end

@kwdef struct SILS_sinfo
    flag::Cint = -1
    stat::Cint = -1
    cond::Cdouble = -1.0
    cond2::Cdouble = -1.0
    berr::Cdouble = -1.0
    berr2::Cdouble = -1.0
    error::Cdouble = -1.0
end

@kwdef mutable struct LANCELOT_inform_type
    status::Cint = 0
    alloc_status::Cint = 0
    iter::Cint = -1
    itercg::Cint = -1
    itcgmx::Cint = -1
    ncalcf::Cint = 0
    ncalcg::Cint = 0
    nvar::Cint = 0
    ngeval::Cint = 0
    iskip::Cint = 0
    ifixed::Cint = 0
    nsemib::Cint = 0
    aug::Cdouble = floatmax(Cdouble)
    obj::Cdouble = floatmax(Cdouble)
    pjgnrm::Cdouble = floatmax(Cdouble)
    cnorm::Cdouble = 0.
    ratio::Cdouble = 0.
    mu::Cdouble = 0.
    radius::Cdouble = 0.
    ciccg::Cdouble = 0.
    newsol::FortranBool = false
    bad_alloc::NTuple{80,Cchar} = ntuple(_ -> ' ', 80)
    SCU_info::SCU_inform_type = SCU_inform_type()
    SILS_infoa::SILS_ainfo = SILS_ainfo()
    SILS_infof::SILS_finfo = SILS_finfo()
    SILS_infos::SILS_sinfo = SILS_sinfo()
end