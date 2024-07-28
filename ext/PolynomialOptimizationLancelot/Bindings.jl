export LANCELOT_simple

include("./FortranSpecials.jl")

_desc(::Nothing) = C_NULL
_desc(x::AbstractArray) = Ref(convert(GFortranArrayDescriptor, x))

#region LANCELOT_simple
@doc raw"""
    LANCELOT_simple(n, X, MY_FUN; MY_GRAD=missing, MY_HESS=missing, BL=nothing, BU=nothing,
        neq=0, nin=0, CX=nothing, Y=nothing, maxit=1000, gradtol=1e-5, feastol=1e-5,
        print_level=1)

# Purpose
A simple and somewhat NAIVE interface to LANCELOT B for solving the nonlinear optimization problem
```math
\min_x f(x)
```
possibly subject to constraints of the one or more of the forms
```math
\begin{aligned}
   b_{\mathrm l} & \leq x \leq b_{\mathrm u}, \\
   c_{\mathrm e}( x ) & = 0, \\
   c_{\mathrm i}( x ) & \leq 0
\end{aligned}
```

where ``f\colon \mathbb R^n \to \mathbb R``, ``c_{\mathrm e}: \mathbb R^n \to \mathbb R^{n_{\mathrm{eq}}}`` and
``c_{\mathrm i}\colon \mathbb R^n \to \mathbb R^{n_{\mathrm{in}}}`` are twice-continuously differentiable functions.

# Why naive?
At variance with more elaborate interfaces for LANCELOT, the present one completely *ignores underlying partial separability or
sparsity structure, restricts the possible forms under which the problem may be presented to the solver, and drastically
limits the range of available algorithmic options*. If simpler to use than its more elaborate counterparts, it therefore
provides a possibly substantially inferior numerical performance, especially for difficult/large problems, where structure
exploitation and/or careful selection of algorithmic variants matter.

!!! warning
    The best performance obtainable with LANCELOT B is probably not with the present interface.

# How to use it?
## Unconstrained problems
The user should provide, at the very minimum, suitable values for the following input arguments:

- `n::Integer`: the number of variables,
- `X::AbstractVector{Float64}` (strided vector of size `n`): the starting point for the minimization
- `MY_FUN::Callable`: a function for computing the objective function value for any `X`, whose interface has the default form
  `MY_FUN(X::AbstractVector{Float64})::Float64` where `X[1:n]` contains the values of the variables on input, and which returns
  a double precision scalar representing the value ``f(X)``.
- If the gradient of ``f`` can be computed, then the (optional) keyword argument `MY_GRAD` must be specified and given a
  function computing the gradient, whose interface must be of the form
  `MY_GRAD(G::AbstractVector{Float64}, X::AbstractVector{Float64})`, where `G` is a double precision vector of size `n` in
  which the function returns the value of the gradient of ``f`` at `X`.
- If, additionally, the second-derivative matrix of ``f`` at `X` can be computed, the (optional) keyword argument `MY_HESS`
  must be specified and given a function computing the Hessian, whose interface must be of the form
  `MY_HESS(H::SPMatrix{Float64}, X::AbstractVector{Float64})`, where `H` is a double precision symmetrix matrix in packed
  storage format (upper triangular by column, see the [`StandardPacked.jl`](https://github.com/projekter/StandardPacked.jl)
  package) of the Hessian of ``f`` at `X`.

In all cases, the best value of ``x`` found by LANCELOT B is returned to the user in the vector `X` and the associated
objective function value is the first return value.
The second return value reports the number of iterations performed by LANCELOT before exiting.
Finally, the last return value contains the exit status of the LANCELOT run, the value `0` indicating a successful run. Other
values indicate errors in the input or unsuccessful runs, and are detailed in the specsheet of LANCELOT B (with the exception
of the value 19, which reports a negative value for one or both input arguments `nin` and `neq`).

### Example
Let us consider the optimization problem
```math
\min_{x_1, x_2} f(x_1, x_2) = 100 ( x_2 - x_1^2 )^2 + ( 1 - x_1 )^2
```
which is the ever-famous Rosenbrock "banana" problem.
The most basic way to solve the problem (but NOT the most efficient) is, assuming the starting point `X = [-1.2, 1.]` known, to
perform the call `PolynomialOptimization.LANCELOT_simple(2, X, FUN)` where the user-provided function `FUN` is given by
```julia
FUN(X) = @inbounds 100 * (X[2] - X[1]^2)^2 + (1 - X[1])^2
```

The solution is returned in 60 iterations with exit code `0`.

If we now wish to use first and second derivatives of the objective function, one should use the call
```julia
PolynomialOptimization.LANCELOT_simple(2, X, FUN, MY_GRAD=GRAD!, MY_HESS=HESS!)
```
and provide the additional routines
```julia
GRAD!(G, X) = @inbounds begin
    G[1] = -400 * (X[2] - X[1]^2) * X[1] - 2 * (1 - X[1])
    G[2] = 200 * (X[2] - X[1]^2)
end

HESS!(H, X) = @inbounds begin
    H[1, 1] = -400 * (X[2] - 3 * X[1]^2) + 2
    H[1, 2] = -400 * X[1]
    H[2, 2] = 200
end
```

Convergence is then obtained in 23 iterations. Note that using exact first-derivatives only is also possible: `MY_HESS` should
then be absent from the calling sequence and providing the subroutine `HESS!` unnecessary.

## Bound constrained problems
Bound on the problem variables may be imposed by specifying one or both of
- `BL::AbstractVector{Float64}` (double precision vector of size `n`): the lower bounds on `X`,
- `BU::AbstractVector{Float64}` (double precision vector of size `n`): the upper bounds on `X`.
Note that infinite bounds (represented by a number larger than `1e20` in absolute value) are acceptable, as well as equal
lower and upper bounds, which amounts to fixing the corresponding variables. Except for the specification of `BL` and/or `BU`,
the interface is identical to that for unconstrained problems.

### Example
If one now wishes to impose zero upper bounds on the variables of our unconstrained problem, one could use the following call
```julia
PolynomialOptimization.LANCELOT_simple(2, X, FUN, MY_GRAD=GRAD!, MY_HESS=HESS!,
    BU=zeros(2))
```
in which case convergence is obtained in 6 iterations.

## Equality constrained problems
If, additionally, general equality constraints are also present in the problem, this must be declared by specifying the
following (optional) input argument:
- `neq::Integer`: the number of equality constraints.
In this case, the equality constraints are numbered from 1 to `neq` and the value of the `i`-th equality constraint must be
computed by a user-supplied routine of the form `FUN(X, i)` (with `i = 1, ..., neq`) where the function now returns the value
of the `i`-th equality constraint evaluated at `X` if `i` is specified. (This extension of the unconstrained case can be
implemented by adding an optional argument `i` to the unconstrained version of `FUN` or by defining a three-parameter method on
its own.)
If derivatives are available, then the `MY_GRAD` and `MY_HESS` subroutines must be adapted as well: `MY_GRAD(G, X, i)` and
`MY_HESS(H, X, i)` for computing the gradient and Hessian of the `i`-th constraint at `X`. Note that, if the gradient of the
objective function is available, so must be the gradients of the equality constraints. The same level of derivative
availability is assumed for all problem functions (objective and constraints). The final values of the constraints and the
values of their associated Lagrange multipliers is optionally returned to the user in the (optional) double precision keyword
arguments `CX` and `Y`, respectively (both being of size `neq`).

## Inequality constrained problems
If inequality constraints are present in the problem, their inclusion is similar to that of equality constraints. One then
needs to specify the (optional) input argument
- `nin::Integer`: the number of inequality constraints.
The inequality constraints are then numbered from `neq+1` to `neq+nin` and their values or that of their derivatives is again
computed by calling, for `i = 1, ..., nin`, `FUN(X, i)`, `MY_GRAD(G, X, i)`, `MY_HESS(H, X, i)`.
The inequality constraints are internally converted in equality ones by the addition of a slack variables, whose names are set
to 'Slack_`i`', where the character `i` in this string takes the integers values `1` to `nin`.
The values of the inequality constraints at the final `X` are finally returned (as for equalities) in the optional double
precision keyword argument `CX` of size `nin`. The values of the Lagrange multipliers are returned in the optional double
precision output argument `Y` of size `nin`.

## Problems with equality and inequality constraints
If they are both equalities and inequalities, `neq` and `nin` must be specified and the values and derivatives of the
constraints are computed by `FUN(X, i)`, `GRAD(G, X, i)`, `HESS(H, X, i)` (`i = 1, ..., neq`) for the equality constraints, and
`FUN(X, i)`, `GRAD(G, X, i)`, `HESS(H, X, i)` (`i = neq+1, ..., neq+nin`) for the inequality constraints. Again, the same level
of derivative availability is assumed for all problem functions (objective and constraints). Finally, the optional arguments
`CX` and/or `Y`, if used, are then of size `neq+nin`.

### Example
If we now wish the add to the unconstrained version the new constraints
```math
\begin{aligned}
    0 & \leq x_1 \\
    x_1 + 3x_2 - 3 & = 0 \\
    x_1^2 + x_2^2 - 4 & \leq 0,
\end{aligned}
```
we may transform our call to
```julia
CX = Vector{Float64}(undef, 2)
Y = Vector{Float64}(undef, 2)
LANCELOT_simple(2, X, FUN; MY_GRAD=GRAD!, MY_HESS=HESS!, BL=[0., -1e20], neq=1, nin=1,
    CX, Y)
```
(assuming we need `CX` and `Y`), and add methods for `FUN`, `GRAD!` and `HESS!` as follows
```julia
FUN(X, i) = @inbounds begin
    if i == 1 # the equality constraint
        return X[1] + 3X[2] - 3
    elseif i == 2 # the inequality constraint
        return X[1]^2 + X[2]^2 - 4
    end
    return NaN # should never happen
end

GRAD!(G, X, i) = @inbounds begin
    if i == 1 # equality constraint's gradient components
        G[1] = 1
        G[2] = 3
    elseif i == 2 # inequality constraint's gradient components
        G[1] = 2X[1]
        G[2] = 2X[2]
    end
    return
end

HESS!(H, X, i) = @inbounds begin
    if i == 1 # equality constraint's Hessian
        fill!(H, 1.)
    elseif i == 2 # inequality constraint's Hessian
        H[1] = 2
        H[2] = 0
        H[3] = 2
    end
    return
end
```

Convergence is then obtained in 8 iterations. Note that, in our example, the objective function or its derivatives is/are
computed if the index `i` is omitted (see above).
Of course, the above examples can easily be modified to represent new minimization problems :-).

# Available algorithmic options
Beyond the choice of derivative level for the problem functions, the following arguments allow a (very limited) control of the
algorithmic choices used in LANCELOT.
- `maxit::Integer`: maximum number of iterations (default: `1000`)
- `gradtol::Real`: the threshold on the infinity norm of the gradient (or of the lagrangian's gradient) for declaring
  convergence  (default: `1.0e-5`)
- `feastol::Real`: the threshold on the infinity norm of the constraint violation for declaring convergence (for constrained
  problems) (default: `1.0e-5`)
- `print_level::Integer`: a positive number proportional to the amount of output by the package: `0` corresponds to the silent
  mode, `1` to a single line of information per iteration (default), while higher values progressively produce more output.

# Other sources
The user is encouraged to consult the specsheet of the (non-naive) interface to LANCELOT within the GALAHAD software library
for a better view of all possibilities offered by an intelligent use of the package. The library is described in the paper
```
N. I. M. Gould, D. Orban, Ph. L. Toint,
GALAHAD, a library of thread-sage Fortran 90 packages for large-scale
nonlinear optimization,
Transactions of the AMS on Mathematical Software, vol 29(4),
pp. 353-372, 2003
```

The book
```
A. R. Conn, N. I. M. Gould, Ph. L. Toint,
LANCELOT, A Fortan Package for Large-Scale Nonlinear Optimization
(Release A),
Springer Verlag, Heidelberg, 1992
```
is also a good source of additional information.

Main author: Ph. Toint, November 2007.
Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
"""
function LANCELOT_simple(n::Integer, X::AbstractVector{Cdouble}, MY_FUN::Base.Callable;
    MY_GRAD::Union{Base.Callable,Missing}=missing, MY_HESS::Union{Base.Callable,Missing}=missing,
    BL::Union{<:AbstractVector{Cdouble},Nothing}=nothing, BU::Union{<:AbstractVector{Cdouble},Nothing}=nothing, neq::Integer=0,
    nin::Integer=0, CX::Union{<:AbstractVector{Cdouble},Nothing}=nothing, Y::Union{<:AbstractVector{Cdouble},Nothing}=nothing,
    maxit::Integer=1000, gradtol::Real=1e-5, feastol::Real=1e-5, print_level::Integer=1)
    MY_FUN_F = @cfunction($((X, fx, i) -> begin
        unsafe_store!(fx, i == C_NULL ? MY_FUN(X[]) : MY_FUN(X[], Int(unsafe_load(i))))
        return
    end), Cvoid, (Ref{GFortranVectorDescriptorRaw{Cdouble}}, Ptr{Cdouble}, Ptr{Cint}))
    fx = Ref{Cdouble}(NaN)
    iters = Ref{Cint}(-1)
    exit_code = Ref{Cint}(-1)
    if !ismissing(MY_GRAD)
        MY_GRAD_F = @cfunction($((X, G, i) -> begin
            if i == C_NULL
                MY_GRAD(G[], X[])
            else
                MY_GRAD(G[], X[], Int(unsafe_load(i)))
            end
            return
        end), Cvoid, (Ref{GFortranVectorDescriptorRaw{Cdouble}}, Ref{GFortranVectorDescriptorRaw{Cdouble}}, Ptr{Cint}))
    end
    if !ismissing(MY_HESS)
        MY_HESS_F = let n=n
            @cfunction($((X, H, i) -> begin
                if i == C_NULL
                    MY_HESS(SPMatrix(n, H[]), X[])
                else
                    MY_HESS(SPMatrix(n, H[]), X[], unsafe_load(i))
                end
                return
            end), Cvoid, (Ref{GFortranVectorDescriptorRaw{Cdouble}}, Ref{GFortranVectorDescriptorRaw{Cdouble}}, Ptr{Cint}))
        end
    end
    GC.@preserve BL BU CX Y (@ccall libgalahad_double.__lancelot_simple_double_MOD_lancelot_simple(
        n::Ref{Cint},                                   # INTEGER ( KIND = ip_ ), INTENT( IN ) :: n
        _desc(X)::Ref{DoubleGVec},                      # REAL ( KIND = rp_ ), INTENT( INOUT ) :: X( : )
        MY_FUN_F::Ptr{Cvoid},
        fx::Ref{Cdouble},                               # REAL ( KIND = rp_ ), INTENT( OUT ) :: fx
        exit_code::Ref{Cint},                           # INTEGER ( KIND = ip_ ), INTENT( OUT ) :: exit_code
        (ismissing(MY_GRAD) ? C_NULL : MY_GRAD_F)::Ptr{Cvoid},
        (ismissing(MY_HESS) ? C_NULL : MY_HESS_F)::Ptr{Cvoid},
        _desc(BL)::Ptr{DoubleGVec},                     # REAL ( KIND = rp_ ), OPTIONAL :: BL ( : )
        _desc(BU)::Ptr{DoubleGVec},                     # REAL ( KIND = rp_ ), OPTIONAL :: BU ( : )
        C_NULL::Ptr{Cvoid},                             # CHARACTER ( LEN = 10 ), OPTIONAL :: VNAMES( : )
        C_NULL::Ptr{Cvoid},                             # CHARACTER ( LEN = 10 ), OPTIONAL :: CNAMES( : )
        neq::Ref{Cint},                                 # INTEGER ( KIND = ip_ ), OPTIONAL :: neq
                                                        # ^ optional is equivalent to 0, so we just require it
        nin::Ref{Cint},                                 # INTEGER ( KIND = ip_ ), OPTIONAL :: nin
                                                        # ^ optional is equivalent to 0, so we just require it
        _desc(CX)::Ptr{DoubleGVec},                     # REAL ( KIND = rp_ ), OPTIONAL :: CX ( : )
        _desc(Y)::Ptr{DoubleGVec},                      # REAL ( KIND = rp_ ), OPTIONAL :: Y ( : )
        iters::Ref{Cint},                               # INTEGER ( KIND = ip_ ), OPTIONAL :: iters
        maxit::Ref{Cint},                               # INTEGER ( KIND = ip_ ), OPTIONAL :: maxit
                                                        # ^ optional is equivalent to 1000, so we just require it
        gradtol::Ref{Cdouble},                          # REAL ( KIND = rp_ ), OPTIONAL :: gradtol
                                                        # ^ optional is equivalent to 1e-5, so we just require it
        feastol::Ref{Cdouble},                          # REAL ( KIND = rp_ ), OPTIONAL :: feastol
                                                        # ^ optional is equivalent to 1e-5, so we just require it
        print_level::Ref{Cint},                         # INTEGER ( KIND = ip_ ), OPTIONAL :: print_level
                                                        # ^ optional is equivalent to 1, so we just require it
        0::Csize_t,                                     # maybe LEN(VNAMES)?
        0::Csize_t,                                     # maybe LEN(CNAMES)?
    )::Cvoid)
    return fx[], iters[], exit_code[]
end
#endregion

#region LANCELOT - full
include("./Structs.jl")

"""
    LANCELOT_initialize()

Returns a `(data, control)` tuple to be used in subsequent optimization with [`LANCELOT_solve`](@ref) and to be cleaned up
using [`LANCELOT_terminate`](@ref).
"""
function LANCELOT_initialize()
    data = LANCELOT_data_type()
    control = LANCELOT_control_type()
    @ccall libgalahad_double.__lancelot_double_MOD_lancelot_initialize(
        data::LANCELOT_data_type,      # TYPE ( LANCELOT_data_type ), INTENT( INOUT ) :: data
        control::LANCELOT_control_type # TYPE ( LANCELOT_control_type ), INTENT( OUT ) :: control
    )::Cvoid
    return data, control
end

@doc raw"""
    LANCELOT_solve(prob; RANGE[, GVALS][, FT][, XT][, FUVALS][, ICALCF][, ICALCG][, IVAR][, Q][,
        DGRAD,] control, inform, data[, ELDERS][, ELFUN][, ELFUN_flexible][, GROUP])

- `prob` is a [`LANCELOT_problem_type`](@ref). It is used to hold data about the problem being solved. Components `IELING`,
  `ISTADG`, `IELVAR`, `ISTAEV`, `ICNA`, `ISTADA`, `A`, `B`, `BL`, `BU`, `GSCALE`, `ESCALE`, `VSCALE`, `GXEQX`, `INTREP`,
  `VNAMES`, and `GNAMES` must all be set on entry, and will thereafter be unaltered.

  The component `INTVAR` must also be set, but will subsequently be reset by `LANCELOT_solve` so that its ``i``-th value gives
  the position in the array `FUVALS` of the first component of the gradient of the ``i``-th nonlinear element function, with
  respect to its internal variables (see `FUVALS`). The component ``X`` must also be set, but will be altered by
  `LANCELOT_solve`. The component `ISTADH` need not be set on initial entry, but will be set by `LANCELOT_solve`.

  If the dummy argument `ELFUN` (see below) is present, the components `ISTEPA` and `EPVALU` must also be set on entry—they
  will not subsequently be altered—but otherwise they need not be set. Likewise if the dummy argument `GROUP` (see below) is
  present, the components `ISTGPA`, `GPVALU` and `ITYPEG` must also be set on entry—they will not subsequently be altered—but
  otherwise they need not be set.

  If the problem involves general constraints, the components `KNDOFG`, `Y` and `C` must be set, and the first two assigned
  values—`KNDOFG` will not subsequently be altered, while `C` will be set by `LANCELOT_solve`.
  If the problem does not involve general constraints, `KNDOFG`, `Y` and `C` need not be set.
- `RANGE` is a user-supplied subroutine whose purpose is to define the linear transformation of variables for those non-linear
  elements which have different elemental and internal variables. See below for details.
- `GVALS::AbstractMatrix{Cdouble}` is a matrix of shape `(prob.ng, 3)` that is used to store function and derivative
  information for the group functions. The user may be asked to provide values for these functions and/or derivatives,
  evaluated at the argument `FT` when control is returned to the calling program with a negative value of the variable
  `inform.status`. This information needs to be stored by the user in specified locations within `GVALS`. Details of the
  required information are given in Section 2.4.7 of the manual.
- `FT::AbstractVector{Cdouble}` is a vector of dimension `prob.ng` that is set within `LANCELOT_solve` to a trial value of the
  argument of the ``i``-th group function at which the user may be required to evaluate the values and/or derivatives of that
  function. Precisely what group function information is required at ``FT`` is under the control of the variable
  `inform.status` and details are given in Section 2.4.7 of the manual.
- `XT::AbstractVector{Cdouble}` is a vector of dimension `prob.n` that is set within `LANCELOT_solve` to a trial value of the
  variables `\vec x` at which the user may be required to evaluate the values and/or derivatives of the nonlinear elements
  functions. Precisely what element function information is required at `XT` is under the control of the variable
  `inform.status` and details are given in Section 2.4.7 of the manual.
- `FUVALS::AbstractVector{Cdouble}` is a vector that is used to store function and derivative information for the nonlinear
  element functions. The user is asked to provide values for these functions and/or derivatives, evaluated at the arguments
  `XT`, at specified locations within `FUVALS`, when control is returned to the calling program with a negative value of the
  variable `inform.status`. Alternatively, the user may have provided a suitable subroutine `ELFUN` to compute the required
  function or derivative values (see below). Details of the required information are given in Sections 2.4.4 and 2.4.7 of the
  manual.

  The first segment of `FUVALS` contains the values of the nonlinear element functions; the next two segments contain their
  gradients and Hessian matrices, taken with respect to their internal variables. The remaining two used segments contain the
  gradient of the objective function and the diagonal elements of the second derivative approximation, respectively. At the
  solution, the components of the gradient of the augmented Lagrangian function corresponding to variables which lie on one of
  their bounds are of particular interest in many applications areas. In particular they are often called shadow prices, and
  are used to assess the sensitivity of the solution to variations in the bounds on the variables.

  The required length of `FUVALS` may be calculated using the following code segment:
  ```julia
    prob.nel + 2prob.n + sum(i -> i * (i +3) ÷ 2, Iterators.take(prob.INTVAR, prob.nel))
  ```
- `ICALCF::AbstractVector{Cint}` is a vector of dimension `prob.nel`. If the user has chosen not to perform internal element
  evaluations, and if the value of `inform.status` on return from `LANCELOT_solve` indicates that further element functions
  values or their derivatives are required prior to a reentry, the first `inform.ncalcf` components of `ICALCF` give the
  indices of the group functions which need to be recalculated at `XT`. Precisely what group function information is required
  is under the control of the variable `inform.status` and details are given in Section 2.4.7 of the manual.
- `ICALCG::AbstractVector{Cint}` is a vector of dimension `prob.ng`. If the user has chosen not to perform internal group
  evaluations, and if the value of `inform.status` on return from `LANCELOT_solve` indicates that further group functions
  values or their derivatives are required prior to a re-entry, the first `inform.ncalcg` components of `ICALCG` give the
  indices of the group functions which need to be recalculated at `FT`. Precisely what group function information is required
  is under the control of the variable `inform.status` and details are given in Section 2.4.7 of the manual.
- `IVAR::AbstractVector{Cint}` in a vector of dimension `prob.n` that is required when the user is providing a special
  preconditioner for the conjugate gradient inner iteration.
- `Q` and `DGRAD` are `AbstractVector{Cdouble}`s of dimension `prob.n` that are required when the user is providing a special
  preconditioner for the conjugate gradient inner iteration.
- `control` is a [`LANCELOT_control_type`](@ref). On exit, `control` contains default values for the components. These values
  should only be changed after calling [`LANCELOT_initialize`](@ref).
- `inform` is a [`LANCELOT_inform_type`](@ref). The component status must be set to `0` on initial entry, and a successful call
  to `LANCELOT_solve` is indicated when the component status has the value `0`. For other return values of status, see Sections
  2.4.7 and 2.5 of the manual.
- `data` is a [`LANCELOT_data_type`](@ref). It is used to hold data about the problem being solved. It must not have been
  altered by the user since the last call to [`LANCELOT_initialize`](@ref).
- `ELDERS::AbstractMatrix{Cint}` is an optional matrix argument of shape `(2, prob.nel)` that may be used to specify what kind
  of first- and second-derivative approximations will be required for each nonlinear element function. If `ELDERS` is not
  present, the first- and second-derivative requirements for every element will be as specified by the variables
  `control.first` derivatives and `control.second` derivatives respectively. For finer control, if `ELDERS` is present, and the
  user is able to provide analytical first derivatives for the ``i``-th nonlinear element function, `ELDERS[1, i]` must be set
  to to be non-positive. If analytical first derivatives are unavailable, they may be estimated by forward differences by
  setting `ELDERS[1, i] = 1` or, more accurately but at additional expense, by central differences by setting
  `ELDERS[1, i] ≥ 2`. Similarly, if `ELDERS` is present, and the user is able to provide analytical second derivatives for the
  ``i``-th nonlinear element function, `ELDERS[2, i]` must be set to to be non-positive. If the user is unable to provide
  second derivatives, these derivatives will be approximated using one of four secant approximation formulae. If `ELDERS[2, i]`
  is set to `1`, the BFGS formula is used; if it is set to `2`, the DFP formula is used; if it is set to `3`, the PSB formula
  is used; and if it is set to `4` or larger, the symmetric rank-one formula is used. The user is strongly advised to use
  analytic first and second derivatives if at all possible as this often significantly improves the convergence of the method.
- `ELFUN` and `ELFUN_flexible` are optional user-supplied functions whose purpose is to evaluate the values and derivatives of
  the nonlinear element functions. See Section 2.4.1 of the manual for background information and below for details. Only one
  of `ELFUN` and `ELFUN_flexible` may be present at once. Which of the two arguments is permitted depends on whether `ELDERS`
  (see above) is present. `ELFUN_flexible` is only permitted when `ELDERS` is present, since `ELFUN_flexible` allows the user
  to provide element-specific levels of derivative information. In the absence of `ELDERS`, `ELFUN` should be used instead of
  `ELFUN_flexible`. If both `ELFUN` and `ELFUN_flexible` are absent, `LANCELOT_solve` will use reverse communication to obtain
  element function values and derivatives.
- `GROUP` is an optional user-supplied subroutine whose purpose is to evaluate the values and derivatives of the group
  functions. See Section 2.4.1 of the manual for background information and below for details. If `GROUP` is absent,
  `LANCELOT_solve` will use reverse communication to obtain group function values and derivatives. `GROUP` need not be present
  if all components of `prob.GXEQX` are `true`.

# Callback arguments
## `RANGE(ielemn, transp, W1, W2, nelvar, ninvar, ieltyp)`
The purpose of the `RANGE` callback is to define the transformation between internal and elemental variables for nonlinear
elements with useful internal representations.
- `ielemn::Int` gives the index of the nonlinear element whose transformation is required by `LANCELOT_solve`.
- `transp::Bool`. If `transp` is `false`, the callback must put the result of the transformation ``W \vec v`` in the array
  `W2`, where ``\vec v`` is input in the array `W1`. Otherwise, the callback must supply the result of the transposed
  transformation ``W^\top \vec u`` in the array `W2`, where ``\vec u`` is input in the array `W1`.
- `W1::Vector{Cdouble}` is a vector whose dimension is the number of elemental variables if `transp` is `false` and the number
  of internal variables otherwise.
- `W2::Vector{Cdouble}` is a vector whose dimension is the number of internal variables if `transp` is `false` and the number
  of elemental variables otherwise. The result of the transformation of `W1` or its transpose, as defined by `transp`, must be
  set in `W2`.
- `nelvar::Int` gives the number of elemental variables for the element specified by `ielemn`.
- `ninvar::Int` gives the number of internal variables for the element specified by `ielemn`.
- `ieltyp::Int` defines the type for the element specified by `ielemn`.

Restriction: `length(W1) ≥ ninvar` if `transp` is `true` and `length(W1) ≥ nelvar` if transp is `false`.
Restriction: `length(W2) ≥ nelvar` if `transp` is `true` and `length(W2) ≥ ninvar` if transp is `false`.

The user will already have specified which elements have useful transformations in the array `INTREP`. `RANGE` will only be
called for elements for which the corresponding component of `INTREP` is `true`.

## `ELFUN(FUVALS, XVALUE, EPVALU, ncalcf, ITYPEE, ISTAEV, IELVAR, INTVAR, ISTADH, ISTEPA, ICALCF, ifflag) -> Integer`
If the argument `ELFUN` is present when calling `LANCELOT_solve`, the user is expected to provide a callback to evaluate
function or derivative values (with respect to internal variables) of (a given subset of) the element functions.
- `FUVALS::Vector{Cdouble}` is used to store function and derivative information for the nonlinear element functions. The
  callback is asked to provide values for these functions and/or derivatives, evaluated at the argument `XVALUE`, at specified
  locations within `FUVALS`.

  The first segment of `FUVALS` is used to hold the values of the nonlinear element functions, component ``i`` holding the
  value of the ``i``-th element function. The second segment holds the components of the gradients of the element functions
  taken with respect to their internal variables, as described by `INTVAR` below. The final segment contains the elements'
  Hessian matrices, again taken with respect to their internal variables, as described by `ISTADH` below.
- `XVALUE::Vector{Cdouble}` contains the values of ``x`` at which the callback is required to evaluate the values or
  derivatives of the nonlinear elements functions.
- `EPVALU::Vector{Cdouble}` contains the values of the element parameters ``\var p^e_i``, ``i = 1, \dots, n_e``. The indices
  for the parameters for element ``i`` immediately precede those for element ``i + 1``, and each element's parameters appear in
  a contiguous list.
- `ncalcf::Int` specifies how many of the nonlinear element functions or their derivatives are to be evaluated.
- `ITYPEE::Vector{Cint}`'s ``i``-th component specifies the type of element ``i``.
- `ISTAEV::Vector{Cint}`: see [`LANCELOT_problem_type`](@ref)
- `IELVAR::Vector{Cint}`: see [`LANCELOT_problem_type`](@ref)
- `INTVAR::Vector{Cint}`'s ``i``-th component (``1 \leq i \leq n_e``) gives the position in the array `FUVALS` of the first
  component of the gradient of the ``i``-th nonlinear element function, with respect to its internal variables.
- `ISTADH::Vector{Cint}`'s ``i``-th component (``1 \leq i \leq n_e``) gives the position in the array `FUVALS` of the first
  component of the Hessian matrix of the ``i``-th nonlinear element function ``e_i``, with respect to its internal variables.
  Only the upper triangular part of each Hessian matrix is stored and the storage is by columns. That is to say that the
  component of the Hessian of the ``k``-th nonlinear element with respect to internal variables ``i`` and ``j``, ``i \leq j``,
  ``\frac{\partial^2 e_k}{\partial u_i\partial u_j}`` must be placed in `FUVALS[ISTADH[k]+(j(j-1)÷2)+i-1]`. The element
  `ISTADH[ne + 1]` is space required to finish storing the Hessian of the last nonlinear element in `FUVALS` plus one.
- `ISTEPA::Vector{Cint}`'s ``i``-th component gives the position in `EPVALU` of the first parameter for element function ``i``.
  In addition, `ISTEPA[ne+1]` is the position in `EPVALU` of the last parameter for element function ``n_e`` plus one.
- `ICALCF::Vector{Cint}`'s first `ncalcf` components gives the indices of the nonlinear element functions whose values or
  derivatives are to be evaluated.
- `ifflag::Int` defines whether it is the values of the element functions that are required (`ifflag = 1`) or if it is the
  derivatives (`ifflag > 1`). Possible values and their requirements are:
  - `ifflag = 1`. The values of nonlinear element functions `ICALCF[i]`, `i = 1, ..., ncalcf`, are to be computed and placed in
    `FUVALS[ICALCF[i]]`.
  - `ifflag = 2`. The gradients of nonlinear element functions `ICALCF[i]`, `i = 1, ..., ncalcf`, are to be computed. The
    gradient of the `ICALCF[i]`-th element is to be placed in the segment of `FUVALS` starting at `INTVAR[ICALCF[i]]`.
  - `ifflag = 3`. The gradients and Hessians of nonlinear element functions `ICALCF[i]`, `i = 1, ..., ncalcf`, are to be
    computed. The gradient of the `ICALCF[i]`-th element is to be placed in the segment of `FUVALS` starting at
    `INTVAR[ICALCF[i]]`. The Hessian of this element should be placed in the segment of `FUVALS` starting at
    `ISTADH[ICALCF[i]]`.
  *N.B.* If the user intends to use approximate second derivatives (`control.second_derivatives > 0`), `LANCELOT_solve` will
  never call `ELFUN` with `ifflag = 3`, so the user need not provide Hessian values. Furthermore, if the user intends to use
  approximate first derivatives (`control.first_derivatives > 0`), `LANCELOT_solve` will never call `ELFUN` with `ifflag = 2`
  or `3`, so the user need then not provide gradient or Hessian values.

Restrictions: `length(ITYPEE) ≥ prob.nel`, `length(ISTAEV) ≥ prob.nel +1`, `length(IELVAR) ≥ ISTAEV[end] -1`,
`length(INTVAR) ≥ prob.nel +1`, `length(ISTADH) ≥ prob.nel +1`, `length(ISTEPA) ≥ prob.nel +1`, `length(ICALCF) ≥ ncalcf`,
`length(FUVALS) ≥ ISTADH[end] -1`, `length(XVALUE) ≥ prob.n`, and `length(EPVALU) ≥ ISTEPA(prob.nel +1) -1`.

Return value: The return value should be set to `0` if all of the required values have been found, and to any nonzero value if,
for any reason, one or more of the required values could not be determined. For instance, if the value of a nonlinear element
(or its derivative) was required outside of its domain of definition, a nonzero value should be returned.

## `ELFUN_FLEXIBLE(FUVALS, XVALUE, EPVALU, ncalcf, ITYPEE, ISTAEV, IELVAR, INTVAR, ISTADH, ISTEPA, ICALCF, ifflag, ELDERS) -> Integer`
If the argument `ELFUN_flexible` is present when calling `LANCELOT_solve`, the user is expected to provide a callback to
evaluate function or derivative values (with respect to internal variables) of (a given subset of) the element functions.
- see `ELFUN` for all parameters not documented in the following
- `ifflag::Int` defines whether it is the values of the element functions that are required (`ifflag = 1`) or if it is the
  derivatives (`ifflag > 1`).
  - If `ifflag = 1`, the values of nonlinear element functions `ICALCF[i]`, `i = 1, ..., ncalcf`, are to be computed and placed
    in `FUVALS[ICALCF[i]]`.
  - If `ifflag > 1`, those derivative values specified by ELDERS (see below) are required.
- `ELDERS::Matrix{Cint}` of shape `(2, ?)` specifies what first- and second-derivative information (if any) is required for
  each nonlinear element function when `ifflag > 1`. Specifically, in this case, the gradient of each nonlinear element
  function `ICALCF[i]`, `i = 1, ..., ncalcf`, for which `ELDERS[1, ICALCF[i]] ≤ 0` should be computed. The gradient of the
  `ICALCF[i]`-th element is to be placed in the segment of `FUVALS` starting at `INTVAR[ICALCF[i]]`. Furthermore, if
  additionally `ELDERS[2,ICALCF[i]] ≤ 0`, the Hessian of the nonlinear element function `ICALCF[i]` should be computed and
  placed in the segment of `FUVALS` starting at `ISTADH[ICALCF[i]]`.

Restriction: `size(ELDERS, 2) ≥ prob.nel`.

## `GROUP(GVALUE, FVALUE, GPVALU, ncalcg, ITYPEG, ISTGPA, ICALCG, derivs) -> Integer`
If the argument `GROUP` is present when calling `LANCELOT_solve`, the user is expected to provide a callback to evaluate
function or derivative values of (a given subset of) the group functions.
- `GVALUE::Matrix{Cdouble}` of shape `(?, 3)`. The value and first and second derivative of the ``i``-th group function are
  held in `GVALUE[i, 1]`, `GVALUE[i, 2]` and `GVALUE[i, 3]` respectively.
- `FVALUE::Vector{Cdouble}`'s ``i``-th component contains the value of group variable at which the subroutine is required to
  evaluate the value or derivatives of the ``i``-th group function.
- `GPVALU::Vector{Cdouble}` contains the values of the group parameters ``\vec p^g_i``, ``i = 1, \dots, n_g``. The indices for
  the parameters for group ``i`` immediately precede those for group ``i + 1``, and each group's parameters appear in a
  contiguous list.
- `ncalcg::Int` specifies how many of the group functions or their derivatives are to be evaluated.
- `ITYPEG::Vector{Cint}`'s ``i``-th component specifies the type of group ``i``.
- `ISTGPA::Vector{Cint}`'s ``i``-th component gives the position in `GPVALU` of the first parameter for group function ``i``.
  In addition, `ISTGPA[ng+1]` is the position in `GPVALU` of the last parameter for element function ``n_g`` plus one.
- `ICALCG::Vector{Cint}`'s first `ncalcg` components gives the indices of the nonlinear group functions whose values or
  derivatives are to be evaluated.
- `derivs::Bool`. When `derivs` is `false`, the callback must return the values of group functions `ICALCG[i]`,
  `i = 1, ..., ncalcg` in `GVALUE[ICALCF[i], 1]`. When `derivs` is `true`, the callback must return the first and second
  derivatives of group functions `ICALCG[i]`, `i = 1, ..., ncalcg` in `GVALUE[ICALCF[i], 2]` and `GVALUE[ICALCF[i], 3]`
  respectively.

Restrictions: `size(GVALUE, 2) ≥ prob.ng`, `length(ITYPEG) ≥ prob.n`, `length(ISTGPA) ≥ prob.ng +1`, `length(ICALCG) ≥ ncalcg`,
`length(FVALUE) ≥ prob.ng`, and `length(GPVALU) ≥ ISTGPA[end] -1`.

Return value: The return value should be set to `0` if all of the required values have been found, and to any nonzero value if,
for any reason, one or more of the required values could not be determined. For instance, if the value of a group function (or
its derivative) was required outside of its domain of definition, a nonzero value should be returned.
"""
function LANCELOT_solve(prob::LANCELOT_problem_type; RANGE::Base.Callable,
    GVALS::AbstractMatrix{Cdouble}=Matrix{Cdouble}(undef, prob.ng, 3),
    FT::AbstractVector{Cdouble}=Vector{Cdouble}(undef, prob.ng),
    XT::AbstractVector{Cdouble}=Vector{Cdouble}(undef, prob.n),
    FUVALS::AbstractVector{Cdouble}=Vector{Cdouble}(undef, prob.nel + 2prob.n +
                                                    sum(i -> i * (i +3) ÷ 2, Iterators.take(prob.INTVAR, prob.nel), init=0)),
    ICALCF::AbstractVector{Cint}=Vector{Cint}(undef, prob.nel),
    ICALCG::AbstractVector{Cint}=Vector{Cint}(undef, prob.ng),
    IVAR::AbstractVector{Cint}=Vector{Cint}(undef, prob.n),
    Q::AbstractVector{Cdouble}=Vector{Cdouble}(undef, prob.n),
    DGRAD::AbstractVector{Cdouble}=Vector{Cdouble}(undef, prob.n),
    control::LANCELOT_control_type, inform::LANCELOT_inform_type, data::LANCELOT_data_type,
    ELFUN::Union{Base.Callable,Missing}=missing, GROUP::Union{Base.Callable,Missing}=missing,
    ELFUN_FLEXIBLE::Union{Base.Callable,Missing}=missing, ELDERS::Union{<:AbstractVector{Cint},Missing}=missing)
    LinearAlgebra.chkstride1(GVALS)
    LinearAlgebra.chkstride1(FT)
    LinearAlgebra.chkstride1(XT)
    LinearAlgebra.chkstride1(FUVALS)
    LinearAlgebra.chkstride1(ICALCF)
    LinearAlgebra.chkstride1(DGRAD)
    ismissing(ELDERS) || LinearAlgebra.chkstride1(ELDERS)
    ((ismissing(ELFUN) || ismissing(ELDERS)) &&
        (ismissing(ELFUN_FLEXIBLE) || !ismissing(ELDERS))) ||
        throw(ArgumentError("Wrong combination of ELDERS, ELFUN, and ELFUN_FLEXIBLE"))
    ismissing(ELDERS) && ismissing(ELFUN_FLEXIBLE)
    if size(GVALS) != (prob.ng, 3) || length(FT) != prob.ng || length(XT) != prob.n || length(ICALCF) != prob.nel ||
        length(ICALCG) != prob.ng || length(IVAR) != prob.n || length(Q) != prob.n || length(DGRAD) != prob.n ||
        (!ismissing(ELDERS) && size(ELDERS) != (2, prob.nel))
        throw(ArgumentError("Invalid dimensions"))
    end
    #=
        SUBROUTINE RANGE ( ielemn, transp, W1, W2, nelvar, ninvar, ieltyp, lw1, lw2 )
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: ielemn, nelvar, ninvar, ieltyp
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: lw1, lw2
            LOGICAL, INTENT( IN ) :: transp
            REAL ( KIND = rp_ ), INTENT( IN ), DIMENSION( lw1 ) :: W1
            REAL ( KIND = rp_ ), DIMENSION( lw2 ) :: W2
        END SUBROUTINE RANGE
    =#
    RANGE_F = @cfunction($((ielemn, transp, W1, W2, nelvar, ninvar, ieltyp, lw1, lw2) -> begin
        println("In RANGE")
        flush(stdout)
        RANGE(
            Int(unsafe_load(ielemn)),
            !iszero(unsafe_load(transp)),
            unsafe_wrap(Array, W1, unsafe_load(lw1)),
            unsafe_wrap(Array, W2, unsafe_load(lw2)),
            Int(unsafe_load(nelvar)),
            Int(unsafe_load(ninvar)),
            Int(unsafe_load(ieltyp))
        )
        return
    end), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}))
    #=
        SUBROUTINE ELFUN ( FUVALS, XVALUE, EPVALU, ncalcf, ITYPEE, ISTAEV,      &
                           IELVAR, INTVAR, ISTADH, ISTEPA, ICALCF, ltypee,      &
                           lstaev, lelvar, lntvar, lstadh, lstepa, lcalcf,      &
                           lfuval, lxvalu, lepvlu, ifflag, ifstat )
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: ncalcf, ifflag, ltypee, lstaev
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: lelvar, lntvar, lstadh, lstepa
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: lcalcf, lfuval, lxvalu, lepvlu
            INTEGER ( KIND = ip_ ), INTENT( OUT ) :: ifstat
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: ITYPEE(ltypee), ISTAEV(lstaev)
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: IELVAR(lelvar), INTVAR(lntvar)
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: ISTADH(lstadh), ISTEPA(lstepa)
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: ICALCF(lcalcf)
            REAL ( KIND = rp_ ), INTENT( IN ) :: XVALUE(lxvalu)
            REAL ( KIND = rp_ ), INTENT( IN ) :: EPVALU(lepvlu)
            REAL ( KIND = rp_ ), INTENT( INOUT ) :: FUVALS(lfuval)
       END SUBROUTINE ELFUN
    =#
    if !ismissing(ELFUN)
        ELFUN_F = @cfunction(
            $((FUVALS, XVALUE, EPVALU, ncalcf, ITYPEE, ISTAEV, IELVAR, INTVAR, ISTADH, ISTEPA, ICALCF, ltypee, lstaev, lelvar,
               lntvar, lstadh, lstepa, lcalcf, lfuval, lxvalu, lepvlu, ifflag, ifstat) -> begin
                println("in ELFUN")
                flush(stdout)
                unsafe_store!(ifstat,
                    ELFUN(
                        unsafe_wrap(Array, FUVALS, unsafe_load(lfuval)),
                        unsafe_wrap(Array, XVALUE, unsafe_load(lxvalu)),
                        unsafe_wrap(Array, EPVALU, unsafe_load(lepvlu)),
                        Int(unsafe_load(ncalcf)),
                        unsafe_wrap(Array, ITYPEE, unsafe_load(ltypee)),
                        unsafe_wrap(Array, ISTAEV, unsafe_load(lstaev)),
                        unsafe_wrap(Array, IELVAR, unsafe_load(lelvar)),
                        unsafe_wrap(Array, INTVAR, unsafe_load(lntvar)),
                        unsafe_wrap(Array, ISTADH, unsafe_load(lstadh)),
                        unsafe_wrap(Array, ISTEPA, unsafe_load(lstepa)),
                        unsafe_wrap(Array, ICALCF, unsafe_load(lcalcf)),
                        unsafe_load(ifflag)
                    )
                )
                return
            end), Cvoid,
            (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
             Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
             Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
             Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint})
        )
    end
    #=
        SUBROUTINE GROUP ( GVALUE, lgvalu, FVALUE, GPVALU, ncalcg,              &
                           ITYPEG, ISTGPA, ICALCG, ltypeg, lstgpa,              &
                           lcalcg, lfvalu, lgpvlu, derivs, igstat )
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: lgvalu, ncalcg, ltypeg, lstgpa
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: lcalcg, lfvalu, lgpvlu
            INTEGER ( KIND = ip_ ), INTENT( OUT ) :: igstat
            LOGICAL, INTENT( IN ) :: derivs
            INTEGER ( KIND = ip_ ), INTENT( IN ), DIMENSION( ltypeg ) :: ITYPEG
            INTEGER ( KIND = ip_ ), INTENT( IN ), DIMENSION( lstgpa ) :: ISTGPA
            INTEGER ( KIND = ip_ ), INTENT( IN ), DIMENSION( lcalcg ) :: ICALCG
            REAL ( KIND = rp_ ), INTENT( IN ), DIMENSION( lfvalu ) :: FVALUE
            REAL ( KIND = rp_ ), INTENT( IN ), DIMENSION( lgpvlu ) :: GPVALU
            REAL ( KIND = rp_ ), INTENT( INOUT ), DIMENSION( lgvalu, 3 ) :: GVALUE
        END SUBROUTINE GROUP
    =#
    if !ismissing(GROUP)
        GROUP_F = @cfunction(
            $((GVALUE, lgvalu, FVALUE, GPVALU, ncalcg, ITYPEG, ISTGPA, ICALCG, ltypeg, lstgpa, lcalcg, lfvalu, lgpvlu, derivs,
               igstat) -> begin
                println("In GROUP")
                flush(stdout)
                unsafe_store!(igstat,
                    GROUP(
                        unsafe_wrap(Array, GVALUE, (unsafe_load(lgvalu), 3)),
                        unsafe_wrap(Array, FVALUE, unsafe_load(lfvalu)),
                        unsafe_wrap(Array, GPVALU, unsafe_load(lgpvlu)),
                        Int(unsafe_load(ncalcg)),
                        unsafe_wrap(Array, ITYPEG, unsafe_load(ltypeg)),
                        unsafe_wrap(Array, ISTGPA, unsafe_load(lstgpa)),
                        unsafe_wrap(Array, ICALCG, unsafe_load(lcalcg)),
                        !iszero(unsafe_load(derivs))
                    )
                )
                return
            end), Cvoid,
            (Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint},
             Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
             Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint})
        )
    end
    #=
        SUBROUTINE ELFUN_flexible ( FUVALS, XVALUE, EPVALU, ncalcf, ITYPEE, ISTAEV,      &
                                    IELVAR, INTVAR, ISTADH, ISTEPA, ICALCF, ltypee,      &
                                    lstaev, lelvar, lntvar, lstadh, lstepa, lcalcf,      &
                                    lfuval, lxvalu, lepvlu, llders, ifflag, ELDERS,      &
                                    ifstat )
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: ncalcf, ifflag, ltypee, lstaev
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: lelvar, lntvar, lstadh, lstepa
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: lcalcf, lfuval, lxvalu, lepvlu
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: llders
            INTEGER ( KIND = ip_ ), INTENT( OUT ) :: ifstat
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: ITYPEE(ltypee), ISTAEV(lstaev)
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: IELVAR(lelvar), INTVAR(lntvar)
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: ISTADH(lstadh), ISTEPA(lstepa)
            INTEGER ( KIND = ip_ ), INTENT( IN ) :: ICALCF(lcalcf), ELDERS(2,llders)
            REAL ( KIND = rp_ ), INTENT( IN ) :: XVALUE(lxvalu)
            REAL ( KIND = rp_ ), INTENT( IN ) :: EPVALU(lepvlu)
            REAL ( KIND = rp_ ), INTENT( INOUT ) :: FUVALS(lfuval)
        END SUBROUTINE ELFUN_flexible
    =#
    if !ismissing(ELFUN_FLEXIBLE)
        ELFUN_FLEXIBLE_F = @cfunction(
            $((FUVALS, XVALUE, EPVALU, ncalcf, ITYPEE, ISTAEV, IELVAR, INTVAR, ISTADH, ISTEPA, ICALCF, ltypee, lstaev, lelvar,
               lntvar, lstadh, lstepa, lcalcf, lfuval, lxvalu, lepvlu, llders, ifflag, ELDERS, ifstat) -> begin
                println("in ELFUN_FLEXIBLE")
                flush(stdout)
                unsafe_store!(ifstat,
                    ELFUN(
                        unsafe_wrap(Array, FUVALS, unsafe_load(lfuval)),
                        unsafe_wrap(Array, XVALUE, unsafe_load(lxvalu)),
                        unsafe_wrap(Array, EPVALU, unsafe_load(lepvlu)),
                        Int(unsafe_load(ncalcf)),
                        unsafe_wrap(Array, ITYPEE, unsafe_load(ltypee)),
                        unsafe_wrap(Array, ISTAEV, unsafe_load(lstaev)),
                        unsafe_wrap(Array, IELVAR, unsafe_load(lelvar)),
                        unsafe_wrap(Array, INTVAR, unsafe_load(lntvar)),
                        unsafe_wrap(Array, ISTADH, unsafe_load(lstadh)),
                        unsafe_wrap(Array, ISTEPA, unsafe_load(lstepa)),
                        unsafe_wrap(Array, ICALCF, unsafe_load(lcalcf)),
                        unsafe_load(ifflag),
                        unsafe_wrap(Array, ELDERS, (2, unsafe_load(llders)))
                    )
                )
                return
            end), Cvoid,
            (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
             Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
             Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
             Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
             Ptr{Cint})
        )
    end
    GC.@preserve ELDERS begin
        @ccall libgalahad_double.__lancelot_double_MOD_lancelot_solve(
            prob::LANCELOT_problem_type,    # TYPE ( LANCELOT_problem_type ), INTENT( INOUT ), TARGET :: prob
            RANGE_F::Ptr{Cvoid},
            GVALS::Ref{Cdouble},            # REAL ( KIND = rp_ ), INTENT( INOUT ), DIMENSION( prob%ng, 3 ) :: GVALS
            FT::Ref{Cdouble},               # REAL ( KIND = rp_ ), INTENT( INOUT ), DIMENSION( prob%ng ) :: FT
            XT::Ref{Cdouble},               # REAL ( KIND = rp_ ), INTENT( INOUT ), DIMENSION( prob%n ) :: XT
            FUVALS::Ref{Cdouble},           # REAL ( KIND = rp_ ), INTENT( INOUT ), DIMENSION( lfuval ) :: FUVALS
            length(FUVALS)::Ref{Cint},      # INTEGER ( KIND = ip_ ), INTENT( IN ) :: lfuval
            ICALCF::Ref{Cint},              # INTEGER ( KIND = ip_ ), INTENT( INOUT ), DIMENSION( prob%nel ) :: ICALCF
            ICALCG::Ref{Cint},              # INTEGER ( KIND = ip_ ), INTENT( INOUT ), DIMENSION( prob%ng ) :: ICALCG
            IVAR::Ref{Cint},                # INTEGER ( KIND = ip_ ), INTENT( INOUT ), DIMENSION( prob%n  ) :: IVAR
            Q::Ref{Cdouble},                # REAL ( KIND = rp_ ), INTENT( INOUT ), DIMENSION( prob%n ) :: Q
            DGRAD::Ref{Cdouble},            # REAL ( KIND = rp_ ), INTENT( INOUT ), DIMENSION( prob%n ) :: DGRAD
            control::LANCELOT_control_type, # TYPE ( LANCELOT_control_type ), INTENT( INOUT ) :: control
            inform::LANCELOT_inform_type,   # TYPE ( LANCELOT_inform_type ), INTENT( INOUT ) :: inform
            data::LANCELOT_data_type,       # TYPE ( LANCELOT_data_type ), INTENT( INOUT ) :: data
            (ismissing(ELFUN) ? C_NULL : ELFUN_F)::Ptr{Cvoid},
            (ismissing(GROUP) ? C_NULL : GROUP_F)::Ptr{Cvoid},
            (ismissing(ELFUN_FLEXIBLE) ? C_NULL : ELFUN_FLEXIBLE_F)::Ptr{Cvoid},
            (ismissing(ELDERS) ? Ptr{Cint}(C_NULL) : ELDERS)::Ptr{Cint}
                                        # INTEGER ( KIND = ip_ ), INTENT( INOUT ), OPTIONAL, DIMENSION( 2, prob%nel ) :: ELDERS
        )::Cvoid
    end
end

"""
    LANCELOT_terminate(data, control, inform)

All previously allocated arrays are deallocated with this function.
- `data::LANCELOT_data_type` exactly as for [`LANCELOT_solve`](@ref) that must not have been altered by the user since the
  last call to [`LANCELOT_initialize`](@ref).
- `control::LANCELOT_control_type`
- `inform::LANCELOT_inform_type`. Only the components `status`, `alloc_status` and `bad_alloc` might have been altered on exit,
  and a successful call to [`LANCELOT_terminate`](@ref) is indicated when this component `status` has the value `0`. For other
  return values of `status`, see Section 2.5 of the manual.
"""
function LANCELOT_terminate(data::LANCELOT_data_type, control::LANCELOT_control_type, inform::LANCELOT_inform_type)
    @ccall libgalahad_double.__lancelot_double_MOD_lancelot_terminate(
        data::LANCELOT_data_type,       # TYPE ( LANCELOT_data_type ), INTENT( INOUT ) :: data
        control::LANCELOT_control_type, # TYPE ( LANCELOT_control_type ), INTENT( IN ) :: control
        inform::LANCELOT_inform_type    # TYPE ( LANCELOT_inform_type ), INTENT( INOUT ) :: inform
    )::Cvoid
end
#endregion