export LANCELOT_simple

# We choose the same approach as GALAHAD.jl. We might even want to introduce a weak dependency on it once it provides a JLL.
if haskey(ENV, "JULIA_GALAHAD_LIBRARY_PATH")
    const libgalahad_double = joinpath(ENV["JULIA_GALAHAD_LIBRARY_PATH"], "libgalahad_double.$(Libc.Libdl.dlext)")
else
    @warn("For Lancelot support, enter the complete file path to the libgalahad_double_d shared library")
    const libgalahad_double = readline()
end
if isfile(libgalahad_double)

#region Fortran Array Descriptors
# Note that LANCELOT is the only part of GALAHAD that does not have a C interface yet. Hence, we first need to provide some
# Fortran magic - unfortunately, all the array arguments are passed as assumed-shape, so that we have to pass array
# descriptors. Warning: Here, we assume that GALAHAD was compiled using GFortran with a version ≥ 9. This is a moderately safe
# assumption, as the Intel compiler does not work out-of-the-box... but in case you wish to adapt this to Intel, this is pretty
# easy, as their array descriptor is well-documented:
# https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2023-1/handle-fortran-array-descriptors.html
# In contrast, the GFortran array descriptor is not documented properly, the following structure was extracted directly from
# the source code:
# https://github.com/gcc-mirror/gcc/blob/386df7ce7b38ef00e28080a779ef2dfd6949cf15/libgfortran/libgfortran.h#L433-L440
# Note: This changed with GCC 9 (PR37577 and PR34640), when the step to Fortran 2008 was made. Before, there was no dedicated
# DType struct, but a bitmask instead.
struct GFortranArrayDescriptorDType
    elem_len::Csize_t
    version::Cint
    rank::Cchar
    type::Cchar
    attribute::Cshort
end

struct GFortranArrayDescriptorDimension
    stride::Cptrdiff_t
    lower_bound::Cptrdiff_t
    upper_bound::Cptrdiff_t
end

abstract type GFortranArrayDescriptor{T,N} end

struct GFortranArrayDescriptorRaw{T,N} <: GFortranArrayDescriptor{T,N}
    data::Ptr{T} # cannot be Ref, as this may be C_NULL
    offset::Cptrdiff_t
    dtype::GFortranArrayDescriptorDType
    span::Cptrdiff_t
    dim::NTuple{N,GFortranArrayDescriptorDimension}
end

struct GFortranArrayDescriptorBacked{T,N,A<:AbstractArray{T,N}} <: GFortranArrayDescriptor{T,N}
    data::Ptr{T} # cannot be Ref, as this may be C_NULL
    offset::Cptrdiff_t
    dtype::GFortranArrayDescriptorDType
    span::Cptrdiff_t
    dim::NTuple{N,GFortranArrayDescriptorDimension}
    # here the official struct ends
    dataref::A # so it does not harm to add this field, plus it gives us GC compliance
end

gfortran_array_element_type(::Type{<:Integer}) = Cchar(1)
gfortran_array_element_type(::Type{Bool}) = Cchar(2)
gfortran_array_element_type(::Type{<:Real}) = Cchar(3)
gfortran_array_element_type(::Type{<:Complex}) = Cchar(4)
gfortran_array_element_type(::Any) = Cchar(5)
gfortran_array_element_type(::Type{<:Union{<:AbstractChar,AbstractString}}) = Cchar(6)
# BT_CLASS, BT_PROCEDURE, BT_HOLLERITH, BT_VOID, BT_ASSUMED, BT_UNION, BT_BOZ not supported

Base.convert(G::Type{<:GFortranArrayDescriptor}, X::AbstractArray{T,N}) where {T,N} = convert(G, X)
function Base.convert(::Type{<:GFortranArrayDescriptorBacked{T,N}}, X::AbstractArray{T,N}) where {T,N}
    @assert(N ≤ 15)
    s = strides(X)
    a = axes(X)
    @assert(eltype(a) <: AbstractUnitRange)
    return GFortranArrayDescriptorBacked{T,N,typeof(X)}(
        isempty(X) ? Ptr{T}() : pointer(X),
        -sum(first(aᵢ) * sᵢ for (aᵢ, sᵢ) in zip(a, s); init=0),
        GFortranArrayDescriptorDType(
            sizeof(T),
            0, # version appears to be unused
            N,
            gfortran_array_element_type(T),
            0 # attributes seem to be unused - in particular, GFC_ARRAY_ASSUMED_SHAPE_CONT or GFC_ARRAY_ASSUMED_SHAPE appear to
              # be irrelevant
        ),
        sizeof(T), # span appears to be identical to the elem_len (at least in all considered cases)
        NTuple{N,GFortranArrayDescriptorDimension}(
            GFortranArrayDescriptorDimension(sᵢ, first(aᵢ), last(aᵢ)) for (sᵢ, aᵢ) in zip(s, a)
        ),
        X
    )
end
Base.convert(::Type{Array}, D::GFortranArrayDescriptor{T,N}) where {T,N} = convert(Array{T, N}, D)
function Base.convert(::Type{Array{T,N}}, D::GFortranArrayDescriptor{T,N}) where {T,N}
    @assert(sizeof(T) == D.dtype.elem_len && N == D.dtype.rank)
    # We only do an imperfect translation, as our arrays will always start at 1
    return unsafe_wrap(Array, D.data, NTuple{N,Int}(d.upper_bound - d.lower_bound + 1 for d in D.dim))
end
const GFortranVectorDescriptor{T} = GFortranArrayDescriptor{T,1}
const GFortranMatrixDescriptor{T} = GFortranArrayDescriptor{T,2}
const GFortranVectorDescriptorRaw{T} = GFortranArrayDescriptorRaw{T,1}
const GFortranVectorDescriptorBacked{T} = GFortranArrayDescriptorBacked{T,1}
const GFortranMatrixDescriptorRaw{T} = GFortranArrayDescriptorRaw{T,1}
const GFortranMatrixDescriptorBacked{T} = GFortranArrayDescriptorBacked{T,1}
const DoubleGVec = GFortranArrayDescriptorBacked{Cdouble,1,Vector{Float64}}
#endregion

@doc raw"""
    LANCELOT_simple(n, X, MY_FUN; MY_GRAD=missing, MY_HESS=missing, BL=nothing, BU=nothing, neq=0, nin=0, CX=nothing,
        Y=nothing, maxit=1000, gradtol=1e-5, feastol=1e-5, print_level=1)

# Purpose
A simple and somewhat NAIVE interface to LANCELOT B for solving the nonlinear optimization problem

``\min_x f(x)``

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
- `MY_FUN::Function`: a function for computing the objective function value for any `X`, whose interface has the default form
  `MY_FUN(X::AbstractVector{Float64})::Float64` where `X[1:n]` contains the values of the variables on input, and which returns
  a double precision scalar representing the value ``f(X)``.
- If the gradient of ``f`` can be computed, then the (optional) keyword argument `MY_GRAD` must be specified and given a
  function computing the gradient, whose interface must be of the form
  `MY_GRAD(G::AbstractVector{Float64}, X::AbstractVector{Float64})`, where `G` is a double precision vector of size `n` in
  which the function returns the value of the gradient of ``f`` at `X`.
- If, additionally, the second-derivative matrix of ``f`` at `X` can be computed, the (optional) keyword argument `MY_HESS`
  must be specified and given a function computing the Hessian, whose interface must be of the form
  `MY_HESS(H::PackedMatrix{Float64}, X::AbstractVector{Float64})`, where `H` is a double precision symmetrix matrix in packed
  storage format (upper triangular by column) of the Hessian of ``f`` at `X`.

In all cases, the best value of ``x`` found by LANCELOT B is returned to the user in the vector `X` and the associated
objective function value is the first return value.
The second return value reports the number of iterations performed by LANCELOT before exiting.
Finally, the last return value contains the exit status of the LANCELOT run, the value `0` indicating a successful run. Other
values indicate errors in the input or unsuccessful runs, and are detailed in the specsheet of LANCELOT B (with the exception
of the value 19, which reports a negative value for one or both input arguments `nin` and `neq`).

### Example
Let us consider the optimization problem

``\min_{x_1, x_2} f(x_1, x_2) = 100 ( x_2 - x_1^2 )^2 + ( 1 - x_1 )^2``,

which is the ever-famous Rosenbrock "banana" problem.
The most basic way to solve the problem (but NOT the most efficient) is, assuming the starting point `X = [-1.2, 1.]` known, to
perform the call `LANCELOT_simple(2, X, FUN)` where the user-provided function `FUN` is given by
```julia
function FUN(X)
    @inbounds return 100 * (X[2] - X[1]^2)^2 + (1 - X[1])^2
end
```

The solution is returned in 60 iterations with exit code `0`.

If we now wish to use first and second derivatives of the objective function, one should use the call
`LANCELOT_simple(2, X, FUN, MY_GRAD=GRAD!, MY_HESS=HESS!)` and provide the additional routines

```julia
function GRAD!(G, X)
    @inbounds G[1] = -400 * (X[2] - X[1]^2) * X[1] - 2 * (1 - X[1])
    @inbounds G[2] = 200 * (X[2] - X[1]^2)
end

function HESS!(H, X)
    @inbounds H[1, 1] = -400 * (X[2] - 3 * X[1]^2) + 2
    @inbounds H[1, 2] = -400 * X[1]
    @inbounds H[2, 2] = 200
end
```

Convergence is then obtained in 23 iterations. Note that using exact first-derivatives only is also possible: `MY_HESS` should
then be absent from the calling sequence and providing the subroutine HESS unnecessary.

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
LANCELOT_simple(2, X, FUN, MY_GRAD=GRAD!, MY_HESS=HESS!, BU=zeros(2))
```

in which case convergence is obtained in 6 iterations.

## Equality constrained problems
If, additionally, general equality constraints are also present in the problem, this must be declared by specifying the
following (optional) input argument:
- `neq::Integer`: the number of equality constraints.
In this case, the equality constraints are numbered from 1 to `neq` and the value of the `i`-th equality constraint must be
computed by a user-supplied routine of the form `FUN(X, i)` (`i = 1, ..., neq`) where the function now returns the value of the
`i`-th equality constraint evaluated at `X` if `i` is specified. (This extension of the unconstrained case can be implemented
by adding an optional argument `i` to the unconstrained version of `FUN` or by defining a three-parameter method on its own.)
If derivatives are available, then the `GRAD` and `HESS` subroutines must be adapted as well: `GRAD(G, X, i)` and
`HESS(H, X, i)` (`i = 1, ..., neq`) for computing the gradient and Hessian of the `i`-th constraint at `X`.
Note that, if the gradient of the objective function is available, so must be the gradients of the equality constraints. The
same level of derivative availability is assumed for all problem functions (objective and constraints). The final values of the
constraints and the values of their associated Lagrange multipliers is optionally returned to the user in the (optional) double
precision keyword arguments `CX` and `Y`, respectively (both being of size `neq`).

## Inequality constrained problems
If inequality constraints are present in the problem, their inclusion is similar to that of equality constraints. One then
needs to specify the (optional) input argument
- `nin::Integer`: the number of inequality constraints.
The inequality constraints are then numbered from `neq+1` to `neq+nin` and their values or that of their derivatives is again
computed by calling, for `i = 1, ..., nin`, `FUN(X, i)`, `GRAD(G, X, i)`, `HESS(H, X, i)`.
The inequality constraints are internally converted in equality ones by the addition of a slack variables, whose names are set
to 'Slack_i', where the character `i` in this string takes the integers values 1 to `nin`.)
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
LANCELOT_simple(2, X, FUN; MY_GRAD=GRAD!, MY_HESS=HESS!, BL=[0., -1e20], neq=1, nin=1, CX, Y)
```
(assuming we need `CX` and `Y`), and add methods for `FUN`, `GRAD!` and `HESS!` as follows
```julia
function FUN(X, i)
    if i == 1 # the equality constraint
        @inbounds return X[1] + 3X[2] - 3
    elseif i == 2 # the inequality constraint
        @inbounds return X[1]^2 + X[2]^2 - 4
    else
        @assert(false)
    end
end

function GRAD!(G, X, i)
    if i == 1 # equality constraint's gradient components
        @inbounds G[1] = 1
        @inbounds G[2] = 3
    elseif i == 2 # inequality constraint's gradient components
        @inbounds G[1] = 2X[1]
        @inbounds G[2] = 2X[2]
    else
        @assert(false)
    end
end

function HESS!(H, X, i)
    if i == 1 # equality constraint's Hessian
        fill!(H, 1.)
    elseif i == 2 # inequality constraint's Hessian
        @inbounds H[1] = 2
        @inbounds H[2] = 0
        @inbounds H[3] = 2
    else
        @assert(false)
    end
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
function LANCELOT_simple(n::Integer, X::AbstractVector{Cdouble}, MY_FUN::Function; MY_GRAD::Union{Function,Missing}=missing,
    MY_HESS::Union{Function,Missing}=missing, BL::Union{<:AbstractVector{Cdouble},Nothing}=nothing,
    BU::Union{<:AbstractVector{Cdouble},Nothing}=nothing, neq::Integer=0, nin::Integer=0,
    CX::Union{<:AbstractVector{Cdouble},Nothing}=nothing, Y::Union{<:AbstractVector{Cdouble},Nothing}=nothing,
    maxit::Integer=1000, gradtol::Real=1e-5, feastol::Real=1e-5, print_level::Integer=1)
    MY_FUN_F = @cfunction($((X, fx, i) -> begin
        unsafe_store!(fx, i == C_NULL ? MY_FUN(convert(Array, X)) : MY_FUN(convert(Array, X), Int(unsafe_load(i))))
        return
    end), Cvoid, (Ref{GFortranVectorDescriptorRaw{Cdouble}}, Ptr{Cdouble}, Ptr{Cint}))
    fx = Ref{Cdouble}(NaN)
    iters = Ref{Cint}(-1)
    exit_code = Ref{Cint}(-1)
    preserves::NTuple{5,Union{DoubleGVec,Nothing}} = (X, BL, BU, CX, Y)
    if !ismissing(MY_GRAD)
        MY_GRAD_F = @cfunction($((X, G, i) -> begin
            if i == C_NULL
                MY_GRAD(convert(Array, G), convert(Array, X))
            else
                MY_GRAD(convert(Array, G), convert(Array, X), Int(unsafe_load(i)))
            end
            return
        end), Cvoid, (Ref{GFortranVectorDescriptorRaw{Cdouble}}, Ref{GFortranVectorDescriptorRaw{Cdouble}}, Ptr{Cint}))
    end
    if !ismissing(MY_HESS)
        MY_HESS_F = let n=n
            @cfunction($((X, H, i) -> begin
                if i == C_NULL
                    MY_HESS(PackedMatrix(n, convert(Array, H)), convert(Array, X))
                else
                    MY_HESS(PackedMatrix(n, convert(Array, H)), convert(Array, X), unsafe_load(i))
                end
                return
            end), Cvoid, (Ref{GFortranVectorDescriptorRaw{Cdouble}}, Ref{GFortranVectorDescriptorRaw{Cdouble}}, Ptr{Cint}))
        end
    end
    GC.@preserve preserves (@ccall libgalahad_double.__lancelot_simple_double_MOD_lancelot_simple(
        n::Ref{Cint},                                   # INTEGER ( KIND = ip_ ), INTENT( IN ) :: n
        Ref(preserves[1])::Ref{DoubleGVec},             # REAL ( KIND = rp_ ), INTENT( INOUT ) :: X( : )
        MY_FUN_F::Ptr{Cvoid},
        fx::Ref{Cdouble},                               # REAL ( KIND = rp_ ), INTENT( OUT ) :: fx
        exit_code::Ref{Cint},                           # INTEGER ( KIND = ip_ ), INTENT( OUT ) :: exit_code
        (ismissing(MY_GRAD) ? C_NULL : MY_GRAD_F)::Ptr{Cvoid},
        (ismissing(MY_HESS) ? C_NULL : MY_HESS_F)::Ptr{Cvoid},
        (isnothing(BL) ? C_NULL : Ref(preserves[2]))::Ptr{DoubleGVec},
                                                        # REAL ( KIND = rp_ ), OPTIONAL :: BL ( : )
        (isnothing(BU) ? C_NULL : Ref(preserves[3]))::Ptr{DoubleGVec},
                                                        # REAL ( KIND = rp_ ), OPTIONAL :: BU ( : )
        C_NULL::Ptr{Cvoid},                             # CHARACTER ( LEN = 10 ), OPTIONAL :: VNAMES( : )
        C_NULL::Ptr{Cvoid},                             # CHARACTER ( LEN = 10 ), OPTIONAL :: CNAMES( : )
        neq::Ref{Cint},                                 # INTEGER ( KIND = ip_ ), OPTIONAL :: neq
                                                        # ^ optional is equivalent to 0, so we just require it
        nin::Ref{Cint},                                 # INTEGER ( KIND = ip_ ), OPTIONAL :: nin
                                                        # ^ optional is equivalent to 0, so we just require it
        (isnothing(CX) ? C_NULL : Ref(preserves[4]))::Ptr{DoubleGVec},
                                                        # REAL ( KIND = rp_ ), OPTIONAL :: CX ( : )
        (isnothing(Y) ? C_NULL : Ref(preserves[5]))::Ptr{DoubleGVec},
                                                        # REAL ( KIND = rp_ ), OPTIONAL :: Y ( : )
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

end