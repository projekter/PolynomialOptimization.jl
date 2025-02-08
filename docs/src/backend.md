```@meta
CurrentModule = PolynomialOptimization.Solver
CollapsedDocStrings = true
```

# Backend
`PolynomialOptimization` mostly uses external solvers that have Julia bindings or are implemented in Julia, but it also
provides own solver implementations particularly for the purpose of polynomial optimization. This page is only relevant if you
intend to implement an interface between `PolynomialOptimization` and a new solver or if you want to provide missing
functionality for an existing solver. It is of no relevance if you only want to use the existing solvers.

!!! warning
    The package does not introduce any hard dependencies on the external solvers. Therefore, you may or may not decide to
    install them on your system. Instead, the solvers are addressed as weak dependencies. This means that *you have to load the
    solver Julia package manually* before you are able to use it for optimization.
`PolynomialOptimization` also provides an interface that can be implemented if another solver should be supported. This
consists of just a few methods for the various functions.

!!! info "MathOptInterface"
    Why do we re-invent the wheel and don't just use `MathOptInterface`? This would immediately bring compatibility with a huge
    number of solvers that are available under Julia.

    This would certainly be possible. It might even be a not-too-bad idea to provide an automatic fallback to
    `MathOptInterface` in the future for solvers that are not supported natively.

    However, `MathOptInterface` is enormous in its feature range. Considering that, it is amazingly fast. But in
    `PolynomialOptimization`, only a very small subset of all the features is required, and a lot of additional assumptions on
    how the problem is constructed can be made. This allows to use the solver's API in a most efficient way, which typically
    is not the way in which implementations of `MathOptInterface` would address the solver.

    Additionally, in some situations, a lot of very similar sub-problems need to be solved; in these cases,
    `PolynomialOptimization`'s own interface allows to keep the optimizer task alive and just do the tiny modification instead
    of setting things up from the start again - which for `MathOptInterface` only works if a solver is implemented in a
    particular way.

    This focus on efficiency is very important as (relevant) polynomial optimization problems are huge. It is just not possible
    to waste time and memory in bookkeeping that is not really needed.

## [`poly_optimize`](@ref)
The optimization of polynomial problems requires a solver that understands linear and semidefinite constraints. All functions
in this section are defined (and exported) in the submodule `PolynomialOptimization.Solver`.

### Solver interface
In general, a solver implementation can do whatever it wants; it just needs to implement the [`poly_optimize`](@ref) method
with the appropriate `Val`-wrapped solver method as its first parameter. However, it is very helpful to just do some basic
setup such as creating the solver object in this function and delegate all the work of setting up the actual problem to
[`moment_setup!`](@ref) or [`sos_setup!`](@ref).
In order to do so, a solver implementation should create a new type that contains all the relevant data during setup of the
problem. Usually, a solver falls in one of three categories:
- Problem data has to be supplied in matrix/vector form; in this case, the new type should be a descendant of
  [`AbstractSparseMatrixSolver`](@ref). Usually, open-source solvers fall in this category.
- Problem data is constructed incrementally via various calls to API functions of the solver, which does not provide access to
  its internals. In this case, the new type should be a descendant of [`AbstractAPISolver`](@ref). Usually, commercial solvers
  fall in this category.
However, it is not required that the type is in fact a subtype of either of those; the most general possible supertype is
[`AbstractSolver`](@ref), which does not make any assumptions or provide any but the most skeleton fallback implementations and
a default for [`mindex`](@ref).
```@docs
AbstractSolver
mindex
```
Every implementation of [`poly_optimize`](@ref) should return a tuple that contains some internal state of the solver as well
as the optimal value and the status of the solver. A method for [`issuccess`](@ref issuccess(::Val, ::Any)) should then
translate this status into a simple boolean, where deciding on ambiguities (near success) is up to the solver implementation.
```@docs; canonical=false
poly_optimize(::Val, ::AbstractRelaxation, ::RelaxationGroupings)
issuccess(::Val, ::Any)
```
Once a solver has been implemented, it should add its solver symbol to the vector `solver_methods`, which enables this solver
to be chosen automatically. Apart from the exact specification `:<solvername>Moment` or `:<solvername>SOS`, a short form
`:<solvername>` that chooses the recommended method should also be implemented. For this, the [`@solver_alias`](@ref) macro can
be used. When details on the solution data a requested, the [`extract_moments`](@ref), [`extract_sos`](@ref), or
[`extract_info`](@ref) function is called, where at least the former two have to be implemented for each solver:
```@docs
extract_moments
extract_sos
extract_sos_prepare
extract_info
```
In order to relate aspects of the problem with data in the solver, the cones that are added are counted. This works
automatically, keeping a separate counter for every type of cone and counting vector-valued cones (which are most) with their
actual length. This behavior can be customized:
```@docs
@counter_alias
@counter_atomic
addtocounter!
Counters
```
Using this information, an additional implementation may be provided for a faster re-optimization of the same problem:
```@docs
poly_optimize(::Val, ::Any, ::AbstractRelaxation, ::RelaxationGroupings)
```

While this page documents in detail how a new solver can be implemented, the explanation is far more extensive than an actual
implementation. In order to implement a new solver, it is therefore recommended to first determine the category in which it
falls, then copy and modify an appropriate existing implementation.

#### [`AbstractSparseMatrixSolver`](@ref)
This solver type accumulates data in a COO matrix-like format. The callback functions can use
[`append!`](@ref append!(::SparseMatrixCOO{I,K,V,Offset}, ::IndvalsIterator{K,V}) where {I<:Integer,K<:Integer,V<:Real,Offset})
to quickly add given data to the temporary storage.
However, the format is not yet suitable for passing data to the solver, as all monomials are densely indexed. Therefore, in a
postprocessing step, the COO matrices have to be converted to CSC matrices with continuous monomial indices using
[`coo_to_csc!`](@ref).
After the optimization is done, the optimal moment vector (the decision variables for moment optimization, the
constraint duals for SOS optimization) can be constructed using [`MomentVector`](@ref MomentVector(::AbstractRelaxation, ::Vector{V}, ::SparseMatrixCOO{<:Integer,K,V,Offset}, ::SparseMatrixCOO{<:Integer,K,V,Offset}...) where {K<:Integer,V<:Real,Offset}).
```@docs
AbstractSparseMatrixSolver
SparseMatrixCOO
append!(::SparseMatrixCOO{I,K,V,Offset}, ::Indvals{K,V}) where {I<:Integer,K<:Integer,V<:Real,Offset}
append!(::SparseMatrixCOO{I,K,V,Offset}, ::IndvalsIterator{K,V}) where {I<:Integer,K<:Integer,V<:Real,Offset}
coo_to_csc!(::Union{SparseMatrixCOO{I,K,V,Offset},<:Tuple{AbstractVector{K},AbstractVector{V}}}...) where {I<:Integer,K<:Integer,V<:Real,Offset}
MomentVector(::AbstractRelaxation{<:Problem{<:SimplePolynomial{<:Any,Nr,Nc}}}, ::Vector{V}, ::K, ::SparseMatrixCOO{<:Integer,K,V,Offset}, ::SparseMatrixCOO{<:Integer,K,V,Offset}...) where {Nr,Nc,K<:Integer,V<:Real,Offset}
```

#### [`AbstractAPISolver`](@ref)
This solver type accumulates the data continuously by the means of API calls provided by the solver. The monomials therefore
have to have contiguously defined indices from the beginning. This is done by internal bookkeeping; however, it requires the
implementation of an additional callback function [`append!`](@ref append!(::AbstractAPISolver{K}, ::K) where {K<:Integer}) to
add new monomials to the solver.
After the optimization is done, the optimal moment vector (the decision variables for moment optimization, the
constraint duals for SOS optimization) can be constructed using [`MomentVector`](@ref MomentVector(::AbstractRelaxation, ::Vector{V}, ::AbstractAPISolver{K}) where {K<:Integer,V<:Real}).
```@docs
AbstractAPISolver
append!(::AbstractAPISolver{K}, ::K) where {K<:Integer}
MomentVector(::AbstractRelaxation{<:Problem{<:SimplePolynomial{<:Any,Nr,Nc}}}, ::Vector{V}, ::AbstractAPISolver{K}) where {Nr,Nc,K<:Integer,V<:Real}
```

#### Defining solver capabilities
There are some functions that should be implemented to tell `PolynomialOptimization` what kind of data the solver expects and
which cones are supported; these should return constants.
```@docs
supports_rotated_quadratic
supports_quadratic
supports_psd_complex
supports_dd
supports_dd_complex
supports_lnorm
supports_lnorm_complex
supports_sdd
supports_sdd_complex
PSDIndextype
PSDIndextypeMatrixCartesian
PSDIndextypeVector
PSDIndextypeCOOVectorized
psd_indextype
negate_fix
negate_free
prepend_fix
prepend_free
```

#### Working with data from the interface
The interface functions that need to be implemented will get their data in the form of index-value pairs or a specific
structure related to the desired matrix format.
```@docs
Indvals
PSDMatrixCartesian
IndvalsIterator
```

#### Interface for the moment optimization
The custom implementation of [`poly_optimize`](@ref) first has to set up all the necessary initial data; then, a call to
[`moment_setup!`](@ref) is sufficient to trigger the process.
```@docs
moment_setup!
moment_add_matrix!
moment_add_equality!
```

For this to work, the following methods (or a subset as previously indicated) must be implemented.
```@docs
add_constr_nonnegative!
add_constr_rotated_quadratic!
add_constr_quadratic!
add_constr_psd!
add_constr_psd_complex!
add_constr_dddual!
add_constr_dddual_complex!
add_constr_linf!
add_constr_linf_complex!
add_constr_sdddual!
add_constr_sdddual_complex!
add_constr_fix_prepare!
add_constr_fix!
add_constr_fix_finalize!
fix_objective!
add_var_slack!
```

#### Interface for the SOS optimization
This is in complete analogy to the moment case; the entry point is now [`sos_setup!`](@ref).
```@docs
sos_setup!
sos_add_matrix!
sos_add_equality!
```

The following methods (or a subset as previously indicated) must be implemented.
```@docs
add_var_nonnegative!(::AbstractSolver{T,V}, ::Indvals{T,V}) where {T,V}
add_var_nonnegative!(::AbstractSolver{T,V}, ::IndvalsIterator{T,V}) where {T,V}
add_var_rotated_quadratic!
add_var_quadratic!
add_var_psd!(::AbstractSolver{T,V}, ::Int, ::PSDMatrixCartesian{T,V}) where {T,V}
add_var_psd!(::AbstractSolver{T,V}, ::Int, ::IndvalsIterator{T,V}) where {T,V}
add_var_psd_complex!
add_var_dd!
add_var_dd_complex!
add_var_l1!
add_var_l1_complex!
add_var_sdd!
add_var_sdd_complex!
add_var_free_prepare!
add_var_free!
add_var_free_finalize!
fix_constraints!
add_constr_slack!
```

#### Interface for the moment optimization in primal form
A very particular case is if the moment optimization should be done using semidefinite variables that the solver allows to
define, but only as a single variable (not as a cone in which to put variables); this is the primal form, most suitable for SOS
optimizations. However, there can be a good reason to use the primal form for moment optimizations instead: namely, if the
solver can exploit a low-rank assumptions on this matrix. In this case, `poly_optimize` should call
[`primal_moment_setup!`](@ref) instead:
```@docs
primal_moment_setup!
```

The following methods must be implemented:
```@docs
add_var_nonnegative!(::AbstractSolver{<:Integer,V}, ::Int, ::Int, ::SparseMatrixCOO{I,I,V}, ::Tuple{FastVec{I},FastVec{V}}) where {I,V}
add_var_psd!(::AbstractSolver{<:Integer,V}, ::Int, ::I, ::SparseMatrixCOO{I,I,V}, ::Union{Nothing,Tuple{FastVec{I},FastVec{V}}}) where {I,V}
add_var_psd!(::AbstractSolver{<:Integer,V}, ::Int, ::I, ::SparseMatrixCOO{I,I,V}, ::Union{Nothing,Tuple{Tuple{FastVec{I},FastVec{I}},FastVec{V}}}) where {I,V}
add_var_psd!(::AbstractSolver{<:Integer,V}, ::Int, ::I, ::Tuple{FastVec{I},Tuple{FastVec{I},FastVec{I}},FastVec{V}}, ::Union{Nothing,Tuple{FastVec{I},FastVec{V}}}) where {I,V}
add_var_psd!(::AbstractSolver{<:Integer,V}, ::Int, ::I, ::Tuple{FastVec{I},Tuple{FastVec{I},FastVec{I}},FastVec{V}}, ::Union{Nothing,Tuple{Tuple{FastVec{I},FastVec{I}},FastVec{V}}}) where {I,V}
fix_constraints!(::AbstractSolver{<:Integer,V}, ::Int, ::Indvals{<:Integer,V}) where {V}
objective_indextype
```
Note that these methods work with the COO representation, which can be quickly converted to either CSR or CSC using
[`coo_to_csr!`](@ref) and [`coo_to_csc!`](@ref coo_to_csc!(::Integer, ::SparseMatrixCOO{I,I,V,offset}) where {I,V,offset}),
which respects the offset desired by the solver. Only this interface allows to set [`psd_indextype`](@ref) to a
[`PSDIndextypeCOOVectorized`](@ref); but [`PSDIndextypeVector`](@ref) is now forbidden.

#### Helper functions
The solver module exports a number of helper functions which may be of use in implementations:
```@docs
@solver_alias
monomial_count
trisize
count_uniques
coo_to_csc!(::Integer, ::SparseMatrixCOO{I,I,V,offset}) where {I,V,offset}
coo_to_csr!
```
Additionally, `Solver` reexports a number of useful types and functions for developing the interface (see
`src/optimization/solver/Solver.jl`); therefore, usually only the `Solver` submodule itself has to be used and not
`PolynomialOptimization` itself. However, note that it is highly recommended to say
`using PolynomialOptimization: @assert, @inbounds` in the solver implementation; this will replace the Base implementations of
the macros by ones that, depending on a debugging constant, enable or disable the desired functionality.

```@meta
CurrentModule = PolynomialOptimization.Newton
```
## [`Newton.halfpolytope`](@ref)
Finding the Newton halfpolytope requires a linear solver that supports problem modification for quick reoptimization.  All the
functions here are defined in the submodule `PolynomialOptimization.Newton` and they are not exported.

### [List of supported solvers](@id solvers_newton)
The following list contains all the solvers and the required packages that provide access to the solver. The name of the solver
is identical with the solver method (as a Symbol).

|  Solver |                            Package                          |   License  |  Speed  | Accuracy | Memory  |
| ------: | :---------------------------------------------------------: | :--------: | :-----: | :------: | :-----: |
| COPT    | [COPT.jl](https://github.com/COPT-Public/COPT.jl/tree/main) | commercial | 👍👍👍 | 👍👍👍  | 👍👍👍 |
| Mosek   | [Mosek.jl](https://github.com/MOSEK/Mosek.jl)               | commercial | 👍👍👍 | 👍👍👍  | 👍👍👍 |

### Solver interface
The following functions need to be implemented so that a solver is available via [`Newton.halfpolytope`](@ref). The
preprocessing functions can be omitted if no preprocessing capabilities should be provided.
```@docs
preproc
prepare
alloc_global
alloc_local
clonetask
work
```

Once a solver has been implemented, it should add its solver symbol to the vector `Newton.newton_methods`,
which enables this solver to be chosen automatically.

## Automatic tightening
```@meta
CurrentModule = PolynomialOptimization
```
Automatic tightening of a polynomial optimization problem requires a linear solver that finds a solution to a system of linear
equations that minimizes the ℓ₁ norm (better yet, the ℓ₀-norm, if you can implement this). The solver is only called if the
number of rows is smaller than the number of columns; else, the solution is calculated using SPQR's direct solver.

### [List of supported solvers](@id solvers_tighten)
The following list contains all the solvers and the required packages that provide access to the solver. The name of the solver
is identical with the solver method (as a Symbol).

|  Solver |                            Package                          |   License  |  Speed  | Accuracy | Memory  |
| ------: | :---------------------------------------------------------: | :--------: | :-----: | :------: | :-----: |
| COPT    | [COPT.jl](https://github.com/COPT-Public/COPT.jl/tree/main) | commercial | 👍👍👍 | 👍👍👍  | 👍👍👍 |
| Mosek   | [Mosek.jl](https://github.com/MOSEK/Mosek.jl)               | commercial | 👍👍👍 | 👍👍👍  | 👍👍👍 |

### Solver interface
The following function needs to be implemented so that a solver is available via automatic tightening.
```@docs
tighten_minimize_l1
```

Once a solver has been implemented, it should add its solver symbol to the vector `PolynomialOptimization.tightening_methods`,
which enables this solver to be chosen automatically.