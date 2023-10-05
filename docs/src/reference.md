# Reference

```@meta
CurrentModule = PolynomialOptimization
```

## Problem definition
```@docs
PolyOptProblem
poly_problem
newton_polytope
EqualityMethod
variables
nvariables
degree
length
isreal
```

## Optimization
```@docs
sparse_optimize
poly_optimize
```

## Working with problem solutions
All functions listed here requires that the problem be optimized before.
```@docs
poly_solutions
poly_solutions_heuristic
poly_all_solutions
poly_solution_badness
optimality_certificate
moment_matrix
last_moments
last_objective
```

## Sparsity
```@docs
AbstractSparsity
SparsityNone
SparsityCorrelative
SparsityTerm
SparsityTermBlock
SparsityTermCliques
SparsityCorrelativeTerm
sparse_iterate!
sparse_groupings
sparse_problem
TermMode
```

## Chordal graphs
```@docs
chordal_completion!
chordal_cliques!
```

## SpecBM
While the solver was implemented for the purpose of being used within `PolynomialOptimization`, it also works as a standalone
routine (and could in principle be a separate package). SpecBM is a
[spectral bundle algorithm for primal semidefinite programs](https://doi.org/10.48550/arXiv.2307.07651) and is based on the
assumption that the optimal dual solution has low rank; indeed, in polynomial optimizations, if there is an optimal point for
the problem that can be encoded in the chosen basis, then this automatically gives rise to a rank-one semidefinite encoding of
this point.
The implementation also allows for free variables and multiple semidefinite constraints and contains further improvements
compared to the [reference implementation](https://github.com/soc-ucsd/specBM). It requires either Hypatia or a very recent
version of Mosek as subsolvers.
```@docs
specbm_primal
SpecBMResult
```

## FastVector
To improve the speed in some implementation details, `PolynomialOptimization` provides a "fast" vector type. This is basically
just a wrapper around the stdlib `Vector`; however, it actually takes proper advantage of sizehints. The fact that Julia does
this badly has been known for quite some time ([#24909](https://github.com/JuliaLang/julia/issues/24909)), but the default
implementation has not changed. Our own `FastVec` is a bit more specific than the
[PushVector](https://github.com/tpapp/PushVectors.jl), but also allows for more aggressive optimization.
```@docs
FastVec
Base.sizehint!(::FastVec, ::Integer)
prepare_push!
Base.push!(::FastVec{V}, el) where {V}
unsafe_push!
Base.append!(::FastVec{V}, ::AbstractVector) where {V}
unsafe_append!
Base.prepend!(::FastVec{V}, ::AbstractVector) where {V}
unsafe_prepend!
Base.similar(::FastVec{V}) where {V}
finish!
```

## Lancelot
`PolynomialOptimization` provides an interface to LANCELOT. This is only temporary, as the GALAHAD developers are developing
Julia interfaces themselves (in fact, this is already there - but LANCELOT is the only part that not even has a C interface).
```@docs
LANCELOT_simple
```

## PackedMatrices
The SpecBM solver works with matrices in packed storage format. There are a lot of LAPACK routines that natively support this
format, which offers superior storage (at the cost of worse runtime and precision); however, Julia does not natively provide
a convenient way to work with packed matrices or even just wrapper functions. The submodule `PackedMatrices` provides this
functionality (though not every LAPACK function was exported; the development was guided by what was needed for
`PolynomialOptimization`).

```@meta
CurrentModule = PolynomialOptimization.PackedMatrices
```
### The PackedMatrix
```@docs
PackedMatrix
PackedMatrix(::Integer, ::AbstractVector{R}, ::Symbol) where {R}
PackedMatrix{R}(::UndefInitializer, ::Integer, ::Symbol) where {R}
PackedMatrix(::Symmetric{R,<:AbstractMatrix{R}}) where {R}
PackedMatrixUpper
PackedMatrixLower
PackedMatrixUnscaled
PackedMatrixScaled
packedsize
packed_isupper
packed_islower
packed_isscaled
PackedDiagonalIterator
rmul_diags!
rmul_offdiags!
lmul_diags!
lmul_offdiags!
packed_scale!
packed_unscale!
getindex
setindex!
vec(::PackedMatrix)
Matrix{R}(::PackedMatrixUnscaled{R}) where {R}
```
### LAPACK wrappers
Note that unless noted, these functions are only wrappers for real-valued LAPACK functions. No other data types than native single and double precision are therefore supported.
```@docs
spev!
spevd!
spevx!
pptrf!
spmv!
spr!
tpttr!
trttp!
gemmt!
```
### Extensions in LinearAlgebra
```@docs
axpy!
mul!
spr!(::Any, ::AbstractVector, ::PackedMatrix)
dot
eigen!
eigvals
eigvals!
eigmin!
eigmax!
cholesky!
isposdef
```