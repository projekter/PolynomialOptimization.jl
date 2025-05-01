# Reference

```@meta
CurrentModule = PolynomialOptimization
```

# Optimization reference
This reference page lists all functions that are relevant for polynomial optimization.

## Problem definition
```@docs
Problem
poly_problem(::P) where {P<:AbstractPolynomialLike}
variables
nvariables
isreal
```

## Relaxations
Types and functions related to relaxations of polynomial optimization problems are found in the submodule `Relaxation`. The
types in this module are mostly not exported, so that a qualified name is required.
```@meta
CurrentModule = PolynomialOptimization.Relaxation
```
```@docs
AbstractRelaxation
poly_problem(::AbstractRelaxation)
basis
MultivariatePolynomials.degree(::AbstractRelaxation)
groupings
iterate!(::AbstractRelaxation)
Core.Type(::Problem, ::Tuple{Vararg{Any}})
RelaxationGroupings
```

### Relaxations based on a global basis
```@docs
AbstractRelaxationBasis
Dense
Custom
```

### Relaxations based on individual sparsity
```@docs
AbstractRelaxationSparse
```
Two exact sparsity methods are available: a highly-optimized Newton polytope algorithm which will often reduce the size of the
basis for the moment matrix quite well, but cannot help for the constraints; and a less efficient facial reduction which is the
best currently known exact sparsity method (also being able to reduce constraints).
A reverse-sparsity method complements the set of exact methods.
```@docs
Newton
FacialReduction
CliqueMerged
```
Several inexact methods are available which will potentially worsen the quality of the bounds.
```@docs
SparsityCorrelative
SparsityTerm
SparsityTermBlock
SparsityTermChordal
SparsityCorrelativeTerm
TermMode
iterate!(::SparsityTerm)
```

## Optimization and problem solutions
```@meta
CurrentModule = PolynomialOptimization
```
```@docs
poly_optimize(::Val, ::AbstractRelaxation)
poly_optimize(::Val, ::Problem, ::Vararg{Any})
poly_optimize(::Result)
Solver.RepresentationMethod
RepresentationPSD
RepresentationSDD
RepresentationDD
RepresentationIAs
Result
issuccess(::Result)
poly_problem(::Result)
optimality_certificate
poly_all_solutions
poly_solutions
poly_solution_badness
moment_matrix
MomentVector
MomentAssociation
SOSCertificate
sos_matrix
IterateRepresentation
```

## Newton polytope construction (manually)
Note that using these functions is usually not necessary; construct a [`Newton`](@ref Relaxation.Newton) relaxation instead.
```@docs
Newton.halfpolytope
Newton.halfpolytope_from_file
```

## Facial reduction (manually)
Note that using this function is usually not necessary; construct a [`FacialReduction`](@ref Relaxation.FacialReduction)
relaxation instead.
```@docs
FacialReduction.facial_reduction!!
```