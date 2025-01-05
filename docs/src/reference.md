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
Newton
Custom
```

### Relaxations based on individual sparsity
```@docs
AbstractRelaxationSparse
SparsityCorrelative
SparsityTerm
SparsityTermBlock
SparsityTermChordal
SparsityCorrelativeTerm
TermMode
CliqueMerged
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