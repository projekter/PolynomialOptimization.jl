# Introduction

`PolynomialOptimization` is a Julia package that allows to easily optimize large-scale polynomial optimization problems.
It builds on `MultivariatePolynomials` to provide a concise interface for the specification of the problem and allows to
directly control the problem's sparsity using correlative sparsity, (chordal) term sparsity, and a combination of both. It also
supports complex-valued problems and positive semidefinite constraints, and allows to extract solutions even for sparse
problems.
It _directly_ interfaces the solvers [Mosek](https://www.mosek.com/) (version 10 strongly preferred, less features available
with version 9), [COSMO](https://oxfordcontrol.github.io/COSMO.jl/stable/),
and [Hypatia](https://github.com/chriscoey/Hypatia.jl), without using `JuMP`. Despite `JuMP` being very performant for a
modelling framework, it introduces a significant overhead that is omitted in this way.

## About this branch
This branch contains an implementation of the [SketchyCGAL](https://doi.org/10.1137/19M1305045) solver in its original form.
It also contains an extension with multiple semidefinite matrices (not thouroughly tested). However, as experiments did not
show much promise for solving polynomial optimization problems, no integration with `PolynomialOptimization` is currently
provided.

## Compatibility notice
Currently, the required complex-valued functionality depends on a not-yet-merged request in
[`MultivariatePolynomials`](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl/pull/213). For this reason, the package
is not available on the registry at this moment. You have to manually install the necessary parts:
```
julia> ]
pkg> add https://github.com/projekter/MultivariatePolynomials.jl
pkg> add https://github.com/projekter/DynamicPolynomials.jl.git
pkg> add https://github.com/projekter/PolynomialOptimization.jl.git
```

Note that with regard to solvers, by far the most superior performance is obtained using Mosek. All features are available only
with the latest version.