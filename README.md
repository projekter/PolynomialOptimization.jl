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

This branch provides an interface to the [LANCELOT sovler](https://github.com/ralna/GALAHAD). LANCELOT is a sophisticated
non-convex solver which can now be invoked using the solver `:LANCELOT`. Note that contrary to all other methods, this will
*not* result in an optimization that is based on convex relaxations; no global guarantees are made. Futher notes:

- We call `LANCELOT_simple`, which does not make use of all the potential sophistication LANCELOT provides as regards the use
  of variables.
- You will have to download and compile the GALAHAD library first using the GFortran compiler (at least version 9). Since
  LANCELOT currently only has a Fortran interface, we need to rely on compiler internals to call it. Memory corruption will
  definitely occur if the library was compiled with a version of GFortran smaller than 9 or e.g., Intel Fortran. Both use
  different formats for the array descriptors.
- You will have to set up an environment variable `JULIA_GALAHAD_LIBRARY_PATH` containing the path to the `libgalahad_double_d`
  shared library whenever you import `PolynomialOptimization`. If this variable is not present, you'll be asked for the file
  path every you import the package.

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