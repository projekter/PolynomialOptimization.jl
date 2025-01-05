[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://projekter.github.io/PolynomialOptimization.jl/dev)
[![Build status (Github Actions)](https://github.com/projekter/PolynomialOptimization.jl/workflows/CI/badge.svg)](https://github.com/projekter/PolynomialOptimization.jl/actions)
[![codecov.io](http://codecov.io/github/projekter/PolynomialOptimization.jl/coverage.svg?branch=main)](http://codecov.io/github/projekter/PolynomialOptimization.jl?branch=main)

# PolynomialOptimization.jl

`PolynomialOptimization` is a Julia package that allows to easily optimize large-scale polynomial optimization problems.
It builds on `MultivariatePolynomials` to provide a concise interface for the specification of the problem and allows to apply
many kinds of sparsity methods. It also fully supports complex-valued problems and positive semidefinite constraints and allows
to extract solutions even for sparse problems.
It _directly_ interfaces the solvers [Clarabel](https://github.com/oxfordcontrol/Clarabel.jl),
[COPT](https://www.shanshu.ai/copt), [Hypatia](https://github.com/jump-dev/Hypatia.jl), [Mosek](https://www.mosek.com/), and
[SCS](https://github.com/cvxgrp/scs) without using `JuMP`, avoiding this bottleneck so that indeed almost all the time is spent
in the solver, not in the problem formulation.

## About this branch
This branch implements the [STRIDE solver](https://doi.org/10.1007/s10107-022-01912-6) and provides the `:STRIDE` algorithm
as an interface to `PolynomialOptimization`. Note that STRIDE requires local optimizer. If none is specified, by default
LANCELOT is invoked (see also the notes on the `lancelot` branch).
We also do not perform the BFGS step that is detailed in the solver. Experiments (including their own code) have not shown any
significant advantage of this step (see also the outdated `stride-experiments` branch).
In general, STRIDE does not appear to perform particularly well.