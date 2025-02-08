[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://projekter.github.io/PolynomialOptimization.jl/dev)
[![Build status (Github Actions)](https://github.com/projekter/PolynomialOptimization.jl/workflows/CI/badge.svg)](https://github.com/projekter/PolynomialOptimization.jl/actions)
[![codecov.io](http://codecov.io/github/projekter/PolynomialOptimization.jl/coverage.svg?branch=main)](http://codecov.io/github/projekter/PolynomialOptimization.jl?branch=main)

# PolynomialOptimization.jl

`PolynomialOptimization` is a Julia package that allows to easily optimize large-scale polynomial optimization problems.
It builds on `MultivariatePolynomials` to provide a concise interface for the specification of the problem and allows to apply
many kinds of sparsity methods. It also fully supports complex-valued problems and positive semidefinite constraints and allows
to extract solutions even for sparse problems.
It _directly_ interfaces the solvers
[Clarabel](https://github.com/oxfordcontrol/Clarabel.jl),
[COPT](https://www.shanshu.ai/copt),
[Hypatia](https://github.com/jump-dev/Hypatia.jl),
[Mosek](https://www.mosek.com/),
[ProxSDP](https://github.com/mariohsouto/ProxSDP.jl), and
[SCS](https://github.com/cvxgrp/scs) without using `JuMP`, avoiding this bottleneck so that indeed almost all the time is spent
in the solver, not in the problem formulation. It also features new research-level solvers: an own implementation of the
[primal spectral bundling](https://arxiv.org/abs/2307.07651v1) method;
an efficiency-oriented refactoring of the low-rank solver [Loraine](https://github.com/kocvara/Loraine.jl),
and an interface to the experimental [LoRADS](https://github.com/projekter/LoRADS) solver.