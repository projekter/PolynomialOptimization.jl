# Introduction

`PolynomialOptimization` is a Julia package that allows to easily optimize large-scale polynomial optimization problems.
It builds on `MultivariatePolynomials` to provide a concise interface for the specification of the problem and allows to
directly control the problem's sparsity using correlative sparsity, (chordal) term sparsity, and a combination of both. It also
supports complex-valued problems and positive semidefinite constraints, and allows to extract solutions even for sparse
problems.
It provides a solver interface specifically designed for the optimization problems arising in polynomial optimization. This
interface makes it easy to implement new solvers. On purpose, `MathOptInterface`/`JuMP` is not employed; despite `JuMP` being
very performant for a modelling framework, it introduces a significant overhead that is omitted in this way.
The following solvers are supported:
- [Clarabel](https://github.com/oxfordcontrol/Clarabel.jl)
- [COPT](https://www.shanshu.ai/copt)
- [Hypatia](https://github.com/jump-dev/Hypatia.jl)
- [Loraine](https://github.com/kocvara/Loraine.jl) (own implementation)
- [LoRADS](https://github.com/projekter/LoRADS)
- [LANCELOT](https://github.com/ralna/GALAHAD) (not based on relaxations)
- [Mosek](https://www.mosek.com/)
- [ProxSDP](https://github.com/mariohsouto/ProxSDP.jl)
- [SCS](https://github.com/cvxgrp/scs)
- [Sketchy CGAL](https://doi.org/10.1137/19M1305045) (own implementation)
- [SpecBM Primal](https://doi.org/10.48550/arXiv.2307.07651) (own implementation)

## Overview
```@contents
Depth=3
```