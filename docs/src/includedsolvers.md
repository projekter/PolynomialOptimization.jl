```@meta
CurrentModule = PolynomialOptimization
```
# [Supported solvers](@id solvers_poly_optimize)
The following list contains all the solvers and the required packages that provide access to the solver.
A solver of name `X` will always provide at least one of the two methods `:XSOS` or `:XMoment`. The former models the
sum-of-squares formulation of the problem - each monomial corresponds to one constraint, and the solution comes from the dual
vector. The latter models the Lasserre moment hierarchy of the problem - each monomial corresponds to a primal variable.
Which method (or both) is offered by the solver depends on the particular interface that is supported by the solver and which
form fits more naturally to this interface. Every solver will also provide an alias `:X` method that defaults to the
recommended method for this solver.
The given performance indicators are merely based on experience and may not be accurate for your particular problem. In
particular the maximally recommended basis size depends heavily on the structure of the final problem, which can easily put the
number up or down by 100 or more. All solvers may expose options that can influence the runtime behavior.

|  Solver     |                            Package                          |   License  | Methods       | Speed    | Accuracy | Memory  | max. recomm. basis size |
| ------:     | :---------------------------------------------------------: | :--------: | :---------:   | :-----:  | :------: | :-----: | :---------------------- |
| Clarabel    | [Clarabel.jl](https://github.com/oxfordcontrol/Clarabel.jl) | Apache     | moment        | 👍👍👍 | 👍👍👍 | 👍👍    | ~200                    |
| COPT        | [COPT.jl](https://github.com/COPT-Public/COPT.jl/tree/main) | commercial | moment        | 👍👍👍 | 👍👍👍 | 👍👍👍 | ~700                    |
| Hypatia[^1] | [Hypatia.jl](https://github.com/jump-dev/Hypatia.jl)        | MIT        | moment        | 👍👍    | 👍👍   | 👍      | ~100                    |
| LANCELOT[^2]| [GALAHAD.jl](https://github.com/ralna/GALAHAD/tree/master/GALAHAD.jl) | BSD | nonlinear  | n.a.     | n.a.    | 👍👍👍 | n.a.                    |
| Loraine     | ∅[^3]                                                       | MIT        | primal moment | 👍👍👍 | 👍👍   | 👍👍👍 | moderately large         |
| LoRADS      | [LoRADS](https://github.com/COPT-Public/LoRADS)[^4]         | MIT        | primal moment | 👍👍👍 | 👍      | 👍👍👍 | very large              |
| Mosek[^5]   | [Mosek.jl](https://github.com/MOSEK/Mosek.jl)               | commercial | SOS, moment   | 👍👍👍 | 👍👍👍 | 👍👍    | ~300 - 500             |
| ProxSDP     | [ProxSDP.jl](https://github.com/mariohsouto/ProxSDP.jl)     | MIT        | primal moment | 👍👍👍 | 👍👍    | 👍👍👍 | very large             |
| SCS         | [SCS.jl](https://github.com/jump-dev/SCS.jl)                | MIT        | moment        | 👍      | 👍      | 👍👍👍 |                         |
| SketchyCGAL | ∅                                                           |            | primal moment | 👎      | 👍      | 👍👍👍 |                        |
| SpecBM      | ∅[^6]                                                       |            | SOS           | n.a.     | n.a.    | 👍👍👍 |                        |

[^1]: Note that by default, a sparse solver is used (unless the problem was constructed with a `factor_coercive` different from
      one). This is typically a good idea for large systems with not too much monomials. However, if you have a very dense
      system, the sparse solver will take forever; better pass `dense=true` to the optimization routine. This will then be much
      faster (and always much more accurate).
[^2]: LANCELOT is a nonlinear solver that directly works on the problem itself. It does not use a relaxation. Therefore, it
      cannot provide lower-bound guarantees on the objective value; however, there is no problem with the extraction of a
      solution, as the solver directly works on the decision variables. When invoking the LANCELOT solver, a function is
      returned which performs the optimization and which requires a vector of initial values as parameter. This function will
      then return a 2-tuple with the (locally) optimal objective value and the point of the local optimum.
      Currently, the LANCELOT interface does not support complex-valued problems.
[^3]: There is a separate [Julia package for Loraine](https://github.com/kocvara/Loraine.jl). However, the implementation is so
      bad that it is not only unnecessarily inefficient, but would also not allow to solve large-scale systems. Therefore,
      `PolynomialOptimization` provides a rewritten implementation, which is heavily based on the original source code, but
      tries to use the available memory more efficiently (though there is still room for improvement). Since only a subset of
      the features of the original package has been implemented, this is not yet ready to be contributed to the solver; but it
      can be expected that in the future, an external package will be required to use Loraine.
[^4]: There is no Julia package for LoRADS available. You first have to clone the linked Git repositiory and compile the solver
      for your system; then, call [`Solvers.LoRADS.set_solverlib`](@ref) in order to tell `PolynomialOptimization` where to
      look for the binary. After restarting Julia, the method is available.
[^5]: `:MosekMoment` requires at least version 10, `:MosekSOS` already works with version 9.
      The moment variant is more prone to failure in case of close-to-illposed problems; sometimes, this is an issue of the
      presolver, which can be turned off by passing `MSK_IPAR_PRESOLVE_USE="MSK_PRESOLVE_MODE_OFF"` to [`poly_optimize`](@ref).
      The performance indicators in the table are valid for `:MosekSOS`. The new PSD cone interface of Mosek 10 that is used by
      the moment-based variant proves to be much slower than the old one; therefore, using `:MosekMoment` is not recommended.
[^6]: `SpecBM` is provided by `PolynomialOptimization`; however, it requires a subsolver for the quadratic master problem.
      Currently, `Mosek` and `Hypatia` are implemented and must therefore be loaded to make `SpecBM` work.

# Packaged solvers
During the development of this package, several interesting solvers were proposed in research. The ones that were implemented
are documented on this page. They can be accessed from the `PolynomialOptimization.Solvers` submodule.

```@meta
CurrentModule = PolynomialOptimization.Solvers.Loraine
```
## Loraine
The Loraine solver is suiteable for large-scale low-rank semidefinite programming. Note that the implementation provided by
`PolynomialOptimization` supports only a subset of Loraine's features (but this much more efficiently): the direct solver is
not available, data matrices are assumed not to be rank-one and always sparse.
```@docs
Model
Preconditioner
Solver
solve!
```

```@meta
CurrentModule = PolynomialOptimization.Solvers.LoRADS
```
## LoRADS
The experimental LoRADS solver is suitable for large-scale low-rank semidefinite programming. Note that it is currently in a
very early stage of development and not yet very customizable. In particular, at the time of writing, it is not possible to
configure the tolerance of the termination criterion, which is always set to ``10^{-5}``. There is also no interface defined to
extract the solution data; `PolynomialOptimization` relies on the internal representation of the state to access this
information. Note that this is likely to change in a new version; the supported solver version is a patched version of `1.0.0`.
You will also see a lot of solver output, regardless of the `verbose` flag, as this cannot be turned off.
LoRADS has to be compiled from the source into a shared library; its path must then be provided to the package using
[`set_solverlib`](@ref). This setting will take effect after the Julia session is restarted.
```@docs
set_solverlib
ASDP
ConeType
Strategy
init_solver
set_dual_objective
conedata_to_userdata
set_cone
set_lp_cone
init_cone_data
preprocess
load_sdpa
solve
get_X
get_Xlin
```

```@meta
CurrentModule = PolynomialOptimization.Solvers.SketchyCGAL
```
## SketchyCGAL
While the solver was implemented for the purpose of being used within `PolynomialOptimization`, it also works as a standalone
routine (and could in principle be a separate package). [SketchyCGAL](https://doi.org/10.1137/19M1305045) is a solver that
scales very well for large problem sizes by maintaining a sketch of the semidefinite matrices based on the assumption that the
optimal solution has low rank; indeed, in polynomial optimizations, if there is an optimal point for the problem that can be
encoded in the chosen basis, then this automatically gives rise to a rank-one semidefinite encoding of this point.
```@docs
sketchy_cgal
Status
```

```@meta
CurrentModule = PolynomialOptimization.Solvers.SpecBM
```
## SpecBM
While the solver was implemented for the purpose of being used within `PolynomialOptimization`, it also works as a standalone
routine (and could in principle be a separate package). SpecBM is a
[spectral bundle algorithm for primal semidefinite programs](https://doi.org/10.48550/arXiv.2307.07651) and is based on the
assumption that the optimal dual solution has low rank; indeed, in polynomial optimizations, if there is an optimal point for
the problem that can be encoded in the chosen basis, then this automatically gives rise to a rank-one semidefinite moment
matrix this point.
The implementation also allows for free variables and multiple semidefinite constraints and contains further improvements
compared to the [reference implementation](https://github.com/soc-ucsd/specBM). It requires either Hypatia or a rather recent
version of Mosek (at least 10.1.13) as subsolvers.
```@docs
specbm_primal
Result
```

```@meta
CurrentModule = PolynomialOptimization.Solvers.LANCELOT
```
## LANCELOT
`PolynomialOptimization` provides an interface to LANCELOT. While LANCELOT is part of the
[GALAHAD](https://github.com/ralna/GALAHAD) library which has a quite recent Julia interface, the LANCELOT part is still pure
Fortran without even a C interface. Therefore, here, we exploit that the pre-packaged binaries are compiled with GFortran,
version at least 9, so that we know the binary layout of the parameters and can pretend that we like Fortran. Currently, only
`LANCELOT_simple` is supported, which is of course not quite ideal[^7]. Since `Galahad.jl` is a weak dependency, the package
has to be loaded first before the `Solvers.LANCELOT` module becomes available:
```@docs
LANCELOT_simple
```

[^7]: The branch `lancelot` goes further and defines an interface for the full version of LANCELOT, which is a lot more
      sophisticated. Unfortunately, it also seems to be broken at the moment and bugfixing will require some debugging of the
      disassembly. This is not a priority at the moment (which is whenever you read this statement).