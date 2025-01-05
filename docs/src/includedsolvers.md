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

|  Solver     |                            Package                          |   License  | Methods     | Speed    | Accuracy | Memory  | max. recomm. basis size |
| ------:     | :---------------------------------------------------------: | :--------: | :---------: | :-----:  | :------: | :-----: | :---------------------- |
| Clarabel    | [Clarabel.jl](https://github.com/oxfordcontrol/Clarabel.jl) | Apache     | moment      | 👍👍👍  | 👍👍👍 | 👍👍    | ~200                   |
| COPT        | [COPT.jl](https://github.com/COPT-Public/COPT.jl/tree/main) | commercial | moment      | 👍👍👍  | 👍👍👍 | 👍👍👍 | ~700                   |
| Hypatia[^1] | [Hypatia.jl](https://github.com/jump-dev/Hypatia.jl)        | MIT        | moment      | 👍👍    | 👍👍    | 👍      | ~100                   |
| LANCELOT[^2]| [GALAHAD.jl](https://github.com/ralna/GALAHAD/tree/master/GALAHAD.jl) | BSD | nonlinear | n.a.    | n.a.     | 👍👍👍 | n.a.                   |
| Mosek[^3]   | [Mosek.jl](https://github.com/MOSEK/Mosek.jl)               | commercial | SOS, moment | 👍👍👍  | 👍👍👍 | 👍👍    | ~300 - 500             |
| SCS         | [SCS.jl](https://github.com/jump-dev/SCS.jl)                | MIT        | moment      | 👍       | 👍      | 👍👍👍 |                        |
| SpecBM      | ∅[^4]                                                       |            | SOS         | n.a.     | n.a.     | 👍👍👍 |                        |

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
[^3]: `:MosekMoment` requires at least version 10, `:MosekSOS` already works with version 9.
      The moment variant is more prone to failure in case of close-to-illposed problems; sometimes, this is an issue of the
      presolver, which can be turned off by passing `MSK_IPAR_PRESOLVE_USE="MSK_PRESOLVE_MODE_OFF"` to [`poly_optimize`](@ref).
      The performance indicators in the table are valid for `:MosekSOS`. The new PSD cone interface of Mosek 10 that is used by
      the moment-based variant proves to be much slower than the old one; therefore, using `:MosekMoment` is not recommended.
[^4]: `SpecBM` is provided by `PolynomialOptimization`; however, it requires a subsolver for the quadratic master problem.
      Currently, `Mosek` and `Hypatia` are implemented and must therefore be loaded to make `SpecBM` work.

# Packaged solvers
During the development of this package, several interesting solvers were proposed in research. The ones that were implemented
are documented on this page. They can be accessed from the `PolynomialOptimization.Solvers` submodule.

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
`LANCELOT_simple` is supported, which is of course not quite ideal[^5]. Since `Galahad.jl` is a weak dependency, the package
has to be loaded first before the `Solvers.LANCELOT` module becomes available:
```@docs
LANCELOT_simple
```

[^5]: The branch `lancelot` goes further and defines an interface for the full version of LANCELOT, which is a lot more
      sophisticated. Unfortunately, it also seems to be broken at the moment and bugfixing will require some debugging of the
      disassembly. This is not a priority at the moment (which is whenever you read this statement).