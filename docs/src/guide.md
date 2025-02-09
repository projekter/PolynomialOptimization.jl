```@meta
DocTestFilters = [
  r"\d+\.\d+ seconds" => "seconds",
  r"\d\.\d+e-\d+" => s"~0~",
  r"(\d*)\.(\d{6})\d+" => s"\1.\2~"
]
DocTestSetup = :(import Random; Random.seed!(1234))
CurrentModule = PolynomialOptimization
```
# Walkthrough

We start the Julia session by including the required packages.
```jldoctest walkthrough
julia> using PolynomialOptimization, DynamicPolynomials
```

## A simple unconstrained problem
### Constructing the problem
Next, we define some simple optimization problem.
```jldoctest walkthrough
julia> @polyvar x[1:3];

julia> prob = poly_problem(1 + x[1]^4 + x[2]^4 + x[3]^4 + x[1]^2*x[2]^2 + x[1]^2*x[3]^2 + x[2]^2*x[3]^2 + x[2]*x[3])
Real-valued polynomial optimization problem in 3 variables
Objective: 1.0 + x₂x₃ + x₃⁴ + x₂²x₃² + x₂⁴ + x₁²x₃² + x₁²x₂² + x₁⁴
```
This is a very simple problem: We have three variables and want to minimize an unconstrained objective function.
Currently, `prob` is just an instance of a [`Problem`](@ref): some elementary checks and conversions have been done, but the
heavy machinery of polynomial optimization was not applied yet.
During the process of constructing the problem, it is possible to automatically perform modifications. These are available via
keyword parameters of [`poly_problem`](@ref).

### Densely solving the problem
Since this problem is so small, we can solve it directly without any sparsity consideration.
Note that `PolynomialOptimization` works with a variety of solvers; however, they are included only as weak dependencies. You
have to load the appropriate solver Julia package first to make the solvers available. For a list of supported solvers, see the
documentation for [`poly_optimize`](@ref) or the [details](@ref) section.
```jldoctest walkthrough
julia> import Clarabel

julia> res = poly_optimize(:Clarabel, prob)
[ Info: Automatically selecting minimal degree cutoff 2
Polynomial optimization result
Relaxation method: Dense
Used optimization method: ClarabelMoment
Status of the solver: SOLVED
Lower bound to optimum (in case of good status): 0.9166666672624658
Time required for optimization: 0.9356329 seconds
```
The solution that we found was indeed optimal and the value is `0.9166...`.
Note that "optimal" here means that the solver converged for the given problem. However, a polynomial optimization problem is
difficult to solve in general; therefore, it cannot be optimized directly. Instead, a _relaxation_ of the problem has to be
constructed. Normally, this must be done explicitly by instantiating a decendant of [`AbstractRelaxation`](@ref); but when
[`poly_optimize`](@ref) is called with a problem instead of a relaxation, it will by default construct a dense relaxation of
minimal degree.
Therefore, "optimal" in fact only means that the _relaxation_ was solved to global optimality, which in general will only yield
an underestimator to the original problem.
Note that while Clarabel has a very clear return code - `SOLVED` says that things went well - this is not necessarily the case
for other solvers. Use [`issuccess`](@ref) on the result object to check whether the reported solver status is a good one:
```jldoctest walkthrough
julia> issuccess(res)
true
```

Further note that the optimization time seems to be pretty high for such a small problem. However, this is purely due to the
compilation time. Running the optimization again will give a time of the order of a millisecond. Finally, it is not necessary
to specify the solver explicitly if only one solver package is loaded; `poly_optimize(prob)` would work as well. However, if
multiple solvers are available, which one is then chosen may depend on the loading order of the packages.

### Checking optimality
There are two different ways to check whether the given bound is optimal for the original problem. Obviously, if we find a
point ``x`` such that the objective evaluated at ``x`` gives our bound, then this bound must have been optimal. Of course, the
difficulty now lies in finding the point. `PolynomialOptimization` implements a state of the art
[solution extraction algorithm](https://doi.org/10.1016/j.laa.2017.04.015), which can relatively quickly (the cost is
essentially that of performing an SVD on the moment matrix) obtain solutions. This will only be guaranteed to work if the
problem was indeed optimal, there were finitely many solutions in the first place, and a "good" moment matrix is obtained
(i.e., a dense matrix and no low-rank solver was employed) - but there is an alternative which might work well in case these
conditions are not satisfied (apart from optimality, obviously), more on this below.
The function [`poly_solutions`](@ref) gives an iterator that delivers all the (potential) solutions one at a time in an
arbitrary order.
Alternatively, [`poly_all_solutions`](@ref) directly calculates all the solutions and grades them according to how much they
violate the bound or constraints, if any were given. The solutions are then returned in a best-to-worst order.
```jldoctest walkthrough; filter=r"^ .+$"m
julia> poly_all_solutions(res)
2-element Vector{Tuple{Vector{Float64}, Float64}}:
 ([-1.1700613807653743e-18, 0.4082426580485429, -0.408242660645461], 5.266275193704928e-10)
 ([-5.352664236308434e-20, -0.4082426580485437, 0.40824266064546183], 5.266276303927953e-10)
```
Every element in the vector is a tuple, where the first entry corresponds to the optimal variables, and the second term is the
badness of this solution (which can also be calculated manually using [`poly_solution_badness`](@ref)). Since here, the badness
is of the order of ``10^{-8}``, i.e., numerically zero, the points are indeed valid global minima.

Note that the solution extraction functions just give a vector of numbers; to assign these numbers to variables, the function
[`variables`](@ref) can be applied to the problem. This is particularly useful if there are multiple variables created at
different times, occurring differently in the constraints, such that the order is not clear beforehand.
```jldoctest walkthrough
julia> variables(prob)
3-element Vector{Variable{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}}:
 x₁
 x₂
 x₃
```

!!! warning "Variables in the problem"
    A polynomial optimization problem is always constructed using [`poly_problem`](@ref), and this function supports any input
    that implements the `MultivariatePolynomials` interface. In the example here, we used `DynamicPolynomials`, which is
    probably the most common choice.
    However, these inputs are not necessarily well-suited to deliver high performance when constructing the actual
    optimizations. Therefore, `poly_problem` will convert all inputs to an internal polynomial format. This internal format
    does not know about the name of variables - real-valued variables will always be printed as ``x_i``, complex-valued
    variables as ``z_j`` with continuous indices. The [`variables`](@ref) function now becomes even more important: It contains
    the _original_ variables of the polynomials that were given to [`poly_problem`](@ref).

A second way to check for optimality is to use the flat extension/truncation criterion originally due to
[Curto and Fialkow](https://doi.org/10.1090/S0002-9947-00-02472-7) and improved by [Nie](https://arxiv.org/abs/1106.2384v2).
This is a sufficient criterion for optimality, and it can be manually checked by calling [`optimality_certificate`](@ref) on
the problem (it will only work with non-sparse problems). The function will return `:Optimal` if the given minimum value can
be certified to be optimal; else it will return `:Unknown`:
```jldoctest walkthrough
julia> optimality_certificate(res)
:Optimal
```
Note that this is just a sufficient criterion, and the solution might be optimal even if it is violated. As calculating the
certificate will involve calculating the ranks of several matrices (and is more complicated in the complex case), it is not
necessarily cheaper than trying to extract solutions; as the latter is more informative, it should usually be the way to go.

### Extracting a SOS certificate
Whenever the optimization was successful, a valid sums-of-squares certificate will be available, i.e., a decomposition of the
objective (in this simple, unconstrained, case). Here, the minimum value of the objective was found to be `0.9166...`.
We can therefore obtain a certificate for the positivity of the original objective minus this global minimum:
```jldoctest walkthrough
julia> cert = SOSCertificate(res)
Sum-of-squares certificate for polynomial optimization problem
1.0 + x₂x₃ + x₃⁴ + x₂²x₃² + x₂⁴ + x₁²x₃² + x₁²x₂² + x₁⁴ - 0.9166666672624658
= (-0.2492464642498061 + 0.0 + 0.0 + 0.0 + 0.6811832979888123x₃² - 0.13316036665354078x₂x₃ + 0.6812065898234094x₂² + 0.0 + 0.0 + 0.7152479291533441x₁²)²
+ (-0.00015840159683466377 + 0.0 + 0.0 + 0.0 + 0.6360091467299869x₃² + 0.0011902230323755974x₂x₃ - 0.6338687294141168x₂² + 0.0 + 0.0 - 0.0018514462038562963x₁²)²
+ (0.050599062148224884 + 0.0 + 0.0 + 0.0 - 0.35322230868481075x₃² - 0.40669620906370346x₂x₃ - 0.3569973936050725x₂² + 0.0 + 0.0 + 0.6183225665112776x₁²)²
+ (0.0 - 0.6307886589586179x₃ - 0.6307886578526084x₂ + 0.0 + 0.0 + 0.0 + 0.0 + 0.0 + 0.0 + 0.0)²
+ (0.0 + 0.0 + 0.0 + 0.0 + 0.0 + 0.0 + 0.0 - 0.5707311041008931x₁x₃ - 0.5706316416095126x₁x₂ + 0.0)²
+ (0.1365622366906534 + 0.0 + 0.0 + 0.0 - 0.08194926810609472x₃² + 0.6553389692882523x₂x₃ - 0.08198096663343958x₂² + 0.0 + 0.0 + 0.32572101226607164x₁²)²
+ (0.0 + 0.0 + 0.0 + 0.0 + 0.0 + 0.0 + 0.0 - 0.43861947273884583x₁x₃ + 0.4386959251861786x₁x₂ + 0.0)²
+ (0.0 + 0.0 + 0.0 + 0.4527802849918676x₁ + 0.0 + 0.0 + 0.0 + 0.0 + 0.0 + 0.0)²
```
This certificate consists of a number of polynomials that, when squared and added, should give rise to the original objective.
Note that when printing the certificate, values that are below a certain threshold will be set to zero by default.
We can also explicitly iterate through all the polynomials and sum them up, although we have to be careful to map them back to
their original representation for this:
```jldoctest walkthrough
julia> p = zero(polynomial_type(x, Float64));

julia> for pᵢ in cert[:objective, 1] # no sparsity, so there is just a single grouping
           p += PolynomialOptimization.change_backend(pᵢ, x)^2
       end

julia> map_coefficients!(x -> round(x, digits=8), p)
0.08333334 + x₂x₃ + x₃⁴ + x₂²x₃² + x₂⁴ + x₁²x₃² + x₁²x₂² + x₁⁴
```
Note how this is precisely the objective minus the global minimum.

### Using the Newton polytope
The current example is an unconstrained optimization problem; hence, the size of the full basis, which is 10, may be larger
than actually necessary. It is not a simple problem to determine the relevant basis elements in general; but unconstrained
problems allow for the Newton polytope technique. To use it, we first need to load a supported solver for the Newton polytope,
then we simply explicitly construct the [`Newton`](@ref Relaxation.Newton) relaxation object:
```jldoctest walkthrough
julia> import Mosek

julia> Relaxation.Newton(prob)
[ Info: Automatically selecting minimal degree cutoff 2
Relaxation.Newton of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
PSD block sizes:
  [10 => 1]
Relaxation degree: 2
```
In this case, no basis reduction was possible. However, in other cases, this can work. For example, if you want to inspect the
Newton polytope of the polynomials whose squares might make up a certain objective, you can call [`Newton.halfpolytope`](@ref)
directly (the name comes from the fact that you pass the objective to the function [`Newton.halfpolytope`](@ref), and by a
[theorem by Reznick](https://doi.org/10.1215/S0012-7094-78-04519-2), the Newton polytope of the decomposition functions that
have to be squared to give the objective will be contained in half the Newton polytope of the objective itself):
```jldoctest walkthrough
julia> @polyvar x y;

julia> Newton.halfpolytope(x^4*y^2 + x^2*y^4 - 3x^2*y^2 +1)
4-element MonomialVector{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}:
 1
 xy
 xy²
 x²y
```
This reveals that, were the Motzkin representable by a sum of squares, the equality
```math
x^4 y^2 + x^2 y^4 - 3 x^2 y^2 + 1 = \sum_i (\alpha_i + \beta_i x y + \gamma_i x y^2 + \delta_i x^2 y)^2
```
would have to hold; but expanding the right-hand side will lead to the coefficient ``\sum_i \beta_i^2`` in front of the
monomial ``x^2 y^2``, which cannot be negative; hence, the Motzkin polynomial is not a sum of squares.

Note that the calculation of the Newton polytope currently requires Mosek or COPT. There are some preprocessing options that
may be able to speed up the calculation, although it is already extremely fast by itself and can calculate the correct basis
for objectives with hundreds of terms in a decent time (which can be further reduced by multithreading or distributed
computing). Check out the documentation for [`Newton.halfpolytope`](@ref) for more information.

In case you already happen to know a (better) choice of basis, you may opt for [`Relaxation.Custom`](@ref). Note that
relaxations are built incrementally, where the only relaxation that can be constructed directly from a problem is the dense
one. So what actually happened in calling `Relaxation.Newton(prob)` is that first a dense relaxation of the problem was
constructed, which was then passed on to the Newton relaxation constructor: `Relaxation.Newton(Relaxation.Dense(prob))`. The
info message about the minimal degree cutoff was generated by the dense relaxation. "Construction" here does not mean that the
full dense basis was actually built in memory; a lazy representation is used for the dense basis.
This has the consequence that if you use a custom basis, you can then decide to refine this further by passing the custom
relaxation to the Newton relaxation constructor. Similarly, every relaxation can serve as the starting point for another one.

### Applying inexact sparsity
There are four kinds of inexact sparsity method implemented in `PolynomialOptimization`:
- [`Relaxation.SparsityCorrelative`](@ref)
- [`Relaxation.SparsityTermBlock`](@ref)
- [`Relaxation.SparsityTermChordal`](@ref)
- [`Relaxation.SparsityCorrelativeTerm`](@ref)
Applying the analysis as simple as passing the problem to the respective type. For this particular problem, there is no
correlative sparsity:
```jldoctest walkthrough
julia> Relaxation.SparsityCorrelative(prob)
[ Info: Automatically selecting minimal degree cutoff 2
Relaxation.SparsityCorrelative of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
PSD block sizes:
  [10 => 1]
```
So there is only a single clique, leading to a basis of size `10`. However, there is term sparsity:
```jldoctest walkthrough
julia> tbs = Relaxation.SparsityTermBlock(prob)
[ Info: Automatically selecting minimal degree cutoff 2
Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
PSD block sizes:
  [5 => 1, 2 => 1, 1 => 3]
```
We get a basis of size `5`, one of size `2`, and three bases of size `1` (here, by _basis_ we mean a set of monomials that
indexes the moment/SOS matrices). `PolynomialOptimization` will model these by a `5 × 5` semidefinite matrix, a rotated
second-order cone, as well as three linear constraints. This is much cheaper than a `10 × 10` semidefinite matrix.
Let's optimize the sparse problem:
```jldoctest walkthrough
julia> poly_optimize(:Clarabel, tbs)
Polynomial optimization result
Relaxation method: SparsityTerm
Used optimization method: ClarabelMoment
Status of the solver: SOLVED
Lower bound to optimum (in case of good status): 0.9166666718972408
Time required for optimization: 0.00306 seconds
```
Again, we get the same optimal value, so introducing the sparsity did not make our relaxation worse (which is _per se_ not
guaranteed), and we are still able to get the same optimal solutions:
```jldoctest walkthrough; filter=r"^ .+$"m
julia> poly_all_solutions(ans)
2-element Vector{Tuple{Vector{Float64}, Float64}}:
 ([0.0, 0.4082661947158492, -0.408266194715714], 4.589421509493263e-9)
 ([0.0, -0.4082661947158492, 0.408266194715714], 4.589421509493263e-9)
```
Note that perhaps surprisingly, `PolynomialOptimization` can still deliver good optimal points despite the fact that term
sparsity was in effect. The usual extraction algorithms will fail, as for every moment ``m``, they also require the moment
``m x`` to be present for all variables ``x`` - term sparsity cannot to provide this (just check that the moment matrix
contains lots of `NaN` values).
Hence, [`poly_all_solutions`](@ref) will automatically switch to a different heuristic solution extraction algorithm that is
always successful in the simple case of a rank-1 moment matrix, but can also often also give good results in the more general
case such as here. As a rule of thumb, if all the solutions encoded in the moment matrix differ only by the signs or phases of
individual components, the heuristic will be successful. Still, the fact that the moments may encode multiple solutions may be
an issue that can prevent successfully obtaining a solution vector. We will introduce a way to bypass this problem below.

Assume that our term sparsity gave a worse bound than the dense case (which in general we would not know, since the dense
problem is typically far too large to be solved; but we just don't get a proper optimal point, though this could also be due to
the Lasserre hierarchy level being insufficient). Then, we could try to iterate the term sparsity hierarchy, keeping the same
level in the Lasserre hierarchy.
```jldoctest walkthrough
julia> Relaxation.iterate!(tbs)
Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
PSD block sizes:
  [5 => 1, 2 => 2, 1 => 1]

julia> Relaxation.iterate!(tbs)
```
In general, we simply use [`Relaxation.iterate!`](@ref) to move to the next higher level (note how two of the linear
constraints were merged into a quadratic constraint). If `iterate!` returns the new sparsity object, something changed and we
might try to optimize the new sparse problem, getting a potentially better bound. If `iterate!` returns `nothing`, the hierachy
terminated and nothing more can (for term sparsity: must) be done (as the last level for term sparsity is as good as the dense
problem).

We can also try what happens if we use term sparsity with chordal cliques instead of connected components:
```jldoctest walkthrough
julia> tcs = Relaxation.SparsityTermChordal(prob)
[ Info: Automatically selecting minimal degree cutoff 2
Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
PSD block sizes:
  [4 => 1, 2 => 2, 1 => 3]

julia> res = poly_optimize(:Clarabel, tcs)
Polynomial optimization result
Relaxation method: SparsityTerm
Used optimization method: ClarabelMoment
Status of the solver: SOLVED
Lower bound to optimum (in case of good status): 0.9166666672685418
Time required for optimization: 0.0027146 seconds
```
So again, we get the same optimal result. We could even extract a solution point with better accuary, and the problem was
smaller, as we now have just a basis of size `4` instead of `5` (at the cost of another quadratic constraint, which is much
cheaper than larger semidefinite matrices).

### [Details on the optimization process](@id details)
The first parameter for [`poly_optimize`](@ref) is the solver/method that is used to optimize the problem. For a list of
supported methods, see [the solver reference](@ref solvers_poly_optimize).

Note that by passing the keyword argument `verbose=true` to the optimization function, we get some more insight into what
happens behind the hood. Let's redo the last optimization.
```julia-repl
julia> poly_optimize(:Clarabel, tcs, verbose=true)
Beginning optimization...
Clique merging disabled.
PSD block sizes:
  [4 => 1, 2 => 2, 1 => 3]
Starting solver...
Setup complete in 0.000266 seconds
-------------------------------------------------------------
           Clarabel.jl v0.9.0  -  Clever Acronym
                   (c) Paul Goulart
                University of Oxford, 2022
-------------------------------------------------------------

problem:
  variables     = 11
  constraints   = 20
  nnz(P)        = 0
  nnz(A)        = 24
  cones (total) = 7
    : Zero        = 1,  numel = 1
    : Nonnegative = 3,  numel = (1,1,1)
    : SecondOrder = 2,  numel = (3,3)
    : PSDTriangle = 1,  numel = 10

settings:
  linear algebra: direct / qdldl, precision: Float64
  max iter = 200, time limit = Inf,  max step = 0.990
  tol_feas = 1.0e-08, tol_gap_abs = 1.0e-08, tol_gap_rel = 1.0e-08,
  static reg : on, ϵ1 = 1.0e-08, ϵ2 = 4.9e-32
  dynamic reg: on, ϵ = 1.0e-13, δ = 2.0e-07
  iter refine: on, reltol = 1.0e-13, abstol = 1.0e-12,
               max iter = 10, stop ratio = 5.0
  equilibrate: on, min_scale = 1.0e-04, max_scale = 1.0e+04
               max iter = 10

iter    pcost        dcost       gap       pres      dres      k/t        μ       step
---------------------------------------------------------------------------------------------
  0   1.0000e+00   1.0000e+00  0.00e+00  4.99e-01  5.83e-01  1.00e+00  1.91e+00   ------
  1   1.0240e+00   1.0238e+00  2.08e-04  6.91e-02  8.37e-02  9.27e-02  3.10e-01  9.09e-01
  2   9.2888e-01   9.2859e-01  2.93e-04  5.10e-03  6.37e-03  6.34e-03  2.54e-02  9.21e-01
  3   9.1842e-01   9.1835e-01  7.08e-05  1.12e-03  1.42e-03  1.40e-03  5.85e-03  8.71e-01
  4   9.1762e-01   9.1759e-01  3.15e-05  4.04e-04  5.09e-04  4.96e-04  2.15e-03  7.53e-01
  5   9.1673e-01   9.1672e-01  3.94e-06  5.76e-05  7.26e-05  7.15e-05  3.06e-04  8.64e-01
  6   9.1668e-01   9.1668e-01  4.44e-07  7.88e-06  9.94e-06  9.86e-06  4.16e-05  9.30e-01
  7   9.1667e-01   9.1667e-01  1.05e-07  1.48e-06  1.86e-06  1.83e-06  7.83e-06  8.23e-01
  8   9.1667e-01   9.1667e-01  3.09e-08  3.49e-07  4.40e-07  4.26e-07  1.85e-06  8.84e-01
  9   9.1667e-01   9.1667e-01  6.54e-09  6.68e-08  8.41e-08  8.08e-08  3.54e-07  8.23e-01
 10   9.1667e-01   9.1667e-01  1.35e-09  1.32e-08  1.66e-08  1.59e-08  6.99e-08  9.04e-01
 11   9.1667e-01   9.1667e-01  2.79e-10  2.51e-09  3.16e-09  3.01e-09  1.34e-08  8.23e-01
---------------------------------------------------------------------------------------------
Terminated with status = solved
solve time = 5.69ms
Optimization complete, retrieving moments
Polynomial optimization result
Relaxation method: SparsityTerm
Used optimization method: ClarabelMoment
Status of the solver: SOLVED
Lower bound to optimum (in case of good status): 0.9166666672685418
Time required for optimization: 0.1599203 seconds
```
So first, `PolynomialOptimization` will determine the bases for the matrices according to the sparsity pattern.
Note that after this step, the if the resulting relaxation is wrapped into a [`Relaxation.CliqueMerged`](@ref), an attempt will
be made to merge bases if their heuristic cost for treating them separately would be worse than joining them (this concept is
nicely explained in the [COSMO documenation](https://oxfordcontrol.github.io/COSMO.jl/stable/decomposition/#Clique-merging)).
In general, doing clique merging will lead to faster optimizations; however, the merging process itself can be quite costly and
in fact for large problems might cost much more time than it gains - hence, it must be enabled by explicitly constructing the
merged relaxation.
After this step is done, the Clarabel data (or any other optimizer structure, which we all address directly without `JuMP`) is
constructed; then the solver runs.

Indeed, due to sparsity, the moment matrix is full of unknowns:
```julia-repl
julia> show(stdout, "text/plain", moment_matrix(res))
10×10 LinearAlgebra.Symmetric{Float64, Matrix{Float64}}:
   1.0         NaN         NaN         NaN             0.166666    -0.166666     0.166666   NaN          NaN            1.60942e-8
 NaN             0.166666   -0.166666  NaN           NaN          NaN          NaN          NaN          NaN          NaN
 NaN            -0.166666    0.166666  NaN           NaN          NaN          NaN          NaN          NaN          NaN
 NaN           NaN         NaN           1.60942e-8  NaN          NaN          NaN          NaN          NaN          NaN
   0.166666    NaN         NaN         NaN             0.0277777  NaN            0.0277777  NaN          NaN            2.5981e-9
  -0.166666    NaN         NaN         NaN           NaN            0.0277777  NaN          NaN          NaN          NaN
   0.166666    NaN         NaN         NaN             0.0277777  NaN            0.0277777  NaN          NaN            2.5981e-9
 NaN           NaN         NaN         NaN           NaN          NaN          NaN            2.5981e-9  NaN          NaN
 NaN           NaN         NaN         NaN           NaN          NaN          NaN          NaN            2.5981e-9  NaN
   1.60942e-8  NaN         NaN         NaN             2.5981e-9  NaN            2.5981e-9  NaN          NaN           -3.42786e-9
```
The rows and columns of the matrix are indexed by the basis of the relaxation:
```jldoctest walkthrough
julia> Relaxation.basis(tcs)
10-element PolynomialOptimization.SimplePolynomials.SimpleMonomialVector{3, 0, UInt64, PolynomialOptimization.SimplePolynomials.MultivariateExponents.ExponentsDegree{3, UInt64}, PolynomialOptimization.SimplePolynomials.SimpleMonomial{3, 0, UInt64, PolynomialOptimization.SimplePolynomials.MultivariateExponents.ExponentsDegree{3, UInt64}}}:
 1
 x₃
 x₂
 x₁
 x₃²
 x₂x₃
 x₂²
 x₁x₃
 x₁x₂
 x₁²
```
Combining the basis information with the moment matrix, we can see how the package is able to return solutions without having
access to the full moment matrix. There are values for the squares of the variables available, so we can deduce two possible
candidates for the original variables - at least, if the values assigned to the moments are consistent. Choosing among the
signs becomes possible by looking for mixed terms.

### Always extracting a solution
The fact that `PolynomialOptimization` was unable to extract a (valid) solution can either mean that the relaxation was
insufficient and did not converge to the optimum of the _actual_ problem, or that there are multiple solutions which are too
difficult for the heuristic to grasp.
There is a simple remedy of this problem: By introducing a small, linear perturbation to the objective, the solution will
almost surely be unique; so if the extracted solution is bad, this means that the relaxation was insufficient.
Of course, now the solution is not a solution to the original problem, but the perturbed one - but assuming a robust problem,
the returned optimal point will be close to the actual global optimum. It can therefore then be used as an initial point to
another nonlinear solver that will deliver the true global optimum.
`PolynomialOptimization` makes adding a perturbation easy: Just call [`poly_problem`](@ref) with the keyword parameter
`perturbation=...`, where the magnitude of the perturbation should be specified (typically between `1e-3` and `1e-6` is a good
guess).
Note that adding a perturbation may degrade sparsity. For this, you may also give a vector of the same length as the number of
variables, specifying a different perturbation magnitude for each variable (or just disabling the perturbation by passing `0`).

### Changing the internal representation
Usually, the moment matrix is modeled as a semidefinite constraint in the solver. However, semidefinite programs scale much
less favorably than other types of convex optimization programs such as linear or quadratic ones. It is possible to change
between the internal representations that are used for the moment matrix: apart from semidefinite, also the
[diagonally dominant and scaled diagonally dominant cones](https://doi.org/10.1137/18M118935X) are supported. While they scale
better, they usually provide worse bounds:
```jldoctest walkthrough
julia> res_dd = poly_optimize(:Clarabel, prob, representation=RepresentationDD())
Polynomial optimization result
Relaxation method: Dense
Used optimization method: ClarabelMoment
Status of the solver: SOLVED
Lower bound to optimum (in case of good status): 0.5000000105896211
Time required for optimization: 0.5579216 seconds
```
Again, the long time is purely due to precompilation of the new methods.

This time, no semidefinite constraint was employed, only linear ones. As is clearly visible, the bound is quite bad. But this
process can be iterated and the diagonally dominant cone can be [rotated](https://doi.org/10.48550/arXiv.1510.01597) based on
the data from the previous optimization. This is call re-optimization, and while a lot of details can be customized, if none
are specified, `PolynomialOptimization` will take the Cholesky decomposition of the matrix underlying the SOS certificate of
the previous iteration as a new rotation basis.
!!! info
    Different representations lead to different problems from the perspective of the solver - in particular, different with
    respect to which solver functions are used or how the data is aligned in memory. Therefore, changing the _type_ of the
    representation requires a completely new problem to be set up. However, if only the _data within_ the representatin is
    changed, the old problem can be re-used, reducing the setup time. This is irrelevant if a solver is used whose interface
    in `PolynomialOptimization` does not support this faster reoptimization.
    Note that by default, the default rotations used for the DD and SDD representation are identities of the type
    `UniformScaling`. This is only compatible with other diagonal rotations. Use [`RepresentationIAs`](@ref) to "fake" an
    upper triangular identity at the beginning.

```jldoctest walkthrough
julia> res_dd_rotated = poly_optimize(res_dd)
# output truncated
Lower bound to optimum (in case of good status): 0.7882579368640297
Time required for optimization: 1.5101965 seconds

julia> res_dd_rotated = poly_optimize(res_dd_rotated)
# output truncated
Lower bound to optimum (in case of good status): 0.8744697853614662
Time required for optimization: 0.0355882 seconds

julia> res_dd_rotated = poly_optimize(res_dd_rotated)
# output truncated
Lower bound to optimum (in case of good status): 0.9122827868379025
Time required for optimization: 0.0097173 seconds

julia> res_dd_rotated = poly_optimize(res_dd_rotated)
# output truncated
Lower bound to optimum (in case of good status): 0.9160797411325013
Time required for optimization: 0.0102907 seconds

julia> res_dd_rotated = poly_optimize(res_dd_rotated)
# output truncated
Lower bound to optimum (in case of good status): 0.9165662778912137
Time required for optimization: 0.0123682 seconds
```
So indeed, after a couple of iterations, the optimum is approached pretty well. We could have used the scaled diagonally
dominant representation instead, which relies on quadratic instead of linear programs, which for this particular example would
have been exact without any iteration.

## Constraints

### Equality constraints
Equality constraints are accessible by passing the keyword argument `zero` to [`poly_problem`](@ref), which constrains those
polynomials to be zero. They are relatively cheap to realize in the solver, as they don't require another semidefinite matrix,
just linear constraints or free scalar variables depending on the approach.
```jldoctest walkthrough; filter=r"^ .+$"m
julia> @polyvar x[1:2];

julia> poly_optimize(:Clarabel, poly_problem(-(x[1] -1)^2 - (x[1] - x[2])^2 - (x[2] -3)^2,
                                             zero=[(x[1] - 1.5)^2 + (x[2] - 2.5)^2 - .5]), 1)
Polynomial optimization result
Relaxation method: Dense
Used optimization method: ClarabelMoment
Status of the solver: SOLVED
Lower bound to optimum (in case of good status): -3.999999965663831
Time required for optimization: 0.0014009 seconds

julia> poly_all_solutions(ans)
1-element Vector{Tuple{Vector{Float64}, Float64}}:
 ([1.0000155981763819, 3.0000155844516696], 2.0076497353471723e-8)

```
Note that when grading the quality of a solution, the package will determine the violation of the constraints as well as how
far the actual value is away from what it should be, and return the worst of all violations.
Note that in principle, Gröbner basis methods would allow to incorporate equality constraints with a potentially even higher
reduction in the number of variables. While an early version of `PolynomialOptimization` supported this, any Gröbner basis
method has been removed from the package. Experience showed that the cost of calculating a Gröbner basis can easily be many
times larger than working with the original problem; furthermore, then taking everything modulo this basis prevents some
optimizing assumptions to be made during the problem construction. Lastly, removing a variable or constraint does not help a
lot with respect to scaling, as the main issue is the size of the semidefinite cones - which would be given by a basis of
standard monomials with respect to the Gröbner basis, and the savings there are often minuscule.

### Inequality constraints
Inequality constraints are implemented using Putinar's Positivstellensatz or localizing matrices. They can be specified by
passing the keyword argument `nonneg` to [`poly_problem`](@ref), which constraints those polynomials to be greater or
equal to zero.
```jldoctest walkthrough
julia> @polyvar x[1:2];

julia> prob = poly_problem(-(x[1]-1)^2 - (x[1]-x[2])^2 - (x[2]-3)^2,
                           nonneg=[1-(x[1]-1)^2, 1-(x[1]-x[2])^2, 1-(x[2]-3)^2])
Real-valued polynomial optimization problem in 2 variables
Objective: -10.0 + 6.0x₂ + 2.0x₁ - 2.0x₂² + 2.0x₁x₂ - 2.0x₁²
3 nonnegative constraints
1: 2.0x₁ - x₁² ≥ 0
2: 1.0 - x₂² + 2.0x₁x₂ - x₁² ≥ 0
3: -8.0 + 6.0x₂ - x₂² ≥ 0

julia> poly_optimize(:Clarabel, prob, 1)
Polynomial optimization result
Relaxation method: Dense
Used optimization method: ClarabelMoment
Status of the solver: SOLVED
Lower bound to optimum (in case of good status): -2.9999999986040407
Time required for optimization: 0.002556 seconds

julia> poly_optimize(:Clarabel, prob, 2)
Polynomial optimization result
Relaxation method: Dense
Used optimization method: ClarabelMoment
Status of the solver: ALMOST_SOLVED
Lower bound to optimum (in case of good status): -2.000000014367033
Time required for optimization: 0.0057509 seconds
```
This is an example where the first relaxation level is not optimal, but the second is, as inspecting the solutions will show
(which also allows us to ignore the somewhat uncertain status of the solver).

### PSD constraints
`PolynomialOptimization` also supports conditions that constrain a matrix that is made up of polynomials to be positive
semidefinite. They can be specified by passing the keyword argument `psd` to [`poly_problem`](@ref); note that the
matrices must be symmetric/hermitian.
```jldoctest walkthrough
julia> @polyvar x[1:2];

julia> prob = poly_problem(-x[1]^2 - x[2]^2, zero=[x[1]+x[2]-1],
                           psd=[[1-4x[1]*x[2]  x[1]; x[1]  4-x[1]^2-x[2]^2]])
Real-valued polynomial optimization problem in 2 variables
Objective: -x₂² - x₁²
1 equality constraint
1: -1.0 + x₂ + x₁ = 0
1 semidefinite constraint
2: [1.0 - 4.0x₁x₂  x₁
    x₁             4.0 - x₂² - x₁²] ⪰ 0

julia> poly_optimize(:Clarabel, prob, 1)
Polynomial optimization result
Relaxation method: Dense
Used optimization method: ClarabelMoment
Status of the solver: SOLVED
Lower bound to optimum (in case of good status): -3.999999994240309
Time required for optimization: 0.1618452 seconds

julia> poly_optimize(:Clarabel, prob, 2)
Polynomial optimization result
Relaxation method: Dense
Used optimization method: ClarabelMoment
Status of the solver: SOLVED
Lower bound to optimum (in case of good status): -3.904891539034092
Time required for optimization: 0.0041964 seconds

julia> optimality_certificate(ans)
:Optimal
```
At second level, we get the optimal solution.

### Improving the optimization without changing the level
The problem can be further tightened by a careful analysis, as [Nie](https://doi.org/10.1007/s10107-018-1276-2) noted, by
rewriting the Lagrange multipliers as polynomials - which will not modify the problem if the minimum is attained at a critical
point (but not that non-critical global minima will be missed).
`PolynomialOptimization` is able to automatically analyze the problem and add the tightening constraints (Mosek or COPT are
required at the moment). For this, simply pass `tighter=true` (or `tighter=:Mosek` resp. `tighter=:COPT`) to `poly_problem`.
This will result in a preprocessing that adds constraints, so expect the problem to grow. To see the progress during the
preprocessing stage, use `verbose=true`.
It may be the case that the required tightening polynomials cannot be determined since their degree always turns out to be
insufficient to satisfy the conditions. Since `PolynomialOptimization` cannot distinguish this from the case where the degree
is just quite high, the procedure may run into an infinite(ly-seeming) loop.
Complex-valued problems are not supported at the moment; and PSD constraints will be skipped during the tightening.
```jldoctest walkthrough; filter=r"-1\.50\d+"=>"-1.50"
julia> @polyvar x y;

julia> poly_optimize(:Clarabel, poly_problem(x^4*y^2 + x^2*y^4 - 3x^2*y^2 +1), 5)
Polynomial optimization result
Relaxation method: Dense
Used optimization method: ClarabelMoment
Status of the solver: ALMOST_SOLVED
Lower bound to optimum (in case of good status): -1.5097034199113162
Time required for optimization: 0.4214607 seconds
```
The given problem is quite hard, as it leads to ill-posed optimization problems with which most solvers expectedly struggle (in fact, while Clarabel gives a solution, the result will differ significantly depending on the operating system).
Adding the tightening equalities (here, as there are no additional constraints, this just means to add the condition
``\nabla\mathrm{objective} = 0``), the fifth order is already sufficient:
```jldoctest walkthrough
julia> prob = poly_problem(x^4*y^2 + x^2*y^4 - 3x^2*y^2 +1, tighter=true)
Real-valued polynomial optimization problem in 2 variables
Objective: 1.0 - 3.0x₁²x₂² + x₁²x₂⁴ + x₁⁴x₂²
2 equality constraints
1: -6.0x₁x₂² + 2.0x₁x₂⁴ + 4.0x₁³x₂² = 0
2: -6.0x₁²x₂ + 4.0x₁²x₂³ + 2.0x₁⁴x₂ = 0

julia> res = poly_optimize(:Clarabel, prob, 5)
Polynomial optimization result
Relaxation method: Dense
Used optimization method: ClarabelMoment
Status of the solver: SOLVED
Lower bound to optimum (in case of good status): 8.741532679530553e-8
Time required for optimization: 0.1160416 seconds
```
Here, it appears that the default solution extraction mechanism does not work well (in fact, since the algorithm is randomized,
you'll get a vastly different result whenever the extraction is performced), so let's try to get the solution via the heuristic
method:
```jldoctest walkthrough; filter=r"^ .+$"m
julia> poly_all_solutions(:heuristic, res)
4-element Vector{Tuple{Vector{Float64}, Float64}}:
 ([1.000000651012571, 1.0000006504036543], 7.809732969654704e-6)
 ([-1.000000651012571, 1.0000006504036543], 7.809732969654704e-6)
 ([1.000000651012571, -1.0000006504036543], 7.809732969654704e-6)
 ([-1.000000651012571, -1.0000006504036543], 7.809732969654704e-6)
```
This was successful in delivering multiple solutions.

### Helping convergence
Another way to modify the problem is to exploit a prefactor in the objective.
[Mai et al.](https://doi.org/10.1007/s10107-021-01634-1) showed that by changing the objective from ``f(x)`` to
``\theta^k(x) \bigl(f(x) + \epsilon \theta^d(x)\bigr)``, where ``\theta(x) = 1 + \Vert x\Vert^2`` there is a bound on ``k``
that guarantees membership in the degree-``(k + d)`` sums-of-squares cone. In contrast to other known bounds, this one is very
easy to calculate, and it is not exponential (at least for unconstrained problems...). It holds even for noncompact feasible
sets, in contrast to Putinar's result.
The price to pay is that the objective itself is of course modified and therefore the optimal value of the problem is only in a
neighborhood of the original problem.

By using the `noncompact=(ϵ, k)` when constructing the problem using [`poly_problem`](@ref), this is done automatically.
Let us apply this to the Motzkin case:
```jldoctest walkthrough
julia> prob = poly_problem(x^4*y^2 + x^2*y^4 - 3x^2*y^2 +1, noncompact=(1e-5, 1))
Real-valued polynomial optimization problem in 2 variables
Objective: 1.00001 + 1.00004x₂² + 1.00004x₁² + 6.000000000000001e-5x₂⁴ - 2.99988x₁²x₂² + 6.000000000000001e-5x₁⁴ + 4.0e-5x₂⁶ - 1.9998799999999999x₁²x₂⁴ - 1.9998799999999999x₁⁴x₂² + 4.0e-5x₁⁶ + 1.0e-5x₂⁸ + 1.00004x₁²x₂⁶ + 2.00006x₁⁴x₂⁴ + 1.00004x₁⁶x₂² + 1.0e-5x₁⁸
Objective was scaled by the prefactor 1.0 + x₂² + x₁²

julia> poly_optimize(:Clarabel, prob)
[ Info: Automatically selecting minimal degree cutoff 4
Polynomial optimization result
Relaxation method: Dense
Used optimization method: ClarabelMoment
Status of the solver: SOLVED
Lower bound to optimum (in case of good status): 0.0002699763854160192
Time required for optimization: 0.0046539 seconds
```
Indeed, now a basis of degree 4 was sufficient to find that the minimum value looks pretty nonnegative. However, this is hard
to quantify, as for this, we'd have to extract a solution from the perturbed `prob`. The algorithm to do this is not
implemented at the moment, as it would require the successive construction and solution of multiple polynomial optimization
problems, which is not very efficient.

## Complex-valued problems
`PolynomialOptimization` fully supports the [complex-valued Lasserre hierarchy](https://doi.org/10.1137/15M1034386), including
its [sparse analysis](https://doi.org/10.1007/s10957-021-01975-z). For this, simply use `@complex_polyvar` instead of
`@polyvar` to declare your variables as complex. Note that feature of `DynamicPolynomials` requires at least version `0.6`.
Use `conj` at your discretion, but note that `real` and `imag` should not be used in the problem description! Instead, use
`(z + conj(z))/2` for the real and `im*(conj(z) - z)/2` for the imaginary part, as well as `z*conj(z)` for the absolute value
square.

As soon as [`poly_problem`](@ref) detects complex variables, it switches to the complex-valued hierarchy.
```jldoctest walkthrough; filter=r"^ .+$"m
julia> @complex_polyvar z;

julia> prob = poly_problem(z + conj(z), zero=[z*conj(z)-1])
Complex-valued polynomial optimization problem in 1 variable
Objective: z̅₁ + z₁
1 equality constraint
1: (-1.0 + 0.0im) + z₁z̅₁ = 0

julia> poly_optimize(:Clarabel, prob)
[ Info: Automatically selecting minimal degree cutoff 1
Polynomial optimization result
Relaxation method: Dense
Used optimization method: ClarabelMoment
Status of the solver: SOLVED
Lower bound to optimum (in case of good status): -1.9999999928826857
Time required for optimization: 1.1383336 seconds

julia> poly_all_solutions(ans)
1-element Vector{Tuple{Vector{ComplexF64}, Float64}}:
 ([-1.0000000000000042 + 0.0im], 7.117322731176046e-9)
```
The dense solution extraction mechanism also works in the complex case.

Let's try a more complicated example from the paper on the complex-valued Lasserre hierarchy (example 4.1):
```jldoctest walkthrough; filter=r"^ .+$"m
julia> @complex_polyvar z[1:2];

julia> prob = poly_problem(3 - z[1]*conj(z[1]) - .5im*z[1]*conj(z[2])^2 + .5im*z[2]^2*conj(z[1]),
                           zero=[z[1]*conj(z[1])-.25z[1]^2-.25conj(z[1])^2-1, # abs(z₁)^2 - z₁^2/4 - conj(z₁)^2/4 = 1
                                 z[1]*conj(z[1])+z[2]*conj(z[2])-3, # abs(z₁)^2 + abs(z₂)^2 = 3
                                 im*z[2]-im*conj(z[2])], # i z₂ - i conj(z₂) = 0
                           nonneg=[z[2]+conj(z[2])]); # z₂ + conj(z₂) ≥ 0

julia> poly_optimize(:Clarabel, prob, 3)
Polynomial optimization result
Relaxation method: Dense
Used optimization method: ClarabelMoment
Status of the solver: SOLVED
Lower bound to optimum (in case of good status): 0.42817470663218404
Time required for optimization: 1.2524744 seconds

julia> poly_all_solutions(ans)
1-element Vector{Tuple{Vector{ComplexF64}, Float64}}:
 ([2.2522776434581143e-16 - 0.8164965543178686im, 1.527525200200192 + 3.4269901696444247e-22im], 1.3954041122588023e-7)
```
Indeed, this solution gives the same objective value and satisfies the constraints, so we found the optimum! Again, note that
the longer optimization times are due to compilation times, as some of the methods for handling complex-valued variables needed
to be compiled first.

And finally something with matrices:
```jldoctest walkthrough; filter=r"^ .+$"m
julia> res = poly_optimize(:Clarabel, poly_problem(-z[1]*conj(z[1]) - z[2]*conj(z[2]),
                                                   psd=[[1-2*(z[1]*z[2]+conj(z[1]*z[2]))  z[1]
                                                         conj(z[1])  4-z[1]*conj(z[1])-z[2]*conj(z[2])]]), 3)
Polynomial optimization result
Relaxation method: Dense
Used optimization method: ClarabelMoment
Status of the solver: SOLVED
Lower bound to optimum (in case of good status): -3.9999999661436303
Time required for optimization: 0.3807631 seconds

julia> poly_all_solutions(res)
5-element Vector{Tuple{Vector{ComplexF64}, Float64}}:
 ([0.26030598038884634 - 0.00026736336860256715im, 0.3407935395847438 + 0.00035003309799230153im], 3.816100332088393)
 ([-0.2603059803888315 + 0.0002673633686024347im, -0.34079353958474984 - 0.00035003309799188205im], 3.8161003320883964)
 ([-0.06036535235750472 + 6.200197140731072e-5im, -0.18480144151354266 - 0.00018981175865375325im], 3.962204377720153)
 ([0.06036535235750221 - 6.200197140723483e-5im, 0.1848014415135279 + 0.00018981175865374902im], 3.962204377720159)
 ([-4.503345903954037e-15 + 5.126149114067468e-17im, 8.070569843533215e-15 - 6.681107940528965e-18im], 3.9999999661436303)

julia> poly_all_solutions(:heuristic, res)
1-element Vector{Tuple{Vector{ComplexF64}, Float64}}:
 ([8.786670321838077e-19 - 0.0im, -6.275532682330681e-16 - 0.0im], 3.9999999661436303)

julia> optimality_certificate(res)
:Unknown
```
Note that the solution extraction algorithm in principle also works in the complex case even though the moment matrix is no
longer of Hankel form; the theory is powerful enough to handle this "minor detail." The built-in heuristic will still try
to find good solutions and can sometimes do so even in the case of multiple solutions if they only differ in the phase of
variables. However, as in the real case, there is no guarantee that the solutions can be decomposed in atomic measures, and
therefore the extraction or certification may also fail, which is shown here.