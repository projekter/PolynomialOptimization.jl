export poly_problem, poly_optimize, EqualityMethod, emCalculateGröbner, emCalculateGroebner, emAssumeGröbner,
    emAssumeGroebner, emInequalities, emSimple

# Internal representation of the type of a constraint.
@enum PolyConstraintType pctEqualitySimple pctEqualityGröbner pctEqualityNonneg pctNonneg pctPSD
struct PolyOptConstraint{P,M}
    type::PolyConstraintType
    constraint::Union{P,<:AbstractMatrix{P}}
    basis::AbstractVector{M}
end

"""
    PolyOptProblem

The basic structure for a polynomial optimization problem.

Generate this type using [`poly_problem`](@ref); perform optimizations by constructing an appropriate
[`SparseAnalysisState`](@ref).

See also [`poly_problem`](@ref), [`poly_optimize`](@ref), [`SparsityNone`](@ref), [`SparsityCorrelative`](@ref),
[`SparsityTermBlock`](@ref), [`SparsityTermCliques`](@ref), [`SparsityCorrelativeTerm`](@ref).
"""
struct PolyOptProblem{P,M,V,GB,MV<:AbstractVector{M}}
    objective::P
    prefactor::P
    variables::Vector{V}
    var_map::Dict{V,Int}
    degree::Int
    basis::MV
    constraints::Vector{PolyOptConstraint{P,M}}
    gröbner_basis::GB
    complex::Bool
    last_moments::Dict{<:Union{MonomialComplexContainer{M},M},Float64}

    function PolyOptProblem{P,M,V,GB,B}(objective::P, prefactor::P, variables::AbstractVector{V}, var_map::Dict{V,Int},
        degree::Int, basis::B, constraints::Vector{PolyOptConstraint{P,M}}, gröbner_basis::GB,
        complex::Bool) where {P,M,V,GB,B<:AbstractVector{M}}
        @assert(V <: variable_union_type(P) && M <: monomial_type(P))
        return new{P,M,V,GB,B}(objective, prefactor, variables, var_map, degree, basis, constraints, gröbner_basis,
            complex, Dict{complex ? MonomialComplexContainer{M} : M,Float64}())
    end
end

function Base.show(io::IO, m::MIME"text/plain", p::PolyOptProblem)
    print(io, p.complex ? "Complex" : "Real", "-valued polynomial optimization hierarchy of degree ", p.degree, " in ",
        length(p.variables), " variable(s)\nObjective: ")
    show(io, m, p.objective)
    if !isone(p.prefactor)
        print(io, "\nObjective was scaled by the prefactor ")
        show(io, m, p.prefactor)
    end
    (p.gröbner_basis isa EmptyGröbnerBasis) ||
        print(io, "\nEquality constraints are modeled using Gröbner basis methods (basis length: ", length(p.gröbner_basis),
            ")")
    if !isempty(p.constraints)
        print(io, "\n", length(p.constraints), " constraints")
        len = ceil(Int, log10(length(p.constraints)))
        for (i, constr) in enumerate(p.constraints)
            print(io, "\n", lpad(i, len, "0"), ": ")
            if constr.type == pctEqualitySimple || constr.type == pctEqualityGröbner
                print(io, "0 = ")
            elseif constr.type == pctEqualityNonneg
                print(io, "(0 ≤ x) ∧ (0 ≤ -x) for x = ")
            elseif constr.type == pctNonneg
                print(io, "0 ≤ ")
            elseif constr.type == pctPSD
                print(io, "0 ⪯ ")
            else
                @assert(false)
            end
            show(io, m, constr.constraint)
        end
    end
    print(io, "\nSize of full basis: ", length(p.basis))
end

MultivariatePolynomials.polynomial_type(::Union{PolyOptProblem{P},Type{PolyOptProblem{P}}}) where {P} = P
MultivariatePolynomials.polynomial_type(::Union{PolyOptProblem{P},Type{PolyOptProblem{P}}}, T) where {P} = polynomial_type(P, T)
MultivariatePolynomials.monomial_type(::Union{PolyOptProblem{P,M},Type{PolyOptProblem{P,M}}}) where {P,M} = M
MultivariatePolynomials.variable_union_type(::Union{PolyOptProblem{P,M,V},Type{PolyOptProblem{P,M,V}}}) where {P,M,V} = V
"""
    variables(prob::Union{PolyOptProblem,SparseAnalysisState})

Returns the variables associated to a given polynomial optimization problem. This defines the order in which solutions are
returned.

See also [`poly_optimize`](@ref), [`poly_solutions`](@ref), [`poly_all_solutions`](@ref).
"""
MultivariatePolynomials.variables(prob::PolyOptProblem) = prob.variables
"""
    nvariables(prob::Union{PolyOptProblem,SparseAnalysisState})

Returns the number of variables associated to a given polynomial optimization problem. This defines the order in which
solutions are returned.

See also [`poly_optimize`](@ref), [`poly_solutions`](@ref), [`poly_all_solutions`](@ref).
"""
MultivariatePolynomials.nvariables(prob::PolyOptProblem) = length(prob.variables)
"""
    degree(prob::Union{PolyOptProblem,SparseAnalysisState})

Returns the degree associated with a polynomial optimization problem.

See also [`poly_problem`](@ref).
"""
MultivariatePolynomials.degree(prob::PolyOptProblem) = prob.degree
"""
    isreal(prob::Union{PolyOptProblem,SparseAnalysisState})

Returns whether a given polynomial optimization problem contains only real-valued variables or also complex ones.
"""
Base.isreal(prob::PolyOptProblem) = !prob.complex
"""
    length(prob::PolyOptProblem)

Return the length of the full basis associated with a polynomial optimization problem.
"""
Base.length(prob::PolyOptProblem) = length(prob.basis)

struct EmptyGröbnerBasis{T} end
Base.length(::EmptyGröbnerBasis) = 0
Base.eachindex(::EmptyGröbnerBasis) = Int[]
Base.isempty(::EmptyGröbnerBasis) = true
Base.iterate(::EmptyGröbnerBasis) = nothing
Base.getindex(::EmptyGröbnerBasis{T}, a::AbstractVector) where {T} = (@assert(isempty(a)); T[])
Base.eltype(::EmptyGröbnerBasis{T}) where {T} = T

Base.rem(p::AbstractPolynomialLike, ::EmptyGröbnerBasis) = polynomial(p)

"""
    EqualityMethod

Defines how equality constraints are internally handled.
- `emCalculateGröbner` (or `emCalculateGroebner`) assumes that some arbitrary equality constraints are passed. A Gröbner basis
  is calculated and used, and all calculations will be done modulo the ideal; the basis will be chosen as the standard
  monomials. This increases preprocessing time, but can reduce the size of the optimization program itself.
- `emAssumeGröbner` (or `emAssumeGroebner`) works as `emCalculateGröber`, but assumes that the given equality constraints
  already constitute a Gröbner basis. No check of this property is performed. If constraints of this type are mixed with
  `emCalculateGröbner`, a new Gröbner basis has to be calculated anyway; it is recommended to put all known Gröbner basis
  elements before the calculated ones to speed up calculations.
- `emSimple` just adds the equality constraints to the problem without apply any Gröbner-based methods and does not exploit
  ideal reduction. The preprocessing time can be greatly simplified.
- `emInequalities` transparently rewrites the equality constraint in terms of two inequality constraints. No ideal methods are
  used whatsoever.

See also [`poly_problem`](@ref).
"""
@enum EqualityMethod emCalculateGröbner emAssumeGröbner emSimple emInequalities

const emCalculateGroebner = emCalculateGröbner
const emAssumeGroebner = emAssumeGröbner

function effective_variables_in_real(p::AbstractPolynomialLike, in)
    for v in variables(p)
        if !iszero(maxdegree(p, v)) && !(v ∈ in)
            return false
        end
    end
    return true
end

function effective_variables_in_complex(p::AbstractPolynomialLike, in)
    for v in variables(p)
        if !iszero(maxdegree(p, v)) && !(ordinary_variable(v) ∈ in)
            return false
        end
    end
    return true
end

subbasis(b::EmptyGröbnerBasis, _) = b
function subbasis(basis::AbstractVector{P}, variables) where {P}
    # a basis is always real (or contains only the ordinary_variables of the complex monomials, so that we always use the real
    # variant for performance reasons)
    result = FastVec{P}(buffer=length(basis))
    for b in basis
        effective_variables_in_real(b, variables) && unsafe_push!(result, b)
    end
    return finish!(result)
end

# provide a faster function in the case of DynamicPolynomials monomials, as they are stored in a more efficient way
# we need a fast check whether a vector of exponents is nonzero on some subset of indices. This generates pretty much optimal
# assembler code, and it is based on the assumption that ind is not empty.
@inline function unsafe_any_positive(vec::AbstractVector{N}, ind::AbstractVector{<:Integer}) where {N<:Number}
    len = length(ind)
    @inbounds while true
        vec[ind[len]] > zero(N) && return true
        iszero(len -= 1) && return false
    end
end

function subbasis(basis::DynamicPolynomials.MonomialVector{VV,VM}, variables::Vector{V}) where {VV,VM,V<:DynamicPolynomials.Variable{VV,VM}}
    isempty(variables) && return typeof(basis)(basis.vars, Vector{Int}[])
    @assert issorted(variables, rev=true)
    isempty(basis) && return basis

    # Efficiently get the indices of all the variables of basis that are _not_ present within variables. Note that both vectors
    # are sorted reversely, but the need not be any subset relation. For efficiency reasons, we populate indices in reverse
    # order.
    indices = FastVec{Int}(buffer=length(variables))
    i, j = length(basis.vars), length(variables)
    @inbounds bv, v = basis.vars[i], variables[j]
    @inbounds while true
        if bv == v
            iszero(i -= 1) && break
            if iszero(j -= 1)
                append!(indices, i:-1:1)
                break
            end
            bv, v = basis.vars[i], variables[j]
        elseif bv > v
            if iszero(j -= 1)
                append!(indices, i:-1:1)
                break
            end
            v = variables[j]
        else
            push!(indices, i)
            iszero(i -= 1) && break
            bv = basis.vars[i]
        end
    end
    isempty(indices) && return basis

    # Indeed, we do not reduce the number of the variables or the size of the coefficient vectors - in this way, we can just
    # re-use their reference instead of allocating new ones.
    length(indices) == length(basis.vars) && return typeof(basis)(basis.vars, Vector{Int}[])
    exps = FastVec{Vector{Int}}(buffer=length(basis.Z))
    for b in basis.Z
        unsafe_any_positive(b, indices) || unsafe_push!(exps, b)
    end
    return typeof(basis)(basis.vars, finish!(exps))
end

squarebasis(basis::AbstractVector{M}, ::EmptyGröbnerBasis) where {M<:AbstractMonomialLike} = basis .^ 2
squarebasis(basis::DynamicPolynomials.MonomialVector, ::EmptyGröbnerBasis) = typeof(basis)(basis.vars, 2 .* basis.Z)
squarebasis(basis::AbstractVector{M}, gröbner_basis) where {M<:AbstractMonomialLike} =
    merge_monomial_vectors(monomials.(rem.(squarebasis(basis, EmptyGröbnerBasis{polynomial_type(M)}()), (gröbner_basis,))))

@doc raw"""
    poly_problem(objective, variables, degree; zero=[], nonneg=[], psd=[], custom_basis=[], perturbation=0.,
        equality_method=emSimple, add_gröbner=false, factor_coercive=1, perturbation_coefficient=0., perturbation_form=0,
        noncompact=(0., 0), tighter=false, verbose=False)

Analyze a polynomial optimization problem and return a [`PolyOptProblem`](@ref) that can be used for sparse analysis and
optimization.

# Arguments
## Problem formulation
- `objective::AbstractPolynomial`: the objective that is to be minimized. Note that all other polynomials will be cast to the
  same type as the objective; so make sure that the coefficient type is correct!
- `degree::Int`: the degree of the Lasserre relaxation, which must be larger or equal to half of the (complex) degree of all
  polynomials that are involved. Set this value to `0` in order to automatically determine the minimum required degree (this
  may fail when `tighter` is set to `true`, and it will also not detect degree reductions due to Gröbner basis methods).
- `zero::AbstractVector{<:AbstractPolynomialLike}`: a vector with all polynomials that should be constrained to be zero. In the
  real-valued case, this is implemented via Gröbner basis methods. In the complex case, these constraints are converted into
  two inequality constraints; hence, the values of the polynomials must always be real-valued.
- `nonneg::AbstractVector{<:AbstractPolynomialLike}`: a vector with all polynomials that should be constrained to be
  nonnegative. The values of the polynomials must always be real-valued.
- `psd::AbstractVector{<:AbstractMatrix{<:AbstractPolynomialLike}}`: a vector of matrices that should be constrainted to be
  positive semidefinite. The matrices must be symmetric/hermitian.
## Problem representation
- `custom_basis::AbstractVector{<:AbstractTermLike}`: this allows to overwrite the default full basis of (standard) monomials
  by another choice, if it is known in advance that certain monomials will not be needed (for example, for unconstrained
  optimizations, the Newton polytope may be used - see the second version of this method to do this automatically).
- `equality_method::Union{EqualityMethod, <:AbstractVector{EqualityMethod}}`: this allows to overwrite the way how equality
  constraints are handled. It is not possible to use Gröbner-basis related methods for equality constraints containing complex
  variables. It is supported to use a different method for every equality constraint.
  Note that all constraints that use Gröbner basis methods will be removed; and in the end, every element of the Gröbner basis
  will be appended to the list of constraints only if `add_gröbner=true`. Be aware of this when inspecting the task or
  selectively iterating single constraints.
- `add_gröbner::Bool`: if set to true, this adds all elements of the Gröbner basis as equality constraints to the problem. Not
  doing so will potentially allow for better sparsity and faster solutions. The optimum value will not change, regardless of
  this parameter; however, the extracted solutions might not respect the equality constraints if the Gröbner basis elements are
  not added and there are solutions with the same optimum value that do not obey the constraints.
  This parameter can also be called `add_groebner`; if any of those two is `true`, it is assumed that this was the intended
  value. Default is `false`.
## Problem modification
### For unique solution extraction
- `perturbation::Union{Float64, <:AbstractVector{Float64}}`: adds a random linear perturbation with an absolute value not
  greater than this value to the objective for every variable (or, in the vector case, with different magnitudes for each
  variable). This will ensure that the result is unique and hence solution extraction will always work, at the cost of
  potentially reducing sparsity and slightly changing the actual result.
### For noncompact sets
The following four parameters allow to automatically modify the problem according to a strategy by
[Mai, Lasserre, and Magron](https://doi.org/10.1007/s10107-021-01634-1) that was mainly developed for noncompact semialgebraic
sets. It will modify the objective to

```math
\mathrm{factor\_coercive} \bigl(
    \mathrm{objective} + \mathrm{perturbation\_coefficient} \cdot \mathrm{perturbation\_form}
  \bigr)\text.
```

Usually, `perturbation_form` and `factor_coercive` are both given by ``1 + \lVert\mathrm{variables}\rVert^2``.
If `perturbation_coefficient` is strictly positive, then for almost all degrees, the optimal value of the modified problem is
then in ``\bigl[f_{\mathrm{opt}}, f_{\mathrm{opt}} + \mathrm{perturbation\_coefficient} \cdot \mathrm{perturbation\_form}^{d_\mathrm o}(x_{\mathrm{opt}})\bigr]``.
Often, this even works for a strictly zero coefficient (relatively generic conditions were found by
[Huang, Nie, and Yuan](https://doi.org/10.1007/s10107-022-01878-5)).
Note that when modifying a problem in such a way, all sparsity methods provided by this package will be useless.
- `factor_coercive::AbstractPolynomial`: Let ``k`` be divisible by ``2r``, then this must be the dehomogenization of a coercive
  positive form in ``n+1`` variables of degree ``2r`` to the power ``k/(2r)``.
  Be aware that using this parameter will multiply the objective by another polynomial and therefore require a higher total
  relaxation degree to model the problem!
- `perturbation_coefficient::Float64`: a nonnegative prefactor that determines the strength of the perturbation and whose
  inverse dictates the scaling of a sufficient lower degree bound that guarantees optimality.
- `perturbation_form::AbstractPolynomial`: must be a the dehomogenization of a positive form in `n+1` variables of degree
  `2(1+degree(objective)÷2)`.
- `noncompact::Tuple{Real,Int}`, now called ``(\epsilon, k)``: this is a shorthand that will set the previous three
  parameters to their standard values, `factor_coercive=(1 + sum(variables.^2))^k`, `perturbation_coefficient` to the value
  passed to this parameter, and `perturbation_form=(1 + sum(variables.^2))^maxhalfdegree(objective)`.
  Be aware that using this parameter will multiply the objective by another polynomial and therefore require a higher total
  relaxation degree to model the problem!
### For convergence at earlier levels
[Nie](https://doi.org/10.1007/s10107-018-1276-2) provides a way to add additional constraints based on optimality conditions to
the problems. This can speed up or make possible convergence at all. However, not every problem can be tightened in such a way,
and sometimes, tightening might also increase the minimal degree required to optimize the problem.
Note that the problem will end up with more equality and inequality constraints than originally entered. The augmentation not
change the solution of the original problem in case the minimum is attained at a critical point; if it is not, tightening will
lead to missing this minimum.
- `tighter::Bool`: if set to `true`, tries to automatically construct constraints using Nie's method. Note that the algorithm
  internally needs to create lots of dense polynomials of appropriate degrees before solving for the coefficients. It is
  therefore possible that for larger problems, this can take a very long time.
  Automatic tightening requires the Mosek solver to be installed.
  The automatic detection of the minimal degree may not work together with this procedure. It will only consider the degrees
  required before tightning, not the increase in degree that may occur during the process.
  This parameter can also be called `tighten`; if any of those two is `true`, it is assumed that this was the intended value.
## Progress monitoring
- `verbose::Bool`: if set to true, information about the current state of the method is printed; this may be useful for large
  and complicated problems whose construction can take some time (e.g., if Gröbner bases are calculated or a tightening process
  is requested).

# AbstractAlgebra
In case `AbstractAlgebra` is loaded, the polynomials may instead be of the time `MPolyElem`, if they all belong to the same
ring. In this case, all polynomials in `custom_basis` must be monomials, and `equality_method` cannot be `emAssumeGröbner`.


    poly_problem(objective; perturbation=0., verbose=false, preprocess=true)

Analyze an unconstrained polynomial optimization problem and return a [`PolyOptProblem`](@ref) that can be used for sparse
analysis and optimization. This differs from using the full version of this method in that the Newton polytope is calculated in
order to automatically determine a suitable basis. The keyword arguments are passed through to [`newton_polytope`](@ref).

See also [`PolyOptProblem`](@ref), [`poly_optimize`](@ref), [`SparsityNone`](@ref), [`SparsityCorrelative`](@ref),
[`SparsityTermBlock`](@ref), [`SparsityTermCliques`](@ref), [`SparsityCorrelativeTerm`](@ref), [`newton_polytope`](@ref),
[`EqualityMethod`](@ref).
"""
function poly_problem(objective::P, degree::Int;
    zero::AbstractVector{<:AbstractPolynomialLike}=P[],
    nonneg::AbstractVector{<:AbstractPolynomialLike}=P[],
    psd::AbstractVector{<:AbstractMatrix{<:AbstractPolynomialLike}}=AbstractMatrix{P}[],
    custom_basis::AbstractVector{<:AbstractTermLike}=monomial_type(P)[],
    perturbation::Union{Float64,<:AbstractVector{Float64}}=0.,
    equality_method::Union{EqualityMethod,<:AbstractVector{EqualityMethod}}=emSimple,
    add_gröbner::Bool=false, add_groebner::Bool=false,
    factor_coercive::AbstractPolynomialLike=one(P), perturbation_coefficient::Float64=0.,
    perturbation_form::AbstractPolynomialLike=Base.zero(P), noncompact::Tuple{Real,Integer}=(0.,0),
    tighter::Bool=false, tighten::Bool=false, verbose::Bool=false) where {P<:AbstractPolynomialLike}
    T = polynomial_type(objective)
    @assert all(coefficient.(custom_basis) .== 1)
    return poly_problem(convert(T, objective), degree, convert(Vector{T}, zero),
        convert(Vector{T}, nonneg), convert(Vector{Matrix{T}}, psd),
        convert(Vector{monomial_type(T)}, custom_basis), perturbation, equality_method, add_gröbner || add_groebner,
        convert(T, factor_coercive), perturbation_coefficient, convert(T, perturbation_form), noncompact, tighter || tighten;
        verbose)
end

function poly_problem(objective::P; perturbation::Union{Float64,<:AbstractVector{Float64}}=0.0,
    verbose::Bool=false, preprocess::Bool=true) where {P}
    vars = variables(objective)
    @assert(all(isreal, vars))
    # no constraints, so we use the Newton polytope as the basis
    deg = maxdegree(objective) ÷ 2
    basis = newton_polytope(:Mosek, objective; verbose, preprocess)
    return poly_problem(objective, deg, P[], P[], AbstractMatrix{P}[], basis, perturbation; verbose)
end

function poly_problem(objective::P, degree::Int, zero::AbstractVector{P}, nonneg::AbstractVector{P},
    psd::AbstractVector{<:AbstractMatrix{P}}, custom_basis::AbstractVector{M},
    perturbation::Union{Float64,<:AbstractVector{Float64}}=0.,
    equality_method::Union{EqualityMethod,<:AbstractVector{EqualityMethod}}=emSimple, add_gröbner::Bool=false,
    factor_coercive::P=one(P), perturbation_coefficient::Float64=0., perturbation_form::P=Base.zero(P),
    noncompact::Tuple{Real,Integer}=(0.,0), tighter::Bool=false; verbose::Bool=false) where {P,M}
    @assert(M <: monomial_type(P))
    @assert(degree ≥ 0 && perturbation_coefficient ≥ 0 && noncompact[2] ≥ 0)

    variables = union(MultivariatePolynomials.variables(objective), MultivariatePolynomials.variables.(zero)...,
        MultivariatePolynomials.variables.(nonneg)..., MultivariatePolynomials.variables.(psd)...)
    complex = !all(isreal, variables)
    complex && filter!(!isconj, variables) # we only consider the "true" variables, not conjugates
    oscar = !(P <: AbstractPolynomialLike)
    complex && @assert(!(oscar || tighter ||
                         any(isconj, variables) || any(isrealpart, variables) || any(isimagpart, variables)))
    if equality_method isa EqualityMethod
        equality_method = fill(equality_method, length(zero))
    else
        @assert(length(equality_method) == length(zero))
    end
    if !iszero(noncompact[1])
        perturbation_coefficient = noncompact[1]
        # complex case? which kind of degree?
        if isempty(zero) && isempty(nonneg) && !tighter
            perturbation_form = (1 + variables' * variables)^ceil(Int, maxdegree(objective)/2)
        else
            perturbation_form = (1 + variables' * variables)^(maxdegree(objective)÷2 +1)
        end
    end
    if !iszero(noncompact[2])
        factor_coercive = (1 + variables' * variables)^noncompact[2]
    end
    if !iszero(perturbation_coefficient)
        objective += perturbation_coefficient * perturbation_form
    end
    if !isone(factor_coercive)
        objective *= factor_coercive
    end
    if complex
        @assert(imag(objective) == 0)
        @assert(imag(factor_coercive) == 0)
        @assert(all(constr -> imag(constr) == 0, nonneg))
        @assert(all(constr -> constr' == constr, psd))
        @inbounds for i in 1:length(zero)
            if equality_method[i] == emCalculateGröbner
                @assert(all(isreal, effective_variables(zero[i])))
            elseif equality_method[i] == emAssumeGröbner
                @assert(!oscar)
                @assert(all(isreal, effective_variables(zero[i])))
            end
        end
    else
        @assert(all(constr -> transpose(constr) == constr, psd))
    end
    if perturbation isa AbstractVector || !iszero(perturbation)
        @verbose_info("Constructing perturbation")
        # To extract the solution, we add small random linear terms to the objective. This lifts potential degeneracies and
        # allows for an easy extraction of the solution.
        if complex
            # the objective must not be complex-valued. We also don't want to "just" introduce an absolute-value
            # dependency, which would not lift phase degeneracies. So we introduce real and imaginary part
            # perturbations, but in the language of the conjugates.
            objective += sum(((rand(Float64, length(variables)) .- 0.5) .* perturbation) .* (variables .+ conj(variables))) +
                sum((1im .* (rand(Float64, length(variables)) .- 0.5) .* perturbation) .* (variables .- conj(variables)))
        else
            objective += sum(((rand(Float64, length(variables)) .- 0.5) .* perturbation) .* variables)
        end
    end
    new_p = promote_type(P, polynomial_type(objective)) # int polynomial may now be float polynomial due to perturbation
    if tighter
        # Tightening uses numerical solvers, requires machine precision
        new_p = polynomial_type(new_p, Float64)
        objective = convert(new_p, objective)
    end
    # all these may be empty, leading to Any[]
    zero = AbstractVector{new_p}(convert.((new_p,), zero))
    nonneg = AbstractVector{new_p}(convert.((new_p,), nonneg))
    psd = AbstractVector{Matrix{new_p}}(convert.((Matrix{new_p},), psd))
    factor_coercive = convert(new_p, factor_coercive)

    if isempty(zero)
        gröbner_basis = EmptyGröbnerBasis{P}()
    else
        zero_gröbner_fv = FastVec{Int}(buffer=length(zero))
        # we need to find the indices of all the constraints for which we want to calculate the Gröbner basis. However, if some
        # constraints are already known to be a GB (= are part of the GB), then those must come first in the list.
        have_calc = false
        for (i, method) in enumerate(equality_method)
            if method == emCalculateGröbner
                have_calc = true
                unsafe_push!(zero_gröbner_fv, i)
            elseif method == emAssumeGröbner
                @assert(!have_calc, "If some parts of the Gröbner basis are already known while others must still be calculated, the known parts must precede the unknown ones.")
                unsafe_push!(zero_gröbner_fv, i)
            end
        end
        zero_gröbner = finish!(zero_gröbner_fv)
        if !have_calc
            # everything is known
            @inbounds gröbner_basis = zero[zero_gröbner]
        elseif !oscar
            # mix known with unknown parts. Everything must be part of the calculation, even the known ones.
            @verbose_info("Starting Gröbner basis calculation (SemialgebraicSets)")
            @inbounds gröbner_basis = SemialgebraicSets.gröbner_basis(zero[zero_gröbner])
            @verbose_info("Gröbner basis calculation completed")
        else
            @verbose_info("Starting Gröbner basis calculation (Oscar)")
            @inbounds gröbner_basis = Oscar.ideal(zero[zero_gröbner]) # as SemialgebraicSets is the standard, the ideal is what we call gröbner_basis
            gb = Oscar.groebner_basis(gröbner_basis) # should really be called groebner_basis!, this is just a precomputation
            @verbose_info("Gröbner basis calculation completed")
        end
        if isempty(gröbner_basis)
            gröbner_basis = EmptyGröbnerBasis{P}()
        else
            @verbose_info("Reducing all polynomials modulo the ideal")
            gröbner_basis = AbstractVector{new_p}(convert.((new_p,), gröbner_basis))
            objective = convert(new_p, rem(objective, gröbner_basis))
            @inbounds for i in eachindex(zero)
                zero[i] = rem(zero[i], gröbner_basis)
            end
            @inbounds for i in eachindex(nonneg)
                nonneg[i] = rem(nonneg[i], gröbner_basis)
            end
            @inbounds for i in eachindex(psd)
                psd[i] = rem.(psd[i], (gröbner_basis,))
            end
            factor_coercive = convert(new_p, rem(factor_coercive, gröbner_basis))
        end
    end
    # We can only now check whether the specified degree was actually valid, as all the Gröbner reductions might have reduced
    # the degree (and the tightening might have increased).
    degrees_eqs = maxhalfdegree.(zero)
    degrees_nonnegs = maxhalfdegree.(nonneg)
    degrees_psds = convert.(Int, maxhalfdegree.(psd)) # somehow, this is the only expression that gives Any[] when empty
    degrees_gröbner = add_gröbner ? maxhalfdegree.(gröbner_basis) : Int[]
    mindeg = max(maxhalfdegree(objective), maximum(degrees_eqs, init=0), maximum(degrees_nonnegs, init=0),
        maximum(degrees_psds, init=0), maximum(degrees_gröbner, init=0))
    if iszero(degree)
        degree = mindeg
        @info("Automatically selecting minimal degree $degree for the relaxation")
    end
    if degree < mindeg
        @warn("The minimum required degree for the relaxation hierarchy is $mindeg, but $degree was given. Should we use the minimum degree [Y], abort [N], or use different degree altogether [number]?")
        while true
            action = readline()
            if action == "Y" || action == "y"
                degree = mindeg
                break
            elseif action == "N" || action == "n"
                error("The problem construction failed due to a too small degree.")
            else
                degree = tryparse(Int, action)
                if isnothing(degree)
                    @warn("The specified value was neither Y, N, nor any number. Please retry.")
                elseif degree < mindeg
                    @warn("The specified degree must be at least $mindeg. Please retry.")
                else
                    break
                end
            end
        end
    end
    if tighter
        @assert(!complex, "Tightening is currently only implemented for real-valued problems.")
        @verbose_info("Beginning tightening process: Constructing the matrix C out of the constraints")
        zero_len = length(zero)
        nonneg_len = length(nonneg)
        tighten!(objective, variables, degree, zero, nonneg, equality_method; verbose)
        # tightening will lead to new zero constraints (of type simple) and new inequality constraints. We should take care of
        # our Gröbner stuff with the new constraints as well.
        if !(gröbner_basis isa EmptyGröbnerBasis)
            @verbose_info("Reducing new polynomials modulo the ideal")
            @inbounds for i in zero_len+1:length(zero)
                zero[i] = rem(zero[i], gröbner_basis)
            end
            @inbounds for i in nonneg_len+1:length(nonneg)
                nonneg[i] = rem(nonneg[i], gröbner_basis)
            end
        end
        @verbose_info("Tightening completed")
        @inbounds append!(degrees_eqs, maxhalfdegree.(@view(zero[zero_len+1:end])))
        @inbounds append!(degrees_nonnegs, maxhalfdegree.(@view(nonneg[nonneg_len+1:end])))
        mindeg = max(mindeg, maximum(@view(degrees_eqs[zero_len+1:end]), init=0),
            maximum(@view(degrees_nonnegs[nonneg_len+1:end]), init=0))
        if degree < mindeg
            @warn("The mimimum required degree for the relaxation hierarchy has increased due to tightening and is now $mindeg, but $degree was given. Should we use the minimum degree [Y], abort [N], or use different degree altogether [number]?")
            while true
                action = readline()
                if action == "Y" || action == "y"
                    degree = mindeg
                    break
                elseif action == "N" || action == "n"
                    error("The problem construction failed due to a too small degree.")
                else
                    degree = tryparse(Int, action)
                    if isnothing(degree)
                        @warn("The specified value was neither Y, N, nor any number. Please retry.")
                    elseif degree < mindeg
                        @warn("The specified degree must be at least $mindeg. Please retry.")
                    else
                        break
                    end
                end
            end
        end
    end

    @verbose_info("Constructing basis")
    if isempty(custom_basis)
        # in our simple case, the ideal can only contain real-valued polynomial variables; hence, it does not matter whether we
        # take the original monomial or its conjugate to check ideal membership
        if gröbner_basis isa EmptyGröbnerBasis
            basis = monomials(oscar ? Oscar.parent(objective) : sort(variables, rev=true), 0:degree)
        elseif !oscar
            leading_terms_ideal = leading_monomial.(gröbner_basis)
            basis = monomials(sort(variables, rev=true), 0:degree, b -> !divides(leading_terms_ideal, b))
        else
            # We could construct the leading ideal and do an ideal membership check, but this appears to be much slower than
            # simple division checks.
            leading_terms_ideal = leading_monomial.(gb)
            basis = monomials(Oscar.parent(objective), 0:degree, b -> !divides(leading_terms_ideal, b))
        end
    else
        basis = monomial_vector(convert.(M, custom_basis))
    end
    @verbose_info("Constructing constraint bases")
    constr_bases_unique = Dict{Int,AbstractVector{M}}(
        deg => basis[maxdegree.(basis).≤degree-deg] # we can use the ordinary maxdegree, as all monomials are un-conjugated
        for deg in union!(Set(degrees_eqs), degrees_nonnegs, degrees_psds, degrees_gröbner)
    )
    @verbose_info("Assembling problem")
    @inbounds constraints = PolyOptConstraint{new_p,M}[
        (PolyOptConstraint{new_p,M}(
            method == emInequalities ? pctEqualityNonneg : pctEqualitySimple, constr,
            constr_bases_unique[deg]
        ) for (constr, method, deg) in zip(zero, equality_method, degrees_eqs) if method ∈ (emSimple, emInequalities))...,
        (PolyOptConstraint{new_p,M}(
            pctEqualityGröbner, constr, constr_bases_unique[deg]
        ) for (constr, deg) in zip(gröbner_basis, degrees_gröbner))...,
        (PolyOptConstraint{new_p,M}(
            pctNonneg, constr, constr_bases_unique[deg]
        ) for (constr, deg) in zip(nonneg, degrees_nonnegs))...,
        (PolyOptConstraint{new_p,M}(
            pctPSD, constr, constr_bases_unique[deg]
        ) for (constr, deg) in zip(psd, degrees_psds))...
    ]
    @inbounds return PolyOptProblem{new_p,M,variable_union_type(new_p),typeof(gröbner_basis),typeof(basis)}(objective,
        factor_coercive, variables, Dict(variables[i] => i for i in 1:length(variables)), degree, basis, constraints,
        gröbner_basis, complex)
end

"""
    poly_optimize(method, prob::PolyOptProblem, args...; kwargs...)

Directly optimize a polynomial optimization problem without using any sparsity analysis.
This is a shorthand for calling [`sparse_optimize`](@ref) on the [`SparsityNone`](@ref) object of the given problem.

All additional arguments are passed to [`sparse_optimize`](@ref).

See also [`poly_problem`](@ref), [`SparsityNone`](@ref), [`sparse_optimize`](@ref).


    poly_optimize(:LANCELOT, prob::PolyOptProblem; verbose, feastol, gradtol)

Construct a local optimizer out of the polynomial problem using LANCELOT. This returns a function which then, given a suitable
vector with an initial point (which will be mutated), finds a local optimium. The result of the function is given as a tuple
which contains the optimal value and the vector of the optimal point. Contrary to all other optimization routines, this one
does not give any global guarantees and no moments can be extracted.
"""
poly_optimize(method::Union{<:Val,Symbol}, prob::PolyOptProblem, rest...; kwrest...) =
    sparse_optimize(method, SparsityNone(prob), rest...; kwrest...)

"""
    last_moments(state::PolyOptProblem)
    last_moments(state::SparseAnalysisState)

Returns the moments dictionary that was the result of the last optimization.
Note that the results are associated with a _problem_, not with the sparse states. Therefore, the second form is merely a
convenience function: calling `last_moments` on any sparse state will always give the unique dictionary of the problem.
This function is not exported; its interface or particular return type may change without notice.

See also [`moment_matrix`](@ref).
"""
last_moments(prob::PolyOptProblem) = prob.last_moments

function poly_structure_indices(prob::PolyOptProblem, objective::Bool, zero::Union{Bool,<:AbstractSet{<:Integer}},
    nonneg::Union{Bool,<:AbstractSet{<:Integer}}, psd::Union{Bool,<:AbstractSet{<:Integer}})
    len = objective ? 1 : 0
    if zero === true
        len += count(c -> c.type ∈ (pctEqualitySimple, pctEqualityGröbner, pctEqualityNonneg), prob.constraints)
    elseif zero !== false
        len += length(zero)
    end
    if nonneg === true
        len += count(c -> c.type == pctNonneg, prob.constraints)
    elseif nonneg !== false
        len += length(nonneg)
    end
    if psd === true
        len += count(c -> c.type == pctPSD, prob.constraints)
    elseif psd !== false
        len += length(psd)
    end

    indices = FastVec{Int}(buffer=len)
    objective && unsafe_push!(indices, 1)
    eq_idx = 1
    nonneg_idx = 1
    psd_idx = 1
    for (i, constr) in enumerate(prob.constraints)
        if constr.type == pctNonneg
            (nonneg === true || nonneg_idx ∈ nonneg) && unsafe_push!(indices, i +1)
            nonneg_idx += 1
        elseif constr.type == pctPSD
            (psd === true || psd_idx ∈ psd) && unsafe_push!(indices, i +1)
            psd_idx += 1
        else
            (zero === true || eq_idx ∈ zero) && unsafe_push!(indices, i +1)
            eq_idx += 1
        end
    end
    @assert(len == length(indices))
    return finish!(indices)
end