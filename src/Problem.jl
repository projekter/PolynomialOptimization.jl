export poly_problem, poly_optimize

"""
    Problem

The basic structure that describes a polynomial optimization problem. In order to perform optimizations on this problem,
construct [`AbstractRelaxation`](@ref)s from it.
Note that the variables in a `Problem` are rewritten to internal data types, i.e., they will probably not display in the same
way as the original variables (they are simply numbered consecutively).

This type is not exported.

See also [`poly_problem`](@ref), [`poly_optimize`](@ref), [`AbstractRelaxation`](@ref).
"""
struct Problem{P<:SimplePolynomial,OV}
    objective::P
    prefactor::P
    mindegree::Int
    constr_zero::Vector{P}
    constr_nonneg::Vector{P}
    constr_psd::Vector{Matrix{P}}
    original_variables::Vector{OV}
end

"""
    variables(problem::Union{Problem,<:AbstractRelaxation})

Returns the original variables (not their internal rewrites) associated to a given polynomial optimization problem. This
defines the order in which solutions are returned. In the complex case, they do not contain conjugates.

See also [`poly_optimize`](@ref), [`poly_solutions`](@ref), [`poly_all_solutions`](@ref).
"""
MultivariatePolynomials.variables(problem::Problem) = problem.original_variables
"""
    nvariables(problem::Union{Problem,<:AbstractRelaxation})

Returns the number of variables associated to a given polynomial optimization problem. In the complex case, conjugates are not
counted.

See also [`variables`](@ref variables(::Problem)).
"""
MultivariatePolynomials.nvariables(::Problem{<:SimplePolynomial{<:Any,Nr,Nc}}) where {Nr,Nc} = Nr + Nc

const RealProblem = Problem{<:SimplePolynomial{<:Any,<:Any,0}}
const ComplexProblem = Problem{<:SimplePolynomial{<:Any,0}}

"""
    isreal(problem::Union{Problem,<:AbstractRelaxation})

Returns whether a given polynomial optimization problem contains only real-valued variables or also complex ones.
"""
Base.isreal(::RealProblem) = true
Base.isreal(::Problem) = false


function Base.show(io::IO, m::MIME"text/plain", p::Problem)
    nv = nvariables(p)
    type = isreal(p) ? "Real" : (p isa ComplexProblem ? "Complex" : "Real- and complex")
    print(io, type, "-valued polynomial optimization problem in ", nv,
        " variable", isone(nv) ? "" : "s", "\nObjective: ")
    show(io, m, p.objective)
    if !isone(p.prefactor)
        print(io, "\nObjective was scaled by the prefactor ")
        show(io, m, p.prefactor)
    end
    ∑length = length(p.constr_zero) + length(p.constr_nonneg) + length(p.constr_psd)
    matrix_io = io
    if !haskey(matrix_io, :compact)
        matrix_io = IOContext(matrix_io, :compact => true)
    end
    if !iszero(∑length)
        len = floor(Int, log10(∑length)) +1
        i = 1
        for (c, t, s) in ((p.constr_zero, "equality", " = 0"),
                          (p.constr_nonneg, "nonnegative", " ≥ 0"),
                          (p.constr_psd, "semidefinite", "")) # we'll do it differently for the semidefinite case
            if !isempty(c)
                print(io, "\n", length(c), " ", t, " constraint", isone(length(c)) ? "" : "s")
                for constr in c
                    if constr isa AbstractArray
                        # Do not simply show the array. This will print the whole eltype, which makes this completely
                        # unreadable. Besides, we can do a nice indentation in this way.
                        println(io)
                        Base.print_matrix(matrix_io, constr, lpad(i, len, "0") * ": [", "  ", "] ⪰ 0")
                    else
                        print(io, "\n", lpad(i, len, "0"), ": ")
                        show(io, m, constr)
                        print(io, s)
                    end
                    i += 1
                end
            end
        end
    end
end

# We need to provide "tolerant" versions of the complex degree-related stuff. They will not raise an error when real/imaginary
# parts arise, but are functionally equivalent else. We must then make sure that real/imaginary parts are never mixed with the
# complex representation of the same variable.
function halfdegree_tolerant(t::AbstractTermLike)
    realdeg = 0
    cpdeg = 0
    conjdeg = 0
    for (var, exp) in powers(t)
        if isreal(var)
            realdeg += exp
        else
            if isconj(var)
                conjdeg += exp
            else
                cpdeg += exp
            end
        end
    end
    return div(realdeg, 2, RoundUp) + max(cpdeg, conjdeg)
end
function maxhalfdegree_tolerant(X::AbstractArray{<:AbstractTermLike})
    return isempty(X) ? 0 : maximum(halfdegree_tolerant, X, init=0)
end
maxhalfdegree_tolerant(p::AbstractPolynomialLike) = maxhalfdegree_tolerant(terms(p))
maxhalfdegree_tolerant(m::AbstractMatrix{<:AbstractPolynomialLike}, args...) =
    maximum((maxhalfdegree_tolerant(p, args...) for p in m), init=0)::Int

@doc raw"""
    poly_problem(objective; zero=[], nonneg=[], psd=[], perturbation=0.,
        factor_coercive=1, perturbation_coefficient=0., perturbation_form=0,
        noncompact=(0., 0), tighter=false, verbose=false, monomial_index_type=UInt)

Analyze a polynomial optimization problem and return a [`Problem`](@ref) that can be used for sparse analysis and optimization.

# Arguments
## Problem formulation
- `objective::AbstractPolynomial`: the objective that is to be minimized
- `zero::AbstractVector{<:AbstractPolynomialLike}`: a vector with all polynomials that should be constrained to be zero.
- `nonneg::AbstractVector{<:AbstractPolynomialLike}`: a vector with all polynomials that should be constrained to be
  nonnegative. The values of the polynomials must always be effectively real-valued (in the sense that their imaginary parts
  evaluate to zero even if no values are plugged in).
- `psd::AbstractVector{<:AbstractMatrix{<:AbstractPolynomialLike}}`: a vector of matrices that should be constrainted to be
  positive semidefinite. The matrices must be symmetric/hermitian.
## Problem representation
- `monomial_index_type::Type{<:Integer}`: internally, whatever interface of `MultivariatePolynomials` is used, the data is
  converted to the efficient [`SimplePolynomial`](@ref) representation. Every monomial is represented by a single number of the
  type given for this keyword argument. The default is usually a good choice, allowing quite large problems. For very small
  problems, the index type might be reduced (however, note that the index must be large enough to capture the monomial for
  every desired relaxation, and there will be no warning on overflow!); if the problem is extremely large, it might also be
  enlarged to `UInt128` or `BigInt`, the latter in particular with potentially severe performance and memory consumption
  issues.
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
- `perturbation_form::AbstractPolynomial`: must be the dehomogenization of a positive form in `n+1` variables of degree
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
Note that the problem will end up with more equality and inequality constraints than originally entered. The augmentation does
not change the solution of the original problem in case the minimum is attained at a critical point; if it is not, tightening
will lead to missing this minimum.
- `tighter::Union{Bool,Symbol}`: if set to a valid solver or `true` (= choose default), tries to automatically construct
  constraints using Nie's method. Note that the algorithm internally needs to create lots of dense polynomials of appropriate
  degrees before solving for the coefficients. It is therefore possible that for larger problems, this can take a very long
  time.
  For a list of supported solvers, see [the solver reference](@ref solvers_tighten).
  This parameter can also be called `tighten`; if any of those two is `true`, it is assumed that this was the intended value.
## Progress monitoring
- `verbose::Bool`: if set to true, information about the current state of the method is printed; this may be useful for large
  and complicated problems whose construction can take some time.

See also [`Problem`](@ref), [`poly_optimize`](@ref), [`Relaxation.AbstractRelaxation`](@ref).
"""
function poly_problem(objective::P;
    zero::AbstractVector{<:AbstractPolynomialLike}=P[],
    nonneg::AbstractVector{<:AbstractPolynomialLike}=P[],
    psd::AbstractVector{<:AbstractMatrix{<:AbstractPolynomialLike}}=AbstractMatrix{P}[],
    perturbation::Union{Float64,<:AbstractVector{Float64}}=0.,
    factor_coercive::AbstractPolynomialLike=one(P), perturbation_coefficient::Float64=0.,
    perturbation_form::AbstractPolynomialLike=Base.zero(P), noncompact::Tuple{Real,Integer}=(0.,0),
    tighter::Union{Bool,Symbol}=false, tighten::Union{Bool,Symbol}=false, verbose::Bool=false,
    monomial_index_type::Type{<:Integer}=UInt) where {P<:AbstractPolynomialLike}
    if tighten !== false
        tighter = tighten
    end

    #region List of variables
    vars = sort!(collect(let vars = Set(variables(objective))
        for c in zero
            union!(vars, variables(c))
        end
        for c in nonneg
            union!(vars, variables(c))
        end
        for c in psd
            union!(vars, variables(c))
        end
        # Simple polynomial does not allow for conjugates, but we need to make sure all the base variables are included
        filter!(∘(!, isconj), union!(vars, conj.(vars)))
    end), rev=true) # for better compatibility with DynamicPolynomials, we reverse-sort the variables.
    #endregion

    #region Type of the problem
    nreal = count(isreal, vars)
    ncomplex = length(vars) - nreal
    complex = !iszero(ncomplex) # this only captures the presence of complex-valued variables, not of complex-valued
                                # coefficients
    if complex
        for v in vars
            if isrealpart(v) || isimagpart(v)
                ov = ordinary_variable(v)
                if insorted(ov, vars, rev=true) || insorted(conj(ov), vars, rev=true)
                    error("A complex-valued variable must not be used simultaneously with its real/imaginary decomposition and its complex value")
                end
            end
        end
    end

    T = let checkcomplex=p -> coefficient_type(p)<:Complex
        # our solvers and other operations will require machine float. For now, this is hard-coded here.
        if complex || checkcomplex(objective) || any(checkcomplex, zero) || any(checkcomplex, nonneg) ||
            any(m -> any(checkcomplex, m), psd)
            Complex{Float64}
        else
            Float64
        end
    end
    PT = polynomial_type(objective, T)
    #endregion

    #region Validation: no imaginary parts, correct form of inputs
    if T <: Complex
        tighter !== false && error("Complex-valued problems cannot be tightened")
        imag(objective) == 0 || error("Nonvanishing imaginary part in objective")
        imag(factor_coercive) == 0 || error("Nonvanishing imaginary party in coercive factor")
        all(∘(iszero, imag), nonneg) || error("Nonvanishing imaginary in nonnegative constraint")
        all(constr -> constr' == constr, psd) || error("Nonhermitian PSD constraint")
    else
        all(constr -> transpose(constr) == constr, psd) || error("Nonsymmetric PSD constraint")
    end
    #endregion

    #region Preprocessing: perturbation, coercive factor
    objective = convert(PT, objective)
    zero = AbstractVector{PT}(convert.((PT,), zero))
    nonneg = AbstractVector{PT}(convert.((PT,), nonneg))
    psd = AbstractVector{Matrix{PT}}(convert.((Matrix{PT},), psd))
    factor_coercive = convert(PT, factor_coercive)
    if !iszero(noncompact[1])
        perturbation_coefficient = noncompact[1]
        # complex case? which kind of degree?
        if isempty(zero) && isempty(nonneg) && tighter === false
            perturbation_form = (1 + vars' * vars)^ceil(Int, maxdegree(objective)/2)
        else
            perturbation_form = (1 + vars' * vars)^(maxdegree(objective)÷2 +1)
        end
    end
    if !iszero(noncompact[2])
        factor_coercive = (1 + vars' * vars)^noncompact[2]
    end
    if !iszero(perturbation_coefficient)
        objective += perturbation_coefficient * perturbation_form
    end
    if !isone(factor_coercive)
        objective *= factor_coercive
    end
    if perturbation isa AbstractVector || !iszero(perturbation)
        @verbose_info("Constructing perturbation")
        # To extract the solution, we add small random linear terms to the objective. This lifts potential degeneracies and
        # allows for an easy extraction of the solution.
        if complex
            # the objective must not be complex-valued. We also don't want to "just" introduce an absolute-value
            # dependency, which would not lift phase degeneracies. So we introduce real and imaginary part
            # perturbations, but in the language of the conjugates.
            objective += sum(((rand(Float64, length(vars)) .- 0.5) .* perturbation) .* (vars .+ conj(vars))) +
                sum((im .* (rand(Float64, length(vars)) .- 0.5) .* perturbation) .* (vars .- conj(vars)))
        else
            objective += sum(((rand(Float64, length(vars)) .- 0.5) .* perturbation) .* vars)
        end
    end
    if tighter === true
        tighter = default_tightening_method()
    end
    #endregion

    #region Obtaining information on the minimal degree
    degrees_eqs = maxhalfdegree_tolerant.(zero)
    degrees_nonnegs = maxhalfdegree_tolerant.(nonneg)
    degrees_psds = convert.(Int, maxhalfdegree_tolerant.(psd)) # somehow, this is the only expression that gives Any[] when empty
    mindeg = max(maxhalfdegree_tolerant(objective), maximum(degrees_eqs, init=0), maximum(degrees_nonnegs, init=0),
                 maximum(degrees_psds, init=0))
    #endregion

    #region Tightening and degree adjustment
    if tighter !== false
        @verbose_info("Beginning tightening process: Constructing the matrix C out of the constraints")
        zero_len = length(zero)
        nonneg_len = length(nonneg)
        tighten!(tighter, objective, vars, zero, nonneg; verbose)
        # tightening will lead to new zero constraints (of type simple) and new inequality constraints.
        @verbose_info("Tightening completed")
        @inbounds append!(degrees_eqs, maxhalfdegree_tolerant.(@view(zero[zero_len+1:end])))
        @inbounds append!(degrees_nonnegs, maxhalfdegree_tolerant.(@view(nonneg[nonneg_len+1:end])))
        mindeg = max(mindeg, maximum(@view(degrees_eqs[zero_len+1:end]), init=0),
                     maximum(@view(degrees_nonnegs[nonneg_len+1:end]), init=0))
    end
    #endregion

    #region SimplePolynomial conversion
    @verbose_info("Converting data to simple polynomials")
    sobj = SimplePolynomial{monomial_index_type}(objective, T; vars)
    sprefactor = SimplePolynomial{monomial_index_type}(factor_coercive, T; vars)
    szero = FastVec{typeof(sobj)}(buffer=length(zero))
    for zeroᵢ in zero
        unsafe_push!(szero, SimplePolynomial{monomial_index_type}(zeroᵢ, T; vars))
    end
    snonneg = FastVec{typeof(sobj)}(buffer=length(nonneg))
    for nonnegᵢ in nonneg
        unsafe_push!(snonneg, SimplePolynomial{monomial_index_type}(nonnegᵢ, T; vars))
    end
    spsd = FastVec{Matrix{typeof(sobj)}}(buffer=length(psd))
    for psdᵢ in psd
        m = Matrix{typeof(sobj)}(undef, size(psdᵢ)...)
        for j in eachindex(psdᵢ, m)
            @inbounds m[j] = SimplePolynomial{monomial_index_type}(psdᵢ[j], T; vars)
        end
        unsafe_push!(spsd, m)
    end
    #endregion

    return Problem{typeof(sobj),eltype(vars)}(sobj, sprefactor, mindeg, finish!(szero), finish!(snonneg), finish!(spsd), vars)
end

poly_problem(problem::Problem) = problem