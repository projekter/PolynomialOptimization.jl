struct Newton{P<:Problem,G<:RelaxationGroupings} <: AbstractRelaxationSparse{P}
    parent::AbstractRelaxation{P} # no specialization
    groupings::G

    @doc """
        Newton(relaxation::AbstractRelaxation; [method,] parameters...)

    Constructs a relaxation based on the Newton halfpolytope applied to another relaxation of a polynomial optimization
    problem. This will be a superset of the largest possible representation for a given degree bound, with no negative
    consequences for finding the optimum. It can be much smaller than a dense basis, but solution reconstruction may be harder.

    When constraints are present, the polytope is determined based on the Putinar reformulation, where all constraints and the
    objective are moved to one side (comprising a new virtual objective). The groupings for the constraints are determined by
    the previous relaxation method and not affected by the Newton polytope.

    Note that for the complex-valued hierarchy, strictly speaking there is no "Newton polytope"; as the representation of
    complex-valued polynomials is unique, the process is much simpler there; still, the size reduction is accomplished by using
    `Newton`.

    The `method` determines which solver to use for determining the Newton polytope. If omitted, this will be the default
    solver (in the complex case, it must be `:complex`).
    The `parameters` are passed on to [`Newton.halfpolytope`](@ref PolynomialOptimization.Newton.halfpolytope).
    """
    function Newton(relaxation::AbstractRelaxation{P};
        method::Symbol=iszero(Nc) ? PolynomialOptimization.Newton.default_newton_method() : :complex, verbose::Bool=false,
        parameters...) where {Nr,Nc,I<:Integer,Poly<:IntPolynomial{<:Any,Nr,Nc,<:IntMonomialVector{Nr,Nc,I}},P<:Problem{Poly}}
        if !iszero(Nr) && !iszero(Nc)
            # Well, we could do this. For a polynomial ∑ᵢⱼₖ αᵢⱼₖ xⁱ zʲ z̄ᵏ, we could factor the complex valued part and then
            # apply the Newton polytope to the real factor: ∑ⱼₖ NP(∑ᵢ αᵢⱼₖ xⁱ) zʲ z̄ᵏ; and then simplify the complex valued part.
            # Not implemented at the moment.
            error("Mixing real- and complex-valued variables prevents Newton polytope methods")
        end
        iszero(Nc) || newton_method === :complex || throw(ArgumentError("Complex-valued problems require the :complex method"))
        problem = poly_problem(relaxation)
        parent = groupings(relaxation)
        basis = PolynomialOptimization.Newton.halfpolytope(method, problem.objective; zero=problem.constr_zero,
            nonneg=problem.constr_nonneg, psd=problem.constr_psd, prefactor=problem.prefactor, groupings=parent, verbose,
            parameters...)
        basis isa IntMonomialVector ||
            error("Newton polytope calculation did not give results. Were the results written to a file?")
        @verbose_info("Embedding new groupings in old")
        gentime = @elapsed(gr = embed(RelaxationGroupings(
            IntMonomialVector{Nr,Nc,I}[basis],
            parent.zeros,
            parent.nonnegs,
            parent.psds,
            parent.var_cliques
        ), parent, true))
        @verbose_info("Obtained embedding in ", gentime, " seconds")

        new{P,typeof(gr)}(relaxation, gr)
    end
end

# All the parameters are only for choosing a different method, which might impact runtime behavior; but all methods valid for
# a certain problem must always give the same result. So there's no use in doing the work all over, even if the parameters
# have changed.
Newton(relaxation::Newton; kwargs...) = relaxation

default_solution_method(::Newton) = :heuristic