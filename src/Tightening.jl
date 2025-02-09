const tightening_methods = Symbol[]

function default_tightening_method()
    isempty(tightening_methods) &&
        error("No tightening method is available. Load a solver package that provides such a method (e.g., Mosek)")
    @inbounds return tightening_methods[begin]
end

"""
    tighten_minimize_l1(V, spmat::SparseMatrixCSC, rhs::Vector)

Computes a solution to the underdetermined linear system `spmat * x = rhs` with minimal ℓ₁-norm.
"""
function tighten_minimize_l1 end

function tighten!(method::Val, objective::P, variables::AbstractVector{V}, zero::AbstractVector{P}, nonneg::AbstractVector{P};
    verbose::Bool=false) where {P,V}
    # We apply Nie's method of Lagrange multipliers. This means that we add two additional functions φ and ψ that are made up
    # of particular polynomials. We first need to calculate those polynomials.
    coefficient_type(P) == Float64 || error("Tightening requires Float64 coefficient type")
    ∇f = differentiate.((objective,), variables)
    if !isempty(zero) || !isempty(nonneg)
        ∇zero = [differentiate.((z,), variables) for z in zero]
        ∇nonneg = [differentiate.((n,), variables) for n in nonneg]
        Cᵀ = hcat(transpose(hcat(∇zero..., ∇nonneg...)), Diagonal([constr for zn in (zero, nonneg) for constr in zn]))
        if verbose
            println("C:")
            show(stdout, "text/plain", transpose(Cᵀ))
            println("\nFound C. Starting with ansatzes for the polynomials in L.")
            flush(stdout)
        end
        # we must now find a matrix L(x) such that L(x) C(x) = id. If all the constraints are linear, the degree bound on L
        # is size(C, 2) - rank(A) with coefficient matrix A. There is no a priori bound for higher-degree constraints, so
        # we just start with the degree bound given by the appropriate maximal degree in C and increase if necessary.
        # First check that we are not inconsistent.
        rank((variables => rand(length(variables)),) .|> Cᵀ) == size(Cᵀ, 1) || error("Tightening is not possible")
        # For efficiency reasons, as everything is col-major, we do Cᵀ * Lᵀ instead.
        # (CᵀLᵀ)ᵢⱼ = ∑ₖ Cᵀᵢₖ Lᵀₖⱼ
        n = length(variables)
        Ldeg = maximum(maxdegree.(Cᵀ))
        # define the coefficients for the polynomials in Lᵀ (we don't care about the name, but it might have to be unique
        # depending on the MP implementation).
        coeffprefix = gensym()
        @inbounds coeffsᵀ = [[similar_variable(V, Symbol(coeffprefix, "[", k, ",", j, ",", i, "]"))
                              for i in 1:binomial(n + Ldeg, n)] for k in axes(Cᵀ, 2), j in axes(Cᵀ, 1)]
        mons = [monomials(variables, [deg]) for deg in 0:Ldeg]
        @inbounds Lᵀ = [polynomial(polynomial.(coeffsᵀ[k, j]), vcat(@view(mons[1:Ldeg+1])...))
                        for k in axes(coeffsᵀ, 1), j in axes(coeffsᵀ, 2)]
        CᵀLᵀid = Cᵀ * Lᵀ - I
        # now we have our candidate for L with completely unknown coefficients. Compare the coefficients of L*C with those
        # of the identity. The matrix CᵀLᵀid must be identically zero. It is a matrix of polynomials whose coefficients are
        # linear polynomials. Therefore, we look at the list of coefficients of each entry in CᵀLᵀid and set it to zero.
        coeff_polysᵀ = coefficients.(CᵀLᵀid) # the coefficients themselves are polynomials
        # Every column in coeff_polysᵀ shares the same variables (which are listed in the corresponding column of coeffsᵀ),
        # so every column constitutes an independent system of linear equations.
        L₁ᵀ = Matrix{polynomial_type(first(Lᵀ),Float64)}(undef, length(variables), size(Lᵀ, 2))
        @verbose_info("First raw set of L is constructed. Trying to match coefficients in ", size(coeffsᵀ, 2), " columns")
        for (j, (sys, vars)) in enumerate(zip(eachcol(coeff_polysᵀ), eachcol(coeffsᵀ)))
            @verbose_info("Column ", j)
            flatvars = collect(Iterators.flatten(vars))
            deg = Ldeg
            mons2 = Dict{typeof(first(flatvars)),Int}(var => i for (i, var) in enumerate(flatvars))
            m = length(mons2)
            while true
                flsys = Iterators.flatten(sys)
                flsys_len = sum(length, sys; init=0)
                cols = Vector{SparseVector{Float64,Int}}(undef, flsys_len)
                rhs = zeros(Float64, flsys_len) # SuiteSparse can only work with machine precision
                @inbounds for (i, poly) in enumerate(flsys)
                    @assert(maxdegree(poly) ≤ 1)
                    idxs = FastVec{Int}(buffer=nterms(poly))
                    vals = FastVec{Float64}(buffer=nterms(poly))
                    for t in terms(poly)
                        mon = monomial(t)
                        if isconstant(mon)
                            rhs[i] -= coefficient(t)
                        else
                            unsafe_push!(idxs, mons2[variable(mon)])
                            unsafe_push!(vals, coefficient(t))
                        end
                    end
                    cols[i] = sparsevec(finish!(idxs), finish!(vals), m)
                end
                spmat = copy(transpose(hcat(cols...))) # todo: can we instead directly create the transposed version?
                try
                    # solution = spmat \ rhs - but the system may be underdetermined and we want the solution with the most
                    # zeros
                    local solution
                    if <(size(spmat)...)
                        # We don't just want to solve the equation, but actually find the solution with minimal cardinality
                        # (the basic solution that SPQR returns does not necessarily satisfy this criterion!). As an
                        # approximation, we minimize the ℓ₁ norm.
                        solution = tighten_minimize_l1(method, spmat, rhs)
                    else
                        solution = spmat \ rhs
                    end
                    @inbounds for i in eachindex(solution)
                        abs(solution[i]) < 1e-10 && (solution[i] = Base.zero(eltype(solution)))
                    end
                    # since spmat is typically nonsquare, \ uses a least-squares algorithm which will not raise an error if
                    # there is no proper solution. We must check this manually.
                    if norm(spmat * solution - rhs) < 1e-6
                        @verbose_info("Found proper coefficients in column ", j)
                        # now we know the values of all the coefficients, which allows us to fix the polynomial itself
                        # note that map_coefficients! would not be better, as the type of the coefficients changes
                        for i in 1:size(L₁ᵀ, 1)
                            @inbounds L₁ᵀ[i, j] = map_coefficients(coeff_poly -> coeff_poly(flatvars => solution), Lᵀ[i, j])
                        end
                        break
                    end
                catch e
                    e isa SingularException || rethrow(e)
                end
                # we did not have a sufficient degree. However, this only means that we have to up our game for this
                # particular column.
                deg += 1
                @verbose_info("Degree was insufficient: Creating ansatz with maxdegree ", deg)
                if length(mons) ≤ deg # ensure the new monomials exist
                    @assert(length(mons) == deg)
                    push!(mons, monomials(variables, [deg]))
                end
                num = binomial(n + deg -1, n -1)
                @inbounds for k in axes(Lᵀ, 1)
                    # we don't waste time in computing the index of the last existing coefficient as an offset just to get
                    # a nice name
                    new_coeffsprefix = gensym()
                    new_coeffs = [similar_variable(V, Symbol(new_coeffsprefix, "[", i, "]")) for i in 1:num]
                    append!(flatvars, new_coeffs)
                    push!(mons2, (new_coeffs .=> (m+1):(m+num))...)
                    m += num
                    # update sys
                    new_poly = polynomial(polynomial.(new_coeffs), mons[deg+1])
                    MutableArithmetics.operate!(+, Lᵀ[k, j], new_poly)
                    for (i, expr) in enumerate(@view(CᵀLᵀid[:, j]))
                        MutableArithmetics.operate!(+, expr, Cᵀ[i, k] * new_poly)
                    end
                end
                for i in axes(CᵀLᵀid, 1)
                    @inbounds sys[i] = coefficients(CᵀLᵀid[i, j])
                end
            end
        end
        if verbose
            println("Found L₁:")
            show(stdout, "text/plain", transpose(L₁ᵀ))
            println("\nNow constructing additional constraints")
            flush(stdout)
        end
        λ = transpose(L₁ᵀ) * ∇f
        peqs = @view(λ[1:length(zero)])
        pnonneg = @view(λ[length(zero)+1:end])
        # we don't know what to do with PSD constraints
        φ = [∇f[k] - sum(p * ∂zero[k] for (p, ∂zero) in zip(peqs, ∇zero); init=0.) -
            sum(p * ∂nonneg[k] for (p, ∂nonneg) in zip(pnonneg, ∇nonneg); init=0.) for k in 1:length(∇f)]
        append!(φ, pnonneg .* nonneg)
        append!(nonneg, filter(p -> !isconstant(p), pnonneg))
    else
        φ = ∇f
    end
    map_coefficients!.(c -> abs(c) < 1e-10 ? Base.zero(c) : c, φ)
    # there might be useless constraints which we can immediately drop
    filter!(p -> !isconstant(p), φ) # we could assert on p ≃ 0
    append!(zero, φ)
    return
end

tighten!(method::Symbol, args...; kwargs...) = tighten!(Val(method), args...; kwargs...)