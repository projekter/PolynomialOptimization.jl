# TODO: we completely ignore PSD constraints for the moment
# This is different from the Newton case, where we just added all the possible monomials due to the positivstellensatz. In
# facial reduction, the coefficients actually matter; since we don't know the coefficients of the multipliers, they enter by
# means of decision variables.
# Permenter's form for facial reduction is
# minᵤ dot(w, u) s.t. gⱼ₀(x) + dot(u, gⱼ(x)) ∈ Σ ∀ j, where w is a fixed weight vector, gⱼ₀ are fixed polynomials, and gⱼ is a
#                                                     vector of fixed polynomials.
# Without constraints, our problem is
# maxₗ l s.t. obj(x) - prefactor(x) * l ∈ Σ = -minᵤ u₁ s.t. obj(x) + prefactor(x) * u₁ ∈ Σ.
# We can easily identify w = [1], g₀₀(x) = obj(x), g₀(x) = [prefactor(x)].
# With constraints, our problem is
# -minᵤ u s.t. obj(x) + prefactor(x) * u₁ = σ₀ + ∑ᵢ pᵢ(x) zeroᵢ(x) + ∑ⱼ σⱼ(x) nonnegⱼ(x) + ∑ₖ ⟨Sₖ(x), psdₖ(x)⟩
#              σ₀, σⱼ ∈ Σ, Sₖ ∈ Σmat, pᵢ ∈ ℝ[x]
# This is again re-cast as obj(x) + prefactor(x) * u₁ + ∑ᵢ ... + ∑ⱼ ... + ∑ₖ ... ∈ Σ
# Now, the decision variables u must contain the coefficients of the pᵢ, σⱼ, Sₖ, so there are many:
# w = [0, 1, 0, ...], g₀₀(x) = obj(x), g₀(x) = [prefactor(x), zeroᵢ(x) * monomials in pᵢ...,
#                                               nonnegⱼ(x) * monomials in σⱼ...,
#                                               psdₖ(x)ₘₙ * monomials in Sₖₘₙ for every entry...],
# g₁(x) = [0, 0..., monomials in σⱼ...] TODO: matrix?
struct ReductionGset{I<:Integer,P<:(IntPolynomial{<:Any,Nr,Nc,<:IntMonomialVector{Nr,Nc,I}} where {Nr,Nc}),Pr<:Problem{P},
                     R<:AbstractRelaxation{Pr}} <: AbstractMatrix{P}
    relaxation::R
    prefactor::Bool
    monsets::Vector{Set{FastKey{I}}}
    startpos::Vector{Int}

    function ReductionGset(relaxation::R) where {I<:Integer,Nr,Nc,P<:IntPolynomial{<:Any,Nr,Nc,<:IntMonomialVector{Nr,Nc,I}},
                                                 Pr<:Problem{P},R<:AbstractRelaxation{Pr}}
        g = Relaxation.groupings(relaxation)
        prob = poly_problem(relaxation)
        monsets = Vector{Set{FastKey{I}}}(undef, length(prob.constr_zero) + sum(length, g.nonnegs, init=0)) # TODO psd
        startpos = Vector{Int}(undef, length(monsets) +1)
        # For every monomial that occurs in a SOS decomposition, we need exactly one variable u. So determine all the unique
        # monomials.
        i = 0
        l = 0
        @inbounds for (constr, groupings) in zip(prob.constr_zero, g.zeros)
            i += 1
            unique_groupings = monsets[i] = Set{FastKey{I}}()
            totalsize = 0
            for groupingⱼ in groupings
                _, _, totalsizeⱼ = unique_outer_groupings(groupingⱼ, unique_groupings)
                totalsize += totalsizeⱼ
            end
            if !isreal(constr)
                totalsize *= 2
            end
            startpos[i] = l +1
            l += totalsize
        end
        @inbounds for (constr, groupings) in zip(prob.constr_nonneg, g.nonnegs)
            for groupingsⱼ in groupings
                i += 1
                monsets[i], _, totalsize = unique_outer_groupings(groupingsⱼ)
                startpos[i] = l +1
                l += totalsize
            end
        end
        startpos[end] = l +1
        new{I,P,Pr,R}(relaxation, !iszero(poly_problem(relaxation).prefactor), monsets, startpos)
    end
end

Base.size(r::ReductionGset) = ((r.prefactor ? 1 : 0) + r.startpos[end], # because startpos[end] is already one more
                               1 + sum(length, Relaxation.groupings(r.relaxation).nonnegs, init=0))

@generated _cached_zeropolynomial(p::IntPolynomial) = :(@inline; return $(zero(p)))

@generated _cached_onepolynomial(::Type{C}, m::IntMonomial{Nr,Nc}) where {C,Nr,Nc} =
    :(@inline; return IntPolynomial($([one(C)]), IntMonomialVector{Nr,Nc}(unsafe, m.e, m.index:m.index)))

# This method is extremely inefficient if a slice of the matrix should be obtained, but it serves as a proof-of-concept.
function Base.getindex(r::ReductionGset{I,<:IntPolynomial{C,Nr,Nc}}, mon::Integer, sosⱼ::Integer) where {I<:Integer,C,Nr,Nc}
    prob = poly_problem(r.relaxation)
    g = Relaxation.groupings(r.relaxation)
    e = ExponentsAll{Nr+2Nc,I}()
    if isone(sosⱼ)
        isone(mon) && return prob.objective
        mon -= 1
        if r.prefactor
            isone(mon) && return prob.prefactor
            mon -= 1
        end
        # There's a single entry in startpos/monsets for a zero constraint
        whichmonset = searchsortedlast(r.startpos, mon)
        if whichmonset ≤ length(prob.constr_zero)
            constr = prob.constr_zero[whichmonset]
            constr_real = isreal(constr)
            realgroupingsub, complexgroupingsub = constr_real ? (1, 2) : (2, 4)
            mon -= r.startpos[whichmonset] -1
            @assert(mon ≥ 1)
            monset = r.monsets[whichmonset]
            local grouping
            it = iterate(monset)
            while !isnothing(it)
                if iszero(Nc)
                    mon -= realgroupingsub
                else
                    grouping = IntMonomial{Nr,Nc}(unsafe, e, convert(I, it[1]))
                    mon -= isreal(grouping) ? realgroupingsub : complexgroupingsub
                end
                mon ≥ 1 || break
                it = iterate(monset, it[2])
            end
            if iszero(Nc)
                grouping = IntMonomial{Nr,Nc}(unsafe, e, convert(I, it[1]))
                return IntPolynomial(constr.coeffs, constr.monomials * grouping)
            end
            if iszero(mon)
                return IntPolynomial(constr.coeffs, constr.monomials * grouping)
            elseif mon == -1
                return IntPolynomial(constr.coeffs, conj(constr.monomials * grouping))
            elseif mon == -2
                return IntPolynomial(constr.coeffs, constr.monomials * IntConjMonomial(grouping))
            elseif mon == -3
                return IntPolynomial(constr.coeffs, conj(constr.monomials * IntConjMonomial(grouping)))
            else
                error("Internal error")
            end
        end
        searchconstr = whichmonset - length(prob.constr_zero)
        # There's a single entry in startpos/monsets for every grouping of a nonneg constraint
        @inbounds for (constr, groupings) in zip(prob.constr_nonneg, g.nonnegs)
            if searchconstr > length(groupings)
                searchconstr -= length(groupings)
                continue
            end
            mon -= r.startpos[whichmonset] -1
            @assert(mon ≥ 1)
            monset = r.monsets[whichmonset]
            local grouping
            it = iterate(monset)
            while !isnothing(it)
                if iszero(Nc)
                    mon -= 1
                else
                    grouping = IntMonomial{Nr,Nc}(unsafe, e, convert(I, it[1]))
                    mon -= isreal(grouping) ? 1 : 2
                end
                mon ≥ 1 || break
                it = iterate(monset, it[2])
            end
            if iszero(Nc)
                grouping = IntMonomial{Nr,Nc}(unsafe, e, convert(I, it[1]))
                return IntPolynomial(constr.coeffs, constr.monomials * grouping)
            end
            if iszero(mon)
                return IntPolynomial(constr.coeffs, constr.monomials * grouping)
            elseif mon == -1
                return IntPolynomial(constr.coeffs, conj(constr.monomials * grouping))
            else
                error("Internal error")
            end
        end
        throw(BoundsError(r, (mon + 1 + r.prefactor, sosⱼ)))
    else
        # we could do a bounds check here if sosⱼ is valid...
        isone(mon) && return _cached_zeropolynomial(prob.objective)
        mon -= 1
        if r.prefactor
            isone(mon) && return _cached_zeropolynomial(prob.objective)
            mon -= 1
        end
        whichmonset = searchsortedlast(r.startpos, mon)
        whichmonset ≤ length(prob.constr_zero) && return _cached_zeropolynomial(prob.objective)
        searchconstr = whichmonset - length(prob.constr_zero)
        sosⱼ -= 1
        @inbounds for groupings in g.nonnegs
            if searchconstr > length(groupings)
                searchconstr -= length(groupings)
                iszero(sosⱼ -= 1) && return _cached_zeropolynomial(prob.objective)
                continue
            end
            isone(sosⱼ) || return _cached_zeropolynomial(prob.objective)
            mon -= r.startpos[whichmonset] -1
            @assert(mon ≥ 1)
            monset = r.monsets[whichmonset]
            local grouping
            it = iterate(monset)
            while !isnothing(it)
                if iszero(Nc)
                    mon -= 1
                else
                    grouping = IntMonomial{Nr,Nc}(unsafe, e, convert(I, it[1]))
                    mon -= isreal(grouping) ? 1 : 2
                end
                mon ≥ 1 || break
                it = iterate(monset, it[2])
            end
            if iszero(Nc)
                grouping = IntMonomial{Nr,Nc}(unsafe, e, convert(I, it[1]))
                return _cached_onepolynomial(C, grouping)
            end
            if iszero(mon)
                return _cached_onepolynomial(C, grouping)
            elseif mon == -1
                return _cached_onepolynomial(C, conj(grouping))
            else
                error("Internal error")
            end
        end
        throw(BoundsError(r, (mon + 1 + r.prefactor, sosⱼ)))
    end
end