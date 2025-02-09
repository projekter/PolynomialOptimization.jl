# We only accept ExponentsAll; else, we cannot guarantee that all the resulting exponents will in fact be possible. In this
# way, we don't have to check for equality of the exponents sets (or form the largest cover), which is also beneficial.
# And converting polynomials from the MP interface will yield ExponentsAll anyway.
function merge_constraints(objective::SimplePolynomial{<:Any,Nr,Nc}, zero::AbstractVector{<:SimplePolynomial{<:Any,Nr,Nc}},
    nonneg::AbstractVector{<:SimplePolynomial{<:Any,Nr,Nc}},
    psd::AbstractVector{<:AbstractMatrix{<:SimplePolynomial{<:Any,Nr,Nc}}}, groupings::RelaxationGroupings{Nr,Nc},
    verbose::Bool, need_copy::Bool) where {Nr,Nc}
    @verbose_info("Incorporating constraints into set of exponents")
    # Note that in the complex-valued case there's no mixing - i.e., no real variables. And every monomial appears once in its
    # original form, once in its conjugate. We are only interested in the "original" (whichever it is), so we discard the conjugate part to avoid even generating duplicates.
    e = ExponentsAll{Nr+2Nc,UInt}()

    if isempty(zero) && isempty(nonneg) && isempty(psd)
        # short path: just directly add the stuff from the objective. While we can get duplicates (due to the elimination of
        # conjugates), we have a clear and manageable upper bound to the size, so no need to determine duplicates on-the-fly.
        @verbose_info("No constraints found, adding all objective monomials")
        mons = monomials(objective)
        if iszero(Nc)
            if monomials(objective).e == e
                @verbose_info("Aliasing monomial indices")
                return need_copy ? copy(mons) : mons
            else
                @verbose_info("Converting monomial indices")
                return SimpleMonomialVector{Nr,0}(e, monomial_index.((e,), mons))
            end
        else
            @verbose_info("Converting monomial indices")
            candidates = FastVec{UInt}(buffer=length(mons))
            for mon in mons
                @inbounds unsafe_push!(candidates, exponents_to_index(e, KillConjugates{Nr}(exponents(mon))))
            end
            @verbose_info("Sorting and removing duplicates")
            return SimpleMonomialVector{Nr,Nc}(e, Base._groupedunique(sort!(finish!(candidates))))
        end
    end

    mons_idx_set = sizehint!(Set{FastKey{UInt}}(), length(objective))
    @verbose_info("├ objective")
    for mon in monomials(objective)
        push!(mons_idx_set, FastKey(iszero(Nc) ? monomial_index(e, mon) :
                                    exponents_to_index(e, KillConjugates{Nr}(exponents(mon)))))
    end
    # If there are constraints present, things are not so simple. We assume a Putinar certificate:
    # f ≥ 0 on {zero == 0, nonneg ≥ 0, psd ⪰ 0} ⇐ f = σ₀ + ∑ᵢ nonnegᵢ σᵢ + ∑ⱼ ⟨psdⱼ, Mⱼ⟩ + ∑ₖ zeroₖ pₖ
    #                                              where σ₀, σᵢ ∈ SOS, Mⱼ ∈ SOSmatrix, pₖ ∈ poly
    # This can simply be reformulated into f - ∑ᵢ nonnegᵢ σᵢ - ∑ⱼ ⟨psdⱼ, Mⱼ⟩ - ∑ₖ zeroₖ pₖ ∈ SOS, i.e., we can now apply Newton
    # methods to the polynomial with subtracted constraint certifiers. We get our exact form for the multipliers from
    # `groupings` - if sparsity methods have already been applied, we can exploit them.
    # Note that since the coefficients of the σᵢ, Mⱼ, and pₖ are unknowns, we don't need to ask ourselves whether some
    # cancellation may occur - we don't know. So every monomial that is present in any of the constraints, multiplied
    # by any monomial of allowed degree for the multiplier, will give rise to an additional entry in the coeffs array.
    # This can quickly become disastrous if `degree` is high but the degree of the constraints is low (as then, the
    # prefactors have lots of entries), but it is not so harmful in the other regime.

    # Out of necessity, our groupings contain abstract types. Given that the number of groupings will be manageable, but the
    # groupings themselves can be quite large, it pays off to introduce a function barrier specializing on the actual type.
    grouping_loop = @capture(((grouping, monₜ, ::Val{constr_is_real}) where {constr_is_real}) -> begin
        if !iszero(Nc)
            monₜe = KillConjugates{Nr}(exponents(monₜ))
            if !constr_is_real
                monₜec = KillConjugates{Nr}(exponents(SimpleConjMonomial(monₜ)))
            end
        end
        for (i, g₁) in enumerate(grouping)
            for g₂ in Iterators.take(grouping, i)
                if iszero(Nc)
                    push!($mons_idx_set,
                        FastKey(monomial_index($e, g₁, monₜ, g₂)),
                        FastKey(monomial_index(e, g₂, monₜ, g₁))
                    )
                elseif constr_is_real
                    # While the prefactor is g₁*conj(g₂), we drop the conjugates anyway.
                    # Note that no prefactor has a single conjugate component; and while the monomials might have them, the
                    # conjugates will come up later anyway!
                    push!(mons_idx_set,
                        FastKey(exponents_sum(e, g₁, monₜe)), # KillConjugates(g₁*monₜ*ḡ₂)
                        FastKey(exponents_sum(e, g₂, monₜe)), # KillConjugates(g₂*monₜ*ḡ₁)
                    )
                else
                    # Here, it is not guaranteed that we will see the conjugates later on.
                    push!(mons_idx_set,
                        FastKey(exponents_sum(e, g₁, monₜe)), # KillConjugates(g₁*monₜ*ḡ₂)
                        FastKey(exponents_sum(e, g₁, monₜec)), # KillConjugates(g₁*conj(monₜ)*ḡ₂)
                        FastKey(exponents_sum(e, g₂, monₜe)), # KillConjugates(g₂*monₜ*ḡ₁)
                        FastKey(exponents_sum(e, g₂, monₜec)), # KillConjugates(g₂*conj(monₜ)*ḡ₁)
                    )
                end
            end
        end
        nothing
    end)

    # 1./2. zero/nonneg constraints
    # Note that we can still exploit that the σⱼ must be made up of products of elements in the grouping: despite being
    # polynomials of degree deg(σⱼ) with unknown coefficients, _some_ of these must for sure be zero, namely those that cannot
    # be reached by combining any two of the possible coefficients that are in the valid multidegree range of σⱼ.
    for (s, constr_groupings, constrs) in (("zero", groupings.zeros, zero), ("nonnegative", groupings.nonnegs, nonneg))
        isempty(constrs) || @verbose_info("├ ", s, " constraints")
        isz = constrs === zero
        for (groupings, constrᵢ) in zip(constr_groupings, constrs)
            newbound = length(constrᵢ) * sum(∘(trisize, length), groupings)
            sizehint!(mons_idx_set, length(mons_idx_set) + (iszero(Nc) ? newbound : (isz ? 4newbound : 2newbound)))
            for t in constrᵢ
                monₜ = monomial(t)
                for grouping in groupings
                    grouping_loop(grouping, monₜ, Val(!isz))
                end
            end
        end
    end

    # 3. psd constraints
    # Those are modeled in terms of SOS matrices. Given the basis u of `degree`, an m×m-matrix M is a SOS matrix iff
    # M(x) = (u ⊗ 1ₘ)ᵀ Z (u ⊗ 1ₘ) with Z ⪰ 0. Still, there is no duplication of the Z-coefficients (apart from
    # symmetry), so there is no way in which these unknown coefficients could potentially cancel: every entry in Z will
    # appear in exactly one cell in a triangle of M. So basically, all that we have to do is to apply the SOS
    # decomposition cell-wise.
    isempty(psd) || @verbose_info("├ PSD constraints")
    for (groupings, psdᵢ) in zip(groupings.psds, psd)
        dim = size(psdᵢ, 1)
        newbound = sum(∘(trisize, length), groupings) * sum(@capture(length($psdᵢ[i, j]) for j in 1:dim for i in 1:j), init=0)
        sizehint!(mons_idx_set, length(mons_idx_set) + (iszero(Nc) ? newbound : 2newbound))
        @inbounds for j in 1:dim, i in 1:j
            for t in psdᵢ[i, j]
                monₜ = monomial(t)
                for grouping in groupings
                    grouping_loop(grouping, monₜ, Val(i == j))
                end
            end
        end
    end

    # now we need to re-cast the indices into the exponent-representations
    @verbose_info("├ Total number of coefficients: ", length(mons_idx_set))
    return SimpleMonomialVector{Nr,Nc}(e, sort!(convert.(UInt, mons_idx_set)))
end