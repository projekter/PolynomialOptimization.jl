function grouping_loop(::Val, ::Val{0}, ::Val, mons_idx_set, e, range₁, grouping, constr, sync)
    this_set = isnothing(sync) ? mons_idx_set : Set{FastKey{UInt}}()
    # TODO: veciter? range₁ are at most two units ranges
    @inbounds for t in constr
        monₜ = monomial(t)
        for i in range₁
            g₁ = grouping[i]
            for g₂ in Iterators.take(grouping, i)
                push!(this_set, FastKey(monomial_index(e, g₁, monₜ, g₂)))
            end
        end
    end
    if !isnothing(sync)
        lock(sync)
        try
            union!(mons_idx_set, this_set)
        finally
            unlock(sync)
        end
    end
    return
end

function grouping_loop(::Val{0}, ::Val, ::Val{true}, mons_idx_set, e, range, grouping, constr, sync)
    this_set = isnothing(sync) ? mons_idx_set : Set{FastKey{UInt}}()
    gview = @inbounds @view grouping[range]
    for t in constr
        monₜe = KillConjugates{0}(exponents(monomial(t)))
        for g in gview
            # While the prefactor is g₁*conj(g₂), we drop the conjugates anyway.
            # Note that no prefactor has a single conjugate component; and while the monomials might have them, the
            # conjugates will come up later anyway!
            push!(this_set, FastKey(exponents_sum(e, exponents(g), monₜe)[1])) # KillConjugates(g*monₜ*ḡ₂)
        end
    end
    if !isnothing(sync)
        lock(sync)
        try
            union!(gl.mons_idx_set, this_set)
        finally
            unlock(sync)
        end
    end
    return
end

function grouping_loop(::Val{0}, ::Val, ::Val{false}, mons_idx_set, e, range, grouping, constr, sync)
    this_set = isnothing(sync) ? mons_idx_set : Set{FastKey{UInt}}()
    gview = @inbounds @view grouping[range]
    for t in constr
        monₜ = monomial(t)
        monₜe = KillConjugates{0}(exponents(monₜ))
        monₜec = KillConjugates{0}(exponents(IntConjMonomial(monₜ)))
        for g in gview
            # Here, it is not guaranteed that we will see the conjugates later on.
            ge = exponents(g)
            push!(this_set,
                FastKey(exponents_sum(e, ge, monₜe)[1]), # KillConjugates(g₁*monₜ*ḡ₂)
                FastKey(exponents_sum(e, ge, monₜec)[1]) # KillConjugates(g₁*conj(monₜ)*ḡ₂)
            )
        end
    end
    if !isnothing(sync)
        lock(sync)
        try
            union!(mons_idx_set, this_set)
        finally
            unlock(sync)
        end
    end
    return
end

# Out of necessity, our groupings contain abstract types. Given that the number of groupings will be manageable, but the
# groupings themselves can be quite large, it pays off to introduce a function barrier specializing on the actual type.
function grouping_loop(Nr::Val, ::Val{Nc}, constr_is_real::Val, mons_idx_set, e, grouping, constr) where {Nc}
    nthreads = Threads.nthreads()
    if isone(nthreads) || length(grouping) < 100
        grouping_loop(Nr, Val(Nc), constr_is_real, mons_idx_set, e, 1:length(grouping), grouping, constr, nothing)
    elseif iszero(Nc)
        # The first entry in grouping has length 1, the second length 2, ...; consider this when dividing the search space.
        # A single task will always do n items from the beginning plus n items from the end, giving a constant number of
        # iterations.
        iterations = div(length(grouping), 2, RoundUp)
        iterationsperthread = div(iterations, nthreads, RoundUp)
        startindex = 1
        ccall(:jl_enter_threaded_region, Cvoid, ())
        try
            threads = Vector{Task}(undef, nthreads)
            sync = Threads.SpinLock()
            @inbounds for tid in 1:nthreads
                ranges₁ = Iterators.flatten(
                    (startindex:min(startindex + iterationsperthread -1, iterations),
                     range(max(length(grouping) - startindex +1 - iterationsperthread, iterations) +1,
                           length(grouping) - startindex +1))
                )
                startindex += iterationsperthread
                threads[tid] = Threads.@spawn grouping_loop($Nr, $(Val(Nc)), $constr_is_real, $mons_idx_set, $e, $ranges₁,
                                                            $grouping, $constr, $sync)
            end
            for thread in threads
                wait(thread)
            end
        finally
            ccall(:jl_exit_threaded_region, Cvoid, ())
        end
    else
        iterations = length(grouping)
        iterationsperthread = div(iterations, nthreads, RoundUp)
        startindex = 1
        try
            threads = Vector{Task}(undef, iszero(Nc) ? nthreads : 2nthreads)
            sync = Threads.SpinLock()
            @inbounds for tid in 1:nthreads
                ranges = startindex:min(startindex + iterationsperthread -1, iterations)
                startindex += iterationsperthread
                threads[tid] = Threads.@spawn grouping_loop($Nr, $(Val(Nc)), $constr_is_real, $mons_idx_set, $e, $ranges,
                                                            $grouping, $constr, $sync)
            end
            for thread in threads
                wait(thread)
            end
        finally
            ccall(:jl_exit_threaded_region, Cvoid, ())
        end
    end
    return
end

# We only accept ExponentsAll; else, we cannot guarantee that all the resulting exponents will in fact be possible. In this
# way, we don't have to check for equality of the exponents sets (or form the largest cover), which is also beneficial.
# And converting polynomials from the MP interface will yield ExponentsAll anyway.
function merge_constraints(objective::IntPolynomial{<:Any,Nr,Nc}, zero::AbstractVector{<:IntPolynomial{<:Any,Nr,Nc}},
    nonneg::AbstractVector{<:IntPolynomial{<:Any,Nr,Nc}},
    psd::AbstractVector{<:AbstractMatrix{<:IntPolynomial{<:Any,Nr,Nc}}}, prefactor::IntPolynomial{<:Any,Nr,Nc},
    groupings::RelaxationGroupings{Nr,Nc}, verbose::Bool, need_copy::Bool) where {Nr,Nc}
    @verbose_info("Incorporating constraints into set of exponents")
    # Note that in the complex-valued case there's no mixing - i.e., no real variables. And every monomial appears once in its
    # original form, once in its conjugate. We are only interested in the "original" (whichever it is), so we discard the
    # conjugate part to avoid even generating duplicates.
    e = ExponentsAll{Nr+2Nc,UInt}()

    if isempty(zero) && isempty(nonneg) && isempty(psd)
        # short path: just directly add the stuff from the objective and prefactor. While we can get duplicates (due to the
        # elimination of conjugates), we have a clear and manageable upper bound to the size, so no need to determine
        # duplicates on-the-fly.
        @verbose_info("No constraints found, adding all objective and prefactor monomials")
        mons = monomials(objective)
        pmons = monomials(prefactor)
        if iszero(Nc)
            if mons.e == e && pmons.e == e && issubset(pmons, mons)
                @verbose_info("Aliasing monomial indices")
                return need_copy ? copy(mons) : mons
            else
                @verbose_info("Merging monomial indices")
                return merge_monomial_vectors(Val(Nr), Val(Nc), e, [pmons, mons])
            end
        else
            @verbose_info("Converting monomial indices")
            candidates = FastVec{UInt}(buffer=length(mons) + length(pmons))
            for mon in mons
                @inbounds unsafe_push!(candidates, exponents_to_index(e, KillConjugates{Nr}(exponents(mon))))
            end
            for mon in pmons
                @inbounds unsafe_push!(candidates, exponents_to_index(e, KillConjugates{Nr}(exponents(mon))))
            end
            @verbose_info("Sorting and removing duplicates")
            return IntMonomialVector{Nr,Nc}(e, Base._groupedunique!(sort!(finish!(candidates))))
        end
    end

    mons_idx_set = sizehint!(Set{FastKey{UInt}}(), length(objective) + length(prefactor))
    @verbose_info("├ objective")
    for mon in monomials(objective)
        push!(mons_idx_set, FastKey(iszero(Nc) ? monomial_index(e, mon) :
                                    exponents_to_index(e, KillConjugates{Nr}(exponents(mon)))))
    end
    # The prefactor is already multiplied into the objective. We only need it explicitly because of the lower bound which
    # ensures that we always have the constant (or, more precisely, the prefactor) as a part of the polynomial.
    @verbose_info("├ prefactor")
    for mon in monomials(prefactor)
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
            for grouping in groupings
                grouping_loop(Val(Nr), Val(Nc), Val(!isz), mons_idx_set, e, grouping, constrᵢ)
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
        @inbounds for j in 1:dim
            for i in 1:j-1
                for grouping in groupings
                    grouping_loop(Val(Nr), Val(Nc), Val(false), mons_idx_set, e, grouping, psdᵢ[i, j])
                end
            end
            for grouping in groupings
                grouping_loop(Val(Nr), Val(Nc), Val(true), mons_idx_set, e, grouping, psdᵢ[j, j])
            end
        end
    end

    # now we need to re-cast the indices into the exponent-representations
    @verbose_info("└ Total number of coefficients: ", length(mons_idx_set))
    return IntMonomialVector{Nr,Nc}(e, sort!(convert.(UInt, mons_idx_set)))
end