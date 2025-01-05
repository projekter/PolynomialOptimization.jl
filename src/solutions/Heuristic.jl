"""
    poly_solutions(:heuristic, result::Result; verbose=false)

Apply a heuristic solution extraction mechanism to a solved polynomial optimization problem. If there is only a single
solution, the heuristic will always work (but is much faster than a rank factorization). It may in certain cases also work if
there are multiple solutions if the individual components only differ by their signs/phases. This method can also be applied in
the presence of term sparsity, where [`poly_solutions`](@ref) fails.
This function gives a vector of solutions (not an iterator, as [`poly_solutions`](@ref) does).

See also [`poly_optimize`](@ref).
"""
function poly_solutions(::Val{:heuristic}, result::Result; verbose::Bool=false)
    # This whole procedure works under the assumption that there be a single unique solution, or that are multiple solutions
    # which may differ in the signs of arguments which only show up with even exponents. This may not necessarily be the case;
    # then, when analyzing the problem, set the perturbation parameter to some small value, which should lift all degeneracies.
    # However, this will also slightly modify the problem, changing the optimal value (and its point) a bit. For small
    # perturbations, the solution is close enough to provide ideal starting points for another optimizer.
    solutions = poly_solutions(Val(Symbol("heuristic-postprocess")), result.moments,
        poly_solutions(Val(Symbol("heuristic-magnitudes")), result.relaxation, result.moments)...)
    @verbose_info("Extracted all signs/phases, found ", length(solutions), " possible solution",
        isone(length(solutions)) ? "" : "s")
    return solutions
end

function poly_solutions(::Val{Symbol("heuristic-magnitudes")},
    relaxation::AbstractRelaxation{<:Problem{<:SimplePolynomial{<:Any,Nr,Nc,<:SimpleMonomialVector{Nr,Nc,MI}}}},
    moments::MomentVector{R,V,Nr,Nc}) where {Nr,Nc,MI<:Integer,R<:Real,V<:Union{R,Complex{R}}}
    # In this part, we extract the parts of the solution by looking at the monomials that are powers of variables.
    solution = fill(iszero(Nc) ? R(NaN) : Complex(R(NaN), R(NaN)), Nr + Nc)
    I = SimplePolynomials.smallest_unsigned(Nr + 2Nc)
    zero_checks = Set{I}()
    unknown_signs = Set{I}()
    unknown_phases = Dict{I,Set{R}}()
    e = ExponentsAll{Nr+2Nc,MI}()
    deg = degree(relaxation)
    for i in 1:Nr
        mon = SimpleMonomial(e, SimpleRealVariable{Nr,Nc}(i))
        # The monomial is present in the original basis, so we try to reconstruct it by searching for its powers. We will
        # favor even powers (even over the variable itself), as they are not affected by multiple sign-symmetric solutions
        for pow in 2:2:2deg
            val = moments[mon^pow]
            if !isnan(val)
                if abs(val) ≤ R(1e-7) # almost zero
                    solution[i] = zero(V)
                    @goto real_done
                elseif val ≥ zero(R)
                    solution[i] = val^inv(pow)
                    push!(unknown_signs, i)
                    @goto real_done
                end
                # else something bad happened - probably the solution was not unique. We will proceed and try to get it by
                # using remaining powers, but the result will probably be useless.
            end
        end
        for pow in 1:2:2deg
            val = moments[mon^pow]
            if !isnan(val)
                solution[i] = abs(val)^inv(pow)
                if val < 0
                    solution[i] = -solution[i]
                end
                @goto real_done
            end
        end
        push!(zero_checks, i)
        @label real_done
    end
    for i in Nr+1:Nr+Nc
        mon = SimpleMonomial(e, SimpleComplexVariable{Nr,Nc}(i))
        # For complex-valued variables, favor the variable itself, as we cannot say which power would gobble the phases.
        val = moments[mon]
        if !isnan(val)
            solution[i] = moments[mon]
        else
            # The monomial itself does not occur in our sparse optimization list, but it is present in the original basis, so
            # we try to reconstruct it by searching for its powers. The phase reconstruction is ambiguous.
            abs_val = R(NaN)
            phase_candidates = Set{R}()
            for pow1 in 1:deg, pow2 in 0:deg
                pow1 == pow2 && continue # we cannot get any phase information
                val = moments[mon^pow1 * SimpleConjMonomial(mon)^pow2]
                if !isnan(val)
                    Σpow, δpow = pow1 + pow2, pow1 - pow2
                    # we have candidates for the phase
                    if isnan(abs_val)
                        abs_val = abs(val)^inv(Σpow)
                        this_phase = atan(imag(val), real(val)) / δpow
                        ω = R(2) * π / δpow
                        append!(empty!(phase_candidates), (this_phase + j * ω for j in 0:δpow-1))
                    else
                        # consistency check: all powers should give rise to the same magnitude
                        if abs(abs_val - abs(val)^inv(Σpow)) > R(1e-7)
                            abs_val = R(NaN)
                            break
                        end
                        # Now check which phase candidates still survive. This is tricky given the floating-point inaccuracy.
                        new_candidates = Set{R}()
                        pow_phase = atan(imag(val), real(val))
                        for old_candidate in phase_candidates
                            if abs(mod2pi(old_candidate * δpow) - pow_phase) < R(1e-7)
                                new_candiates |= old_candiate
                            end
                        end
                        phase_candidates = new_candidates
                    end
                    # even if we end up with just one phase candidate, we continue the consistency check with higher powers
                end
            end
            if isnan(abs_val) || isempty(phase_candidates)
                push!(zero_checks, i)
            elseif isone(length(phase_candidates))
                solution[i] = abs_val * cis(@inbounds phase_candidates[begin])
            else
                solution[i] = abs_val
                unknown_phases[i] = phase_candidates
            end
        end
    end
    return solution, zero_checks, unknown_signs, unknown_phases
end

function poly_solutions(::Val{Symbol("heuristic-postprocess")}, moments::MomentVector{R,V,Nr,Nc}, candidate::AbstractVector{V},
    zero_checks::Set{I}, unknown_signs::Set{I}, unknown_phases::Dict{I,Set{R}}) where {Nr,Nc,I<:Integer,R<:Real,V<:Union{R,Complex{R}}}
    isempty(zero_checks) && isempty(unknown_signs) && isempty(unknown_phases) && return [candidate]
    # Candidate contains the possible absolute values of the solutions. However,
    # - some unknown entries may in fact be zero: check for those indices in zero_checks
    # - some real-valued variables have an unknown sign: check for those in unknown_signs
    # - some complex-valued variables have multiple phases that could be assigned: check for those in unknown_phases
    dependencies = Set{I}() # We might not be able to resolve everything in the first iteration; keep track of things that
                            # would allow us to do better in the next.
    while true
        retry = false
        for (mon, mon_val) in MomentAssociation(moments)
            isnan(mon_val) && continue
            if abs(mon_val) < R(1e-7)
                isempty(zero_checks) && continue
                zerovar = zero(I)
                # Check whether this monomial contains exactly one of the variables that we want to check, and no other ones
                # that we already know to be zero.
                for (mon_var, var_exp) in mon
                    @assert(!iszero(var_exp))
                    var_idx = variable_index(mon_var)
                    if var_idx ∈ zero_checks
                        if iszero(zerovar)
                            zerovar = var_idx
                        else
                            zerovar = zero(I)
                            break
                        end
                    else
                        val = solution[var_idx]
                        if abs(val) < R(1e-7)
                            # does not help, there's another zero variable
                            zerovar = zero(I)
                            break
                        end
                    end
                end
                if !iszero(zerovar)
                    solution[zerovar] = zero(V)
                    setdiff!(zero_checks, (zerovar,))
                    retry |= zerovar ∈ dependencies
                end
            else
                (iszero(Nr) || isempty(unknown_signs)) && (iszero(Nc) || isempty(unknown_phases)) && continue
                # For this monomial, check: does it contain variables with unknown signs/phases? If so, can we fix the
                # signs/phases of all variables but one to a concrete value?
                unknownvar = zero(I)
                if iszero(Nc)
                    valpos = mon_val > 0
                else
                    valangle = angle(mon_val)
                    unknownvarexp = 0
                end
                for (mon_var, var_exp) in mon
                    @assert(!iszero(var_exp))
                    var_idx = variable_index(mon_var)
                    if isnan(candidate[var_idx])
                        # everything about this variable is unknown, making the monomial useless
                        unknownval = zero(I)
                        push!(dependencies, var_idx)
                        break
                    end
                    if !iszero(Nr) && var_idx ∈ unknown_signs
                        # The sign is unknown. Does this have any effect on the monomial?
                        if isodd(var_exp)
                            # The sign is unknown. Do we already have something else which we don't know anything about?
                            if !iszero(unknownvar)
                                # We do, so we cannot make local statements (and to keep track of the correlations this tells
                                # us about is just too laborsome).
                                push!(dependencies, unknownvar, var_idx)
                                unknownvar = zero(I)
                                break
                            end
                            unknownvar = var_idx
                        end # else, no effect
                    elseif isreal(mon_var)
                        # The sign is known; the variable is real-valued. We can calculate the precise effect removing this
                        # variable from the monomial has: namely, flip the sign if the exponent is odd and the value negative.
                        if isodd(var_exp)
                            var_val = R(candidate[var_idx]) # candidate might be in total a complex-valued vector
                            if var_val < zero(R)
                                iszero(Nc) ? (valpos = !valpos) : (valangle = mod2pi(valangle - pi))
                            end
                        end
                    elseif !iszero(Nc) && var_idx ∈ unknown_phases
                        # The phase is unknown. We have several candidate phases that we could assign, but we don't know which.
                        # However, maybe for this monomial it doesn't play a role because all are give rise to the same value
                        # when exponentiated...
                        exp_phases = unknown_phases[var_idx]
                        exp_phase, rest = Iterators.peel(exp_phases)
                        if any(@capture(ϕ -> mod2pi((ϕ - $exp_phase) * $var_exp) ≥ 1e-7), rest)
                            if !iszero(unknownvar)
                                push!(dependencies, unknownvar, var_idx)
                                unknownvar = zero(I)
                                break
                            end
                            unknownvar = var_idx
                            unknownvarexp = var_exp
                        else
                            # The phase may be unknown, but for this monomial, we don't care.
                            valangle = mod2pi(valangle - exp_phase * var_exp)
                        end
                    else
                        # The phase is known; the variable is complex-valued. We can calculate the precise effect removing this
                        # variable from the monomial has: namely, subtract the phase of the variable to its power.
                        valangle = mod2pi(valangle - angle(candidate[var_idx]) * var_exp)
                    end
                end
                if !iszero(unknownvar)
                    # We have exactly one variable for which we don't know its sign.
                    if iszero(Nc)
                        valpos || (candidate[unknownvar] = -candidate[unknownvar])
                        isempty(setdiff!(unknown_signs, (unknownvar,))) && break
                    elseif unknownvar ≤ Nr
                        # Although we always carry out the modulo for every step to keep the accuracy high, we cannot guarantee
                        # that the phase is exactly ±π. But it should be, so let's enforce this in debug mode.
                        @assert(valangle < R(1e-7) || π - R(1e-7) ≤ valangle ≤ π + R(1e-7))
                        valangle > inv(R(2)) * π && (candidate[unknownvar] = -candidate[unknownvar]) # and use a simple threshold
                        isempty(setdiff!(unknown_signs, (unknownvar,))) && isempty(unknown_phases) && break
                    else
                        candidate[unknownvar] *= cis(valangle / unknownvarexp)
                        isempty(setdiff!(unknown_phases, (unknownvar,))) && isempty(unknown_signs) && break
                    end
                    retry |= unknownvar ∈ dependencies
                end
            end
        end
        retry || break # something changed that was a blocker before - let's do it again!
    end
    # Let the user do something on their own. It might for example be the case that powers of var could be calculated from the
    # the ideal or enough other information is available to somehow get its value.
    isempty(zero_checks) || @warn("Could not discern value for " * string(zero_checks))
    if !isempty(unknown_signs) || !isempty(unknown_phases)
        empty!(zero_checks)
        # There are still issues because we cannot assign some signs/phases definitely. We will try all possible combinations.
        # However, any single choice may lead to eliminating lots of other possibilities, so we don't check for all tuple.
        # First, check the signs, for which there are only two options.
        if !isempty(unknown_signs)
            signpos = pop!(unknown_signs)
            first_sol = poly_solutions(Val(Symbol("heuristic-postprocess")), moments, copy(candidate), zero_checks,
                copy(unknown_signs), copy(unknown_phases))
            candidate[signpos] = -candidate[signpos]
            second_sol = poly_solutions(Val(Symbol("heuristic-postprocess")), moments, candidate, zero_checks, unknown_signs,
                unknown_phases)
            return [first_sol; second_sol]
        else
            # Find the variable with the last possible combinations
            candidate_index, candidate_phases = @inbounds unknown_phases[begin]
            for (k, v) in unknown_phases
                if length(v) < candidate_phases
                    candidate_index, candidate_phases = k, v
                end
            end
            delete!(unknown_phases, candidate_index)
            solutions = typeof(candidate)[]
            for candidate_phase in candidate_phases
                candidateᵢ = copy(candidate)
                candidateᵢ[candidate_index] *= cis(candidate_phase)
                append!(solutions, poly_solutions(Val(Symbol("heuristic-postprocess")), moments, candidateᵢ, zero_checks,
                    copy(unknown_signs), copy(unknown_phases)))
            end
            return solutions
        end
    end
    return [candidate]
end