export optimality_certificate

"""
    optimality_certificate(result::Result, ϵ=1e-6)

This function applies the flat extension/truncation criterion to determine whether the optimality of the given problem can be
certified, in which case it returns `:Optimal`. If no such certificate is found, the function returns `:Unknown`. The criterion
is meaningless for sparse problems or if a full basis is not available.
The parameter `ϵ` controls the bound below which singular values are considered to be zero, and its negative below which
eigenvalues are considered to be negative.

See also [`poly_optimize`](@ref).
"""
function optimality_certificate(result::Result{<:AbstractRelaxation{<:Problem{<:SimplePolynomial{<:Any,Nr,Nc}}}}, ϵ::R=1e-6) where {Nr,Nc,R<:Real}
    relaxation = result.relaxation
    deg = degree(relaxation)
    deg < 1 && return :CertificateUnavailable
    d₀ = maxhalfdegree(relaxation.objective)
    dconstr = max(maximum(maxhalfdegree, relaxation.constr_zero, init=0),
        maximum(maxhalfdegree, relaxation.constr_nonneg, init=0), maximum(maxhalfdegree, relaxation.constr_psd, init=0))
    if isone(d₀) && isone(dconstr)
        # this is the one case that is excluded in the loop, but the rank-1 condition is also valid here
        isone(rank(moment_matrix(result, max_deg=1), atol=ϵ)) && return :Optimal
    end
    dₖ = max(iszero(Nc) ? 1 : 2, dconstr)
    n = Nr + Nc
    vars = variables(relaxation.objective)
    for t in deg:-1:max(d₀, dₖ)
        Mₜ = moment_matrix(result, max_deg=t)
        rkMₜ = rank(Mₜ, atol=ϵ)
        isone(rkMₜ) && return :Optimal
        max_deg = t - dₖ
        m₁₁ = moment_matrix(result; max_deg)
        if rkMₜ == rank(m₁₁, atol=ϵ)
            fail = false
            if n > 1
                @inbounds for i in 1:n-1, j in max(i, Nr +1):n # we don't have to check real*real
                    # Note: there's no theory on mixing real and complex-valued variables yet. Let's take the pragmatic
                    # approach here: As soon as a complex-valued variable is involved, we have the more restrictive variant
                    # that required dₖ ≥ 2, and we check the additional criteria for every monomial that contains a complex
                    # variable, as it will not be the same as its conjugate.
                    zᵢ, zⱼ = vars[i], vars[j]
                    m₁₂ = moment_matrix(result; max_deg, prefix=zᵢ)
                    m₁₃ = moment_matrix(result; max_deg, prefix=zⱼ)
                    m₂₂ = moment_matrix(result; max_deg, prefix=zᵢ*SimpleConjMonomial(zᵢ))
                    m₂₃ = moment_matrix(result; max_deg, prefix=zⱼ*SimpleConjMonomial(zᵢ))
                    m₃₃ = moment_matrix(result; max_deg, prefix=zⱼ*SimpleConjMonomial(zⱼ))
                    if min(eigvals!(Hermitian([m₁₁       m₁₂       m₁₃
                                               conj(m₁₂) m₂₂       m₂₃
                                               conj(m₁₃) conj(m₂₃) m₃₃]))) < -ϵ
                        fail = true
                        break
                    end
                end
            elseif !iszero(Nc)
                @inbounds z = vars[1]
                m₁₂ = moment_matrix(result; max_deg, prefix=z)
                m₂₂ = moment_matrix(result; max_deg, prefix=z*SimpleConjMonomial(z))
                if min(eigvals!(Hermitian([m₁₁ m₁₂; conj(m₁₂) m₂₂]))) < -ϵ
                    fail = true
                end
            end
            fail || return :Optimal
        end
    end
    return :Unknown
end