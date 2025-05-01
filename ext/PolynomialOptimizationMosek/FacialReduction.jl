# this is the multiple-SOS version (III.D); for the single-SOS algorithm, simply use [M] and [gs;;] as arguments
function FacialReduction.facial_reduction!!(::Val{:Mosek}, M::AbstractVector{<:IntMonomialVector{Nr,Nc,I}},
    gs::AbstractMatrix{<:IntPolynomial{<:Any,Nr,Nc,<:IntMonomialVector{Nr,Nc,I}}},
    M⁺::AbstractVector{<:IntMonomialVector{Nr,Nc,I}}; verbose::Bool) where {Nr,Nc,I<:Integer}
    @assert length(M) == size(gs, 2)
    all(x -> length(x) ≤ 1, M) && return
    task = Mosek.Task(Mosek.msk_global_env::Mosek.Env)
    newExponents = Vector{Int}(undef, Nr + 2Nc)

    numM⁺ = sum(length, M⁺, init=0)
    Mosek.appendvars(task, numM⁺) # these are the λ variables
    Mosek.putvarboundsliceconst(task, 1, numM⁺ +1, Mosek.MSK_BK_LO, 0., Inf)
    numvars_avail = numvars = numM⁺
    unusedλ = trues(numM⁺)

    numgs = length(gs)
    Mosek.appendcons(task, numgs +1)
    Mosek.putconboundsliceconst(task, 1, numgs +1, Mosek.MSK_BK_FX, 0., 0.)
    Mosek.@MSK_putarow(task.task, numgs, numM⁺, collect(zero(Int32):Int32(numM⁺ -1)), fill(one(Float64), numM⁺))
    Mosek.putconbound(task, numgs +1, Mosek.MSK_BK_FX, 1.0, 1.0)

    λidxshift = 0
    @inbounds for (gⱼ, Mⱼ, M⁺ⱼ) in zip(eachcol(gs), M, M⁺)
        extdegMⱼ = extdegree(Mⱼ)
        pcoeffs = Dict{I,Int}()
        for (k, gⱼₖ) in enumerate(gⱼ), t in terms(gⱼₖ) # k ∈ {0, ..., m}
            # sᵢ = (p + ∑_{a ∈ M⁺} λ₂ₐ e₂ₐ) × ...
            # ⟨sᵢⱼ, gⱼₖ(x)⟩ = 0, so we only need to check nonzero terms in g
            # In principle, this should be the weaker conditions ∑ⱼ ⟨sᵢⱼ, gⱼₖ₍ⱼ₎⟩ = 0 ∀ k(j); but the number of k(j) is
            # exponential in size(gs, 2). So demanding equality on the individual parts is more restrictive (potentially
            # leading to fewer solutions). However, consider that the g come from a constrained Putinar certificate, which
            # means that the non-zero parts in every column in gs are all disjoint. Therefore, the sums in fact exactly
            # collapse to checking just the individual parts.
            m = monomial(t)
            if FacialReduction.inΣhatofMperp(m, Mⱼ, extdegMⱼ, newExponents)
                if haskey(pcoeffs, m.index)
                    pidx = pcoeffs[m.index]
                else
                    pidx = pcoeffs[m.index] = (numvars += 1)
                    if numvars_avail < numvars
                        Mosek.appendvars(task, 32)
                        numvars_avail += 32
                    end
                end
                Mosek.putaij(task, k, pidx, coefficient(t))
            end
            λidx = FacialReduction.λindex(m, M⁺ⱼ, newExponents)
            if !isnothing(λidx)
                Mosek.putaij(task, k, λidx + λidxshift, coefficient(t))
                unusedλ[λidx+λidxshift] = false
            end
        end
        λidxshift += length(M⁺ⱼ)
    end
    Mosek.putvarboundsliceconst(task, numM⁺ +1, numvars +1, Mosek.MSK_BK_FR, -Inf, Inf)
    # We want to maximize the number of nonzeros (in the used λs). This is not possible in a linear program; however, we should
    # use an interior point algorithm instead of the simplex to get a most non-sparse solution.
    Mosek.putintparam(task, Mosek.MSK_IPAR_OPTIMIZER, Mosek.MSK_OPTIMIZER_INTPNT)
    # Maybe this helps to make the result more dense? It usually does so when uniformly distributed λ are possible - but what
    # if due to constraints, maximizing the sum would actually mean maximizing a single one and putting the others to zero?
    #Mosek.putcslice(task, 1, numM⁺ +1, Float64.(.!unusedλ))
    #Mosek.putobjsense(task, Mosek.MSK_OBJECTIVE_SENSE_MAXIMIZE)
    Mosek.optimize(task)
    if Mosek.getsolsta(task, Mosek.MSK_SOL_ITR) == Mosek.MSK_SOL_STA_OPTIMAL
        λ = Mosek.getxxslice(task, Mosek.MSK_SOL_ITR, 1, numM⁺ +1)
        return FacialReduction.λtoresult(λ, M, M⁺, unusedλ)
    end
    return # no more reduction possible
end