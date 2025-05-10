mutable struct FRState
    const task::Mosek.Task
    const newExponents::Vector{Int}
    const numM⁺::Int
    numvars_avail::Int
    numvars::Int
    const usedλ::BitVector
    const numgs::Int
    λidxshift::Int

    function FRState(data::FacialReduction.FacialReductionData{<:Real,NrPlusNy}) where {NrPlusNy}
        task = Mosek.Task(Mosek.msk_global_env::Mosek.Env)
        newExponents = Vector{Int}(undef, NrPlusNy)

        numM⁺ = data.numM⁺
        Mosek.appendvars(task, numM⁺) # these are the λ variables
        Mosek.putvarboundsliceconst(task, 1, numM⁺ +1, Mosek.MSK_BK_LO, 0., Inf)
        numvars_avail = numvars = numM⁺
        usedλ = falses(numM⁺)

        numgs = data.numMons
        Mosek.appendcons(task, numgs +1)
        Mosek.putconboundsliceconst(task, 1, numgs +1, Mosek.MSK_BK_FX, 0., 0.)
        Mosek.@MSK_putarow(task.task, numgs, numM⁺, collect(zero(Int32):Int32(numM⁺ -1)), fill(one(Float64), numM⁺))
        Mosek.putconbound(task, numgs +1, Mosek.MSK_BK_FX, 1.0, 1.0)

        new(
            task,
            Vector{Int}(undef, NrPlusNy),
            numM⁺, numM⁺, numM⁺, falses(numM⁺),
            numgs, 0
        )
    end
end

function FacialReduction.facial_reduction!(state::FRState,
    dataslice::FacialReduction.AbstractFacialReductionDataSlice{<:(IntPolynomial{<:Any,Nr,0,<:IntMonomialVector{Nr,0,I}} where {Nr})}) where {I<:Integer}
    pcoeffs = Dict{I,Int}()
    @inbounds for (k, gⱼₖ) in dataslice, t in terms(gⱼₖ) # k ∈ {0, ..., m}
        # sᵢ = (p + ∑_{a ∈ M⁺} λ₂ₐ e₂ₐ) × ...
        # ⟨sᵢⱼ, gⱼₖ(x)⟩ = 0, so we only need to check nonzero terms in g
        m = monomial(t)
        if FacialReduction.inΣhatofMperp(dataslice.M, m, dataslice.extdegM, state.newExponents)
            if haskey(pcoeffs, m.index)
                pidx = pcoeffs[m.index]
            else
                pidx = pcoeffs[m.index] = (state.numvars += 1)
                if state.numvars_avail < state.numvars
                    Mosek.appendvars(state.task, 32)
                    state.numvars_avail += 32
                end
            end
            Mosek.putaij(state.task, k, pidx, coefficient(t))
        end
        λidx = FacialReduction.λindex(dataslice.M⁺, m, state.newExponents)
        if !iszero(λidx)
            Mosek.putaij(state.task, k, λidx + state.λidxshift, coefficient(t))
            state.usedλ[λidx+state.λidxshift] = true
        end
    end
    state.λidxshift += length(dataslice.M⁺)
    return
end

function FacialReduction.facial_reduction!(::Val{:Mosek}, data::FacialReduction.FacialReductionData; verbose::Bool)
    state = FRState(data)
    for dataslice in data
        FacialReduction.facial_reduction!(state, dataslice) # function barrier
    end
    Mosek.putvarboundsliceconst(state.task, state.numM⁺ +1, state.numvars +1, Mosek.MSK_BK_FR, -Inf, Inf)
    # We want to maximize the number of nonzeros (in the used λs). This is not possible in a linear program; however, we should
    # use an interior point algorithm instead of the simplex to get a most non-sparse solution.
    Mosek.putintparam(state.task, Mosek.MSK_IPAR_OPTIMIZER, Mosek.MSK_OPTIMIZER_INTPNT)
    # This is a good heuristic to make the results more dense:
    Mosek.putcslice(state.task, 1, state.numM⁺ +1, Float64.(state.usedλ))
    Mosek.putobjsense(state.task, Mosek.MSK_OBJECTIVE_SENSE_MAXIMIZE)
    Mosek.optimize(state.task)
    if Mosek.getsolsta(state.task, Mosek.MSK_SOL_ITR) == Mosek.MSK_SOL_STA_OPTIMAL
        λ = Mosek.getxxslice(state.task, Mosek.MSK_SOL_ITR, 1, state.numM⁺ +1)
        return FacialReduction.truncate!(data, λ, state.usedλ)
    end
    return false # no more reduction possible
end