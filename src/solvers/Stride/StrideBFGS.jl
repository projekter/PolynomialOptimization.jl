bfgsmemα = Vector{Float64}(undef, mem_lbfgs)
bfgsmemρ = Vector{Float64}(undef, mem_lbfgs)
bfgsmemS = Matrix{Float64}(undef, d, mem_lbfgs)
bfgsmemY = Matrix{Float64}(undef, d, mem_lbfgs)
mlbfgs!(ssd, bfgsmemα, bfgsmemρ, bfgsmemS, bfgsmemY, Inf, tol_lbfgs, ls_bfgs, 1/2, one(R), one(R), maxiter_lbfgs, verbose_lbfgs)


function mlbfgs!(ssd::StrideData{R}, α::Vector{R}, ρs::Vector{R}, S::Matrix{R}, Y::Matrix{R}, K::R, tol::R, linesearch, ρ::R,
    τ₁::R, τ₂::R, maxiter::Integer, verbose::Bool) where {R<:Real}
    @assert(K > 0 && tol > 0 && 0 < ρ < 1 && τ₁ > 0 && τ₂ > 0)
    @verbose_info("Entering modified limited-memory BFGS method\nIteration |     Residue")
    mem = length(α)
    A, Astar, b, Z, W, ξ = ssd.A, ssd.Astar, ssd.b, ssd.X, ssd.S, ssd.y
    @myinbounds tmpm, ∇ϕ, tmpv = ssd.tmpm1, ssd.tmpv, view(ssd.tmpvlong, 1:ssd.d)
    # Iterate the following steps for k = 1, ...
    k = 1
    for (tmpmᵢ, Zᵢ, Astarᵢ) in zip(tmpm, Z, Astar)
        # X(W, ξ) = Π(A* ξ + Z)
        copyto!(tmpmᵢ, Zᵢ)
        mul!(tmpmᵢ, Astarᵢ, ξ, true, true)
        psdproject!(tmpmᵢ)
    end
    local idx
    @myinbounds while true
        # ϕ(ξ) := ‖Π(A* ξ + Z)‖^2/2 - ⟨b, ξ⟩
        # ∇ϕ(ξ) = AΠ(A* ξ + Z) - b
        # tmpm will already contain the up-to-date Π(A* ξ + Z)
        ϕ = sum(normsq, tmpm, init=zero(R)) / R(2) - dot(b, ξ)
        ∇ϕ .= .-b
        for (Aᵢ, tmpmᵢ) in zip(A, tmpm)
            mul!(∇ϕ, Aᵢ, tmpmᵢ, true, true)
        end
        @views if k != 1
            # leftover from the last iteration (Y[:, idx] already contains last ∇ϕ)
            # Step 3 (Update memory):
            #     Compute and save uₖ = ξₖ₊₁ - ξₖ and wₖ = ∇ϕ(ξₖ₊₁) - ∇ϕ(ξₖ).
            Y[:, idx] .= ∇ϕ .- Y[:, idx]
            ρs[idx] = one(R) / dot(S[:, idx], Y[:, idx]) # required for the two-loop recursion
        end
        # Step 1 (Search direction):
        #     Choose Qₖ₀ ≻ 0 (we always choose Iₘ - actually, due to Nocedal, Wright, we choose a proportionality)
        #            βₖ := τ₁ ‖∇ϕ(ξₖ)‖^(τ₂)
        β = τ₁*norm(∇ϕ)^τ₂
        #     and compute dₖ = -βₖ ∇ϕ(ξₖ) - gₖ
        #     where gₖ := Qₖ ∇ϕ(ξₖ) with Qₖ ⪰ 0 is obtained via the two-loop recursion as in [60, Algorithm 7.4]
        copyto!(tmpv, ∇ϕ)                               # q ← ∇fₖ
        @views for i in k-1:-1:max(1, k - mem)          # for i = k-1, k-2, ..., k-m
            idx = mod1(i, mem)
            α[idx] = ρs[idx] * dot(S[:, idx], tmpv)     # αᵢ ← ρᵢ sᵢᵀq
            tmpv .-= α[idx] .* Y[:, idx]                # q ← q - αᵢyᵢ
        end
        if k > 1
            idx = mod1(k -1, mem)                       # set Hₖ⁰ = γₖI where γₖ = sₖ₋₁ᵀyₖ₋₁ / (yₖ₋₁ᵀyₖ₋₁)
            @views lmul!(dot(S[:, idx], Y[:, idx]) / norm(Y[:, idx])^2, tmpv) # r ← Hₖ⁰q
        end
        @views for i = max(1, k - mem):k-1              # for i = k-m, k-m+1, ..., k-1
            idx = mod1(i, mem)
            tmpv .+= S[:, idx] .* (α[idx] - ρs[idx] * dot(Y[:, idx], tmpv)) # β ← ρᵢyᵢᵀr; r ← r + sᵢ(αᵢ - β)
        end
        tmpv .= .-tmpv .- β .* ∇ϕ # tmpv now corresponds to dₖ
        #     If ‖dₖ‖ ≥ K, then choose dₖ = -βₖ ∇ϕ(ξₖ) (i.e., set Qₖ = 0)
        if norm(tmpv) ≥ K
            tmpv .= (-β) .* ∇ϕ
        end
        # [Step 3]: store ∇ϕ in Y, so that we can re-use ∇ϕ as a temporary in this iteration.
        idx = mod1(k, mem)
        copyto!(@view(Y[:, idx]), ∇ϕ) # just store, the difference is calculated in the next iteration
        # Step 2 (Line search):
        #     Set αₖ = ρ^(mₖ) where mₖ is the smallest nonnegative integer m such that
        #     ϕ(ξₖ + ρ^k dₖ) ≤ ϕ(ξₖ) + μ ρ^m ⟨∇ϕ(ξₖ), dₖ⟩.
        #              ^ does not make sense, should be ρ^(mₖ).
        #     Note that this is just the Armijo condition, so we simply employ a line search algorithm here.
        function eval_ϕ(α)
            local ϕ = zero(R)
            for (tmpmᵢ, Zᵢ, Astarᵢ) in zip(tmpm, Z, Astar)
                copyto!(tmpmᵢ, Zᵢ)
                mul!(tmpmᵢ, Astarᵢ, ξ, true, true)
                mul!(tmpmᵢ, Astarᵢ, tmpv, α, true)
                psdproject!(tmpmᵢ)
                ϕ += norm(tmpmᵢ)^2
            end
            return ϕ/R(2) - sum(bᵢ * (ξᵢ + α * tmpvᵢ) for (bᵢ, ξᵢ, tmpvᵢ) in zip(b, ξ, tmpv), init=zero(R))
        end
        function eval_dϕ(α)
            ∇ϕ .= .-b
            for (tmpmᵢ, Zᵢ, Aᵢ, Astarᵢ) in zip(tmpm, Z, A, Astar)
                copyto!(tmpmᵢ, Zᵢ)
                mul!(tmpmᵢ, Astarᵢ, ξ, true, true)
                mul!(tmpmᵢ, Astarᵢ, tmpv, α, true)
                psdproject!(tmpmᵢ)
                mul!(∇ϕ, Aᵢ, tmpmᵢ, true, true)
            end
            return dot(∇ϕ, tmpv)
        end
        function eval_ϕdϕ(α)
            local ϕ = zero(R)
            ∇ϕ .= .-b
            for (tmpmᵢ, Zᵢ, Aᵢ, Astarᵢ) in zip(tmpm, Z, A, Astar)
                copyto!(tmpmᵢ, Zᵢ)
                mul!(tmpmᵢ, Astarᵢ, ξ, true, true)
                mul!(tmpmᵢ, Astarᵢ, tmpv, α, true)
                psdproject!(tmpmᵢ)
                ϕ += norm(tmpmᵢ)^2
                mul!(∇ϕ, Aᵢ, tmpmᵢ, true, true)
            end
            return ϕ/R(2) - sum(bᵢ * (ξᵢ + α * tmpvᵢ) for (bᵢ, ξᵢ, tmpvᵢ) in zip(b, ξ, tmpv), init=zero(R)), dot(∇ϕ, tmpv)
        end
        αₖ, _ = linesearch(eval_ϕ, eval_dϕ, eval_ϕdϕ, ρ, ϕ, dot(∇ϕ, tmpv))
        # Step 3 (Update memory):
        #     Compute and save uₖ = ξₖ₊₁ - ξₖ = αₖ dₖ and wₖ = ∇ϕ(ξₖ₊₁) - ∇ϕ(ξₖ).
        S[:, idx] .= αₖ .* tmpv
        # We already stored the true value of ∇ϕ in Y, the actual content will be calculated at the next iteration.
        # [Step2]    Compute ξₖ₊₁ = ξₖ + αₖ dₖ
        ξ .+= @view(S[:, idx])
        #     If k > mem, discard the vector {uₖ₋ₘₑₘ, wₖ₋ₘₑₘ} from storage.
        # this is a no-op, as we cycle through our storage
        # Until: ηproj(Wₖ₊₁, ξₖ₊₁) ≤ tol with Wₖ₊₁ = Π(-A* ξₖ₊₁ - Z)  (but the positive part is expected to be of low rank)
        for (Wᵢ, tmpmᵢ, Zᵢ, Astarᵢ) in zip(W, tmpm, Z, Astar)
            copyto!(Wᵢ, Zᵢ)
            mul!(Wᵢ, Astarᵢ, ξ, true, true)
            copyto!(tmpmᵢ, Wᵢ)
            psdproject!(tmpmᵢ)
            Wᵢ .= tmpmᵢ .- Wᵢ
        end
        ηp = ηproj(ssd, W, ξ, false) # ηproj will use tmpm2, tmpvlong, and tmpv as temporaries
        success = ηp ≤ tol
        (success || k == maxiter || isone(k % 50)) && @verbose_info(@sprintf("%9d | %11g", k, ηp))
        success && return
        if !iszero(maxiter) && k == maxiter
            @verbose_info("Maximum iteration count reached")
            return
        end
        # Output Wₖ₊₁, ξₖ₊₁
        k += 1
    end
end