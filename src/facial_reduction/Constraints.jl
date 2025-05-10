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
# This is again re-cast as obj(x) + prefactor(x) * u₁ + ∑ᵢ ... - ∑ⱼ ... - ∑ₖ ... ∈ Σ
# Now, the decision variables u must contain the coefficients of the pᵢ, σⱼ, Sₖ, so there are many:
# w = [1, 0, ...], g₀₀(x) = obj(x), g₀(x) = [prefactor(x), zeroᵢ(x) * monomials in pᵢ...,
#                                            nonnegⱼ(x) * monomials in σⱼ...,
#                                            psdₖ(x)ₘₙ * monomials in Sₖₘₙ for every entry...],
# gⱼ(x) = [0, 0..., -monomials in σⱼ...]
# For matrix constraints, we exploit Sₖ ∈ Σmat ⇔ yᵀ Sₖ y ∈ Σ, adding dim Sₖ variables. We then automatically know by the Newton
# polytope (see also Miller, Wang, Guo (2025)) that we never need to choose a basis larger than degree 1 with respect to y.
# ↪ gₖ(x) = [0, 0..., -monomials in yᵀ Sₖ y...]
# If we allow the objective to have more than just a single grouping, bring the second, ..., last to the left side.
struct FacialReductionData{C<:Real,Nr,NrPlusNy,I<:Integer,
                           PScalar<:IntPolynomial{C,Nr,0,<:IntMonomialVector{Nr,0,I}},
                           GScalar<:IntMonomialVector{Nr,0,I},
                           PMatrix<:IntPolynomial{C,NrPlusNy,0,<:IntMonomialVector{NrPlusNy,0,I}},
                           GMatrix<:IntMonomialVector{Nr,0,I}}
    objective::PScalar
    M_obj::Vector{IntMonomialVector{Nr,0,I}}
    M⁺_obj::Vector{IntMonomialVector{Nr,0,I}}
    M²_obj::Vector{GScalar} # one item less - first grouping is not squared!
    prefactor::PScalar

    constr_zero::Vector{PScalar}
    M²_zero::Vector{GScalar}

    constr_nonneg::Vector{PScalar}
    M_nonneg::Vector{Vector{IntMonomialVector{Nr,0,I}}}
    M⁺_nonneg::Vector{Vector{IntMonomialVector{Nr,0,I}}}
    M²_nonneg::Vector{Vector{GScalar}}

    constr_psd::Vector{Matrix{PMatrix}}
    M_psd::Vector{Vector{IntMonomialVector{NrPlusNy,0,I}}}
    M⁺_psd::Vector{Vector{IntMonomialVector{NrPlusNy,0,I}}}
    M²_psd::Vector{Vector{GMatrix}}
end

function FacialReductionData(relaxation::AbstractRelaxation{<:Problem{<:IntPolynomial{C,Nr,0,<:IntMonomialVector{Nr,0,I}}}};
    kwargs...) where {C<:Real,Nr,I<:Integer}
    # First figure out which maximum number of variables we need
    maxsidedim = 0
    for constr in poly_problem(relaxation).constr_psd
        maxsidedim = max(maxsidedim, size(constr, 1))
    end
    return FacialReductionData{C,Nr,Nr+maxsidedim,I}(relaxation; kwargs...)
end

function _change_var_numbers(mv::IntMonomialVector{Nr,0,I}, ::Val{newNr}) where {Nr,I<:Integer,newNr}
    # TODO: we might instead go for ExponentsMultideg, which is more economical in terms of index usage, as we can limit the
    # y variables to maxdeg 2.
    newE = ExponentsAll{newNr,I}()
    @assert(typeof(newNr) == typeof(Nr)) # we should always have Ints here, but we cannot enforce it
    if newNr === Nr
        if mv isa IntMonomialVector{Nr,0,I,Tuple{ExponentsAll{Nr,I},Vector{I}}}
            return mv
        elseif mv isa IntMonomialVector{Nr,0,I,<:Tuple{ExponentsAll{Nr,I},AbstractUnitRange{I}}}
            return IntMonomialVector{Nr,0}(unsafe, newE, collect(mv.indices))
        elseif mv isa IntMonomialVector{Nr,0,I,<:Union{ExponentsDegree{Nr,I},<:Tuple{ExponentsDegree{Nr,I},AbstractUnitRange{I}}}}
            return IntMonomialVector{Nr,0}(unsafe, newE, collect(range(convert_index(newE, mv.e, first(mv.indices)),
                                                                       convert_index(newE, mv.e, last(mv.indices)))))
        else
            return IntMonomialVector{Nr,0}(unsafe, newE, convert_index.(newE, mv.e, mv.indices))
        end
    end
    result = FastVec{I}(buffer=length(mv))
    expBuf = Vector{Int}(undef, max(Nr, newNr))
    fill!(@view(expBuf[Nr+1:newNr]), 0)
    oldBuf = @view(expBuf[1:Nr])
    newBuf = @view(expBuf[1:newNr])
    for m in IntPolynomials.veciter(mv, oldBuf)
        unsafe_push!(result, exponents_to_index(newE, newBuf))
    end
    if newNr < Nr
        Base._groupedunique!(sort!(finish!(result)))
    end
    return IntMonomialVector{newNr,0}(unsafe, newE, finish!(result))
end
_change_var_numbers(p::IntPolynomial, newNr::Val) =
    IntPolynomial(p.coeffs, _change_var_numbers(p.monomials, newNr))

function FacialReductionData{C,Nr,NrPlusNy,I}(relaxation::AbstractRelaxation{<:Problem{<:IntPolynomial{C,Nr,0,<:IntMonomialVector{Nr,0,I}}}};
    verbose::Bool=false) where {C<:Real,Nr,NrPlusNy,I<:Integer}
    prob = poly_problem(relaxation)
    g = Relaxation.groupings(relaxation)
    @verbose_info("Assembling data for facial reduction")

    GScalar = IntMonomialVector{Nr,0,I,Tuple{ExponentsAll{Nr,I},Vector{I}},IntMonomial{Nr,0,I,ExponentsAll{Nr,I}}}
    GMatrix = IntMonomialVector{NrPlusNy,0,I,Tuple{ExponentsAll{NrPlusNy,I},Vector{I}},
                                IntMonomial{NrPlusNy,0,I,ExponentsAll{NrPlusNy,I}}}
    eScalar = ExponentsAll{Nr,I}()
    eMatrix = ExponentsAll{NrPlusNy,I}()
    conv_time = @elapsed begin
        M_obj = copy(g.obj) # no deep copy necessary, it mutated, the first level will be overwritten entirely
        M⁺_obj = similar(M_obj)
        M²_obj = Vector{GScalar}(undef, length(g.obj) -1)
    end
    @verbose_info("├ objective and prefactor (", conv_time, " seconds)")

    conv_time = @elapsed begin
        # these will never be updated, so generate them once
        M²_zero = Vector{GScalar}(undef, length(prob.constr_zero))
        @inbounds for (i, groupings) in enumerate(g.zeros)
            unique_groupings = Set{FastKey{I}}()
            for groupingⱼ in groupings
                unique_outer_groupings(groupingⱼ, result=unique_groupings, e=eScalar)
            end
            M²_zero[i] = IntMonomialVector{Nr,0}(IntPolynomials.unsafe, eScalar, sort!(convert.(I, unique_groupings)))
        end
    end
    isempty(prob.constr_zero) || @verbose_info("├ zero constraints (", conv_time, " seconds)")

    conv_time = @elapsed begin
        M_nonneg = Vector{Vector{IntMonomialVector{Nr,0,I}}}(undef, length(prob.constr_nonneg))
        M⁺_nonneg = similar(M_nonneg)
        M²_nonneg = Vector{Vector{GScalar}}(undef, length(prob.constr_nonneg))
        @inbounds for (i, groupings) in enumerate(g.nonnegs)
            M_nonneg[i] = copy(groupings)
            M⁺_nonneg[i] = similar(groupings)
            M²_nonneg[i] = Vector{GScalar}(undef, length(groupings))
        end
    end
    isempty(prob.constr_nonneg) || @verbose_info("├ nonnegative constraints (", conv_time, " seconds)")

    conv_time = @elapsed begin
        # we don't need to assign the lower triangle at all, so let's do it explicitly
        constr_psd = Vector{Matrix{IntPolynomial{C,NrPlusNy,0,GMatrix}}}(undef, length(prob.constr_psd))
        @inbounds for (i, psdᵢ) in enumerate(prob.constr_psd)
            constr_psd[i] = constr_psdᵢ = similar(psdᵢ, IntPolynomial{C,NrPlusNy,0,GMatrix})
            for col in 1:size(psdᵢ, 2), row in 1:col
                constr_psdᵢ[row, col] = _change_var_numbers(psdᵢ[row, col], Val(NrPlusNy))
            end
        end
        M⁺_psd = similar(M_psd)
        M²_psd = Vector{Vector{GMatrix}}(undef, length(prob.constr_psd))
        tmp = Vector{Int}(undef, NrPlusNy)
        @inbounds fill!(@view(tmp[Nr+1:NrPlusNy-1]), 0)
        @inbounds tmp[end] = 1
        @inbounds for (i, (groupings, constr)) in enumerate(zip(g.psds, prob.constr_psd))
            M_psd[i] = dataᵢ = Vector{IntMonomialVector{NrPlusNy,0,I}}(undef, length(groupings))
            M⁺_psd[i] = similar(dataᵢ)
            M²_psd[i] = Vector{GMatrix}(undef, length(groupings))
            dim = size(constr, 1)
            for (j, groupingⱼ) in enumerate(groupings)
                # groupingⱼ is unique groupings, but with the old variable number. It must be combined with the new auxilliary
                # variables - linearly with everyone. The index order is degrevlex; therefore, if we say y = [yₙ, yₙ₋₁, ..., y₁]
                # and add the length-dim range starting from the first index, we get all of them, calculating only one.
                # Note that when squaring, this will automatically yield a continuous range in the columnwise upper triangle
                # for every grouping.
                new_groupings = FastVec{I}(buffer=length(groupingⱼ) * dim)
                @inbounds for mon in IntPolynomials.veciter(groupingⱼ, @view(tmp[1:Nr]))
                    startidx = exponents_to_index(eMatrix, mon)
                    for idx in startidx:startidx+dim-1
                        unsafe_push!(new_groupings, idx)
                    end
                end
                @assert(issorted(new_groupings)) # if we create it in this way, it is already sorted
                dataᵢ[j] = IntMonomialVector{NrPlusNy,0}(IntPolynomials.unsafe, eMatrix, finish!(new_groupings))
            end
        end
    end
    isempty(prob.constr_psd) || @verbose_info("├ PSD constraints (", conv_time, " seconds)")

    @verbose_info("└ Conversion finished")
    return FacialReductionData(
        objective, M_obj, M⁺_obj, M²_obj, prefactor,
        constr_zero, M²_zero,
        constr_nonneg, M_nonneg, M⁺_nonneg, M²_nonneg,
        constr_psd, M_psd, M⁺_psd
    )
end

function Base.show(io::IO, m::MIME"text/plain", data::FacialReductionData)
    size_psd = Dict{Int,Int}()
    for grouping in data.M_obj
        s = length(grouping)
        size_psd[s] = get(size_psd, s, 0) +1
    end
    for groupings in data.M_nonneg, grouping in groupings
        s = length(grouping)
        size_psd[s] = get(size_psd, s, 0) +1
    end
    for (constr, groupings) in zip(data.constr_psd, data.M_psd)
        dim = size(constr, 1)
        for grouping in groupings
            s = length(grouping) * dim
            size_psd[s] = get(size_psd, s, 0) +1
        end
    end
    println(io, "Facial reduction step with PSD block sizes ", sort!(collect(size_psd), rev=true))
end

_sumlen(x) = intsum(length, x)

function Base.getproperty(data::FacialReductionData, field::Symbol)
    field === :numM⁺ && return _sumlen(data.M⁺_obj) + intsum(_sumlen, data.M⁺_nonneg) + intsum(_sumlen, data.M⁺_psd)
    field === :numSOS && return 1 + _sumlen(data.constr_nonneg) + _sumlen(data.constr_psd)
    field === :numMons && return (iszero(data.prefactor) ? 1 : 2) + _sumlen(data.M²_obj) +
        _sumlen(data.M²_zero) + intsum(_sumlen, data.M²_nonneg) + intsum(_sumlen, data.M²_psd)
    return getfield(data, field)
end

function updateM²!(data::FacialReductionData{<:Any,Nr,NrPlusNy,I}; verbose::Bool=false) where {Nr,NrPlusNy,I<:Integer}
    @verbose_info("Populating squared monomials")
    eScalar = ExponentsAll{Nr,I}()
    eMatrix = ExponentsAll{NrPlusNy,I}()

    conv_time = @elapsed begin
        @inbounds for (j, groupingⱼ) in enumerate(Iterators.drop(data.M_obj, 1))
            data.M²_obj[j] = IntMonomialVector{Nr,0}(
                IntPolynomials.unsafe, eScalar, sort!(convert.(I, unique_outer_groupings(groupingⱼ, e=eScalar)[1]))
            )
        end
        # we don't assign the first one
    end
    length(data.M_obj) > 1 && @verbose_info("├ objective (", conv_time, " seconds)")

    conv_time = @elapsed begin
        @inbounds for (i, groupings) in enumerate(data.M_nonneg)
            for (j, groupingⱼ) in enumerate(groupings)
                data.M²_nonneg[i][j] = IntMonomialVector{Nr,0}(
                    IntPolynomials.unsafe, eScalar, sort!(convert.(I, unique_outer_groupings(groupingⱼ, e=eScalar)[1]))
                )
            end
        end
    end
    isempty(data.M_nonneg) || @verbose_info("├ nonnegative constraints (", conv_time, " seconds)")

    conv_time = @elapsed begin
        @inbounds for (i, groupings) in enumerate(data.M_psd)
            for (j, groupingⱼ) in enumerate(groupings)
                data.M²_psd[i][j] = IntMonomialVector{NrPlusNy,0}(
                    IntPolynomials.unsafe, eNatrix, sort!(convert.(I, unique_outer_groupings(groupingⱼ, e=eMatrix)[1]))
                )
            end
        end
    end
    isempty(data.M_psd) || @verbose_info("├ PSD constraints (", conv_time, " seconds)")

    @verbose_info("└ End populating squared monomials")
    return data
end

function Relaxation.embed(data::FacialReductionData, old::Relaxation.RelaxationGroupings{Nr,0,I}; verbose::Bool) where {Nr,I<:Integer}
    @verbose_info("Embedding facial reduction data")
    conv_time = @elapsed begin
        psds = Vector{Vector{IntMonomialVector{Nr,0,I}}}(undef, length(data.M_psd))
        @inbounds for (i, groupings) in enumerate(data.M_psd)
            psds[i] = psdᵢ = Vector{IntMonomialVector{Nr,0,I}}(undef, length(groupings))
            for (j, grouping) in enumerate(groupings)
                psdᵢ[j] = _change_var_numbers(grouping, Val(Nr))
            end
        end
    end
    @verbose_info("Embedding finished in ", conv_time, " seconds")
    return Relaxation.RelaxationGroupings(data.M_obj, old.zeros, data.M_nonneg, psds, old.var_cliques)
end

Base.IteratorSize(::Type{<:FacialReductionData}) = Base.HasLength()
Base.IteratorEltype(::Type{<:FacialReductionData}) = Base.HasEltype()
Base.eltype(::Type{<:FacialReductionData}) = AbstractFacialReductionDataSlice
Base.length(data::FacialReductionData) = data.numMons

Base.iterate(data::FacialReductionData) = FacialReductionDataSliceFirst(data), (:obj, 1, 2, iszero(data.prefactor) ? 2 : 3)
function Base.iterate(data::FacialReductionData{C}, (pos, idxouter, idxinner, k)::Tuple{Symbol,Int,Int,Int}) where {C}
    if pos === :obj
        @assert(isone(idxouter))
        if idxinner ≤ length(data.M_obj)
            mons = @inbounds data.M²_obj[idxinner-1]
            @inbounds return FacialReductionDataSliceRest{C,false}(data.M_obj[idxinner], data.M⁺_obj[idxinner], mons, k),
                             (:obj, 1, idxinner +1, k + length(mons))
        end
        k += _sumlen(data.M²_zero)
        pos = :nonneg
        idxinner = 1
    end
    if pos === :nonneg
        while idxouter ≤ length(data.M_nonneg)
            grouping = @inbounds data.M_nonneg[idxouter]
            if idxinner ≤ length(grouping)
                mons = @inbounds data.M²_nonneg[idxouter][idxinner]
                @inbounds return FacialReductionDataSliceRest{C,false}(grouping[idxinner], data.M⁺_nonneg[idxouter][idxinner],
                                                                       mons, k),
                                 (:nonneg, idxouter, idxinner +1, k + length(mons))
            else
                idxouter += 1
                idxinner = 1
            end
        end
        pos = :psd
        idxouter = 1
        idxinner = 1
    end
    if pos === :psd
        while idxouter ≤ length(data.M_psd)
            grouping = @inbounds data.M_psd[idxouter]
            if idxinner ≤ length(grouping)
                mons = @inbounds data.M²_psd[idxouter][idxinner]
                @inbounds return FacialReductionDataSliceRest{C,true}(grouping[idxinner], data.M⁺_psd[idxouter][idxinner],
                                                                      mons, k, size(data.constr_psd[idxouter], 1)),
                                 (:psd, idxouter, idxinner +1, k + length(mons))
            else
                idxouter += 1
                idxinner = 1
            end
        end
        return nothing
    end
    @assert(false)
    return nothing # to avoid unreachable exceptions
end

abstract type AbstractFacialReductionDataSlice{M<:IntMonomialVector,M⁺<:IntMonomialVector} end

Base.IteratorSize(::Type{<:AbstractFacialReductionDataSlice}) = Base.HasLength()
Base.IteratorEltype(::Type{<:AbstractFacialReductionDataSlice}) = Base.HasEltype()

struct FacialReductionDataSliceFirst{M<:IntMonomialVector,M⁺<:IntMonomialVector,D<:FacialReductionData} <: AbstractFacialReductionDataSlice{M,M⁺}
    data::D
    extdegM::Tuple{Int,Int}

    FacialReductionDataSliceFirst(data::D) where {D<:FacialReductionData} =
        new{typeof(data.M_obj[1]),typeof(data.M⁺_obj[1]),D}(data, extdegree(data.M_obj[1]))
end

function Base.getproperty(frs::FacialReductionDataSliceFirst{M,M⁺}, field::Symbol) where {M<:IntMonomialVector,M⁺<:IntMonomialVector}
    field === :M && @inbounds return frs.data.M_obj[1]::M
    field === :M⁺ && @inbounds return frs.data.M⁺_obj[1]::M⁺
    return getfield(frs, field)
end

Base.eltype(::Type{FacialReductionDataSliceFirst{<:IntMonomialVector,<:IntMonomialVector,
                                                 FacialReductionData{C,Nr,NrPlusNy,I,PScalar,GScalar,PMatrix,GMatrix}}}) where
    {C<:Real,Nr,NrPlusNy,I<:Integer,
     PScalar<:IntPolynomial{C,Nr,0,<:IntMonomialVector{Nr,0,I}},GScalar<:IntMonomialVector{Nr,0,I},
     PMatrix<:IntPolynomial{C,NrPlusNy,0,<:IntMonomialVector{NrPlusNy,0,I}},GMatrix<:IntMonomialVector{Nr,0,I}} =
    Tuple{Int,Union{PScalar,PMatrix}}
Base.length(frs::FacialReductionDataSliceFirst) = frs.data.numM⁺

@generated _cached_onepolynomial(::FacialReductionDataSliceFirst{<:FacialReductionData{C,Nr}}, m::IntMonomial{Nr,0}) where {C,Nr} =
    :(@inline; return IntPolynomial($([one(C)]), IntMonomialVector{Nr,0}(IntPolynomials.unsafe, m.e, m.index:m.index)))

Base.iterate(frs::FacialReductionDataSliceFirst) = (1, frs.data.objective), (:obj, 1, 0, 1, 1, 1, 1)

function Base.iterate(frs::FacialReductionDataSliceFirst, (pos, idxouter, idxinner, idxmon, idxrow, idxcol, k)::Tuple{Symbol,Int,Int,Int,Int,Int,Int})
    # 1. obj(x)
    # 2. + prefactor(x) * u₁
    # 3. - (other SOS groupings from objective) * u...  [→ -⟨u..., mons⟩ ∈ Σ]
    # 4. + ∑ᵢ zeroᵢ(x) * u...
    # 5. - ∑ⱼ nonnegⱼ(x) * u...                         [→ -⟨u..., mons⟩ ∈ Σ]
    # 6. - ∑ₖ ⟨U..., psdₖ(x)⟩                             [→ -⟨U..., mons⟩ ∈ Σᵈ]
    # ∈ Σ
    if pos === :obj
        @assert(isone(idxouter))
        if iszero(idxinner)
            iszero(frs.data.prefactor) || return (k += 1, frs.data.prefactor), (:obj, idxouter, idxinner +1, 1, 1, 1, k)
            idxinner += 1
        end
        while idxinner ≤ length(frs.data.M²_obj)
            mons = @inbounds frs.data.M²_obj[idxinner]
            if idxmon ≤ length(mons)
                return (k += 1, _cached_onepolynomial(frs, @inbounds(mons[idxmon]))),
                       (:obj, idxouter, idxinner, idxmon +1, 1, 1, k)
            else
                idxinner += 1
                idxmon = 1
            end
        end
        pos = :zero
        idxinner = 1
        @assert(isone(idxmon))
    end
    if pos === :zero
        @assert(isone(idxinner))
        while idxouter ≤ length(frs.data.M²_zero)
            mons = @inbounds frs.data.M²_zero[idxouter]
            if idxmon ≤ length(mons)
                constr = @inbounds frs.data.constr_zero[idxouter]
                return (k += 1, IntPolynomial(constr.coeffs, constr.monomials * @inbounds(mons[idxmon]))),
                       (:zero, idxouter, idxinner, idxmon +1, 1, 1, k)
            else
                idxouter += 1
                idxmon = 1
            end
        end
        pos = :nonneg
        idxouter = 1
        @assert(isone(idxmon))
    end
    if pos === :nonneg
        while idxouter ≤ length(frs.data.M²_nonneg)
            grouping = @inbounds frs.data.M²_nonneg[idxouter]
            while idxinner ≤ length(grouping)
                mons = @inbounds grouping[idxinner]
                if idxmon ≤ length(mons)
                    constr = @inbounds frs.data.constr_nonneg[idxouter]
                    return (k += 1, IntPolynomial(constr.coeffs, constr.monomials * @inbounds(mons[idxmon]))),
                           (:nonneg, idxouter, idxinner, idxmon +1, 1, 1, k)
                else
                    idxinner += 1
                    idxmon = 1
                end
            end
            idxouter += 1
            idxinner = 1
            @assert(isone(idxmon))
        end
        pos = :psd
        idxouter = 1
        @assert(isone(idxinner) && isone(idxmon))
    end
    if pos === :psd
        while idxouter ≤ length(frs.data.M²_psd)
            grouping = @inbounds frs.data.M²_psd[idxouter]
            while idxinner ≤ length(grouping)
                mons = @inbounds grouping[idxinner]
                while idxmon ≤ length(mons)
                    constrs = @inbounds frs.data.constr_psd[idxouter]
                    # monomial order follows upper triangle
                    while idxcol ≤ size(constrs, 1)
                        if idxrow ≤ idxcol
                            constr = @inbounds constrs[idxrow, idxcol]
                            return (k += 1, IntPolynomial(constr.coeffs, constr.monomials * @inbounds(mons[idxmon]))),
                                   (:psd, idxouter, idxinner, idxmon, idxrow +1, idxcol, k)
                        end
                        idxcol += 1
                        idxrow = 1
                    end
                    idxcol = 1
                    idxmon += 1
                    @assert(isone(idxcol) && isone(idxrow))
                end
                idxinner += 1
                idxmon = 1
                @assert(isone(idxcol) && isone(idxrow))
            end
            idxouter += 1
            idxinner = 1
            @assert(isone(idxmon) && isone(idxcol) && isone(idxrow))
        end
        return nothing
    end
    @assert(false)
    return nothing # to avoid unreachable exceptions
end

struct FacialReductionDataSliceRest{C,Matrix,M_<:IntMonomialVector,M⁺_<:IntMonomialVector,Monomials<:IntMonomialVector,IN<:Union{Int,Nothing}} <: AbstractFacialReductionDataSlice{M_,M⁺_}
    M::M_
    M⁺::M⁺_
    M²::Monomials
    k::Int
    dim::IN
    extdegM::Tuple{Int,Int}

    FacialReductionDataSliceRest{C,false}(M::M_, M⁺::M⁺_, M²::Monomials, k::Int) where
        {C,Matrix,M_<:IntMonomialVector,M⁺_<:IntMonomialVector,Monomials<:IntMonomialVector} =
        new{C,Matrix,M_,M⁺_,Monomials,Nothing}(M, M⁺, M², k, nothing, extdegree(M))
    FacialReductionDataSliceRest{C,true}(M::M_, M⁺::M⁺_, M²::Monomials, k::Int, dim::Int) where
        {C,Matrix,M_<:IntMonomialVector,M⁺_<:IntMonomialVector,Monomials<:IntMonomialVector} =
        new{C,Matrix,M_,M⁺_,Monomials,Int}(M, M⁺, M², k, extdegree(M), dim)
end

Base.eltype(::Type{<:FacialReductionDataSliceRest{C,<:Any,<:IntMonomialVector,<:IntMonomialVector,Monomials}}) where
    {C,Nr,Monomials<:IntMonomialVector{Nr,0}} = Tuple{Int,IntPolynomial{C,Nr,0,Monomials}}
Base.length(frs::FacialReductionDataSliceRest) = length(frs.monomials)

@generated _cached_minusonepolynomial(::FacialReductionDataSliceRest{C}, m::IntMonomial{Nr,0}) where {C,Nr} =
    :(@inline; return IntPolynomial($([-one(C)]), IntMonomialVector{Nr,0}(IntPolynomials.unsafe, m.e, m.index:m.index)))
@generated _cached_minusonehalfpolynomial(::FacialReductionDataSliceRest{C,true}, m::IntMonomial{Nr,0}) where {C,Nr} =
    :(@inline; return IntPolynomial($([-inv(C(2))]), IntMonomialVector{Nr,0}(IntPolynomials.unsafe, m.e, m.index:m.index)))

function Base.iterate(frs::FacialReductionDataSliceRest{<:Any,false}, i::Int=0)
    i ≥ length(frs.M²) && return nothing
    return (frs.k + i, _cached_minusonepolynomial(frs, @inbounds(frs.M²[i+=1]))), i
end

function Base.iterate(frs::FacialReductionDataSliceRest{<:Any,true}, (i, row, col)::Tuple{Int,Int,Int}=(0, 1, 1))
    i ≥ length(frs.M²) && return nothing
    return (frs.k + i, (row == col ? _cached_minusonepolynomial : _cached_minusonehalfpolynomial)(
                           frs, @inbounds(frs.M²[i+=1])
                       )),
           (row == frs.dim ? (i, 1, 1) : ((row == col ? (i, 1, col +1) : (i, row +1, col))))
end