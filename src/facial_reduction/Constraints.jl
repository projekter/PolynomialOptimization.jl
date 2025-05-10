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
struct FacialReductionData{C<:Real,NrPlusNy,I<:Integer,MV<:IntMonomialVector{NrPlusNy,0,I}}
    objective::IntPolynomial{C,NrPlusNy,0,MV}
    M_obj::Vector{MV}
    M⁺_obj::Vector{MV}
    M²_obj::Vector{MV} # one item less - first grouping is not squared!
    prefactor::IntPolynomial{C,NrPlusNy,0,MV}

    constr_zero::Vector{IntPolynomial{C,NrPlusNy,0,MV}}
    M²_zero::Vector{MV}

    constr_nonneg::Vector{IntPolynomial{C,NrPlusNy,0,MV}}
    M_nonneg::Vector{Vector{MV}}
    M⁺_nonneg::Vector{Vector{MV}}
    M²_nonneg::Vector{Vector{MV}}

    constr_psd::Vector{Matrix{IntPolynomial{C,NrPlusNy,0,MV}}}
    M_psd::Vector{Vector{MV}} # these three contain the products with y, i.e., they go in the gₖ that describe the constraints
    M⁺_psd::Vector{Vector{MV}}
    M²_psd::Vector{Vector{MV}}
    M_psd_noy::Vector{Vector{MV}} # these three correspond to the original groupings; they go in the transformed objective
    M⁺_psd_noy::Vector{Vector{MV}}
    M²_psd_noy::Vector{Vector{MV}}
end

function FacialReductionData(relaxation::AbstractRelaxation{<:Problem{<:IntPolynomial{C,Nr,0,<:IntMonomialVector{Nr,0,I}}}};
    kwargs...) where {C<:Real,Nr,I<:Integer}
    # First figure out which maximum number of variables we need
    maxsidedim = 0
    for constr in poly_problem(relaxation).constr_psd
        maxsidedim = max(maxsidedim, size(constr, 1))
    end
    return FacialReductionData{C,Nr+maxsidedim,I}(relaxation; kwargs...)
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
            return IntMonomialVector{Nr,0}(unsafe, newE, convert_index.((newE,), (mv.e,), mv.indices))
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

function FacialReductionData{C,NrPlusNy,I}(relaxation::AbstractRelaxation{<:Problem{<:IntPolynomial{C,Nr,0,<:IntMonomialVector{Nr,0,I}}}};
    verbose::Bool=false) where {C<:Real,Nr,NrPlusNy,I<:Integer}
    prob = poly_problem(relaxation)
    g = Relaxation.groupings(relaxation)
    @verbose_info("Assembling data for facial reduction")

    MV = IntMonomialVector{NrPlusNy,0,I,Tuple{ExponentsAll{NrPlusNy,I},Vector{I}},
                           IntMonomial{NrPlusNy,0,I,ExponentsAll{NrPlusNy,I}}}
    E = ExponentsAll{NrPlusNy,I}()
    conv_time = @elapsed begin
        objective = _change_var_numbers(prob.objective, Val(NrPlusNy))
        prefactor = _change_var_numbers(prob.prefactor, Val(NrPlusNy))
        M_obj = _change_var_numbers.(g.obj, Val(NrPlusNy))
        M⁺_obj = similar(M_obj)
        M²_obj = Vector{MV}(undef, length(g.obj) -1)
    end
    @verbose_info("├ objective and prefactor (", conv_time, " seconds)")

    conv_time = @elapsed begin
        constr_zero = _change_var_numbers.(prob.constr_zero, Val(NrPlusNy))
        # these will never be updated, so generate them once
        M²_zero = Vector{MV}(undef, length(prob.constr_zero))
        # We'll need to convert to NrPlusNy variables and we need to apply unique_outer_groupings. Let's avoid creating a list
        # with the first indices by exploiting that ExponentsMultideg with sufficient degree and the extra variables restricted
        # to zero has exactly the same indices as whatever exponents were used before with less variables.
        @inbounds for (i, groupings) in enumerate(g.zeros)
            unique_groupings = Set{FastKey{I}}()
            for groupingⱼ in groupings
                if Nr == NrPlusNy
                    unique_outer_groupings(groupingⱼ, result=unique_groupings, e=E)
                elseif groupingⱼ.e isa ExponentsAll
                    d = maxdegree(groupingⱼ)
                    upper = Vector{Int}(undef, NrPlusNy)
                    fill!(@view(upper[1:Nr]), d)
                    fill!(@view(upper[Nr+1:end]), 0)
                    newE = ExponentsMultideg{NrPlusNy,I}(0, d, ConstantVector(0, NrPlusNy), upper)
                    unique_outer_groupings(
                        IntMonomialVector{NrPlusNy,0}(IntPolynomials.unsafe, newE, groupingⱼ.indices),
                        result=unique_groupings, e=E
                    )
                elseif groupingⱼ.e isa ExponentsDegree
                    upper = Vector{Int}(undef, NrPlusNy)
                    fill!(@view(upper[1:Nr]), groupingⱼ.e.maxdeg)
                    fill!(@view(upper[Nr+1:end]), 0)
                    newE = ExponentsMultideg{NrPlusNy,I}(
                        groupingⱼ.e.mindeg, groupingⱼ.e.maxdeg, ConstantVector(0, NrPlusNy), upper
                    )
                    if groupingⱼ isa IntPolynomials.IntMonomialVectorComplete
                        unique_outer_groupings(
                            IntMonomialVector{NrPlusNy,0}(IntPolynomials.unsafe, newE),
                            result=unique_groupings, e=E
                        )
                    else
                        unique_outer_groupings(
                            IntMonomialVector{NrPlusNy,0}(IntPolynomials.unsafe, newE, groupingⱼ.indices),
                            result=unique_groupings, e=E
                        )
                    end
                else
                    lower = Vector{Int}(undef, NrPlusNy)
                    copyto!(lower, groupingⱼ.e.minmultideg)
                    fill!(@view(lower[Nr+1:end]), 0)
                    upper = Vector{Int}(undef, NrPlusNy)
                    copyto!(upper, groupingⱼ.e.maxmultideg)
                    fill!(@view(upper[Nr+1:end]), 0)
                    newE = ExponentsMultideg{NrPlusNy,I}(grouping.e.mindeg, grouping.e.maxdeg, lower, upper)
                    if groupingⱼ isa IntPolynomials.IntMonomialVectorComplete
                        unique_outer_groupings(
                            IntMonomialVector{NrPlusNy,0}(IntPolynomials.unsafe, newE),
                            result=unique_groupings, e=E
                        )
                    else
                        unique_outer_groupings(
                            IntMonomialVector{NrPlusNy,0}(IntPolynomials.unsafe, newE, groupingⱼ.indices),
                            result=unique_groupings, e=E
                        )
                    end
                end
            end
            M²_zero[i] = IntMonomialVector{NrPlusNy,0}(IntPolynomials.unsafe, E, sort!(convert.(I, unique_groupings)))
        end
    end
    isempty(prob.constr_zero) || @verbose_info("├ zero constraints (", conv_time, " seconds)")

    conv_time = @elapsed begin
        constr_nonneg = _change_var_numbers.(prob.constr_nonneg, Val(NrPlusNy))
        M_nonneg = Vector{Vector{MV}}(undef, length(prob.constr_nonneg))
        M⁺_nonneg = similar(M_nonneg)
        M²_nonneg = similar(M_nonneg)
        @inbounds for (i, groupings) in enumerate(g.nonnegs)
            M_nonneg[i] = _change_var_numbers.(groupings, Val(NrPlusNy))
            M⁺_nonneg[i] = similar(M_nonneg[i])
            M²_nonneg[i] = similar(M_nonneg[i])
        end
    end
    isempty(prob.constr_nonneg) || @verbose_info("├ nonnegative constraints (", conv_time, " seconds)")

    conv_time = @elapsed begin
        # we don't need to assign the lower triangle at all, so let's do it explicitly
        constr_psd = Vector{Matrix{IntPolynomial{C,NrPlusNy,0,MV}}}(undef, length(prob.constr_psd))
        @inbounds for (i, psdᵢ) in enumerate(prob.constr_psd)
            constr_psd[i] = constr_psdᵢ = similar(psdᵢ, IntPolynomial{C,NrPlusNy,0,MV})
            for col in 1:size(psdᵢ, 2), row in 1:col
                constr_psdᵢ[row, col] = _change_var_numbers(psdᵢ[row, col], Val(NrPlusNy))
            end
        end
        # For the matrix constraints, we convert to the scalar representation: MV ∈ Σmat[x] ⇔ yᵀ MV y ∈ Σ[x, y]. So every unique
        # grouping must be combined with the new auxilliary variables - linearly with everyone. The index order is degrevlex;
        # therefore, if we say y = [yₙ, yₙ₋₁, ..., y₁] and add the unit range of length dim(MV) starting from the first index, we
        # get all of them, calculating only one. When squaring, this will automatically yield a continuous range in the
        # columnwise upper triangle for every grouping.
        M_psd = Vector{Vector{MV}}(undef, length(prob.constr_psd))
        M⁺_psd = similar(M_psd)
        M²_psd = similar(M_psd)
        M_psd_noy = similar(M_psd)
        M⁺_psd_noy = similar(M_psd)
        M²_psd_noy = similar(M_psd)

        # function barrier - specialize on groupingⱼ
        function matrix_to_scalar_groupings(tmp, groupingⱼ, E, dim)
            new_groupings = FastVec{I}(buffer=length(groupingⱼ) * dim)
            new_groupings_noy = FastVec{I}(buffer=length(groupingⱼ))
            tmpred = @view(tmp[1:Nr])
            @inbounds for mon in IntPolynomials.veciter(groupingⱼ, tmpred)
                startidx = exponents_to_index(E, tmp)
                for idx in startidx:startidx+dim-1
                    unsafe_push!(new_groupings, idx)
                end
                unsafe_push!(new_groupings_noy, exponents_to_index(E, tmpred))
            end
            @assert(issorted(new_groupings)) # if we create it in this way, it is already sorted
            return finish!(new_groupings), finish!(new_groupings_noy)
        end

        tmp = Vector{Int}(undef, NrPlusNy)
        @inbounds fill!(@view(tmp[Nr+1:NrPlusNy-1]), 0)
        @inbounds tmp[end] = 1
        @inbounds for (i, (groupings, constr)) in enumerate(zip(g.psds, prob.constr_psd))
            M_psd[i] = dataᵢ = Vector{MV}(undef, length(groupings))
            M⁺_psd[i] = similar(dataᵢ)
            M²_psd[i] = similar(dataᵢ)
            M_psd_noy[i] = data_noyᵢ = similar(dataᵢ)
            M⁺_psd_noy[i] = similar(dataᵢ)
            M²_psd_noy[i] = similar(dataᵢ)
            dim = size(constr, 1)
            for (j, groupingⱼ) in enumerate(groupings)
                idxᵢ, idx_noyᵢ = matrix_to_scalar_groupings(tmp, groupingⱼ, E, dim)
                dataᵢ[j] = IntMonomialVector{NrPlusNy,0}(IntPolynomials.unsafe, E, idxᵢ)
                data_noyᵢ[j] = IntMonomialVector{NrPlusNy,0}(IntPolynomials.unsafe, E, idx_noyᵢ)
            end
        end
    end
    isempty(prob.constr_psd) || @verbose_info("├ PSD constraints (", conv_time, " seconds)")

    @verbose_info("└ Conversion finished")
    return FacialReductionData(
        objective, M_obj, M⁺_obj, M²_obj, prefactor,
        constr_zero, M²_zero,
        constr_nonneg, M_nonneg, M⁺_nonneg, M²_nonneg,
        constr_psd, M_psd, M⁺_psd, M²_psd, M_psd_noy, M⁺_psd_noy, M²_psd_noy
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
    field === :numSOS && return 1 + _sumlen(data.M²_nonneg) + _sumlen(data.M²_psd)
    field === :numMons && return (iszero(data.prefactor) ? 1 : 2) + _sumlen(data.M²_obj) +
        _sumlen(data.M²_zero) + intsum(_sumlen, data.M²_nonneg) + intsum(_sumlen, data.M²_psd)
    return getfield(data, field)
end

function updateM²!(data::FacialReductionData{<:Any,NrPlusNy,I}; verbose::Bool=false) where {NrPlusNy,I<:Integer}
    @verbose_info("Populating squared monomials")
    E = ExponentsAll{NrPlusNy,I}()

    conv_time = @elapsed begin
        @inbounds for (j, groupingⱼ) in enumerate(Iterators.drop(data.M_obj, 1))
            data.M²_obj[j] = IntMonomialVector{NrPlusNy,0}(
                IntPolynomials.unsafe, E, sort!(convert.(I, unique_outer_groupings(groupingⱼ, e=E)[1]))
            )
        end
        # we don't assign the first one
    end
    length(data.M_obj) > 1 && @verbose_info("├ objective (", conv_time, " seconds)")

    conv_time = @elapsed begin
        @inbounds for (i, groupings) in enumerate(data.M_nonneg)
            for (j, groupingⱼ) in enumerate(groupings)
                data.M²_nonneg[i][j] = IntMonomialVector{NrPlusNy,0}(
                    IntPolynomials.unsafe, E, sort!(convert.(I, unique_outer_groupings(groupingⱼ, e=E)[1]))
                )
            end
        end
    end
    isempty(data.M_nonneg) || @verbose_info("├ nonnegative constraints (", conv_time, " seconds)")

    conv_time = @elapsed begin
        @inbounds for (i, (groupings, groupings_noy)) in enumerate(zip(data.M_psd, data.M_psd_noy))
            for (j, (groupingⱼ, grouping_noyⱼ)) in enumerate(zip(groupings, groupings_noy))
                data.M²_psd[i][j] = IntMonomialVector{NrPlusNy,0}(
                    IntPolynomials.unsafe, E, sort!(convert.(I, unique_outer_groupings(groupingⱼ, e=E)[1]))
                )
                data.M²_psd_noy[i][j] = IntMonomialVector{NrPlusNy,0}(
                    IntPolynomials.unsafe, E, sort!(convert.(I, unique_outer_groupings(grouping_noyⱼ, e=E)[1]))
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
        obj = Vector{IntMonomialVector{Nr,0,I}}(undef, length(data.M_obj))
        @inbounds for (i, grouping) in enumerate(data.M_obj)
            obj[i] = _change_var_numbers(grouping, Val(Nr))
        end

        nonneg = Vector{Vector{IntMonomialVector{Nr,0,I}}}(undef, length(data.M_nonneg))
        @inbounds for (i, groupings) in enumerate(data.M_nonneg)
            nonneg[i] = nonnegᵢ = Vector{IntMonomialVector{Nr,0,I}}(undef, length(groupings))
            for (j, grouping) in enumerate(groupings)
                nonnegᵢ[j] = _change_var_numbers(grouping, Val(Nr))
            end
        end

        # TODO: here, we could do better. We can define groupings per side index (i.e. groupings is a vector of length dim),
        # and we populate the grouping by selecting only those monomials which have y variables pertaining to this index.
        # This is basically the idea from Miller, Wang, Guo (2025).
        psds = Vector{Vector{IntMonomialVector{Nr,0,I}}}(undef, length(data.M_psd))
        @inbounds for (i, groupings) in enumerate(data.M_psd)
            psds[i] = psdᵢ = Vector{IntMonomialVector{Nr,0,I}}(undef, length(groupings))
            for (j, grouping) in enumerate(groupings)
                psdᵢ[j] = _change_var_numbers(grouping, Val(Nr))
            end
        end
    end
    @verbose_info("Embedding finished in ", conv_time, " seconds")
    return Relaxation.RelaxationGroupings(obj, old.zeros, nonneg, psds, old.var_cliques)
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
            @inbounds return FacialReductionDataSliceRest{false,C}(data.M_obj[idxinner], data.M⁺_obj[idxinner], mons, k),
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
                @inbounds return FacialReductionDataSliceRest{false,C}(grouping[idxinner], data.M⁺_nonneg[idxouter][idxinner],
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
                @inbounds return FacialReductionDataSliceRest{true,C}(grouping[idxinner], data.M⁺_psd[idxouter][idxinner],
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

abstract type AbstractFacialReductionDataSlice{P<:IntPolynomial} end

Base.IteratorSize(::Type{<:AbstractFacialReductionDataSlice}) = Base.HasLength()
Base.IteratorEltype(::Type{<:AbstractFacialReductionDataSlice}) = Base.HasEltype()
Base.eltype(::Type{<:AbstractFacialReductionDataSlice{P}}) where {P<:IntPolynomial} = Tuple{Int,P}

struct FacialReductionDataSliceFirst{P<:IntPolynomial,D<:FacialReductionData} <: AbstractFacialReductionDataSlice{P}
    data::D
    extdegM::Tuple{Int,Int}

    FacialReductionDataSliceFirst(data::D) where {C<:Real,NrPlusNy,I<:Integer,MV<:IntMonomialVector{NrPlusNy,0,I},D<:FacialReductionData{C,NrPlusNy,I,MV}} =
        new{IntPolynomial{C,NrPlusNy,0,MV},D}(data, extdegree(data.M_obj[1]))
end

function Base.getproperty(frs::FacialReductionDataSliceFirst, field::Symbol)
    field === :M && @inbounds return frs.data.M_obj[1]
    field === :M⁺ && @inbounds return frs.data.M⁺_obj[1]
    return getfield(frs, field)
end

Base.length(frs::FacialReductionDataSliceFirst) = frs.data.numM⁺

@generated _cached_onepolynomial(::FacialReductionDataSliceFirst{<:IntPolynomial{C,NrPlusNy}}, m::IntMonomial{NrPlusNy,0}) where {C,NrPlusNy} =
    :(@inline; return IntPolynomial($([one(C)]), IntMonomialVector{NrPlusNy,0}(IntPolynomials.unsafe, m.e, m.index:m.index)))

Base.iterate(frs::FacialReductionDataSliceFirst) = (1, frs.data.objective), (:obj, 1, 0, 1, 1, 1, 1)

function Base.iterate(frs::FacialReductionDataSliceFirst, (pos, idxouter, idxinner, idxmon, idxrow, idxcol, k)::Tuple{Symbol,Int,Int,Int,Int,Int,Int})
    # 1. obj(x)
    # 2. + prefactor(x) * u₁
    # 3. - (other SOS groupings from objective) * u...  [→ -⟨u..., mons⟩ ∈ Σ]
    # 4. + ∑ᵢ zeroᵢ(x) * u...
    # 5. - ∑ⱼ nonnegⱼ(x) * u...                         [→ -⟨u..., mons⟩ ∈ Σ]
    # 6. - ∑ₖ ⟨U..., psdₖ(x)⟩                             [→ -⟨U..., mons⟩ ∈ Σᵈ]
    # ∈ Σ
    # Note we shift all the minus signs into the definitions of u, so that we can simply alias the coefficient vectors of the
    # polynomials here. For PSD, we only take the upper triangle and also shift the factor 2 as 1/2 into the monomial
    # definition later.
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
        while idxouter ≤ length(frs.data.M²_psd_noy)
            grouping = @inbounds frs.data.M²_psd_noy[idxouter]
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

struct FacialReductionDataSliceRest{Matrix,C,NrPlusNy,MV<:IntMonomialVector{NrPlusNy,0},IN<:Union{Int,Nothing}} <: AbstractFacialReductionDataSlice{IntPolynomial{C,NrPlusNy,0,MV}}
    M::MV
    M⁺::MV
    M²::MV
    k::Int
    dim::IN
    extdegM::Tuple{Int,Int}

    FacialReductionDataSliceRest{false,C}(M::MV, M⁺::MV, M²::MV, k::Int) where {C,NrPlusNy,MV<:IntMonomialVector{NrPlusNy,0}} =
        new{false,C,NrPlusNy,MV,Nothing}(M, M⁺, M², k, nothing, extdegree(M))
    FacialReductionDataSliceRest{true,C}(M::MV, M⁺::MV, M²::MV, k::Int, dim::Int) where
        {C,NrPlusNy,MV<:IntMonomialVector{NrPlusNy,0}} =
        new{true,C,NrPlusNy,MV,Int}(M, M⁺, M², k, dim, extdegree(M))
end

Base.length(frs::FacialReductionDataSliceRest) = length(frs.M²)

@generated _cached_minusonepolynomial(::FacialReductionDataSliceRest{<:Any,C}, m::IntMonomial{NrPlusNy,0}) where {C,NrPlusNy} =
    :(@inline; return IntPolynomial($([-one(C)]), IntMonomialVector{NrPlusNy,0}(IntPolynomials.unsafe, m.e, m.index:m.index)))
@generated _cached_minusonehalfpolynomial(::FacialReductionDataSliceRest{true,C}, m::IntMonomial{NrPlusNy,0}) where {C,NrPlusNy} =
    :(@inline; return IntPolynomial($([-inv(C(2))]), IntMonomialVector{NrPlusNy,0}(IntPolynomials.unsafe, m.e, m.index:m.index)))

function Base.iterate(frs::FacialReductionDataSliceRest{false,<:Any}, i::Int=0)
    i ≥ length(frs.M²) && return nothing
    return (frs.k + i, _cached_minusonepolynomial(frs, @inbounds(frs.M²[i+=1]))), i
end

function Base.iterate(frs::FacialReductionDataSliceRest{true,<:Any}, (i, row, col)::Tuple{Int,Int,Int}=(0, 1, 1))
    i ≥ length(frs.M²) && return nothing
    return (frs.k + i, (row == col ? _cached_minusonepolynomial : _cached_minusonehalfpolynomial)(
                           frs, @inbounds(frs.M²[i+=1])
                       )),
           (row == frs.dim ? (i, 1, 1) : ((row == col ? (i, 1, col +1) : (i, row +1, col))))
end