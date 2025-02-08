export moment_add_matrix!, moment_add_equality!, moment_setup!
# This file is for the commuting case; nevertheless, we already write the monomial index/multiplication calculation with a
# possible extension to the noncommuting case in mind.

function to_soc!(indices, values, lens, supports_rotated)
    @assert(length(lens) ≥ 2)
    @inbounds if supports_rotated
        for j in length(lens):-1:3
            iszero(lens[j]) && deleteat!(lens, j)
        end
    else
        range₁ = 1:lens[1]
        range₂ = lens[1]+1:lens[1]+lens[2]
        prevlen = last(range₂)
        inds = @views count_uniques(indices[range₁], indices[range₂])
        total = 2inds
        if total > prevlen
            Base._growat!(indices, prevlen +1, total - prevlen)
            Base._growat!(values, prevlen +1, total - prevlen)
        end
        # we must make a copy to avoid potential overwrites, doesn't matter in which direction we work
        oldinds = indices[1:prevlen]
        oldvals = values[1:prevlen]
        count_uniques(@view(oldinds[range₁]), @view(oldinds[range₂]), let oldinds=oldinds, oldvals=oldvals
            (o₁, i₁, i₂) -> @inbounds begin
                if ismissing(i₁)
                    indices[o₁+inds] = indices[o₁] = oldinds[lens[1]+i₂]
                    values[o₁+inds] = -(values[o₁] = oldvals[lens[1]+i₂])
                elseif ismissing(i₂)
                    indices[o₁+inds] = indices[o₁] = oldinds[i₁]
                    values[o₁+inds] = values[o₁] = oldvals[i₁]
                else
                    indices[o₁+inds] = indices[o₁] = oldinds[i₁]
                    values[o₁] = oldvals[i₁] + oldvals[lens[1]+i₂]
                    values[o₁+inds] = oldvals[i₁] - oldvals[lens[1]+i₂]
                end
            end
        end)
        lens[1] = inds
        lens[2] = inds
        for j in length(lens):-1:2
            iszero(lens[j]) && deleteat!(lens, j)
        end
    end
    return IndvalsIterator(unsafe, indices, values, lens)
end

macro _infoset(name::Symbol, symbols::QuoteNode...)
    esc(quote
        const $(Symbol(:INFO_, name)) = ($(symbols...),)
        const $(Symbol(:VAL_, name)) = Union{$((:(Val{$x}) for x in symbols)...),}
    end)
end

@_infoset(PSD, :psd, :psd_complex, :rotated_quadratic, :quadratic, :nonnegative)
@_infoset(DD, :dd, :dd_complex,
    :dd_lnorm_real_diag,    :dd_lnorm_real_triu,    :dd_lnorm_real_tril,    :dd_lnorm_real,
    :dd_lnorm_complex_diag, :dd_lnorm_complex_triu, :dd_lnorm_complex_tril, :dd_lnorm_complex,
    :dd_nonneg_diag,        :dd_nonneg_triu,        :dd_nonneg_tril,        :dd_nonneg,
    :dd_quad_diag,          :dd_quad_triu,          :dd_quad_tril,          :dd_quad)
@_infoset(SDD,
    :sdd, :sdd_complex,
    :sdd_quad_real_diag,    :sdd_quad_real_triu,    :sdd_quad_real_tril,    :sdd_quad_real,
    :sdd_quad_complex_diag, :sdd_quad_complex_triu, :sdd_quad_complex_tril, :sdd_quad_complex)
@_infoset(COMPLEX,
    :psd_complex, :sdd_complex, :dd_complex,
    :dd_lnorm_complex_diag, :dd_lnorm_complex_triu, :dd_lnorm_complex_tril, :dd_lnorm_complex,
    :dd_quad_diag,          :dd_quad_triu,          :dd_quad_tril,          :dd_quad,
    :sdd_quad_complex_diag, :sdd_quad_complex_triu, :sdd_quad_complex_tril, :sdd_quad_complex)
@_infoset(DIAG,
    :dd_lnorm_real_diag, :dd_lnorm_complex_diag, :dd_nonneg_diag, :dd_quad_diag,
    :sdd_quad_real_diag, :sdd_quad_complex_diag)
@_infoset(TRIU,
    :dd_lnorm_real_triu, :dd_lnorm_complex_triu, :dd_nonneg_triu, :dd_quad_triu,
    :sdd_quad_real_triu, :sdd_quad_complex_triu)
@_infoset(TRIL,
    :dd_lnorm_real_tril, :dd_lnorm_complex_tril, :dd_nonneg_tril, :dd_quad_tril,
    :sdd_quad_real_tril, :sdd_quad_complex_tril)
@_infoset(MATRIX, :psd, :psd_complex, :sdd, :sdd_complex, :dd, :dd_complex)
@_infoset(NOMATRIX_REAL,
    :sdd_quad_real_diag,    :sdd_quad_real_triu,    :sdd_quad_real_tril,    :sdd_quad_real,
    :dd_lnorm_real_diag,    :dd_lnorm_real_triu,    :dd_lnorm_real_tril,    :dd_lnorm_real,
    :dd_nonneg_diag,        :dd_nonneg_triu,        :dd_nonneg_tril,        :dd_nonneg)
@_infoset(NOMATRIX_COMPLEX,
    :dd_lnorm_complex_diag, :dd_lnorm_complex_triu, :dd_lnorm_complex_tril, :dd_lnorm_complex,
    :dd_quad_diag,          :dd_quad_triu,          :dd_quad_tril,          :dd_quad,
    :sdd_quad_complex_diag, :sdd_quad_complex_triu, :sdd_quad_complex_tril, :sdd_quad_complex
)

# generic moment matrix constraint with
# - only real-valued monomials involved in the grouping, and only real-valued polynomials involved in the constraint (so if it
#   contains complex coefficients/monomials, imaginary parts cancel out)
# - or complex-valued monomials involved in the grouping, but the solver supports the complex-valued PSD cone explicitly
# - or DD/SDD representations in the real case and in the complex case if the quadratic cone is requested
function moment_add_matrix_helper!(state::AnySolver{T,V}, grouping::AbstractVector{M} where {M<:SimpleMonomial},
    constraint::AbstractMatrix{<:SimplePolynomial}, indextype::PSDIndextype{Tri},
    type::Union{Tuple{Val{true},Val},Tuple{Val{false},Val{true}}}, representation::RepresentationMethod,
    counters::Counters) where {T,V,Tri}
    # Note: We rely on a lot of checks having already been done, so that the combination of arguments is always correct. This
    # function should not be called from the outside.

    lg = length(grouping)
    block_size = LinearAlgebra.checksquare(constraint)
    matrix_indexing = indextype isa PSDIndextypeMatrixCartesian
    dim = lg * block_size
    complex = type isa Tuple{Val{false},Val{true}}
    if matrix_indexing
        tri = :U # we always create the data in U format, as the PSDMatrixCartesian then has to compute the row/col indices
                 # based on the linear index - and the formula is slightly simpler for U. However, if Tri === :L, we will fill
                 # the values with conjugates.
    else
        tri = Tri
    end
    maxlen = maximum(length, constraint, init=0)
    colcount = (tri === :F ? (complex ? 2dim^2 - dim : dim^2) : (complex ? dim^2 : trisize(dim)))
    indlen = 2maxlen * colcount
    if matrix_indexing
        rows = FastVec{T}(buffer=indlen)
        indices = similar(rows)
        values = similar(rows, complex ? Complex{V} : V)
    else
        lens = FastVec{SimplePolynomials.smallest_unsigned(2maxlen)}(buffer=colcount)
        indices = FastVec{T}(buffer=indlen)
        values = similar(indices, V)
    end

    # introduce a method barrier to fix the potentially unknown eltype of lens due to the dynamic maxlen
    return moment_add_matrix_helper!(state, grouping, constraint, indextype, Val(tri), Val(matrix_indexing), Val(complex), lg,
        block_size, dim, matrix_indexing ? (rows, indices, values) : (lens, indices, values), representation, counters)
end

function moment_add_matrix_helper!(state::AnySolver{T,V}, grouping::AbstractVector{M} where {M<:SimpleMonomial},
    constraint::AbstractMatrix{<:SimplePolynomial}, indextype::PSDIndextype{Tri}, ::Val{tri}, ::Val{matrix_indexing},
    ::Val{complex}, lg, block_size, dim, data, representation::RepresentationMethod,
    counters::Counters) where {T,V,Tri,tri,matrix_indexing,complex}
    if matrix_indexing
        rows, indices, values = data
        row = zero(T)
        scaleoffdiags = false
    else
        lens, indices, values = data
        if representation isa RepresentationPSD
            # Off-diagonals are multiplied by √2 in order to put variables into the vectorized PSD cone. Even if state isa
            # SOSWrapper (then, the variables directly correspond to a vectorized PSD cone), the actual values in the PSD
            # matrix are still multiplied by 1/√2, so we must indeed always multiply the coefficients by √2 to undo this.
            # We also have to account for the unwanted factor of 2 in the rotated quadratic cone.
            # In the case of a (normal) quadratic cone, we canonically take the rotated cone and transform it by multiplying
            # the left-hand side by 1/√2, giving (x₁/√2)² ≥ (x₂/√2)² + ∑ᵢ (√2 xᵢ)² ⇔ x₁² ≥ x₂² + ∑ᵢ (2 xᵢ)².
            scaleoffdiags = tri !== :F
            if dim == 2
                rquad = supports_rotated_quadratic(state)
                quad = supports_quadratic(state)
                if scaleoffdiags
                    if !rquad && quad
                        scaling = V(2)
                    else
                        scaling = sqrt(V(2))
                    end
                else
                    scaling = sqrt(V(2))
                end
            elseif scaleoffdiags
                scaling = indextype.scaling
                if scaling isa Bool
                    # If no scaling is desired, set it to true, which allows us to completely eliminate the multiplication -
                    # because the value false does not make any sense, so this path is statically determined.
                    @assert(scaling === true)
                    scaleoffdiags = false
                end
            end
        else
            rquad = quad = false
            if (representation isa RepresentationDD && ((complex && supports_dd_complex(state)) ||
                                                        (!complex && supports_dd(state)))) ||
               (representation isa RepresentationSDD && ((complex && supports_sdd_complex(state)) ||
                                                         (!complex && supports_sdd(state))))
                scaleoffdiags = tri !== :F
                if scaleoffdiags
                    scaling = indextype.scaling
                    if scaling isa Bool
                        @assert(scaling === true)
                        scaleoffdiags = false
                    end
                end
            else
                scaleoffdiags = true
                scaling = V(2)
            end
        end
    end
    @inbounds for (exp2, g₂) in enumerate(grouping)
        isreal_g₂ = !complex || isreal(g₂)
        for block_j in 1:block_size
            exp1_range = (tri === :F ? (1:lg) : (tri === :U ? (1:exp2) : (exp2:lg)))
            for (exp1, g₁) in zip(exp1_range, @view(grouping[exp1_range]))
                isreal_g₁ = !complex || isreal(g₁)
                for block_i in (tri === :F ? (1:block_size) : (tri === :U ? (1:(exp1 == exp2 ? block_j : block_size)) :
                                                                            ((exp1 == exp2 ? block_j : 1):block_size)))
                    # - !complex: all polynomials are real-valued in total. Only take into account canonical monomials, but
                    #   double the prefactor if their conjugated must also occur.
                    # - block_i == block_j && exp1 == exp2: we are on the total diagonal. Even if not all polynomials were
                    #   real, these are for sure.
                    # - block_i == block_j && isreal(g₁) && isreal(g₂): we are on the inner diagonal, but since the outer basis
                    #   elements are real-valued, canonical and non-canonicals will also occur in pairs.
                    total_real = !complex || (block_i == block_j && (exp1 == exp2 || (isreal_g₁ && isreal_g₂)))
                    ondiag = exp1 == exp2 && block_i == block_j
                    @twice secondtime (!matrix_indexing && !total_real) begin
                        curlen = 0
                        for term_constr in constraint[block_i, block_j]
                            mon_constr = monomial(term_constr)
                            coeff_constr = coefficient(term_constr)
                            if scaleoffdiags && !ondiag
                                coeff_constr *= scaling
                            end
                            recoeff = real(coeff_constr)
                            imcoeff = imag(coeff_constr)
                            repart, impart, canonical = getreim(state, g₁, mon_constr, SimpleConjMonomial(g₂))
                            if total_real
                                # Even if `complex` - so we are on the diagonal and the values should be Complex{V} - we can
                                # just do the push!, the conversion is implicitly done; as the imaginary part should be zero,
                                # we don't need to discriminate against matrix_indexing.
                                if canonical # ≡ iscanonical(mon_constr)
                                    if !iszero(recoeff)
                                        matrix_indexing && unsafe_push!(rows, row)
                                        unsafe_push!(indices, repart)
                                        unsafe_push!(values, repart == impart ? recoeff : V(2) * recoeff)
                                        curlen += 1
                                    end
                                    if repart == impart
                                        @assert(iszero(imcoeff)) # else the polynomial on the diagonal would not be real-valued
                                        @assert(isreal(mon_constr))
                                    elseif !iszero(imcoeff)
                                        matrix_indexing && unsafe_push!(rows, row)
                                        unsafe_push!(indices, impart)
                                        unsafe_push!(values, V(-2) * imcoeff)
                                        curlen += 1
                                    end
                                end
                            else
                                # We are not on the diagonal, so we must work with every entry. Updating might become necessary
                                # as the conjugate can pop up later.
                                if matrix_indexing
                                    # Interpretation: Complex(x, y) means that the real part of the entry is given by x and
                                    #                 the imaginary part by y.
                                    let coeff=Complex(recoeff, Tri === :L ? -imcoeff : imcoeff)
                                        if !iszero(coeff)
                                            @pushorupdate!(repart, row, coeff)
                                            repart == impart ||
                                                @pushorupdate!(impart, row,
                                                    canonical ? Complex(-imcoeff, Tri === :L ? -recoeff : recoeff) :
                                                                Complex(imcoeff, Tri === :L ? recoeff : -recoeff))
                                        end
                                    end
                                elseif !secondtime
                                    iszero(recoeff) || @pushorupdate!(repart, row, recoeff)
                                    repart == impart || iszero(imcoeff) ||
                                        @pushorupdate!(impart, row, canonical ? -imcoeff : imcoeff)
                                else
                                    iszero(imcoeff) || @pushorupdate!(repart, row, imcoeff)
                                    repart == impart || iszero(recoeff) ||
                                        @pushorupdate!(impart, row, canonical ? recoeff : -recoeff)
                                end
                            end
                        end
                        if matrix_indexing
                            row += one(row)
                        else
                            unsafe_push!(lens, Core.Intrinsics.trunc_int(eltype(lens), curlen)) # no check for overflows
                        end
                    end
                    if !matrix_indexing && complex && total_real && !ondiag
                        # We skipped the addition of an imaginary part because it is zero. But we have to tell this to the
                        # solver.
                        unsafe_push!(lens, zero(eltype(lens)))
                    end
                end
            end
        end
    end
    if matrix_indexing # implies dim ≥ 3
        mc = PSDMatrixCartesian{_get_offset(indextype)}(dim, Tri, finish!(rows), finish!(indices), finish!(values))
        if representation isa RepresentationPSD
            if complex
                add_constr_psd_complex!(state, dim, mc)
                return :psd_complex, addtocounter!(state, counters, Val(:psd_complex), dim^2)
            else
                add_constr_psd!(state, dim, mc)
                return :psd, addtocounter!(state, counters, Val(:psd), trisize(dim))
            end
        elseif representation isa RepresentationDD
            if complex
                add_constr_dddual_complex!(state, dim, mc)
                return :dd_complex, addtocounter!(state, counters, Val(:dd_complex), dim^2)
            else
                add_constr_dddual!(state, dim, mc)
                return :dd, addtocounter!(state, counters, Val(:dd), trisize(dim))
            end
        else
            if complex
                add_constr_sdddual_complex!(state, dim, mc)
                return :sdd_complex, addtocounter!(state, counters, Val(:sdd_complex), dim^2)
            else
                add_constr_sdddual!(state, dim, mc)
                return :sdd, addtocounter!(state, counters, Val(:sdd), trisize(dim))
            end
        end
    else
        let indices=finish!(indices), values=finish!(values)
            @inbounds if dim == 1
                @assert(length(lens) == 1)
                add_constr_nonnegative!(state, Indvals(indices, values))
                return :nonnegative, addtocounter!(state, counters, Val(:nonnegative), dim)
            elseif dim == 2 && (rquad || quad)
                @assert(length(lens) == (complex ? 4 : 3))
                # Note: indices/values represent the vectorized PSD cone (lower triangle) with off-diagonals already
                # pre-multiplied by sqrt(2). To make this exactly equivalent to a rotated quadratic cone, we must bring the
                # last element to the second position.
                @views @unroll for x in (indices, values)
                    reverse!(x[lens[1]+1:end-lens[end]])
                    reverse!(x[end-lens[end]+1:end])
                    reverse!(x[lens[1]+1:end])
                end
                if complex
                    lens[2], lens[3], lens[4] = lens[4], lens[2], lens[3]
                else
                    lens[2], lens[3] = lens[3], lens[2]
                end
                if rquad
                    add_constr_rotated_quadratic!(state, to_soc!(indices, values, lens, true))
                    return :rotated_quadratic, addtocounter!(state, counters, Val(:rotated_quadratic), complex ? 4 : 3)
                else
                    add_constr_quadratic!(state, to_soc!(indices, values, lens, false))
                    return :quadratic, addtocounter!(state, counters, Val(:quadratic), complex ? 4 : 3)
                end
            else # implies dim ≥ 3
                ii = IndvalsIterator(indices, values, lens)
                cl = length(lens)
                if representation isa RepresentationPSD
                    if complex
                        add_constr_psd_complex!(state, dim, ii)
                        return :psd_complex, addtocounter!(state, counters, Val(:psd_complex), cl)
                    else
                        add_constr_psd!(state, dim, ii)
                        return :psd, addtocounter!(state, counters, Val(:psd), cl)
                    end
                elseif representation isa RepresentationDD
                    if complex
                        if supports_dd_complex(state)
                            add_constr_dddual_complex!(state, dim, ii, representation.u)
                            return :dd_complex, addtocounter!(state, counters, Val(:dd_complex), cl)
                        else
                            return moment_add_dddual_transform!(state, dim, ii, representation.u, true, counters)
                        end
                    else
                        if supports_dd(state)
                            add_constr_dddual!(state, dim, ii, representation.u)
                            return :dd, addtocounter!(state, counters, Val(:dd), cl)
                        else
                            return moment_add_dddual_transform!(state, dim, ii, representation.u, false, counters)
                        end
                    end
                else
                    if complex
                        if supports_sdd_complex(state)
                            add_constr_sdddual_complex!(state, dim, ii, representation.u)
                            return :sdd_complex, addtocounter!(state, counters, Val(:sdd_complex), cl)
                        else
                            return moment_add_sdddual_transform!(state, dim, ii, representation.u, true, counters)
                        end
                    else
                        if supports_sdd(state)
                            add_constr_sdddual!(state, dim, ii, representation.u)
                            return :sdd, addtocounter!(state, counters, Val(:sdd), cl)
                        else
                            return moment_add_sdddual_transform!(state, dim, ii, representation.u, false, counters)
                        end
                    end
                end
            end
        end
    end
end

# generic moment matrix constraint with complex-valued monomials involved in the grouping, but the solver does not support the
# complex PSD cone explicitly
function moment_add_matrix_helper!(state::AnySolver{T,V}, grouping::AbstractVector{M} where M<:SimpleMonomial,
    constraint::AbstractMatrix{<:SimplePolynomial}, indextype::PSDIndextype{Tri},
    ::Tuple{Val{false},Val{false}}, representation::Union{<:RepresentationDD,<:RepresentationSDD,RepresentationPSD},
    counters::Counters) where {T,V,Tri}
    # Note: We rely on a lot of checks having already been done, so that the combination of arguments is always correct. This
    # function should not be called from the outside.

    lg = length(grouping)
    block_size = LinearAlgebra.checksquare(constraint)
    matrix_indexing = indextype isa PSDIndextypeMatrixCartesian
    dim = lg * block_size
    dim2 = 2dim
    if matrix_indexing
        tri = :U # we always create the data in U format, as the PSDMatrixCartesian then has to compute the row/col indices
                 # based on the linear index - and the formula is slightly simpler for U.
        scaleoffdiags = false
    else
        tri = Tri
        if representation isa RepresentationPSD ||
            (representation isa RepresentationDD ? supports_dd(state) : supports_sdd(state))
            scaleoffdiags = tri !== :F
            if scaleoffdiags
                scaling = indextype.scaling
                if scaling isa Bool
                    @assert(scaling === true)
                    scaleoffdiags = false
                end
            end
        else
            scaleoffdiags = true
            scaling = V(2)
        end
    end
    maxlen = maximum(length, constraint, init=0)
    colcount = tri === :F ? 4dim^2 : trisize(2dim)
    rows = FastVec{T}(buffer=2maxlen * colcount)
    indices = similar(rows)
    values = similar(rows, V)
    # Here we implement the default way: Z ⪰ 0 ⟺ [Re Z  -Im Z; Im Z  Re Z] ⪰ 0.
    # This is relatively simple, as "rows" encodes the position of the element in the matrix (in linear index form of the upper
    # triangle columnwise).
    dimT = T(dim)
    items = tri === :F ? dimT^2 : trisize(dimT)
    i1 = zero(T)
    if tri === :U
        i2order = items
        i3 = items + dimT
    elseif tri === :L
        i2order = dimT
        i3 = (T(3) * dimT^2 + dimT) >> 1
    else
        i2order = dimT
        i2other = T(2) * items
        i3 = i2other + dimT
    end
    col = zero(T)
    @inbounds for (exp2, g₂) in enumerate(grouping)
        isreal_g₂ = isreal(g₂)
        for block_j in 1:block_size
            row = tri === :L ? col : zero(T)
            exp1_range = (tri === :F ? (1:lg) : (tri === :U ? (1:exp2) : (exp2:lg)))
            if tri === :U
                i2other = items + col + (row * (row + T(2) * dimT + one(T))) >> 1
            elseif tri === :L
                i2other = i2order
            end
            for (exp1, g₁) in zip(exp1_range, @view(grouping[exp1_range]))
                isreal_g₁ = isreal(g₁)
                for block_i in (tri === :F ? (1:block_size) : (tri === :U ? (1:(exp1 == exp2 ? block_j : block_size)) :
                                                                            ((exp1 == exp2 ? block_j : 1):block_size)))
                    total_real = block_i == block_j && (exp1 == exp2 || (isreal_g₁ && isreal_g₂))
                    ondiag = exp1 == exp2 && block_i == block_j
                    curlen = 0
                    for term_constr in constraint[block_i, block_j]
                        mon_constr = monomial(term_constr)
                        coeff_constr = coefficient(term_constr)
                        if scaleoffdiags && !ondiag
                            coeff_constr *= scaling
                        end
                        recoeff = real(coeff_constr)
                        imcoeff = imag(coeff_constr)
                        repart, impart, canonical = getreim(state, g₁, mon_constr, SimpleConjMonomial(g₂))
                        if total_real
                            if canonical
                                if !iszero(recoeff)
                                    unsafe_push!(rows, i1, i3)
                                    unsafe_push!(indices, repart, repart)
                                    if repart != impart
                                        recoeff *= V(2)
                                    end
                                    unsafe_push!(values, recoeff, recoeff)
                                    curlen += 2
                                end
                                if repart == impart
                                    @assert(iszero(imcoeff)) # else the polynomial on the diagonal would not be real-valued
                                    @assert(isreal(mon_constr))
                                elseif !iszero(imcoeff)
                                    unsafe_push!(rows, i1, i3)
                                    unsafe_push!(indices, impart, impart)
                                    imcoeff *= V(-2)
                                    unsafe_push!(values, imcoeff, imcoeff)
                                    curlen += 2
                                end
                            end
                        else
                            @assert(i2order != i2other)
                            if !iszero(recoeff)
                                @pushorupdate!(repart, i1, recoeff, i3, recoeff)
                                if repart != impart
                                    care = canonical ? recoeff : -recoeff
                                    # We don't really need the L/U distinction, as it does not matter where to place the
                                    # minus. However, to make testing simpler, we keep our promise.
                                    if tri === :L || tri === :F
                                        @pushorupdate!(impart, i2order, care, i2other, -care)
                                    else
                                        @pushorupdate!(impart, i2other, care, i2order, -care)
                                    end
                                end
                            end
                            if !iszero(imcoeff)
                                if tri === :L || tri === :F
                                    @pushorupdate!(repart, i2order, imcoeff, i2other, -imcoeff)
                                else
                                    @pushorupdate!(repart, i2other, imcoeff, i2order, -imcoeff)
                                end
                                if repart != impart
                                    care = canonical ? -imcoeff : imcoeff
                                    @pushorupdate!(impart, i1, care, i3, care)
                                end
                            end
                        end
                    end
                    row += one(T)
                    i1 += one(T)
                    if tri === :U
                        i2order += one(T)
                        i2other += dimT + row
                    elseif tri === :L
                        i2order += one(T)
                        i2other += T(2) * dimT - row
                    else
                        i2order += one(T)
                        i2other += one(T)
                    end
                    i3 += one(T)
                end
            end
            if tri === :U
                i2order += dimT
                i3 += dimT
            elseif tri === :L
                i1 += dimT
                i2order += dimT
            else
                i1 += dimT
                i2order += dimT
                i2other += dimT
                i3 += dimT
            end
            col += one(T)
        end
    end
    if matrix_indexing
        add_constr_psd!(state, dim2, PSDMatrixCartesian{_get_offset(indextype)}(dim2, Tri, finish!(rows),
            finish!(indices), finish!(values)))
        return :psd, addtocounter!(state, counters, Val(:psd), dim2^2)
    else
        # We constructed everything in matrix form, as it was easier to assign arbitrary positions in this way; but now we need
        # to convert it to the iterative form. We cannot re-use rows for this, unfortunately. It may happen that some rows
        # (i.e., entries for the PSD matrix) are zero and don't arise in rows, although they will have to be present in lens.
        # So to avoid potentially overwriting data, we need a whole new vector.
        @inbounds let rows=finish!(rows), indices=finish!(indices), values=finish!(values)
            lens = Vector{SimplePolynomials.smallest_unsigned(2maxlen)}(undef, (Tri === :F ? 4dim^2 : trisize(2dim)))
            sort_along!(rows, indices, values)
            row = zero(T)
            i = 1
            j = 1
            remaining = length(rows)
            while !iszero(remaining)
                currow = rows[j]
                while currow > row
                    lens[i] = zero(eltype(lens))
                    i += 1
                    row += one(T)
                end
                curlen = one(eltype(lens))
                j += 1
                remaining -= 1
                while !iszero(remaining) && rows[j] == currow
                    curlen += one(eltype(lens))
                    j += 1
                    remaining -= 1
                end
                lens[i] = curlen
                i += 1
                row += one(T)
            end
            fill!(@view(lens[i:end]), zero(eltype(lens)))

            ii = IndvalsIterator(indices, values, lens)
            cl = length(lens)
            if representation isa RepresentationPSD
                add_constr_psd!(state, dim2, ii)
                return :psd, addtocounter!(state, counters, Val(:psd), cl)
            elseif representation isa RepresentationDD
                if supports_dd(state)
                    add_constr_dddual!(state, dim2, ii, representation.u)
                    return :dd, addtocounter!(state, counters, Val(:dd), cl)
                else
                    return moment_add_dddual_transform!(state, dim2, ii, representation.u, false, counters)
                end
            else
                if supports_sdd(state)
                    add_constr_sdddual!(state, dim2, ii, representation.u)
                    return :sdd, addtocounter!(state, counters, Val(:sdd), cl)
                else
                    return moment_add_sdddual_transform!(state, dim2, ii, representation.u, false, counters)
                end
            end
        end
    end
end

#region Rewrite the DD cone in terms of other cones.

# Diagonal-dominant representation: this is a relaxation for the SOS formulation, where we replace M ∈ PSD by
# M ∈ {U† D U, D ∈ DD}. Since this is more restrictive than PSD, the SOS maximization will only decrease, so we still have a
# valid lower bound.
# Vectorized version: vec(M) = vec(U† mat(d) U). In component form, this is
# mᵢ = ∑_(diagonal j) Ū[row(j), row(i)] U[col(j), col(i)] dⱼ +
#      ∑_(offdiag j) (Ū[col(j), row(i)] U[row(j), col(i)] + Ū[row(j), row(i)] U[col(j), col(i)]) dⱼ ⇔ m = Ũ d.
# Note that if U is diagonal, mᵢ = Ū[row(i), row(i)] U[col(i), col(i)] dᵢ.
# So define d ∈ vec(DD), m free, then demand 𝟙*m + (-Ũ)*d = 0. While we could eliminate m fully, it is actually advantageous to
# keep them: in this way, we have direct access to the rotated DD cone in the solution returned by the solver. Additionally,
# using m instead of -Ũ*d in the constraints can give a higher sparsity. Of course, if U is the identity, we would be better
# off eliminating m.
# Then, we need to translate DD into a cone that is supported; let's assume for simplicity that the ℓ₁ cone is available.
# DD = ℓ₁ × ... × ℓ₁ plus equality constraints that enforce symmetry.
# However, here we construct the moment representation; so we now need the dual formulation of diagonal dominance. Due to the
# equality constraints, this is more complicated:
# For side dimension n, there are n ℓ_∞ cones (we just take the columns - also taking into account the rows would be even more
# restrictive). For every data entry, we need a new slack variable which we fix to be equal to this data entry (a variable for
# each the real and the imaginary part). We need additional slacks for every entry in the strict upper triangle (twice in the
# complex-valued case).
# Without the U, all this would look as follows (for a 3x3 real matrix):
# data₁ = slack₁, data₂ = slack₂, data₃ = slack₃, data₄ = slack₄, data₅ = slack₅, data₆ = slack₆
# {slack₁, slack₂ - slack₇, slack₃ - slack₈} ∈ ℓ_∞
# {slack₄, slack₇,          slack₅ - slack₉} ∈ ℓ_∞
# {slack₆, slack₈,          slack₉}          ∈ ℓ_∞
# effectively corresponding to the DDDual matrix
# data₁            slack₇           slack₈
# data₂ - slack₇   data₄            slack₉
# data₃ - slack₈   data₅ - slack₉   data₆

# The assignment of slack variables to the data will always be the same and unrelated to U. Therefore, in the following, for
# better clarity, we will omit the equality constraints and simply write slack(row, col), where it is understood that for
# row ≥ col (lower triangle), this corresponds to the equality-constrained slack variable that is (twice) the data, and for
# row < col (strict upper triangle), this corresponds to a free slack variable required due to symmetry.

# For a general U and complex-valued data, we then have for the column j:
# {∑_{col = 1}^dim ∑_{row = col}^dim (Re(U[j, col] Ū[j, row]) slackᵣ(row, col) + Im(U[j, col] Ū[j, row]) slackᵢ(row, col)),
#  slackᵣ(i, j), slackᵢ(i, j) for i ∈ 1, ..., j -1,
#  ∑_{col = 1}^dim ∑_{row = col}^dim ((Re(U[i, row] Ū[j, col] + U[i, col] Ū[j, row]) slackᵣ(row, col) -
#                                     (Im(U[i, row] Ū[j, col] - U[i, col] Ū[j, row]) slackᵢ(row, col)) - slackᵣ(j, i),
#  ∑_{col = 1}^dim ∑_{row = col}^dim ((Im(U[i, row] Ū[j, col] + U[i, col] Ū[j, row]) slackᵣ(row, col) +
#                                     (Re(U[i, row] Ū[j, col] - U[i, col] Ū[j, row]) slackᵢ(row, col)) + slackᵢ(j, i)
#  for i in j +1, ..., dim
# } ∈ ℓ_∞

# Let's specialize the formula. If U is diagonal:
# {|U[j, j]|² slackᵣ(j, j),
#  slackᵣ(i, j), slackᵢ(i, j) for i ∈ 1, ..., j -1,
#  (Re(U[i, i] Ū[j, j]) slackᵣ(i, j) - Im(U[i, i] Ū[j, j]) slackᵢ(i, j)) - slackᵣ(j, i),
#  (Im(U[i, i] Ū[j, j]) slackᵣ(i, j) + Re(U[i, i] Ū[j, j]) slackᵢ(i, j)) + slackᵢ(j, i)
#  for i in j +1, ..., dim
# } ∈ ℓ_∞

# If everything is instead real-valued:
# {∑_{col = 1}^dim ∑_{row = col}^dim U[j, col] U[j, row] slack(row, col),
#  slack(i, j) for i ∈ 1, ..., j -1,
#  ∑_{col = 1}^dim ∑_{row = col}^dim (U[i, row] U[j, col] + U[i, col] U[j, row]) slack(row, col) - slack(j, i)
#  for i in j +1, ..., dim
# } ∈ ℓ_∞

# If everything is real and U is diagonal:
# {U[j, j]² slack(j, j),
#  slack(i, j) for i ∈ 1, ..., j -1,
#  U[i, i] U[j, j] slack(i, j) - slack(j, i) for i in j +1, ..., dim
# } ∈ ℓ_∞

# If the ℓ_∞ cone is not available, we can instead use the nonnegative cone if we are not complex. The most efficient way to
# do this is in fact to directly rewrite the ℓ_∞ cone formulation in terms of two nonnegative constraints for all but the
# diagonal entry. In turn, if we are complex, the most efficient way is to split the ℓ_∞ cone formulation in terms of quadratic
# cones with the diagonal as the first entry and the corresponding off-diagonals in real and imaginary part as the second and
# third entry.

# If we don't have a full matrix u, but also not diagonal, we'll specialize two other cases: the triangular ones (and the upper
# triangular is the important one for our default Cholesky update strategy). This is not as sophisticated as in the diagonal
# case, we don't adjust the creation of slack variables etc., but we suppress adding unnecessary zeros.
# Therefore, the loops will always be col in 1:dim, row in col:dim, and these methods should return whether the result is known
# to be zero, already knowing these bounds.
# ∑_{col = 1}^dim ∑_{row = col}^dim U[j, col] Ū[j, row]
dddual_transform_inrange(::Any, j, col) = true
dddual_transform_inrange(::Val{:diag}, ::Any, j, row) = true
# ↪ TRIU: j ≤ col ≤ row ≤ dim
dddual_transform_inrange(::UpperOrUnitUpperTriangular, j, col) = col ≥ j
dddual_transform_inrange(::Val{:diag}, ::UpperOrUnitUpperTriangular, j, row) = true
# ↪ TRIL: 1 ≤ col ≤ j, col ≤ row ≤ j
dddual_transform_inrange(::LowerOrUnitLowerTriangular, j, col) = col ≤ j
dddual_transform_inrange(::Val{:diag}, ::LowerOrUnitLowerTriangular, j, row) = row ≤ j

# ∑_{col = 1}^dim ∑_{row = col}^dim (U[i, row] Ū[j, col] + U[i, col] Ū[j, row])
dddual_transform_inrange(::Val{:offdiag}, ::Any, i, row) = true
# ↪ TRIU: j ≤ col ≤ dim, max(col, i) ≤ row ≤ dim
dddual_transform_inrange(::Val{:offdiag}, ::UpperOrUnitUpperTriangular, i, row) = row ≥ i
# ↪ TRIL: 1 ≤ col ≤ j, col ≤ row ≤ i
dddual_transform_inrange(::Val{:offdiag}, ::LowerOrUnitLowerTriangular, i, row) = row ≤ i

# Note: We always add the items and don't check whether the value is zero, potentially reducing the sparsity pattern (only for
# diagonal u, we respect it). This is so that the solver data is always at the same position and re-optimizations with a
# different u will not change the structure of the solver problem (of course, we'd have to write support for re-using the
# problem in the first place...).

function dddual_transform_equalities!(state::AnySolver{T,V}, ::Val{complex}, dim::Int, data::IndvalsIterator{T,V}, slacks,
    counters::Counters) where {T,V,complex}
    # Add the equality constraints first; this is independent of U or the method to model the cone
    @inbounds begin
        negate = negate_fix(state)
        if negate
            rmul!(nonzeros(data), -one(V))
        end
        eqstate = @inline add_constr_fix_prepare!(state, complex ? dim^2 : trisize(dim))
        # When iterating through data, we need to add our slack element; it does not matter whether at the beginning or the
        # end. So under the assertion that data contains more than just a single Indvals (which is always true, else this
        # function would not be called), we append the slack at the end, temporarily overwriting what is there for the first
        # item; and we prepend it for all others (no need to restore the value afterwards), avoiding copying around.
        firstindvals, restindvals = Iterators.peel(data)
        indices, values = firstindvals.indices, firstindvals.values
        # while we could simply do @inbounds indices[end+1], this won't work in the interpreter; and we need the index multiple
        # times anyway, so let's just compute it under the assumption that everything is as it should be.
        @assert(indices isa Base.FastContiguousSubArray)
        pos = indices.offset1 + length(indices) +1
        storedind = indices.parent[pos]
        storedval = values.parent[pos]
        indices.parent[pos] = slacks[1]
        values.parent[pos] = negate ? one(V) : -one(V)
        eqstate = @inline add_constr_fix!(state, eqstate, Indvals(view(indices.parent, indices.offset1+1:pos),
                                                                  view(values.parent, indices.offset1+1:pos)), zero(V))
        indices.parent[pos] = storedind
        values.parent[pos] = storedval
        offdiags = dim -1
        col = 1
        slack = 2
        let indvalsiter=iterate(restindvals)
            while !isnothing(indvalsiter)
                indvals = indvalsiter[1]
                indices, values = indvals.indices, indvals.values
                pos = indices.offset1
                indices.parent[pos] = slacks[slack]
                values.parent[pos] = negate ? one(V) : -one(V)
                if iszero(offdiags)
                    col += 1
                    offdiags = dim - col
                else
                    if complex
                        slack += 1
                        eqstate = @inline add_constr_fix!(state, eqstate,
                                                          Indvals(view(indices.parent, pos:indices.offset1+length(indices)),
                                                                  view(values.parent, pos:indices.offset1+length(indices))),
                                                          zero(V))
                        indvalsiter = iterate(restindvals, indvalsiter[2])::Tuple
                        indvals = indvalsiter[1]
                        indices, values = indvals.indices, indvals.values
                        pos = indices.offset1
                        indices.parent[pos] = slacks[slack]
                        values.parent[pos] = negate ? one(V) : -one(V)
                    end
                    offdiags -= 1
                end
                slack += 1
                eqstate = @inline add_constr_fix!(state, eqstate,
                                                  Indvals(view(indices.parent, pos:indices.offset1+length(indices)),
                                                          view(values.parent, pos:indices.offset1+length(indices))), zero(V))
                indvalsiter = iterate(restindvals, indvalsiter[2])
            end
        end
        @inline add_constr_fix_finalize!(state, eqstate)
        return addtocounter!(state, counters, Val(:fix), complex ? dim^2 : trisize(dim))
    end
end

const DiagonalTransform = Union{<:UniformScaling,<:Diagonal}

# We have the lnorm cone available, the transformation is diagonal
function dddual_transform_cone!(state::AnySolver{T,V}, ::Val{true}, ::Val{complex}, dim::Int, data::IndvalsIterator{T,V},
    u::DiagonalTransform, counters::Counters) where {T,V,complex}
    dsq = dim^2
    indices = FastVec{T}(buffer=if complex
            6dim -5 # 1 + 6(dim -1)
        else
            2dim -1 # 1 + 2(dim -1)
        end # either the first col or the longest data assignment are the largest
    )
    lens = FastVec{Int}(buffer=complex ? 2dim -1 : dim)
    slacks = add_var_slack!(state, complex ? dsq + 2trisize(dim -1) : dsq)
    values = similar(indices, V)

    if prepend_fix(state)
        valat = dddual_transform_equalities!(state, Val(complex), dim, data, slacks, counters)
    end

    upperslack = (complex ? dsq : trisize(dim)) +1
    lowerslack = 1
    @inbounds for j in 1:dim
        #region Diagonal (naturally is the first item in the L order, and must be the first in the ℓ_∞ cone)
        unsafe_push!(indices, slacks[lowerslack])
        unsafe_push!(values, u isa Diagonal ? abs2(u[j, j]) : one(V))
        unsafe_push!(lens, 1)
        lowerslack += 1
        #endregion
        #region Above diagonal (slacks)
        δ = complex ? 2j -2 : j -1
        unsafe_append!(indices, @view(slacks[upperslack:upperslack+δ-1]))
        unsafe_append!(values, Iterators.repeated(one(V), δ))
        unsafe_append!(lens, Iterators.repeated(1, δ))
        upperslack += δ
        #endregion
        #region Below diagonal
        sbelow = upperslack + (complex ? 2j - 2 : j -1)
        for i in j+1:dim
            uval = u isa Diagonal ? u[i, i] * conj(u[j, j]) : one(V)
            if complex
                unsafe_push!(indices, slacks[lowerslack], slacks[lowerslack+1], slacks[sbelow], slacks[lowerslack],
                    slacks[lowerslack+1], slacks[sbelow+1])
                unsafe_push!(values, real(uval), -imag(uval), -one(V), imag(uval), real(uval), one(V))
                unsafe_push!(lens, 3, 3)
                sbelow += 2i -2
                lowerslack += 2
            else
                unsafe_push!(indices, slacks[lowerslack], slacks[sbelow])
                unsafe_push!(values, real(uval), -one(V)) # should be real anyway
                unsafe_push!(lens, 2)
                sbelow += i -1
                lowerslack += 1
            end
        end
        #endregion
        #region Add the whole column
        (complex ? add_constr_linf_complex! : add_constr_linf!)(state, IndvalsIterator(unsafe, indices, values, lens))
        #endregion
        empty!(indices)
        empty!(values)
        empty!(lens)
    end

    if !prepend_fix(state)
        valat = dddual_transform_equalities!(state, Val(complex), dim, data, slacks, counters)
    end

    return complex ? (:dd_lnorm_complex_diag,
                      (valat, addtocounter!(state, counters, Val(:lnorm_complex), dim, complex ? 2dim -1 : dim))) :
                     (:dd_lnorm_real_diag,
                      (valat, addtocounter!(state, counters, Val(:lnorm_complex), dim, complex ? 2dim -1 : dim)))
end

# We have the lnorm cone available, the transformation is not diagonal
function dddual_transform_cone!(state::AnySolver{T,V}, ::Val{true}, ::Val{complex}, dim::Int, data::IndvalsIterator{T,V}, u,
    counters::Counters) where {T,V,complex}
    ts = trisize(dim)
    dsq = dim^2
    indices = FastVec{T}(buffer=if complex
            2 * ((ts +1) * (2dim -1) - dim) # 2trisize(dim) + (2trisize(dim) +1) * 2(dim -1)
        else
            (ts +1) * dim -1 # trisize(dim) + (trisize(dim) +1) * (dim -1)
        end # either the first col or the longest data assignment are the largest
    )
    lens = FastVec{Int}(buffer=complex ? 2dim -1 : dim)
    slacks = add_var_slack!(state, complex ? dsq + 2trisize(dim -1) : dsq)
    values = similar(indices, V)

    if prepend_fix(state)
        valat = dddual_transform_equalities!(state, Val(complex), dim, data, slacks, counters)
    end

    upperslack = (complex ? dsq : ts) +1
    @inbounds for j in 1:dim
        #region Diagonal (naturally is the first item in the L order, and must be the first in the ℓ_∞ cone)
        lowerslack = 1
        for col in 1:dim
            if !dddual_transform_inrange(u, j, col)
                lowerslack += 1 + (complex ? 2 : 1) * (dim - col)
                continue
            end
            if dddual_transform_inrange(Val(:diag), u, j, col)
                unsafe_push!(indices, slacks[lowerslack])
                unsafe_push!(values, abs2(u[j, col]))
            end
            lowerslack += 1
            for row in col+1:dim
                if dddual_transform_inrange(Val(:diag), u, j, row)
                    uval = u[j, col] * conj(u[j, row])
                    @twice imdata complex begin
                        unsafe_push!(indices, slacks[lowerslack+imdata])
                        unsafe_push!(values, imdata ? imag(uval) : real(uval))
                    end
                end
                lowerslack += complex ? 2 : 1
            end
        end
        unsafe_push!(lens, length(indices))
        #endregion
        #region Above diagonal (slacks)
        δ = complex ? 2j -2 : j -1
        unsafe_append!(indices, @view(slacks[upperslack:upperslack+δ-1]))
        unsafe_append!(values, Iterators.repeated(one(V), δ))
        unsafe_append!(lens, Iterators.repeated(1, δ))
        upperslack += δ
        #endregion
        #region Below diagonal
        sbelow = upperslack + (complex ? 2j - 2 : j -1)
        for i in j+1:dim
            @twice impart complex begin
                startind = length(indices)
                slack = 1
                for col in 1:dim
                    if !dddual_transform_inrange(u, j, col)
                        slack += 1 + (complex ? 2 : 1) * (dim - col)
                        continue
                    end
                    if dddual_transform_inrange(Val(:offdiag), u, i, col)
                        uval = 2u[i, col] * conj(u[j, col])
                        unsafe_push!(indices, slacks[slack])
                        unsafe_push!(values, impart ? imag(uval) : real(uval))
                    end
                    slack += 1
                    for row in col+1:dim
                        if dddual_transform_inrange(Val(:offdiag), u, i, row)
                            uval = u[i, row] * conj(u[j, col])
                            @twice imdata complex let uval=uval
                                if imdata
                                    uval -= u[i, col] * conj(u[j, row])
                                    thisuval = impart ? real(uval) : -imag(uval)
                                else
                                    uval += u[i, col] * conj(u[j, row])
                                    thisuval = impart ? imag(uval) : real(uval)
                                end
                                unsafe_push!(indices, slacks[slack+imdata])
                                unsafe_push!(values, thisuval)
                            end
                        end
                        slack += complex ? 2 : 1
                    end
                end
                unsafe_push!(indices, slacks[sbelow])
                unsafe_push!(values, impart ? one(V) : -one(V))
                unsafe_push!(lens, length(indices) - startind)
                if complex
                    sbelow += impart ? 2i -3 : 1
                else
                    sbelow += i -1
                end
            end
        end
        #endregion
        #region Add the whole column
        (complex ? add_constr_linf_complex! : add_constr_linf!)(state, IndvalsIterator(unsafe, indices, values, lens))
        #endregion
        empty!(indices)
        empty!(values)
        empty!(lens)
    end

    if !prepend_fix(state)
        valat = dddual_transform_equalities!(state, Val(complex), dim, data, slacks, counters)
    end

    return complex ? ((u isa UpperOrUnitUpperTriangular ? :dd_lnorm_complex_triu :
                       (u isa LowerOrUnitLowerTriangular ? :dd_lnorm_complex_tril : :dd_lnorm_complex)),
                      (valat, addtocounter!(state, counters, Val(:lnorm_complex), dim, complex ? 2dim -1 : dim))) :
                     ((u isa UpperOrUnitUpperTriangular ? :dd_lnorm_real_triu :
                       (u isa LowerOrUnitLowerTriangular ? :dd_lnorm_real_tril : :dd_lnorm_real)),
                      (valat, addtocounter!(state, counters, Val(:lnorm), dim, complex ? 2dim -1 : dim)))
end

# We don't have the lnorm cone available, the problem is real-valued, and the transform is diagonal
function dddual_transform_cone!(state::AnySolver{T,V}, ::Val{false}, ::Val{false}, dim::Int, data::IndvalsIterator{T,V},
    u::DiagonalTransform, counters::Counters) where {T,V}
    ts = trisize(dim)
    dsq = dim^2
    # If we don't have the ℓ_∞ cone, we must use linear constraints. While we could do the whole matrix in a single large
    # nonnegative constraint vector, we'll do it columnwise as well.
    indices = FastVec{T}(buffer=2 * 3(dim -1))
    lens = FastVec{Int}(buffer=2dim -2)
    slacks = add_var_slack!(state, dsq)
    values = similar(indices, V)

    if prepend_fix(state)
        valat = dddual_transform_equalities!(state, Val(false), dim, data, slacks, counters)
    end

    upperslack = ts +1
    lowerslack = 1
    @inbounds for j in 1:dim
        #region Diagonal
        udiagval = u isa Diagonal ? u[j, j]^2 : one(V)
        diagslack = lowerslack
        lowerslack += 1
        #endregion
        #region Above diagonal (slacks)
        for i in 0:j-2
            unsafe_push!(indices, slacks[diagslack], slacks[upperslack+i], slacks[diagslack], slacks[upperslack+i])
            unsafe_push!(values, udiagval, one(V), udiagval, -one(V))
            unsafe_push!(lens, 2, 2)
        end
        #endregion
        #region Below diagonal
        upperslack += j -1
        sbelow = upperslack + j -1
        for i in j+1:dim
            uval = u isa Diagonal ? u[i, i] * u[j, j] : one(V)

            unsafe_push!(indices, slacks[diagslack], slacks[lowerslack], slacks[sbelow], slacks[diagslack], slacks[lowerslack],
                slacks[sbelow])
            unsafe_push!(values, udiagval, uval, -one(V), udiagval, -uval, one(V))
            unsafe_push!(lens, 3, 3)
            sbelow += i -1
            lowerslack += 1
        end
        #endregion
        #region Add the whole column
        add_constr_nonnegative!(state, IndvalsIterator(unsafe, indices, values, lens))
        #endregion
        empty!(indices)
        empty!(values)
        empty!(lens)
    end

    if !prepend_fix(state)
        valat = dddual_transform_equalities!(state, Val(false), dim, data, slacks, counters)
    end

    return :dd_nonneg_diag, (valat, addtocounter!(state, counters, Val(:nonnegative), dim, 2dim -2))
end

# We don't have the lnorm cone available, the problem is real-valued, and the transform is not diagonal
function dddual_transform_cone!(state::AnySolver{T,V}, ::Val{false}, ::Val{false}, dim::Int, data::IndvalsIterator{T,V}, u,
    counters::Counters) where {T,V}
    ts = trisize(dim)
    dsq = dim^2
    # If we don't have the ℓ_∞ cone, we must use linear constraints. While we could do the whole matrix in a single large
    # nonnegative constraint vector, we'll do it columnwise as well.
    indices = FastVec{T}(buffer=2 * (ts +1) * (dim -1))
    slacks = add_var_slack!(state, dsq)
    values = similar(indices, V)
    lens = FastVec{Int}(buffer=2 * (dim -1))

    if prepend_fix(state)
        valat = dddual_transform_equalities!(state, Val(false), dim, data, slacks, counters)
    end

    upperslack = ts +1
    @inbounds for j in 1:dim
        #region Above diagonal (slacks)
        if !isone(j)
            lowerslack = 1
            for col in 1:dim
                if !dddual_transform_inrange(u, j, col)
                    lowerslack += dim - col +1
                    continue
                end
                if dddual_transform_inrange(Val(:diag), u, j, col)
                    unsafe_push!(indices, slacks[lowerslack])
                    unsafe_push!(values, u[j, col]^2)
                end
                lowerslack += 1
                for row in col+1:dim
                    if dddual_transform_inrange(Val(:diag), u, j, row)
                        unsafe_push!(indices, slacks[lowerslack])
                        unsafe_push!(values, u[j, col] * u[j, row])
                    end
                    lowerslack += 1
                end
            end
            diagindices = @view(indices[:])
            diagvalues = @view(values[:])
            unsafe_push!(indices, slacks[upperslack])
            unsafe_push!(values, one(V))
            unsafe_append!(lens, Iterators.repeated(length(indices), 2 * (j -1)))
            unsafe_append!(indices, diagindices)
            unsafe_append!(values, diagvalues)
            unsafe_push!(indices, slacks[upperslack])
            unsafe_push!(values, -one(V))
            for i in 1:j-2
                unsafe_append!(indices, diagindices)
                unsafe_append!(values, diagvalues)
                unsafe_push!(indices, slacks[upperslack+i])
                unsafe_push!(values, one(V))
                unsafe_append!(indices, diagindices)
                unsafe_append!(values, diagvalues)
                unsafe_push!(indices, slacks[upperslack+i])
                unsafe_push!(values, -one(V))
            end
        end
        #endregion
        #region Below diagonal
        upperslack += j -1
        sbelow = upperslack + j -1
        for i in j+1:dim
            @twice negative true begin
                slack = 1
                startind = length(indices)
                for col in 1:dim
                    if !dddual_transform_inrange(u, j, col)
                        slack += dim - col +1
                        continue
                    end
                    if dddual_transform_inrange(Val(:diag), u, j, col) ||
                        dddual_transform_inrange(Val(:offdiag), u, i, col)
                        unsafe_push!(indices, slacks[slack])
                        unsafe_push!(values, u[j, col]^2 + (negative ? -2 : 2) * u[i, col] * u[j, col])
                    end
                    slack += 1
                    for row in col+1:dim
                        if dddual_transform_inrange(Val(:diag), u, j, row) ||
                            dddual_transform_inrange(Val(:offdiag), u, i, row)
                            uval = u[i, row] * u[j, col] + u[i, col] * u[j, row]
                            unsafe_push!(indices, slacks[slack])
                            unsafe_push!(values, u[j, col] * u[j, row] + (negative ? -uval : uval))
                        end
                        slack += 1
                    end
                end
                unsafe_push!(indices, slacks[sbelow])
                unsafe_push!(values, negative ? one(V) : -one(V))
                unsafe_push!(lens, length(indices) - startind)
            end
            sbelow += i -1
        end
        #endregion
        #region Add the whole column
        add_constr_nonnegative!(state, IndvalsIterator(unsafe, indices, values, lens))
        #endregion
        empty!(indices)
        empty!(values)
        empty!(lens)
    end

    if !prepend_fix(state)
        valat = dddual_transform_equalities!(state, Val(false), dim, data, slacks, counters)
    end

    return (u isa UpperOrUnitUpperTriangular ? :dd_nonneg_triu :
            (u isa LowerOrUnitLowerTriangular ? :dd_nonneg_tril : :dd_nonneg)),
           (valat, addtocounter!(state, counters, Val(:nonnegative), dim, 2 * (ts +1)))
end

# We don't have the lnorm cone available (but the quadratic cone), the problem is complex-valued, and the transform is diagonal
function dddual_transform_cone!(state::AnySolver{T,V}, ::Val{false}, ::Val{true}, dim::Int, data::IndvalsIterator{T,V},
    u::DiagonalTransform, counters::Counters) where {T,V}
    dsq = dim^2
    dd = dsq - dim
    # Here we use the quadratic cone to mimick the ℓ_∞ norm cone: x₁ ≥ ∑ᵢ (Re² xᵢ + Im² xᵢ). So we need to submit lots of
    # cones, but all of them pretty small.
    indices = FastVec{T}(buffer=7)
    slacks = add_var_slack!(state, dsq + dd)
    values = similar(indices, V)
    lens = [1, 3, 3] # let's not make it a StackVec to not compile yet another method

    if prepend_fix(state)
        valat = dddual_transform_equalities!(state, Val(true), dim, data, slacks, counters)
    end

    upperslack = dsq +1
    lowerslack = 1
    @inbounds for j in 1:dim
        resize!(indices, 3)
        resize!(values, 3)
        #region Diagonal
        indices[1] = slacks[lowerslack]
        values[1] = u isa Diagonal ? abs2(u[j, j]) : one(V)
        lowerslack += 1
        #endregion
        #region Above diagonal (slacks)
        values[2] = one(V)
        values[3] = one(V)
        for _ in 1:j-1
            indices[2] = slacks[upperslack]
            indices[3] = slacks[upperslack+1]
            upperslack += 2
            @inline add_constr_quadratic!(state, IndvalsIterator(unsafe, indices, values, 1))
        end
        #endregion
        resize!(indices, 7)
        resize!(values, 7)
        values[4] = -(values[7] = one(V))
        #region Below diagonal
        sbelow = upperslack + 2j -2
        for i in j+1:dim
            uval = u isa Diagonal ? u[i, i] * conj(u[j, j]) : one(V)

            indices[2] = slacks[lowerslack];   values[2] = real(uval)
            indices[3] = slacks[lowerslack+1]; values[3] = -imag(uval)
            indices[4] = slacks[sbelow]
            indices[5] = slacks[lowerslack];   values[5] = imag(uval)
            indices[6] = slacks[lowerslack+1]; values[6] = real(uval)
            indices[7] = slacks[sbelow+1]

            sbelow += 2i -2
            lowerslack += 2
            @inline add_constr_quadratic!(state, IndvalsIterator(unsafe, indices, values, lens))
        end
        #endregion
    end

    if !prepend_fix(state)
        valat = dddual_transform_equalities!(state, Val(true), dim, data, slacks, counters)
    end

    return :dd_quad_diag, (valat, addtocounter!(state, counters, Val(:quadratic), dd, 3))
end

# We don't have the lnorm cone available, the problem is complex-valued, and the transform is not diagonal
function dddual_transform_cone!(state::AnySolver{T,V}, ::Val{false}, ::Val{true}, dim::Int, data::IndvalsIterator{T,V}, u,
    counters::Counters) where {T,V}
    dsq = dim^2
    dd = dsq - dim
    # Here we use the quadratic cone to mimick the ℓ_∞ norm cone: x₁ ≥ ∑ᵢ (Re² xᵢ + Im² xᵢ). So we need to submit lots of
    # cones, but all of them pretty small.
    indices = FastVec{T}(buffer=3dim * (dim +1) +2)
    slacks = add_var_slack!(state, dsq + dd)
    values = similar(indices, V)
    lens = Vector{Int}(undef, 3)

    if prepend_fix(state)
        valat = dddual_transform_equalities!(state, Val(true), dim, data, slacks, counters)
    end

    upperslack = dsq +1
    @inbounds for j in 1:dim
        #region Diagonal
        lowerslack = 1
        for col in 1:dim
            if !dddual_transform_inrange(u, j, col)
                lowerslack += 1 + 2(dim - col)
                continue
            end
            if dddual_transform_inrange(Val(:diag), u, j, col)
                unsafe_push!(indices, slacks[lowerslack])
                unsafe_push!(values, abs2(u[j, col]))
            end
            lowerslack += 1
            for row in col+1:dim
                if dddual_transform_inrange(Val(:diag), u, j, row)
                    uval = u[j, col] * conj(u[j, row])
                    unsafe_push!(indices, slacks[lowerslack], slacks[lowerslack+1])
                    unsafe_push!(values, real(uval), imag(uval))
                end
                lowerslack += 2
            end
        end
        lens[1] = length(indices)
        #endregion
        #region Above diagonal (slacks)
        unsafe_push!(values, one(V), one(V))
        lens[3] = lens[2] = 1
        for _ in 1:j-1
            unsafe_push!(indices, slacks[upperslack], slacks[upperslack+1])
            upperslack += 2
            @inline add_constr_quadratic!(state, IndvalsIterator(unsafe, indices, values, lens))
            Base._deleteend!(indices, 2)
        end
        Base._deleteend!(values, 2)
        #endregion
        #region Below diagonal
        sbelow = upperslack + 2j - 2
        for i in j+1:dim
            @twice impart true begin
                startind = length(indices)
                slack = 1
                for col in 1:dim
                    if !dddual_transform_inrange(u, j, col)
                        slack += 1 + 2(dim - col)
                        continue
                    end
                    if dddual_transform_inrange(Val(:offdiag), u, i, col)
                        uval = 2u[i, col] * conj(u[j, col])
                        unsafe_push!(indices, slacks[slack])
                        unsafe_push!(values, impart ? imag(uval) : real(uval))
                    end
                    slack += 1
                    for row in col+1:dim
                        if dddual_transform_inrange(Val(:offdiag), u, i, row)
                            uval = u[i, row] * conj(u[j, col])
                            @twice imdata true let uval=uval
                                if imdata
                                    uval -= u[i, col] * conj(u[j, row])
                                    thisuval = impart ? real(uval) : -imag(uval)
                                else
                                    uval += u[i, col] * conj(u[j, row])
                                    thisuval = impart ? imag(uval) : real(uval)
                                end
                                unsafe_push!(indices, slacks[slack+imdata])
                                unsafe_push!(values, thisuval)
                            end
                        end
                        slack += 2
                    end
                end
                unsafe_push!(indices, slacks[sbelow])
                unsafe_push!(values, impart ? one(V) : -one(V))
                lens[2+impart] = length(indices) - startind
                sbelow += impart ? 2i -3 : 1
            end
            @inline add_constr_quadratic!(state, IndvalsIterator(unsafe, indices, values, lens))
            resize!(indices, lens[1])
            resize!(values, lens[1])
        end
        #endregion
        empty!(indices)
        empty!(values)
    end

    if !prepend_fix(state)
        valat = dddual_transform_equalities!(state, Val(true), dim, data, slacks, counters)
    end

    return (u isa UpperOrUnitUpperTriangular ? :dd_quad_triu :
            (u isa LowerOrUnitLowerTriangular ? :dd_quad_tril : :dd_quad)),
           (valat, addtocounter!(state, counters, Val(:quadratic), dd, 3))
end

function moment_add_dddual_transform!(state::AnySolver{T,V}, dim::Integer, data::IndvalsIterator{T,V}, u, complex,
    counters::Counters) where {T,V}
    !complex && (Base.IteratorEltype(u) isa Base.HasEltype) && eltype(u) <: Complex &&
        throw(MethodError(moment_add_dddual_transform!, (state, dim, data, u, complex)))
    @assert(dim > 1)

    return dddual_transform_cone!(
        state,
        Val(complex ? supports_lnorm_complex(state) : supports_lnorm(state)),
        Val(complex),
        dim,
        data,
        u,
        counters
    )
end
#endregion

#region Rewrite the SDD cone in terms of other cones.

# See the comment for the diagonally-dominant representation. Here, the fallback implementation is done in terms of rotated
# quadratic cones due to the relationship of SDD matrices with factor-width-2 matrices.
# We must take care of scaling the off-diagonal data; as we didn't know about the rotation, this could not have been done
# before. The rotated quadratic cone is 2x₁ x₂ ≥ ∑ᵢ xᵢ², so we'll scale the x₃ by √2 (i.e. multiply all the coefficients
# that use x₃ by 1/√2) to make this equivalent to [x₁ x₃; x₃ x₂] ⪰ 0. While only one triangle is considered, we the scaling by
# 2 is already done in the rewriting with slacks.

# For a general U and complex-valued data, we have the following rotated quadratic constraints for the column j and the row
# i > j:
# {∑_{col = 1}^dim ∑_{row = col}^dim (Re(U[j, col] Ū[j, row]) slackᵣ(row, col) + Im(U[j, col] Ū[j, row]) slackᵢ(row, col)),
#  ∑_{col = 1}^dim ∑_{row = col}^dim (Re(U[i, col] Ū[i, row]) slackᵣ(row, col) + Im(U[i, col] Ū[i, row]) slackᵢ(row, col)),
#  1/√2 ∑_{col = 1}^dim ∑_{row = col}^dim ((Re(U[i, row] Ū[j, col] + U[i, col] Ū[j, row]) slackᵣ(row, col) -
#                                          (Im(U[i, row] Ū[j, col] - U[i, col] Ū[j, row]) slackᵢ(row, col)),
#  1/√2 ∑_{col = 1}^dim ∑_{row = col}^dim ((Im(U[i, row] Ū[j, col] + U[i, col] Ū[j, row]) slackᵣ(row, col) +
#                                          (Re(U[i, row] Ū[j, col] - U[i, col] Ū[j, row]) slackᵢ(row, col))
# } ∈ ℛ𝒬₄

# Let's specialize the formula. If U is diagonal:
# {|U[j, j]|² slackᵣ(j, j),
#  |U[i, i]|² slackᵣ(i, i),
#  1/√2 (Re(U[i, i] Ū[j, j]) slackᵣ(i, j) - Im(U[i, i] Ū[j, j]) slackᵢ(i, j)),
#  1/√2 (Im(U[i, i] Ū[j, j]) slackᵣ(i, j) + Re(U[i, i] Ū[j, j]) slackᵢ(i, j))
# } ∈ ℛ𝒬₄

# If everything is instead real-valued:
# {∑_{col = 1}^dim ∑_{row = col}^dim U[j, col] U[j, row] slack(row, col),
#  ∑_{col = 1}^dim ∑_{row = col}^dim U[i, col] U[i, row] slack(row, col),
#  1/√2 ∑_{col = 1}^dim ∑_{row = col}^dim (U[i, row] U[j, col] + U[i, col] U[j, row]) slack(row, col)
# } ∈ ℛ𝒬₃

# If everything is real and U is diagonal:
# {U[j, j]² slack(j, j),
#  U[i, i]² slack(i, i),
#  1/√2 U[i, i] U[j, j] slack(i, j)
# } ∈ ℛ𝒬₃

# The transformation is diagonal
function sdddual_transform_cone!(state::AnySolver{T,V}, ::Val{complex}, dim::Integer, data::IndvalsIterator{T,V},
    u::DiagonalTransform, counters::Counters) where {T,V,complex}
    have_rot = supports_rotated_quadratic(state)
    scaling = inv(have_rot ? sqrt(V(2)) : V(2))
    indices = Vector{T}(undef, complex ? (have_rot ? 6 : 8) : (have_rot ? 3 : 5))
    lens = complex ? (have_rot ? [1, 1, 2, 2] : 2) : (have_rot ? 1 : [2, 2, 1])
    slacks = add_var_slack!(state, complex ? dim^2 : trisize(dim))
    values = similar(indices, V)

    if prepend_fix(state)
        valat = dddual_transform_equalities!(state, Val(complex), dim, data, slacks, counters)
    end

    rowdiagslack = 1
    offstart = (have_rot ? 3 : 5)
    @inbounds for j in 1:dim-1
        lowerslack = rowdiagslack +1
        urowval = u isa Diagonal ? u[j, j] : one(V)
        #region First item
        indices[1] = slacks[rowdiagslack]
        if have_rot
            values[1] = abs2(urowval)
        else
            indices[3] = slacks[rowdiagslack]
            urowvalsq = abs2(urowval)
        end
        rowdiagslack += (complex ? 2(dim - j) : dim - j) +1
        #endregion
        srowdiag = rowdiagslack
        for i in j+1:dim
            ucolval = u isa Diagonal ? u[i, i] : one(V)
            #region Second item
            indices[2] = slacks[srowdiag]
            if have_rot
                values[2] = abs2(ucolval)
            else
                indices[4] = slacks[srowdiag]
                ucolvalsq = abs2(ucolval)
                values[3] = values[1] = urowvalsq
                values[4] = -(values[2] = ucolvalsq)
            end
            srowdiag += (complex ? 2(dim - i) : dim - i) +1
            #endregion
            #region Third and fourth item
            uval = ucolval * conj(urowval) * scaling
            if complex
                indices[offstart+2] = indices[offstart+0] = slacks[lowerslack]
                indices[offstart+3] = indices[offstart+1] = slacks[lowerslack+1]
                values[offstart+0] = real(uval)
                values[offstart+1] = -imag(uval)
                values[offstart+2] = imag(uval)
                values[offstart+3] = real(uval)
                lowerslack += 2
            else
                indices[offstart] = slacks[lowerslack]
                values[offstart] = uval
                lowerslack += 1
            end
            #endregion
            @inline (have_rot ? add_constr_rotated_quadratic! : add_constr_quadratic!)(
                state, IndvalsIterator(unsafe, indices, values, lens)
            )
        end
    end

    if !prepend_fix(state)
        valat = dddual_transform_equalities!(state, Val(complex), dim, data, slacks, counters)
    end

    return (complex ? :sdd_quad_complex_diag : :sdd_quad_real_diag),
           (valat, addtocounter!(state, counters, Val(have_rot ? :rotated_quadratic : :quadratic),
                                 complex ? 4 : 3, trisize(dim -1)))
end

# The transform is not diagonal
function sdddual_transform_cone!(state::AnySolver{T,V}, ::Val{complex}, dim::Integer, data::IndvalsIterator{T,V}, u,
    counters::Counters) where {T,V,complex}
    have_rot = supports_rotated_quadratic(state)
    scaling = inv(have_rot ? sqrt(V(2)) : V(2))
    total = complex ? dim^2 : trisize(dim)
    indices = FastVec{T}(buffer=complex ? 4total : 3total)
    slacks = add_var_slack!(state, total)
    values = similar(indices, V)
    lens = Vector{Int}(undef, complex ? 4 : 3)

    if prepend_fix(state)
        valat = dddual_transform_equalities!(state, Val(complex), dim, data, slacks, counters)
    end

    if !have_rot
        first_values = FastVec{V}(buffer=total)
        first_indices = u isa UpperOrUnitUpperTriangular || u isa LowerOrUnitLowerTriangular ? FastVec{T}(buffer=total) :
                                                                                               indices
    else
        first_values = values
        first_indices = indices
    end
    @inbounds for j in 1:dim-1
        #region First item
        slack = 1
        @inbounds for col in 1:dim
            if !dddual_transform_inrange(u, j, col)
                slack += 1 + (complex ? 2 : 1) * (dim - col)
                continue
            end
            if dddual_transform_inrange(Val(:diag), u, j, col)
                unsafe_push!(first_indices, slacks[slack])
                unsafe_push!(first_values, abs2(u[j, col]))
            end
            slack += 1
            for row in col+1:dim
                if dddual_transform_inrange(Val(:diag), u, j, row)
                    uval = u[j, col] * conj(u[j, row])
                    if complex
                        unsafe_push!(first_indices, slacks[slack], slacks[slack+1])
                        unsafe_push!(first_values, real(uval), imag(uval))
                    else
                        unsafe_push!(first_indices, slacks[slack])
                        unsafe_push!(first_values, uval)
                    end
                end
                slack += complex ? 2 : 1
            end
        end
        lens[1] = firstlen = length(first_indices)
        if !have_rot
            resize!(values, length(first_indices))
            if u isa UpperOrUnitUpperTriangular || u isa LowerOrUnitLowerTriangular
                resize!(indices, length(first_indices))
            end
        end
        #endregion
        for i in j+1:dim
            #region Second item
            slack = 1
            for col in 1:dim
                if !dddual_transform_inrange(u, i, col)
                    slack += 1 + (complex ? 2 : 1) * (dim - col)
                    continue
                end
                if dddual_transform_inrange(Val(:diag), u, i, col)
                    unsafe_push!(indices, slacks[slack])
                    unsafe_push!(values, abs2(u[i, col]))
                end
                slack += 1
                for row in col+1:dim
                    if dddual_transform_inrange(Val(:diag), u, i, row)
                        uval = u[i, col] * conj(u[i, row])
                        if complex
                            unsafe_push!(indices, slacks[slack], slacks[slack+1])
                            unsafe_push!(values, real(uval), imag(uval))
                        else
                            unsafe_push!(indices, slacks[slack])
                            unsafe_push!(values, uval)
                        end
                    end
                    slack += complex ? 2 : 1
                end
            end
            startind = length(indices)
            lens[2] = startind - lens[1]
            #endregion
            #region Third and fourth item
            @twice impart complex begin
                slack = 1
                for col in 1:dim
                    if !dddual_transform_inrange(u, j, col)
                        slack += 1 + (complex ? 2 : 1) * (dim - col)
                        continue
                    end
                    if dddual_transform_inrange(Val(:offdiag), u, i, col)
                        uval = 2u[i, col] * conj(u[j, col])
                        unsafe_push!(indices, slacks[slack])
                        unsafe_push!(values, (impart ? imag(uval) : real(uval)) * scaling)
                    end
                    slack += 1
                    for row in col+1:dim
                        if dddual_transform_inrange(Val(:offdiag), u, i, row)
                            uval = u[i, row] * conj(u[j, col])
                            @twice imdata complex let uval=uval
                                if imdata
                                    uval -= u[i, col] * conj(u[j, row])
                                    thisuval = impart ? real(uval) : -imag(uval)
                                else
                                    uval += u[i, col] * conj(u[j, row])
                                    thisuval = impart ? imag(uval) : real(uval)
                                end
                                unsafe_push!(indices, slacks[slack+imdata])
                                unsafe_push!(values, thisuval * scaling)
                            end
                        end
                        slack += complex ? 2 : 1
                    end
                end
                lens[3+impart] = length(indices) - startind
                startind = length(indices)
            end
            #endregion
            if have_rot
                @inline add_constr_rotated_quadratic!(state, IndvalsIterator(unsafe, indices, values, lens))
            else
                if !(u isa UpperOrUnitUpperTriangular) && !(u isa LowerOrUnitLowerTriangular)
                    # This is a more efficient version of to_soc!, as we already know exactly all our duplicates: all entries
                    # have exactly the same indices, in the same order.
                    for i in 1:total
                        values[i], values[i+total] = (first_values[i] + values[i+total], first_values[i] - values[i+total])
                    end
                    @inline add_constr_quadratic!(state, IndvalsIterator(unsafe, indices, values, lens))
                else
                    @inbounds copyto!(indices, first_indices)
                    @inbounds copyto!(values, first_values)
                    @inline add_constr_quadratic!(state, to_soc!(indices, values, lens, false))
                    resize!(lens, complex ? 4 : 3) # just to make sure, to_soc! might trim zero items
                    lens[1] = firstlen
                end
            end
            resize!(values, firstlen)
            resize!(indices, firstlen)
        end
        empty!(first_indices)
        empty!(first_values)
    end

    if !prepend_fix(state)
        valat = dddual_transform_equalities!(state, Val(complex), dim, data, slacks, counters)
    end

    return (complex ? (u isa UpperOrUnitUpperTriangular ? :sdd_quad_complex_triu :
                       (u isa LowerOrUnitLowerTriangular ? :sdd_quad_complex_tril : :sdd_quad_complex)) :
                      (u isa UpperOrUnitUpperTriangular ? :sdd_quad_real_triu :
                       (u isa LowerOrUnitLowerTriangular ? :sdd_quad_real_tril : :sdd_quad_real))),
           (valat, addtocounter!(state, counters, Val(have_rot ? :rotated_quadratic : :quadratic), total, trisize(dim -1)))
end

function moment_add_sdddual_transform!(state::AnySolver{T,V}, dim::Integer, data::IndvalsIterator{T,V}, u, complex,
    counters::Counters) where {T,V}
    !complex && (Base.IteratorEltype(u) isa Base.HasEltype) && eltype(u) <: Complex &&
        throw(MethodError(moment_add_sdddual_transform!, (state, dim, data, u, complex)))
    @assert(dim > 2)

    return sdddual_transform_cone!(
        state,
        Val(complex),
        dim,
        data,
        u,
        counters
    )
end
#endregion

"""
    moment_add_matrix!(state::AbstractSolver, grouping::SimpleMonomialVector,
        constraint::Union{<:SimplePolynomial,<:AbstractMatrix{<:SimplePolynomial}},
        representation::RepresentationMethod=RepresentationPSD())

Parses a constraint in the moment hierarchy with a basis given in `grouping` (this might also be a partial basis due to
sparsity), premultiplied by `constraint` (which may be the unit polynomial for the moment matrix) and calls the appropriate
solver functions to set up the problem structure according to `representation`.

To make this function work for a solver, implement the following low-level primitives:
- [`mindex`](@ref)
- [`add_constr_nonnegative!`](@ref)
- [`add_constr_rotated_quadratic!`](@ref) (optional, then set [`supports_rotated_quadratic`](@ref) to `true`)
- [`add_constr_quadratic!`](@ref) (optional, then set [`supports_quadratic`](@ref) to `true`)
- [`add_constr_psd!`](@ref)
- [`add_constr_psd_complex!`](@ref) (optional, then set [`supports_psd_complex`](@ref) to `true`)
- [`add_constr_dddual!`](@ref) (optional, then set [`supports_dd`](@ref) to `true`)
- [`add_constr_dddual_complex!`](@ref) (optional, then set [`supports_dd_complex`](@ref) to `true`)
- [`add_constr_linf!`](@ref) (optional, then set [`supports_lnorm`](@ref) to `true`)
- [`add_constr_linf_complex!`](@ref) (optional, then set [`supports_lnorm_complex`](@ref) to `true`)
- [`add_constr_sdddual!`](@ref) (optional, then set [`supports_sdd`](@ref) to `true`)
- [`add_constr_sdddual_complex!`](@ref) (optional, then set [`supports_sdd_complex`](@ref) to `true`)
- [`psd_indextype`](@ref)
- [`add_var_slack!`](@ref)

Usually, this function does not have to be called explicitly; use [`moment_setup!`](@ref) instead.

See also [`moment_add_equality!`](@ref), [`RepresentationMethod`](@ref).
"""
function moment_add_matrix!(state::AnySolver{<:Any,V}, grouping::AbstractVector{M} where {M<:SimpleMonomial},
    constraint::Union{P,<:AbstractMatrix{P}}, representation::RepresentationMethod=RepresentationPSD(),
    counters::Counters=Counters()) where {P<:SimplePolynomial,V<:Real}
    dim = length(grouping) * (constraint isa AbstractMatrix ? LinearAlgebra.checksquare(constraint) : 1)
    if (dim == 1 || (dim == 2 && (supports_rotated_quadratic(state) || supports_quadratic(state))))
        if representation isa RepresentationPSD
            indextype = PSDIndextypeVector(:U, zero(V)) # we don't care about the scaling if everything is rewritten
        else
            return moment_add_matrix!(state, grouping, constraint, RepresentationPSD(), counters)
        end
    else
        indextype = psd_indextype(state)
    end

    real_valued = (length(grouping) == 1 || isreal(grouping)) && (!(constraint isa AbstractMatrix) || isreal(constraint))
    complex_cone = !real_valued &&
        (representation isa RepresentationDD{<:Any,true} &&
            (supports_dd_complex(state) || supports_lnorm_complex(state) || supports_quadratic(state))) ||
        # TODO (unlikely): we could still do it if the solver supports only the rotated, but not the standard cone. But
        # are there any solvers in this category?
        representation isa RepresentationSDD{<:Any,true} ||
        (representation isa RepresentationPSD &&
            (dim == 1 || (dim == 2 && (supports_rotated_quadratic(state) || supports_quadratic(state))) ||
             supports_psd_complex(state)))

    representation isa RepresentationSDD &&
        !supports_rotated_quadratic(state) && !supports_quadratic(state) &&
        ((complex_cone && !supports_sdd_complex(state)) || (!complex_cone && !supports_sdd(state))) &&
        error("The solver does not support the required scaled diagonally dominant cone or the fallback (rotated) quadratic cones, so a representation via scaled diagonally-dominant matrices is not possible")

    if representation isa RepresentationDD || representation isa RepresentationSDD
        try
            sizecheck = real_valued || complex_cone ? dim : 2dim
            representation.u[sizecheck, sizecheck] # just for bounds checking - we cannot access size, as u might be anything
                                                   # (e.g., UniformScaling)
        catch e
            e isa BoundsError &&
                throw(ArgumentError("The given matrix for rotating the DD cone was not large enough (required dimension: $dim)"))
            rethrow()
        end
    end
    if (representation isa RepresentationDD && !(complex_cone ? supports_dd_complex(state) : supports_dd(state))) ||
        (representation isa RepresentationSDD && !(complex_cone ? supports_sdd_complex(state) : supports_sdd(state)))
        indextype = PSDIndextypeVector(:L, zero(V))
    end

    indextype isa PSDIndextypeMatrixCartesian{:F} &&
        error("The Cartesian full matrix indextype is currently supported only for the primal moment optimization.")

    return moment_add_matrix_helper!(
        state,
        grouping,
        constraint isa AbstractMatrix ? constraint : ScalarMatrix(constraint),
        indextype,
        (Val(real_valued), Val(complex_cone)),
        representation,
        counters
    )
end

"""
    moment_add_equality!(state::AbstractSolver, grouping::SimpleMonomialVector, constraint::SimplePolynomial)

Parses a polynomial equality constraint for moments and calls the appropriate solver functions to set up the problem structure.
`grouping` contains the basis that will be squared in the process to generate the prefactor.

To make this function work for a solver, implement the following low-level primitives:
- [`add_constr_fix_prepare!`](@ref) (optional)
- [`add_constr_fix!`](@ref)
- [`add_constr_fix_finalize!`](@ref) (optional)

Usually, this function does not have to be called explicitly; use [`moment_setup!`](@ref) instead.

See also [`moment_add_matrix!`](@ref).
"""
function moment_add_equality!(state::AnySolver{T}, grouping::AbstractVector{M} where {M<:SimpleMonomial},
    constraint::P, counters::Counters=Counters()) where {T,Nr,Nc,I<:Integer,P<:SimplePolynomial{<:Any,Nr,Nc,<:SimpleMonomialVector{Nr,Nc,I}}}
    # keep in sync with SOSCertificate -> poly_decomposition

    # We need to traverse all unique elements in groupings * groupings†. For purely complex-valued groupings, this is the full
    # list; as soon as we have a real variable present, it is smaller.
    # To avoid rehashings, get an overestimator of the total grouping size first.
    # TODO (maybe): In the first loop to populate unique_groupings, we determine whether the grouping is real-valued. So we
    # could instead populate two sets, saving isreal and a lot of conditionals in the second loop.
    unique_groupings = sizehint!(Set{FastKey{I}}(), iszero(Nc) ? trisize(length(grouping)) : length(grouping)^2)
    real_grouping = true
    totalsize = 0
    for (i, g₁) in enumerate(grouping)
        if !iszero(Nc)
            g₁real = !iszero(Nr) && isreal(g₁)
            # Consider the g₂ = ḡ₁ case separately in the complex case. Explanations below.
            let g₂=g₁
                prodidx = FastKey(monomial_index(g₁, SimpleConjMonomial(g₂)))
                indexug, sh = Base.ht_keyindex2_shorthash!(unique_groupings.dict, prodidx)
                if indexug ≤ 0
                    @inbounds Base._setindex!(unique_groupings.dict, nothing, prodidx, -indexug, sh)
                    totalsize += 1
                end
            end
        end
        # In the real case, we can skip the entries behind i as they would lead to duplicates.
        # In the complex case, we can also skip them, as they would lead to exact conjugates, which in the end give rise to the
        # same conditions (but note that i is already handled above).
        for g₂ in Iterators.take(grouping, iszero(Nc) ? i : i -1)
            # We don't use mindex, as this can have unintended side-effects on the solver state (such as creating a
            # representation for this monomial, although we don't even know whether we need it - if constraint does not contain
            # a constant term, this function must not automatically add all the squared groupings as monomials, even if they
            # will probably appear at some place).
            prodidx = FastKey(monomial_index(g₁, SimpleConjMonomial(g₂)))
            # We need to add the product to the set if it does not exists; we also need to count the number of conditions that
            # we get out of it.
            indexug, sh = Base.ht_keyindex2_shorthash!(unique_groupings.dict, prodidx)
            if indexug ≤ 0
                # It does not exist.
                @inbounds Base._setindex!(unique_groupings.dict, nothing, prodidx, -indexug, sh)
                # Assume we have a grouping g = (gᵣ + im*gᵢ) and a polynomial p = pᵣ + im*pᵢ, where the individual parts are
                # real-valued. Then, add_equality! means that g*p = 0 and ḡ*p = 0. Of course we can also conjugate everything.
                # We must split each constraint into its real and imaginary parts:
                # (I)   Re(g*p) = gᵣ*pᵣ - gᵢ*pᵢ
                # (II)  Im(g*p) = gᵣ*pᵢ + gᵢ*pᵣ
                # (III) Re(ḡ*p) = gᵣ*pᵣ + gᵢ*pᵢ
                # (IV)  Im(ḡ*p) = gᵣ*pᵢ - gᵢ*pᵣ
                # To analyze this (which would be easier if we added and subtracted the equalities, but in the
                # SimplePolynomials setup, the given form is most easy to handle), let's consider linear dependencies.
                # - If the constraint is real-valued, (III) is equal to (I) and (IV) is -(II), so we only take (I) and (II).
                # - If the grouping is real-valued, (III) is equal to (I) and (IV) is equal to (II); we only take (I) and (II).
                # - If both are real-valued, (III) is equal to (I) while (II) and (IV) are zero, so we only take (I).
                # - If both are complex-valued, all constraints are linearly independent.
                # Rearranging this, we always take (I); if at least one is complex-valued, we also take (II); if both are, we
                # take all. Note that we don't have to consider the conjugates of the groupings separately, as they only yield
                # a global sign in the zero-equality.
                # For this loop, this means that we will only check whether g₁*ḡ₂ belongs to a real-valued monomial, in which
                # case we add 1; or to a complex-valued monomial, in which case we add 2. After the loop, we multiply by 2 if
                # the constraint was also complex-valued.
                if iszero(Nc) || (!iszero(Nr) && g₁real && isreal(g₂)) # note that g₁ ≠ ḡ₂
                    totalsize += 1
                else
                    totalsize += 2
                    real_grouping = false
                end
            end
        end
    end

    real_constr = isreal(constraint)
    if !real_constr
        totalsize *= 2
    end

    negate = @inline negate_fix(state)
    constrstate = @inline add_constr_fix_prepare!(state, totalsize)
    V = real(coefficient_type(P))
    indices₁ = FastVec{T}(buffer=2length(constraint))
    values₁ = similar(indices₁, V)
    # While we could conditionally define those variables only if the requirements are satisfied, the compiler might not be
    # able to infer that we only use them later on if the same conditions (potentially stricter) are met. So define them
    # always, but not using any memory.
    indices₂ = similar(indices₁, 0, buffer=real_grouping && real_constr ? 0 : 2length(constraint))
    values₂ = similar(indices₂, V)
    indices₃ = similar(indices₁, 0, buffer=real_grouping || real_constr ? 0 : 2length(constraint))
    values₃ = similar(indices₃, V)
    indices₄ = similar(indices₃)
    values₄ = similar(values₃)

    e = ExponentsAll{Nr+2Nc,I}()

    for grouping_idx in unique_groupings
        grouping = SimpleMonomial{Nr,Nc}(unsafe, e, convert(I, grouping_idx))
        real_grouping = isreal(grouping)
        skip₂ = real_constr && real_grouping
        for term_constr in constraint
            mon_constr = monomial(term_constr)
            coeff_constr = coefficient(term_constr)
            recoeff = real(coeff_constr)
            imcoeff = imag(coeff_constr)
            if negate
                recoeff = -recoeff
                imcoeff = -imcoeff
            end
            repart, impart, canonical = getreim(state, grouping, mon_constr)
            if !iszero(recoeff)
                pushorupdate!(indices₁, repart, values₁, recoeff)
                skip₂ || repart == impart || pushorupdate!(indices₂, impart, values₂, canonical ? recoeff : -recoeff)
            end
            if !iszero(imcoeff)
                repart == impart || pushorupdate!(indices₁, impart, values₁, canonical ? -imcoeff : imcoeff)
                skip₂ || pushorupdate!(indices₂, repart, values₂, imcoeff)
            end
            if !real_grouping && !real_constr
                repart₂, impart₂, canonical₂ = getreim(state, SimpleConjMonomial(grouping), mon_constr)
                if !iszero(recoeff)
                    pushorupdate!(indices₃, repart₂, values₃, recoeff)
                    repart₂ == impart₂ || pushorupdate!(indices₄, impart₂, values₄, canonical₂ ? recoeff : -recoeff)
                end
                if !iszero(imcoeff)
                    repart₂ == impart₂ || pushorupdate!(indices₃, impart₂, values₃, canonical₂ ? -imcoeff : imcoeff)
                    pushorupdate!(indices₄, repart₂, values₄, imcoeff)
                end
            end
        end
        constrstate = @inline add_constr_fix!(state, constrstate, Indvals(indices₁, values₁), zero(V))
        empty!(indices₁); empty!(values₁)
        if !skip₂
            if !isempty(indices₂)
                constrstate = @inline add_constr_fix!(state, constrstate, Indvals(indices₂, values₂), zero(V))
                empty!(indices₂); empty!(values₂)
            end
            if !real_grouping && !real_constr
                if !isempty(indices₃)
                    constrstate = @inline add_constr_fix!(state, constrstate, Indvals(indices₃, values₃), zero(V))
                    empty!(indices₃); empty!(values₃)
                end
                if !isempty(indices₄)
                    constrstate = @inline add_constr_fix!(state, constrstate, Indvals(indices₄, values₄), zero(V))
                    empty!(indices₄); empty!(values₄)
                end
            end
        end
    end
    @inline add_constr_fix_finalize!(state, constrstate)
    return addtocounter!(state, counters, Val(:fix), totalsize)
end

function _fix_setup!(state::AnySolver{T,V}, problem::Problem{P}, groupings::RelaxationGroupings, counters::Counters,
    info::Vector, infoᵢ::Int) where {T,V,P}
    # fixed items
    negate = @inline negate_fix(state)
    # fix constant term to 1
    if isone(problem.prefactor)
        @inline add_constr_fix_finalize!(
            state,
            add_constr_fix!(
                state,
                add_constr_fix_prepare!(state, 1),
                Indvals(StackVec(mindex(state, constant_monomial(P))), StackVec(negate ? -one(V) : one(V))),
                negate ? -one(V) : one(V)
            )
        )
    else
        let buffer=length(problem.prefactor), indices=FastVec{T}(; buffer), values=FastVec{V}(; buffer)
            for t in problem.prefactor
                # We know that if imaginary parts pop up somewhere, they will cancel out in the end, as the conjugates also
                # appear. So use the same strategy as always and only add canonical ones, but doubled.
                mon = monomial(t)
                coeff = coefficient(t)
                recoeff = real(coeff)
                imcoeff = imag(coeff)
                if negate
                    recoeff = -recoeff
                    imcoeff = -imcoeff
                end
                repart, impart, canonical = getreim(state, mon)
                if repart == impart
                    @assert(iszero(imcoeff)) # else the objective would not be real-valued
                    unsafe_push!(indices, repart)
                    unsafe_push!(values, recoeff)
                elseif canonical
                    if !iszero(recoeff)
                        unsafe_push!(indices, repart)
                        unsafe_push!(values, recoeff + recoeff)
                    end
                    if !iszero(imcoeff)
                        unsafe_push!(indices, impart)
                        unsafe_push!(values, -(imcoeff + imcoeff))
                    end
                end
            end
            @inline add_constr_fix_finalize!(
                state,
                add_constr_fix!(state, add_constr_fix_prepare!(state, 1), Indvals(indices, values), negate ? -one(V) : one(V))
            )
        end
    end
    addtocounter!(state, counters, Val(:fix), 1) # we don't store this info

    for (groupingsᵢ, constrᵢ) in zip(groupings.zeros, problem.constr_zero)
        info[infoᵢ] = this_info = Vector{Tuple{Symbol,Any}}(undef, length(groupingsᵢ))
        for (j, grouping) in enumerate(groupingsᵢ)
            this_info[j] = (:fix, moment_add_equality!(state, collect_grouping(grouping), constrᵢ, counters))
        end
        infoᵢ += 1
    end

    return infoᵢ
end

"""
    moment_setup!(state::AbstractSolver, relaxation::AbstractRelaxation,
        groupings::RelaxationGroupings[; representation])

Sets up all the necessary moment matrices, variables, constraints, and objective of a polynomial optimization problem
`problem` according to the values given in `grouping` (where the first entry corresponds to the basis of the objective, the
second of the equality, the third of the inequality, and the fourth of the PSD constraints).
The function returns a `Vector{<:Vector{<:Tuple{Symbol,Any}}}` that contains internal information on the problem. This
information is required to obtain dual variables and re-optimize the problem and should be stored in the `state`.

The following methods must be implemented by a solver to make this function work:
- [`mindex`](@ref)
- [`add_constr_nonnegative!`](@ref)
- [`add_constr_rotated_quadratic!`](@ref) (optional, then set [`supports_rotated_quadratic`](@ref) to `true`)
- [`add_constr_quadratic!`](@ref) (optional, then set [`supports_quadratic`](@ref) to `true`)
- [`add_constr_psd!`](@ref)
- [`add_constr_psd_complex!`](@ref) (optional, then set [`supports_psd_complex`](@ref) to `true`)
- [`add_constr_dddual!`](@ref) (optional, then set [`supports_dd`](@ref) to `true`)
- [`add_constr_dddual_complex!`](@ref) (optional, then set [`supports_dd_complex`](@ref) to `true`)
- [`add_constr_linf!`](@ref) (optional, then set [`supports_lnorm`](@ref) to `true`)
- [`add_constr_linf_complex!`](@ref) (optional, then set [`supports_lnorm_complex`](@ref) to `true`)
- [`add_constr_sdddual!`](@ref) (optional, then set [`supports_sdd`](@ref) to `true`)
- [`add_constr_sdddual_complex!`](@ref) (optional, then set [`supports_sdd_complex`](@ref) to `true`)
- [`psd_indextype`](@ref)
- [`add_constr_fix_prepare!`](@ref) (optional)
- [`add_constr_fix!`](@ref)
- [`add_constr_fix_finalize!`](@ref) (optional)
- [`fix_objective!`](@ref)
- [`add_var_slack!`](@ref)

!!! warning "Indices"
    The variable indices used in all solver functions directly correspond to the indices given back by [`mindex`](@ref).
    However, in a sparse problem there may be far fewer indices present; therefore, when the problem is finally given to the
    solver, care must be taken to eliminate all unused indices. The functionality provided by [`AbstractAPISolver`](@ref) and
    [`AbstractSparseMatrixSolver`](@ref) already takes care of this.

!!! info "Order"
    The individual constraint types can be added in any order (including interleaved).

!!! info "Representation"
    This function may also be used to describe simplified cones such as the (scaled) diagonally dominant one. The
    `representation` parameter can be used to define a representation that is employed for the individual groupings. This may
    either be an instance of a [`RepresentationMethod`](@ref) - which requires the method to be independent of the dimension of
    the grouping - or a callable. In the latter case, it will be passed as a first parameter an identifier[^3] of the current
    conic variable, and as a second parameter the side dimension of its matrix. The method must then return a
    [`RepresentationMethod`](@ref) instance.

    [^3]: This identifier will be a tuple, where the first element is a symbol - either `:objective`, `:nonneg`, or `:psd` - to
          indicate the general reason why the variable is there. The second element is an `Int` denoting the index of the
          constraint (and will be undefined for the objective, but still present to avoid extra compilation). The last element
          is an `Int` denoting the index of the grouping within the constraint/objective.

See also [`sos_setup!`](@ref), [`moment_add_matrix!`](@ref), [`moment_add_equality!`](@ref),
[`RepresentationMethod`](@ref).
"""
function moment_setup!(state::AnySolver{T,V}, relaxation::AbstractRelaxation{<:Problem{P}}, groupings::RelaxationGroupings;
    representation=RepresentationPSD()) where {T,V,P}
    problem = poly_problem(relaxation)
    (real(coefficient_type(problem.objective)) <: V) ||
        @warn("Expected value type for the solver $V might not be compatible with polynomial coefficient type $(real(coefficient_type(problem.objective)))")

    counters = Counters()
    info = Vector{Vector{<:Tuple{Symbol,Any}}}(undef, 1 + length(problem.constr_zero) + length(problem.constr_nonneg) +
                                                      length(problem.constr_psd))
    infoᵢ = prepend_fix(state) ? _fix_setup!(state, problem, groupings, counters, info, 2) : 2

    # SOS term for objective
    constantP = SimplePolynomial(constant_monomial(P), coefficient_type(problem.objective))
    @inbounds info[1] = this_info = Vector{Tuple{Symbol,Any}}(undef, length(groupings.obj))
    for (i, grouping) in enumerate(groupings.obj)
        g = collect_grouping(grouping)
        @inbounds this_info[i] = moment_add_matrix!(state, g, constantP,
            representation isa RepresentationMethod ? representation : representation((:objective, 0, i), length(g)), counters)
    end

    # localizing matrices
    @inbounds for (i, (groupingsᵢ, constrᵢ)) in enumerate(zip(groupings.nonnegs, problem.constr_nonneg))
        info[infoᵢ] = this_info = Vector{Tuple{Symbol,Any}}(undef, length(groupingsᵢ))
        for (j, grouping) in enumerate(groupingsᵢ)
            g = collect_grouping(grouping)
            this_info[j] = moment_add_matrix!(state, g, constrᵢ,
                representation isa RepresentationMethod ? representation : representation((:nonneg, i, j), length(g)),
                counters)
        end
        infoᵢ += 1
    end
    for (i, (groupingsᵢ, constrᵢ)) in enumerate(zip(groupings.psds, problem.constr_psd))
        info[infoᵢ] = this_info = Vector{Tuple{Symbol,Any}}(undef, length(groupingsᵢ))
        for (j, grouping) in enumerate(groupingsᵢ)
            g = collect_grouping(grouping)
            this_info[j] = moment_add_matrix!(state, g, constrᵢ,
                representation isa RepresentationMethod ? representation :
                                                          representation((:psd, i, j), length(g) * size(constrᵢ, 1)), counters)
        end
        infoᵢ += 1
    end

    prepend_fix(state) || _fix_setup!(state, problem, groupings, counters, info, infoᵢ)

    # Riesz functional in the objective
    let buffer=length(problem.objective), indices=FastVec{T}(; buffer), values=FastVec{V}(; buffer)
        for t in problem.objective
            mon = monomial(t)
            coeff = coefficient(t)
            recoeff = real(coeff)
            imcoeff = imag(coeff)
            repart, impart, canonical = getreim(state, mon)
            # repart is the variable that is associated with the real part of the monomial mon;
            # impart is the variable associated with the imaginary part of the canonicalized monomial mon.
            if repart == impart
                @assert(iszero(imcoeff)) # else the objective would not be real-valued
                unsafe_push!(indices, repart)
                unsafe_push!(values, recoeff)
            elseif canonical # skip over noncanonicals - to have a real-valued objective, we know how they must look like
                if !iszero(recoeff)
                    unsafe_push!(indices, repart)
                    unsafe_push!(values, recoeff + recoeff)
                end
                if !iszero(imcoeff)
                    unsafe_push!(indices, impart)
                    unsafe_push!(values, -(imcoeff + imcoeff))
                end
            end
        end
        fix_objective!(state, Indvals(finish!(indices), finish!(values)))
    end

    return info
end