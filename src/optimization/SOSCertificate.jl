export SOSCertificate

function id_to_index(relaxation::AbstractRelaxation, id::Symbol)
    id === :objective || throw(MethodError(id_to_index, (relaxation, state, constraint)))
    return id_to_index(relaxation, (id, 0))
end

function id_to_index(relaxation::AbstractRelaxation, (type, index)::Tuple{Symbol,Integer})
    type === :objective && !iszero(index) && return id_to_index(relaxation, (type, 0, index))
    groupings = Relaxation.groupings(relaxation)
    if type === :objective
        groupingsᵢ = groupings.obj
    elseif type === :zero
        groupingsᵢ = groupings.zeros[index]
    elseif type === :nonneg
        groupingsᵢ = groupings.nonnegs[index]
    elseif type === :psd
        groupingsᵢ = groupings.psds[index]
    else
        throw(MethodError(id_to_index, (relaxation, (type, index))))
    end
    if isone(length(groupingsᵢ))
        return id_to_index(relaxation, (type, index, 1))
    else
        throw(ArgumentError("There is more than just a single grouping for the given constraint; its index has to be specified"))
    end
end

function id_to_index(relaxation::AbstractRelaxation, (type, index, grouping)::Tuple{Symbol,Integer,Integer})
    groupings = Relaxation.groupings(relaxation)
    if type === :objective
        dim = length(groupings.obj[grouping])
        return 1
    end
    prob = poly_problem(relaxation)
    i = index +1
    if type === :zero
        dim = length(groupings.zeros[index][grouping])
        return i
    end
    i += length(prob.constr_zero)
    if type === :nonneg
        dim = length(groupings.nonnegs[index][grouping])
        return i
    end
    i += length(prob.constr_nonneg)
    if type === :psd
        dim = length(groupings.psds[index][grouping])
        return i
    end

    throw(MethodError(id_to_index, (relaxation, (type, index, grouping))))
end

"""
    sos_matrix(relaxation::AbstractRelaxation, state[, constraint])

Extracts the SOS matrices associated with a solved relaxation.
A `constraint` is identified in the same way as the first argument that is passed to a [`RepresentationMethod`](@ref) (see
[^1]). Short forms are allowed when there is just a single grouping.
Usually, a [`SOSCertificate`](@ref) is more desirable than the construction of individual SOS matrices.
"""
function sos_matrix end

function sos_matrix(relaxation::AbstractRelaxation, state, constraint)
    index = id_to_index(relaxation, constraint)
    solvertype, position = extract_info(state)[index][grouping]
    return sos_matrix(relaxation, state, dim, Val(solvertype), position, Solver.extract_sos_prepare(relaxation, state))
end

sos_matrix(relaxation::AbstractRelaxation, state, ::Int, type::Val{:fix}, position, rawdata) =
    Solver.extract_sos(relaxation, state, type, position, rawdata)

sos_matrix(relaxation::AbstractRelaxation, state, ::Int, type::Val{:nonnegative}, position, rawdata) =
    reshape(Solver.extract_sos(relaxation, state, type, position, rawdata), (1, 1))

function sos_matrix(relaxation::AbstractRelaxation, state, dim::Int, type::Val{:quadratic}, position, rawdata)
    @assert(dim == 2)
    # PSD(psd₁₁, √2 psd₂₁, psd₂₂) ⇔ RQUAD(psd₁₁, psd₂₂, √2 psd₂₁) ⇔ QUAD(psd₁₁ + psd₂₂, psd₁₁ - psd₂₂, 2 psd₂₁)
    # psd₁₁ psd₂₂ ≥ psd₂₁²        ⇔ 2psd₁₁ psd₂₂ ≥ 2psd₂₁²        ⇔ 4psd₁₁ psd₂₂ ≥ 4psd₂₁²
    data = Solver.extract_sos(relaxation, state, type, position, rawdata)
    T = eltype(data)
    half = inv(T(2))
    if length(data) == 3
        @inbounds return SPMatrix(2, [data[1] + data[2], data[3], data[1] - data[2]], :U)
    elseif length(data) == 4
        @inbounds return SPMatrix(2, Complex{T}[data[1] + data[2], Complex(data[3], data[4]),
                                                data[1] - data[2]], :U)
    else
        @assert(false) # SDD
    end
end

function sos_matrix(relaxation::AbstractRelaxation, state, dim::Int, type::Val{:rotated_quadratic}, position, rawdata)
    @assert(dim == 2)
    data = Solver.extract_sos(relaxation, state, type, position, rawdata)
    T = eltype(data)
    if length(data) == 3
        @inbounds return SPMatrix(2, [data[1], data[3] / T(sqrt(2)), data[2]])
    elseif length(data) == 4
        @inbounds return SPMatrix(2, Complex{T}[data[1], Complex(data[3], data[4]) / T(sqrt(2)), data[2]])
    else
        @assert(false)
    end
end

function sos_matrix(relaxation::AbstractRelaxation, state, dim::Int, type::Type, position, rawdata) where {Type<:Solver.VAL_MATRIX}
    data = Solver.extract_sos(relaxation, state, type, position, rawdata)
    complex = type isa Solver.VAL_COMPLEX
    complex_to_real = false
    if data isa AbstractVector
        itype = @inline Solver.psd_indextype(state)
        if !complex && Solver.trisize(2dim) == length(data)
            ddim = 2dim
            complex_to_real = true
        else
            ddim = dim
        end
        T = eltype(data)
        if !complex || T <: Complex
            if itype isa Solver.PSDIndextypeVector{:U}
                if isone(itype.scaling)
                    data = SPMatrix(dim, data, :U)
                elseif itype.scaling^2 ≈ 2
                    data = SPMatrix(ddim, data, :US)
                else
                    error("Unsupported: only off-diagonal scalings with the values 1 or √2 are permitted.")
                end
            elseif itype isa Solver.PSDIndextypeVector{:L}
                if isone(itype.scaling)
                    data = SPMatrix(dim, data, :L)
                elseif itype.scaling^2 ≈ 2
                    data = SPMatrix(ddim, data, :LS)
                else
                    error("Unsupported: only off-diagonal scalings with the values 1 or √2 are permitted.")
                end
            else
                itype::Solver.PSDIndextypeVector{:F}
                data = reshape(data, ddim, ddim)
            end
        else
            # We need to recreate the complex-valued matrix, but the vector does not contain imaginary parts for the
            # diagonal; so we need to re-create our data.
            newdata = Vector{Complex{T}}(undef, Solver.trisize(dim))
            newdataptr = Ptr{T}(pointer(newdata))
            dataptr = pointer(data)
            GC.@preserve data begin # newdata is preserved anyway
                if itype isa Solver.PSDIndextypeVector{:L}
                    @inbounds for j in 1:dim
                        copylen = 2(dim - j) # 2(j -1) complex values below diagonal
                        unsafe_store!(newdataptr, unsafe_load(dataptr)) # set real part of diagonal
                        unsafe_store!(newdataptr + sizeof(T), zero(T))  # set imaginary part of diagonal
                        unsafe_copyto!(newdataptr + 2sizeof(T), dataptr + sizeof(T), copylen)
                        copylen *= sizeof(T)
                        newdataptr += copylen + 2sizeof(T)
                        dataptr += copylen + sizeof(T)
                    end
                    if isone(itype.scaling)
                        return SPMatrix(dim, newdata, :L)
                    elseif itype.scaling^2 ≈ 2
                        return SPMatrix(dim, newdata, :LS)
                    else
                        error("Unsupported: only off-diagonal scalings with the values 1 or √2 are permitted.")
                    end
                else
                    @inbounds for j in 1:dim
                        copylen = 2j -1 # 2(j -1) complex values above diagonal + real part of diagonal
                        unsafe_copyto!(newdataptr, dataptr, copylen) # copy above diagonal
                        copylen *= sizeof(T)
                        unsafe_store!(newdataptr + copylen, zero(T)) # then set the imaginary part to zero
                        newdataptr += copylen + sizeof(T) # and skip to the next column
                        if itype isa Solver.PSDIndextypeVector{:U}
                            dataptr += copylen
                        else
                            dataptr += (2dim -1) * sizeof(T)
                        end
                    end
                    if itype isa Solver.PSDIndextypeVector{:F} || isone(itype.scaling)
                        return SPMatrix(dim, newdata, :U)
                    elseif itype.scaling^2 ≈ 2
                        return SPMatrix(dim, newdata, :US)
                    else
                        error("Unsupported: only off-diagonal scalings with the values 1 or √2 are permitted")
                    end
                end
            end
        end
    else
        if !complex && 2dim == LinearAlgebra.checksquare(data)
            complex_to_real = true
        end
    end
    if complex_to_real
        # For the dual problem, [X₁ X₃ᵀ; X₃ X₂] → (X₁ + X₂) + (X₃ - X₃ᵀ)im
        @assert(!(eltype(data) <: Complex))
        result = Matrix{Complex{eltype(data)}}(undef, dim, dim)
        @inbounds for j in 1:dim
            for i in j:dim
                result[i, j] = Complex(data[i, j] + data[i+dim, j+dim], data[i+dim, j] - data[i, j+dim])
            end
        end
        return Hermitian(result, :L)
    else
        return data
    end
end

function sos_matrix(relaxation::AbstractRelaxation, state, dim::Int, ::Solver.VAL_NOMATRIX_REAL,
                    (position, _)::Tuple{AbstractUnitRange,Any}, rawdata)
    data = Solver.extract_sos(relaxation, state, Val(:fix), position, rawdata)
    if Solver.trisize(2dim) == length(position)
        ddim = 2dim
        complex_to_real = true
    else
        ddim = dim
        complex_to_real = false
    end
    data = SPMatrix(ddim, data, :L)
    if complex_to_real
        @assert(!(eltype(data) <: Complex))
        result = Matrix{Complex{eltype(data)}}(undef, dim, dim)
        @inbounds for j in 1:dim
            for i in j:dim
                result[i, j] = Complex(data[i, j] + data[i+dim, j+dim], data[i+dim, j] - data[i, j+dim])
            end
        end
        return Hermitian(result, :L)
    else
        return data
    end
end

function sos_matrix(relaxation::AbstractRelaxation, state, dim::Int, ::Solver.VAL_NOMATRIX_COMPLEX,
                    pos::Tuple{AbstractUnitRange,Vararg}, rawdata)
    position = pos[1]
    data = Solver.extract_sos(relaxation, state, Val(:fix), position, rawdata)
    T = eltype(data)
    # We need to recreate the complex-valued matrix, but the vector does not contain imaginary parts for the diagonal; so we
    # need to re-create our data.
    newdata = Vector{Complex{T}}(undef, Solver.trisize(dim))
    newdataptr = Ptr{T}(pointer(newdata))
    dataptr = pointer(data)
    GC.@preserve data begin # newdata is preserved anyway
        @inbounds for j in 1:dim
            copylen = sizeof(T) * 2(dim - j) # 2(j -1) complex values below diagonal
            unsafe_store!(newdataptr, unsafe_load(dataptr)) # set real part of diagonal
            unsafe_store!(newdataptr + sizeof(T), zero(T))  # set imaginary part of diagonal
            unsafe_copyto!(newdataptr + 2sizeof(T), dataptr + sizeof(T), copylen)
            newdataptr += copylen + 2sizeof(T)
            dataptr += copylen + sizeof(T)
        end
    end
    return SPMatrix(dim, newdata, :L)
end

"""
    SOSCertificate(result::Result)

Construct a SOS certificate from a given optimization result. The returned object will pretty-print to show the decomposition
of the optimization problem in terms of a positivstellensatz.
To obtain the polynomials for the individual terms, the object can be indexed. The first index is one of `:objective`, `:zero`,
`:nonneg`, `:psd`; the second index is the number of the desired element (omitted for `:objective`); the last index is the
index of the desired grouping due to sparsity. If the last index is omitted, a vector over all groupings is returned.

The returned vectors of polynomials are, for `:objective` and `:nonneg`, to be summed over their squares; for `:psd`, the
returned matrix `m` of polynomials is to be left-multiplied with its adjoint: `m' * m`. For `:zero`, a single polynomial is
returned. If all these operations are carried out while multiplying with the corresponding prefactors (i.e., the constraints
themselves), the resulting polynomial should be equal to the original objective.
Note that this need not be the case; only if the relaxation level was sufficient will a SOS certificate in fact be valid.
"""
struct SOSCertificate{R,V}
    relaxation::R
    objective::V
    data::Vector{Vector{Any}}
end

function SOSCertificate(result::Result)
    relaxation = result.relaxation
    state = result.state
    prob = poly_problem(relaxation)

    grouping = Relaxation.groupings(relaxation)
    info = Solver.extract_info(state)
    rawdata = Solver.extract_sos_prepare(relaxation, state)
    data = Vector{Vector{Any}}(undef, length(info))

    data[1] = dataᵢ = Vector{Any}(undef, length(grouping.obj))
    for (j, (solvertype, position)) in enumerate(info[1])
        dataᵢ[j] = sos_matrix(relaxation, state, length(grouping.obj[j]), Val(solvertype), position, rawdata)
    end

    i = 1
    @inbounds for groupings in (grouping.zeros, grouping.nonnegs), groupingᵢ in groupings
        infoᵢ = info[i += 1]
        data[i] = dataᵢ = Vector{Any}(undef, length(infoᵢ))
        for (j, (solvertype, position)) in enumerate(infoᵢ)
            dataᵢ[j] = sos_matrix(relaxation, state, length(groupingᵢ[j]), Val(solvertype), position, rawdata)
        end
    end

    @inbounds for (constrᵢ, groupingᵢ) in zip(prob.constr_psd, grouping.psds)
        infoᵢ = info[i += 1]
        data[i] = dataᵢ = Vector{Any}(undef, length(infoᵢ))
        for (j, (solvertype, position)) in enumerate(infoᵢ)
            dataᵢ[j] = sos_matrix(relaxation, state, length(groupingᵢ[j]) * size(constrᵢ, 1), Val(solvertype), position,
                rawdata)
        end
    end
    return SOSCertificate(relaxation, result.objective, data)
end

function sos_decomposition(matrix, grouping::SimpleMonomialVector{Nr,Nc,I}, ϵ) where {Nr,Nc,I<:Integer}
    # What to do with negative eigenvalues? We might get very small negative eigenvalues just due to numerical error; this is
    # fine. We could raise an error if the smallest eigenvalue is smaller than -ϵ. But if ϵ = 0, the user wants to keep
    # everything, though obviously negative eigenvalues have to be discarded. So we decide not to raise an error at all.
    # However, this risks not catching an error in the optimization (or, beware, the solver implementation).
    if isone(length(grouping))
        # while we reshaped the vector into a matrix, let's skip the overload here
        @assert(isone(length(matrix)))
        val = first(matrix)
        return [SimplePolynomial([val > 0 ? sqrt(val) : zero(val)], grouping)]
    end
    matrixdim = LinearAlgebra.checksquare(matrix)
    n = length(grouping)
    @assert(n == matrixdim)

    eig = eigen(matrix)
    start = searchsortedfirst(eig.values, ϵ * eig.values[end])
    polynomials = Vector{polynomial_type(grouping, eltype(matrix))}(undef, matrixdim - start +1)
    @inbounds for (k, i) in enumerate(matrixdim:-1:start)
        val = sqrt(eig.values[i]) # we know this to be positive, as the last eigenvalue must for sure be positive and we
                                  # truncate by ϵ
        vec = @view(eig.vectors[:, i])
        threshold = maximum(abs, vec) * ϵ
        for j in 1:matrixdim
            if abs(vec[j]) < threshold
                vec[j] = zero(vec[j])
            else
                vec[j] *= val
            end
        end
        polynomials[k] = SimplePolynomial(collect(vec), grouping)
    end
    return polynomials
end

function sos_decomposition(::Type{Matrix}, matrix, grouping::SimpleMonomialVector{Nr,Nc,I}, ϵ) where {Nr,Nc,I<:Integer}
    matrixdim = LinearAlgebra.checksquare(matrix)
    n = length(grouping)
    d = matrixdim ÷ n
    @assert(d * n == matrixdim)

    eig = eigen(matrix)
    start = searchsortedfirst(eig.values, ϵ * eig.values[end])
    Mᵀ = Matrix{polynomial_type(grouping, eltype(matrix))}(undef, d, matrixdim - start +1) # our resulting SOS matrix is Mᵀ M
    @inbounds for (k, i) in enumerate(matrixdim:-1:start)
        val = sqrt(eig.values[i])
        vec = @view(eig.vectors[:, i])
        threshold = maximum(abs, vec) * ϵ
        for j in 1:matrixdim
            if abs(vec[j]) < threshold
                vec[j] = zero(eltype(vec))
            else
                vec[j] *= val
            end
        end
        coeffs = reshape(vec, d, n) # col-major: the second system is of size d, so here the first
        for j in 1:d
            @inbounds Mᵀ[j, k] = SimplePolynomial(coeffs[j, :], grouping)
        end
    end
    return transpose(Mᵀ)
end

@inline _trunc(x, ϵ) = abs(x) < ϵ ? zero(x) : x

function poly_decomposition(data, grouping::SimpleMonomialVector{Nr,Nc,I},
    constraint::SimplePolynomial{<:Any,Nr,Nc,<:SimpleMonomialVector{Nr,Nc,I}}, ϵ) where {Nr,Nc,I<:Integer}
    # this is a decomposition of the equality constraint prefactors, which are arbitrary polynomials, no longer SOS.
    # keep in sync with MomentHelpers -> moment_add_equality

    unique_groupings = sizehint!(Set{FastKey{I}}(), iszero(Nc) ? Solver.trisize(length(grouping)) : length(grouping)^2)
    real_grouping = true
    for (i, g₁) in enumerate(grouping)
        if !iszero(Nc)
            g₁real = !iszero(Nr) && isreal(g₁)
            let g₂=g₁
                prodidx = FastKey(monomial_index(g₁, SimpleConjMonomial(g₂)))
                indexug, sh = Base.ht_keyindex2_shorthash!(unique_groupings.dict, prodidx)
                indexug ≤ 0 && @inbounds Base._setindex!(unique_groupings.dict, nothing, prodidx, -indexug, sh)
            end
        end
        for g₂ in Iterators.take(grouping, iszero(Nc) ? i : i -1)
            prodidx = FastKey(monomial_index(g₁, SimpleConjMonomial(g₂)))
            indexug, sh = Base.ht_keyindex2_shorthash!(unique_groupings.dict, prodidx)
            if indexug ≤ 0
                @inbounds Base._setindex!(unique_groupings.dict, nothing, prodidx, -indexug, sh)
                if !(iszero(Nc) || (!iszero(Nr) && g₁real && isreal(g₂)))
                    real_grouping = false
                end
            end
        end
    end

    real_constr = isreal(constraint)
    real_constr && real_grouping || error("SOS extraction with complex-valued equality constraints is currently unsupported")
    e = ExponentsAll{Nr+2Nc,I}()
    mons = FastVec{I}(buffer=(real_grouping ? 1 : 2) * length(unique_groupings))
    coeffs = similar(mons, real_constr && real_grouping ? eltype(data) : Complex{eltype(data)})
    i = 1


    for grouping_idx in unique_groupings
        grouping = SimpleMonomial{Nr,Nc}(unsafe, e, convert(I, grouping_idx))
        unsafe_push!(mons, monomial_index(e, grouping))
        unsafe_push!(coeffs, data[i])
        i += 1
    end

    mons_v = finish!(mons)
    coeffs_v = finish!(coeffs)
    sort_along!(mons_v, coeffs_v)
    return SimplePolynomial(coeffs_v, SimpleMonomialVector{Nr,Nc}(e, mons_v))
end

function Base.show(io::IO, m::MIME"text/plain", cert::SOSCertificate, ϵ=1e-6)
    prob = poly_problem(cert.relaxation)
    complex = !isreal(prob)

    println(io, "Sum-of-squares certificate for polynomial optimization problem")
    if !isempty(prob.constr_psd)
        println(io, "Note: ℙ[X] = X' X")
    end
    show(io, m, prob.objective)
    if abs(cert.objective) ≥ ϵ
        print(io, " - ", cert.objective)
    end
    println(io)
    groupings = Relaxation.groupings(cert.relaxation)

    # 1. objective - print SOS decomposition
    beginning = true
    for (dataᵢ, groupingᵢ) in zip(cert.data[1], groupings.obj)
        for objpoly in sos_decomposition(dataᵢ, groupingᵢ, ϵ)
            if !beginning
                print(io, "+ ")
            else
                print(io, "= ")
                beginning = false
            end
            print(io, complex ? "|" : "(")
            show(io, m, objpoly)
            println(io, complex ? "|²" : ")²")
        end
    end

    i = 2
    # 2. zero constraints
    isempty(prob.constr_zero) || println(io, "\n# zero constraints")
    for (constr, grouping) in zip(prob.constr_zero, groupings.zeros)
        beginning = true
        print(io, "+ (")
        show(io, m, constr)
        println(io, ") * (")
        for (dataᵢ, groupingᵢ) in zip(cert.data[i], grouping)
            print(io, "   ")
            if !beginning
                print(io, "+ ")
            else
                beginning = false
            end
            print(io, "(") # Not really necessary, but by indicating the grouping we give a reason why we don't sum up
                           # everything.
            show(io, m, poly_decomposition(dataᵢ, groupingᵢ, constr, ϵ))
            println(io, ")")
        end
        println(io, ")")
        i += 1
    end

    # 3. nonnegative constraints
    isempty(prob.constr_nonneg) || println(io, "\n# nonnegative constraints")
    for (constr, grouping) in zip(prob.constr_nonneg, groupings.nonnegs)
        beginning = true
        print(io, "+ (")
        show(io, m, constr)
        println(io, ") * (")
        for (dataᵢ, groupingᵢ) in zip(cert.data[i], grouping)
            for constrpoly in sos_decomposition(dataᵢ, groupingᵢ, ϵ)
                print(io, "   ")
                if !beginning
                    print(io, "+ ")
                else
                    beginning = false
                end
                print(io, complex ? "|" : "(")
                show(io, m, constrpoly)
                println(io, complex ? "|²" : ")²")
            end
        end
        println(io, ")")
        i += 1
    end

    # 4. PSD constraints
    isempty(prob.constr_psd) || println(io, "\n# PSD constraints")
    for (constr, grouping) in zip(prob.constr_psd, groupings.psds)
        beginning = true
        Base.print_matrix(io, constr, "+ ⟨[", "  ", "],")
        for (dataᵢ, groupingᵢ) in zip(cert.data[i], grouping)
            constrpolymat = sos_decomposition(Matrix, dataᵢ, groupingᵢ, ϵ)
            println(io)
            if beginning
                pre = "   ℙ["
                beginning = false
            else
                pre = "   + ℙ["
            end
            Base.print_matrix(io, constrpolymat, pre, "  ", "]")
        end
        println(io, "⟩")
        i += 1
    end

    return
end

function Base.getindex(s::SOSCertificate, type::Symbol)
    type === :objective || throw(MethodError(getindex, (s, type)))
    @inbounds getindex(s, :objective, 0)
end
Base.@propagate_inbounds function Base.getindex(s::SOSCertificate, type::Symbol, index::Integer; ϵ=0)
    type === :objective && !iszero(index) && return getindex(s, type, 0, index)
    idx = id_to_index(s.relaxation, (type, index, 1))
    groupings = Relaxation.groupings(s.relaxation)
    if type === :objective
        grouping = groupings.obj
    elseif type === :zero
        grouping = groupings.zeros[index]
        return [poly_decomposition(m, g, poly_problem(s.relaxation).constr_zero[index], ϵ)
                for (m, g) in zip(s.data[idx], grouping)]
    elseif type === :nonneg
        grouping = groupings.nonnegs[index]
    elseif type === :psd
        grouping = groupings.psds[index]
        return [sos_decomposition(Matrix, m, g, ϵ) for (m, g) in zip(s.data[idx], grouping)]
    end
    return [sos_decomposition(m, g, ϵ) for (m, g) in zip(s.data[idx], grouping)]
end
Base.@propagate_inbounds function Base.getindex(s::SOSCertificate, type::Symbol, index::Integer, grouping::Integer; ϵ=0)
    idx = id_to_index(s.relaxation, (type, index, grouping))
    groupings = Relaxation.groupings(s.relaxation)
    if type === :objective
        groupingᵢ = groupings.obj
    elseif type === :zero
        groupingᵢ = groupings.zeros[index]
        return poly_decomposition(s.data[idx][grouping], groupingᵢ[grouping], poly_problem(s.relaxation).constr_zero[index], ϵ)
    elseif type === :nonneg
        groupingᵢ = groupings.nonnegs[index]
    elseif type === :psd
        groupingᵢ = groupings.psds[index]
        return sos_decomposition(Matrix, s.data[idx][grouping], groupingᵢ[grouping], ϵ)
    end
    return sos_decomposition(s.data[idx][grouping], groupingᵢ[grouping], ϵ)
end