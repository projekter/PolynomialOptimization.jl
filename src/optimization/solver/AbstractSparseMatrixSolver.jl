export AbstractSparseMatrixSolver, SparseMatrixCOO, coo_to_csc!

"""
    AbstractSparseMatrixSolver{I<:Integer,K<:Integer,V<:Real}

Superclass for a solver that requires its data in sparse matrix form. The data is aggregated in COO form using
[`append!`](@ref append!(::SparseMatrixCOO{I,K,V,Offset}, ::IndvalsIterator{K,V}) where {I<:Integer,K<:Integer,V<:Real,Offset})
and can be converted to CSC form using [`coo_to_csc!`](@ref). The type of the indices in final CSC form is `I`, where the
monomials during construction will be represented by numbers of type `K`.
Any type inheriting from this class is supposed to have a field `slack::K` which is initialized to `-one(K)` if `K` is signed
or `typemax(K)` if it is unsigned.

See also [`SparseMatrixCOO`](@ref).
"""
abstract type AbstractSparseMatrixSolver{I<:Integer,K<:Integer,V<:Real} <: AbstractSolver{K,V} end

Solver.mindex(::AbstractSparseMatrixSolver{<:Integer,K,<:Real}, monomials::SimpleMonomialOrConj{Nr,Nc}...) where {K,Nr,Nc} =
    monomial_index(monomials...)::K

function Solver.add_var_slack!(state::AbstractSparseMatrixSolver{<:Integer,K}, num::Int) where {K}
    stop = state.slack
    state.slack -= num
    return (state.slack + one(K)):stop
end

# necessary as SOS helper, as the state could also be a SOS state
Solver.add_constr_slack!(state::AbstractSparseMatrixSolver{<:Integer}, num::Int) = add_var_slack!(state, num)

"""
    SparseMatrixCOO{I<:Integer,K<:Integer,V<:Real,Offset}

Representation of a sparse matrix in COO form. Fields are `rowinds::FastVec{I}`, `moninds::FastVec{K}` (where `K` is of the
type returned by `monomial_index`), and `nzvals::FastVec{V}`. The first row/column for the solver has index `Offset` (of type
`I`).
"""
struct SparseMatrixCOO{I<:Integer,K<:Integer,V<:Real,Offset}
    rowinds::FastVec{I}
    moninds::FastVec{K}
    nzvals::FastVec{V}

    function SparseMatrixCOO{I,K,V,Offset}() where {I<:Integer,K<:Integer,V<:Real,Offset}
        Offset isa I || throw(MethodError(SparseMatrixCOO{I,K,V,Offset}, ()))
        new{I,K,V,Offset}(FastVec{I}(), FastVec{K}(), FastVec{V}())
    end
end

Base.length(smc::SparseMatrixCOO) = length(smc.rowinds)
@inline function Base.size(smc::SparseMatrixCOO{<:Integer,<:Integer,<:Real,Offset}, dim) where {Offset}
    dim == 1 || error("Not implemented")
    @inbounds return isempty(smc.rowinds) ? 0 : Int(smc.rowinds[end]) + (1 - Int(Offset))
end

function FastVector.prepare_push!(smc::SparseMatrixCOO, new_items::Integer)
    prepare_push!(smc.rowinds, new_items)
    prepare_push!(smc.moninds, new_items)
    prepare_push!(smc.nzvals, new_items)
    return smc
end

"""
    append!(coo::SparseMatrixCOO, indvals::Indvals)

Appends the data given in `indvals` into the next row in `coo`. Returns the index of the row that was added.

See also [`Indvals`](@ref).
"""
@inline function Base.append!(coo::SparseMatrixCOO{I,K,V,Offset}, indvals::Indvals{K,V}) where {I<:Integer,K<:Integer,V<:Real,Offset}
    prepare_push!(coo, length(indvals))
    @inbounds v = isempty(coo.rowinds) ? Offset : coo.rowinds[end] + one(I)
    for (monind, nzval) in indvals
        unsafe_push!(coo.rowinds, v)
        unsafe_push!(coo.moninds, monind)
        unsafe_push!(coo.nzvals, nzval)
    end
    return v
end

"""
    append!(coo::SparseMatrixCOO, indvals::IndvalsIterator)

Appends the data given in `indvals` into successive rows in `coo` (`first(indvals)` to the first row, the next to the second,
...). Returns the range of indices of all rows that were added.

See also [`IndvalsIterator`](@ref).
"""
@inline function Base.append!(coo::SparseMatrixCOO{I,K,V,Offset}, indvals::IndvalsIterator{K,V}) where {I<:Integer,K<:Integer,V<:Real,Offset}
    prepare_push!(coo.rowinds, length(rowvals(indvals)))
    @inbounds v = isempty(coo.rowinds) ? Offset : coo.rowinds[end] + one(I)
    firstv = v
    for l in Base.index_lengths(indvals)
        for _ in 1:l
            unsafe_push!(coo.rowinds, v)
        end
        v += one(I)
    end
    append!(coo.moninds, rowvals(indvals))
    append!(coo.nzvals, nonzeros(indvals))
    return firstv:(v - one(I))
end

struct COO_to_CSC_callback{Offset,CC<:Tuple,CP<:Tuple,CV<:Tuple,DV<:Tuple,I<:Tuple}
    coo_columns::CC
    csc_colptrs::CP
    coo_vecs::CV
    dense_vecs::DV
    ivecs::I

    function COO_to_CSC_callback{Offset}(coo_columns::Tuple, csc_colptrs::Tuple, coo_vecs::Tuple, dense_vecs::Tuple) where {Offset}
        ivecs = ntuple(_ -> Ref(1), length(coo_vecs))
        new{Offset,typeof(coo_columns),typeof(csc_colptrs),typeof(coo_vecs),typeof(dense_vecs),typeof(ivecs)}(
            coo_columns, csc_colptrs, coo_vecs, dense_vecs, ivecs
        )
    end
end

function (ctc::COO_to_CSC_callback{Offset,<:Tuple{Any}})(index, i) where {Offset}
    @inbounds begin
        colidx = ctc.coo_columns[1][i]
        for (ivec, coo_vec, dense_vec) in zip(ctc.ivecs, ctc.coo_vecs, ctc.dense_vecs)
            if ivec[] ≤ length(coo_vec[1]) && coo_vec[1][ivec[]] == colidx
                dense_vec[index] = coo_vec[2][ivec[]]
                ivec[] += 1
            end
        end
        ctc.csc_colptrs[1][index+1] = i + Offset
    end
    return
end

function (ctc::COO_to_CSC_callback{Offset,<:Tuple{Any,Any}})(index, i₁, i₂) where {Offset}
    @inbounds begin
        colidx = ismissing(i₁) ? ctc.coo_columns[2][i₂] : ctc.coo_columns[1][i₁]
        for (ivec, coo_vec, dense_vec) in zip(ctc.ivecs, ctc.coo_vecs, ctc.dense_vecs)
            if ivec[] ≤ length(coo_vec[1]) && coo_vec[1][ivec[]] == colidx
                dense_vec[index] = coo_vec[2][ivec[]]
                ivec[] += 1
            end
        end
        ctc.csc_colptrs[1][index+1] = ismissing(i₁) ? ctc.csc_colptrs[1][index] : i₁ + Offset
        ctc.csc_colptrs[2][index+1] = ismissing(i₂) ? ctc.csc_colptrs[2][index] : i₂ + Offset
    end
    return
end

"""
    coo_to_csc!(coo::Union{SparseMatrixCOO{I,K,V}},Tuple{FastVec{K},FastVec{V}}}...)

Converts sparse COO matrix or vector representations, where the monomial indices of the `coo` matrices or the entries of the
vectors can be arbitrarily sparse, to a CSC-based matrix representation with continuous columns, and the vectors are converted
to dense ones. No more than two matrices may be supplied. The input data may be mutated; and this mutated data must be passed
on to [`MomentVector`](@ref MomentVector(::AbstractRelaxation, ::Vector{V}, ::SparseMatrixCOO{<:Integer,K,V,Offset}, ::SparseMatrixCOO{<:Integer,K,V,Offset}...) where {K<:Integer,V<:Real,Offset}).
The following values are returned:
- number of distinct columns
- for each input, if it is a matrix, a tuple containing the colptr, rowval, nzval vectors
- for each input, if it is a vector, the corresponding dense vector
"""
@generated function coo_to_csc!(coo::Union{SparseMatrixCOO{I,K,V,Offset},<:Tuple{AbstractVector{K},AbstractVector{V}}}...) where {I<:Integer,K<:Integer,V<:Real,Offset}
    result = Expr(:block)
    mats = Int[]
    vecs = Int[]
    for (i, data) in enumerate(coo)
        if data <: SparseMatrixCOO
            length(mats) == 2 && throw(MethodError(coo_to_csc!, coo)) # no more than two such matrices are allowed
            # sort according to col indices, but keep row indices in order (sort_along is in general not order-preserving, so
            # use row indices as secondary option)
            push!(result.args, :(sort_along!(coo[$i].moninds, coo[$i].rowinds, coo[$i].nzvals, relevant=2)))
            push!(mats, i)
        else
            push!(result.args, :(sort_along!(coo[$i][1], coo[$i][2])))
            push!(vecs, i)
        end
    end
    isempty(mats) && throw(MethodError(coo_to_csc!, coo))
    # how many distinct monomials do we have? - we only count the matrices, which contains the conic data
    push!(result.args, :(moncount = count_uniques($((:(coo[$i].moninds) for i in mats)...),)))
    for i in mats
        push!(result.args,
            :($(Symbol(:colptr, i)) = Vector{I}(undef, moncount +1)),
            :(@inbounds $(Symbol(:colptr, i))[1] = $Offset)
        )
    end
    append!(result.args, (:($(Symbol(:vec, i)) = zeros(V, moncount)) for i in vecs))
    push!(result.args, :(count_uniques(
        $((:(coo[$i].moninds) for i in mats)...),
        COO_to_CSC_callback{Offset}(
            ($((:(coo[$i].moninds) for i in mats)...),),
            ($((Symbol(:colptr, i) for i in mats)...),),
            ($((:(coo[$i]) for i in vecs)...),),
            ($((Symbol(:vec, i) for i in vecs)...),)
        )
    )))
    returning = Expr(:tuple, :(moncount))
    for (i, data) in enumerate(coo)
        if data <: SparseMatrixCOO
            push!(returning.args, :(($(Symbol(:colptr, i)), finish!(coo[$i].rowinds), finish!(coo[$i].nzvals))))
        else
            push!(returning.args, :($(Symbol(:vec, i))))
        end
    end
    push!(result.args, :(return $returning))
    result
end

"""
    MomentVector(relaxation::AbstractRelaxation, moments::Vector{<:Real},
        slack::Integer, coo::SparseMatrixCOO...)

Given the moments vector as obtained from a [`AbstractSparseMatrixSolver`](@ref) solver, convert it to a
[`MomentVector`](@ref). Note that this function is not fully type-stable, as the result may be based either on a dense or
sparse vector depending on the relaxation. To establish the mapping between the solver output and the actual moments, all the
column-sorted COO data (i.e., as returned by [`coo_to_csc!`](@ref)) used in the problem construction needs to be passed on.
`slack` must contain the current value of the `slack` field of the `AbstractSparseMatrixSolver`.
"""
function MomentVector(relaxation::AbstractRelaxation{<:Problem{<:SimplePolynomial{<:Any,Nr,Nc}}}, _moments::Vector{V},
    slack::K, coo₁::SparseMatrixCOO{<:Integer,K,V,Offset}, cooₙ::SparseMatrixCOO{<:Integer,K,V,Offset}...) where {Nr,Nc,K<:Integer,V<:Real,Offset}
    # we need at least one coo here for dispatch
    if slack isa Signed
        coo₁moninds = @view(coo₁.moninds[searchsortedfirst(coo₁.moninds, one(slack)):end])
        cooₙmoninds = ((@view(coo.moninds[searchsortedfirst(coo.moninds, one(slack)):end]) for coo in cooₙ)...,)
    else
        coo₁moninds = @view(coo₁.moninds[1:searchsortedlast(coo₁.moninds, slack)])
        cooₙmoninds = ((@view(coo.moninds[1:searchsortedlast(coo.moninds, slack)]) for coo in cooₙ)...,)
    end
    max_mons = max(coo₁moninds[end], (moninds[end] for moninds in cooₙmoninds)...) # coos are sorted according to the columns
    moments = @view(_moments[slack isa Signed ? (-slack:K(length(_moments))) : (one(K):K(length(_moments))-(typemax(slack)-slack))])
    @assert(length(moments) ≤ max_mons)
    if length(moments) == max_mons # (real) dense case
        solution = moments
    else
        # We need to build the vector of monomial indices.
        mon_pos = Vector{K}(undef, length(moments))
        if length(cooₙ) == 0
            count_uniques(coo₁moninds, @capture((index, i) -> @inbounds begin
                $mon_pos[index] = $coo₁moninds[i]
            end))
        elseif length(cooₙ) == 1
            count_uniques(coo₁moninds, cooₙmoninds[1], @capture((index, i₁, i₂) -> @inbounds begin
                $mon_pos[index] = ismissing(i₁) ? $cooₙmoninds[1][i₂] : $coo₁moninds[i₁]
            end))
        else
            throw(MethodError(MomentVector, (relaxation, moments, coo₁, cooₙ...)))
        end
        if 3length(moments) < max_mons
            solution = SparseVector(max_mons, mon_pos, iszero(slack) ? _moments : collect(moments))
        else
            solution = fill(NaN, max_mons)
            copyto!(@view(solution[mon_pos]), moments)
        end
    end
    return MomentVector(relaxation, ExponentsAll{Nr+2Nc,K}(), solution)
end