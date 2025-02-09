export primal_moment_setup!, coo_to_csr!

mutable struct PrimalMomentSolver{I<:Integer,K<:Integer,V<:Real,C,B,P<:AbstractSolver{K,V}} <: AbstractSolver{K,V}
    const parent::P
    max_con::I
    max_nonneg::I
    const psd_dim::FastVec{I}
    const psds::FastVec{Tuple{FastVec{I},C,FastVec{V}}} # COO
    const nonnegs::Tuple{FastVec{I},FastVec{I},FastVec{V}} # COO
    const b::Tuple{FastVec{I},FastVec{V}}
    const mon_eq::Dict{FastKey{K},Tuple{Int,I,B}} # index 1: > 0 - psd index; -1: single nonnegative; -2: two nonnegatives
                                                  # index 2: if vectorized or idx1 < 0 - index in the vector
                                                  #          else - row (with offset)
                                                  # index 3: if vectorized, but not full: is this a diagonal index?
                                                  #          if vectorized and full: index of transpose
                                                  #          if idx1 < 0: unused
                                                  #          else - col (with offset)
                                                  # index (or row/col) will always refer to the entry in the lower triangle
                                                  # (if the lower or full was requested) or to the upper (if the upper was
                                                  # requested)
    c::Indvals{K,V}

    function PrimalMomentSolver{I,K,V}(parent::P) where {I<:Integer,K<:Integer,V<:Real,P<:AbstractSolver{K,V}}
        it = psd_indextype(parent)
        it isa PSDIndextypePrimal ||
            throw(ArgumentError("The indextype of $parent is not suitable for the primal moment solver"))
        C = it isa PSDIndextypeCOOVectorized ? FastVec{I} : Tuple{FastVec{I},FastVec{I}}
        B = it isa PSDIndextypeCOOVectorized && !(it isa PSDIndextypeCOOVectorized{:F}) ? Bool : I
        new{I,K,V,C,B,P}(
            parent,
            _get_offset(it) - one(I), _get_offset(it) - one(I),
            FastVec{I}(),
            FastVec{Tuple{FastVec{I},C,FastVec{V}}}(),
            (FastVec{I}(), FastVec{I}(), FastVec{V}()),
            (FastVec{I}(), FastVec{V}()),
            Dict{FastKey{K},Tuple{Int,I,B}}(),
        )
    end
end

function Base.getproperty(s::PrimalMomentSolver, f::Symbol)
    f === :num_con && return Int(getfield(s, :max_con) - _get_offset(psd_indextype(getfield(s, :parent)))) +1
    f === :num_nonneg && return Int(getfield(s, :max_nonneg) - _get_offset(psd_indextype(getfield(s, :parent)))) +1
    return getfield(s, f)
end

mindex(state::PrimalMomentSolver, monomials::SimpleMonomialOrConj{Nr,Nc}...) where {Nr,Nc} = mindex(state.parent, monomials...)

psd_indextype(state::PrimalMomentSolver) = PSDIndextypeVector(psd_indextype(state.parent) isa PSDIndextypePrimal{:U} ? :U : :L,
                                                              true)

prepend_fix(::PrimalMomentSolver) = false

addtocounter!(state::PrimalMomentSolver, counters::Counters, type::Val, dim::Integer) =
    addtocounter!(state.parent, counters, type, dim)

addtocounter!(state::PrimalMomentSolver, counters::Counters, type::Val, num::Integer, dim::Integer) =
    addtocounter!(state.parent, counters, type, num, dim)

add_var_slack!(::PrimalMomentSolver, ::Int) = error("Slack variables are not supported for the primal moment setup")

@inline function _primal_push_psd!((cooᵢ, cooⱼ, cooᵥ)::Tuple{FastVec{I},FastVec{I},FastVec{V}},
    it::PSDIndextypeCOOVectorized{:F}, max_con::I, idx::I, idxᵀ::I, v::V) where {I,V}
    if idx == idxᵀ
        push!(cooᵢ, max_con)
        push!(cooⱼ, idx)
        push!(cooᵥ, v)
    else
        push!(cooᵢ, max_con, max_con)
        push!(cooⱼ, idx, idxᵀ)
        v *= inv(V(2))
        push!(cooᵥ, v, v)
    end
    return
end

@inline function _primal_push_psd!((cooᵢ, cooⱼ, cooᵥ)::Tuple{FastVec{I},FastVec{I},FastVec{V}},
    it::Union{<:PSDIndextypeCOOVectorized{:L},<:PSDIndextypeCOOVectorized{:U}}, max_con::I, idx::I, diag::Bool, v::V) where {I,V}
    push!(cooᵢ, max_con)
    push!(cooⱼ, idx)
    push!(cooᵥ, diag ? v : v * it.invscaling)
    return
end

@inline function _primal_push_psd!((cooᵢ, cooₘₙ, cooᵥ)::Tuple{FastVec{I},Tuple{FastVec{I},FastVec{I}},FastVec{V}},
    it::PSDIndextypeMatrixCartesian{:F}, max_con::I, i::I, j::I, v::V) where {I,V}
    if i == j
        push!(cooᵢ, max_con)
        push!(cooⱼ[1], i)
        push!(cooⱼ[2], j)
        push!(cooᵥ, v)
    else
        push!(cooᵢ, max_con, max_con)
        push!(cooⱼ[1], i, j)
        push!(cooⱼ[2], j, i)
        v *= inv(V(2)) # hard-coded scaling
        push!(cooᵥ, v)
    end
    return
end

@inline function _primal_push_psd!((cooᵢ, cooₘₙ, cooᵥ)::Tuple{FastVec{I},Tuple{FastVec{I},FastVec{I}},FastVec{V}},
    it::PSDIndextypeMatrixCartesian, max_con::I, i::I, j::I, v::V) where {I,V}
    # the order is already correct
    push!(cooᵢ, max_con)
    push!(cooⱼ[1], i)
    push!(cooⱼ[2], j)
    push!(cooᵥ, i == j ? v : v * inv(V(2))) # hard-coded scaling
    return
end

#region Conversions for the objective
# from :F to x
@inline function _primal_push_psd_obj!((cooⱼ, cooᵥ)::Tuple{FastVec{I},FastVec{V}}, oldit::PSDIndextypeCOOVectorized{:F},
    newit::PSDIndextypeCOOVectorized{:F}, ::I, idx::I, idxᵀ::I, v::V) where {I,V}
    δ = _get_offset(newit) - _get_offset(oldit)
    if idx == idxᵀ
        push!(cooⱼ, idx + δ)
        push!(cooᵥ, v)
    else
        push!(cooⱼ, idx + δ, idxᵀ + δ)
        v *= inv(V(2))
        push!(cooᵥ, v, v)
    end
    return
end

@inline function _primal_push_psd_obj!((cooⱼ, cooᵥ)::Tuple{FastVec{I},FastVec{V}}, oldit::PSDIndextypeCOOVectorized{:F},
    newit::PSDIndextypeCOOVectorized{T}, dim::I, idx::I, idxᵀ::I, v::V) where {I,V,T}
    # col is the number of columns to our left
    if T === :L
        col = (idx - _get_offset(oldit)) ÷ dim
        # the kᵗʰ column has k entries in the strict disregarded triangle, so subtract ∑_{k = 1}^col k = col(col +1)/2
        newidx = idx - col * (col + one(I)) ÷ I(2)
    else
        @assert(T === :U)
        col = (idxᵀ - _get_offset(oldit)) ÷ dim
        # the kᵗʰ column has dim - k entries in the strict disregarded triangle, so subtract ∑ₖ (dim - k) = col(2dim - col -1)/2
        newidx = idxᵀ - col * (I(2) * dim - col - one(I)) ÷ I(2)
    end
    δ = _get_offset(newit) - _get_offset(oldit)
    push!(cooⱼ, newidx + δ)
    push!(cooᵥ, idx == idxᵀ ? v : v * newit.invscaling)
    return
end

@inline function _primal_push_psd_obj!((cooⱼ, cooᵥ)::Tuple{Tuple{FastVec{I},FastVec{I}},FastVec{V}},
    oldit::PSDIndextypeCOOVectorized{:F}, newit::PSDIndextypeMatrixCartesian{T}, dim::I, idx::I, idxᵀ::I, v::V) where {I,V,T}
    col, row = divrem(idx - _get_offset(oldit), dim)
    @assert(col ≤ row) # idx is in the lower triangle
    row += _get_offset(newit)
    col += _get_offset(newit)
    if T === :L || row == col
        push!(cooⱼ[1], row)
        push!(cooⱼ[2], col)
        push!(cooᵥ, row == col ? v : v * inv(V(2)))
    elseif T === :U
        push!(cooⱼ[1], col)
        push!(cooⱼ[2], row)
        push!(cooᵥ, v * inv(V(2)))
    else
        @assert(T === :F)
        push!(cooⱼ[1], row, col)
        push!(cooⱼ[2], col, row)
        v *= inv(V(2))
        push!(cooᵥ, v, v)
    end
    return
end

# from :U and :L to x
@inline function _primal_push_psd_obj!((cooⱼ, cooᵥ)::Tuple{FastVec{I},FastVec{V}}, oldit::PSDIndextypeCOOVectorized{T},
    newit::PSDIndextypeCOOVectorized{T}, ::I, idx::I, diag::Bool, v::V) where {I,V,T}
    δ = _get_offset(newit) - _get_offset(oldit)
    push!(cooⱼ, idx + δ)
    push!(cooᵥ, diag ? v : v * newit.invscaling)
    return
end

@inline function _primal_push_psd_obj!((cooⱼ, cooᵥ)::Tuple{FastVec{I},FastVec{V}}, oldit::PSDIndextypeCOOVectorized{T},
    newit::PSDIndextypeCOOVectorized{:F}, dim::I, idx::I, diag::Bool, v::V) where {I,V,T}
    # We find the column by iteration, not by a binary search - we expect the objective to have few enough terms that the
    # construction of the lookup would make up a significant time
    idx -= _get_offset(oldit)
    local col
    if T === :L
        for outer col in dim:-one(I):one(I)
            idx < col && break
            idx -= col
        end
        col = dim - col
        row = col + idx
    else
        @assert(T === :U)
        for outer col in one(I):dim
            idx < col && break
            idx -= col
        end
        row = idx
        col -= one(I)
    end
    # now row, col are zero-indexed
    if diag
        push!(cooⱼ, dim * col + row + _get_offset(newit))
        push!(cooᵥ, v)
    else
        push!(cooⱼ, dim * col + row + _get_offset(newit), dim * row + col + _get_offset(newit))
        v *= inv(V(2))
        push!(cooᵥ, v)
    end
end

@inline function _primal_push_psd_obj!((cooⱼ, cooᵥ)::Tuple{Tuple{FastVec{I},FastVec{I}},FastVec{V}},
    oldit::PSDIndextypeCOOVectorized{T}, newit::PSDIndextypeMatrixCartesian{Q}, dim::I, idx::I, diag::Bool, v::V) where {I,V,T,Q}
    # We find the column by iteration, not by a binary search - we expect the objective to have few enough terms that the
    # construction of the lookup would make up a significant time
    idx -= _get_offset(oldit)
    local col
    if T === :L
        for outer col in dim:-one(I):one(I)
            idx < col && break
            idx -= col
        end
        col = dim - col
        row = col + idx
    else
        @assert(T === :U)
        for outer col in one(I):dim
            idx < col && break
            idx -= col
        end
        row = idx
        col -= one(I)
    end
    # now row, col are zero-indexed
    row += _get_offset(newit)
    col += _get_offset(newit)
    @assert((row == col) == diag)
    if T === Q || diag
        push!(cooⱼ[1], row)
        push!(cooⱼ[2], col)
        push!(cooᵥ, diag ? v : v * inv(V(2)))
    elseif Q !== :F
        push!(cooⱼ[1], col)
        push!(cooⱼ[2], row)
        push!(cooᵥ, v * inv(V(2)))
    else
        push!(cooⱼ[1], row, col)
        push!(cooⱼ[2], col, row)
        v *= inv(V(2))
        push!(cooᵥ, v, v)
    end
    return
end

# from Cartesian to x
@inline function _primal_push_psd_obj!((cooⱼ, cooᵥ)::Tuple{Tuple{FastVec{I},FastVec{I}},FastVec{V}},
    oldit::PSDIndextypeMatrixCartesian{T}, newit::PSDIndextypeMatrixCartesian{Q}, ::I, i::I, j::I, v::V) where {I,V,T,Q}
    δ = _get_offset(newit) - _get_offset(oldit)
    if i != j
        v *= inv(V(2))
    end
    if Q === :F
        push!(cooⱼ[1], i + δ, j + δ)
        push!(cooⱼ[2], j + δ, i + δ)
        push!(cooᵥ, v, v)
    else
        push!(cooⱼ[1], (T === Q ? i : j) + δ)
        push!(cooⱼ[2], (T === Q ? j : i) + δ)
        push!(cooᵥ, v)
    end
    return
end

@inline function _primal_push_psd_obj!((cooⱼ, cooᵥ)::Tuple{FastVec{I},FastVec{V}}, oldit::PSDIndextypeMatrixCartesian{T},
    newit::PSDIndextypeCOOVectorized{:F}, dim::I, i::I, j::I, v::V) where {I,V,T}
    i -= _get_offset(oldit)
    j -= _get_offset(oldi)
    if i == j
        push!(cooⱼ, dim * j + i + _get_offset(newit))
        push!(cooᵥ, v)
    else
        push!(cooⱼ, dim * j + i + _get_offset(newit), dim * i + j + _get_offset(newit))
        v *= inv((V(2)))
        push!(cooᵥ, v)
    end
    return
end

@inline function _primal_push_psd_obj!((cooⱼ, cooᵥ)::Tuple{FastVec{I},FastVec{V}}, oldit::PSDIndextypeMatrixCartesian{T},
    newit::PSDIndextypeCOOVectorized{Q}, dim::I, i::I, j::I, v::V) where {I,V,T,Q}
    i -= _get_offset(oldit)
    j -= _get_offset(oldi)
    row, col = (T === Q || (T === :F && Q === :L)) ? (i, j) : (j, i)
    if Q === :L
        newidx = col * (I(2) * dim - col - one(I)) ÷ I(2) + row
    else
        newidx = col * (col + one(I)) ÷ I(2) + row
    end
    push!(cooⱼ, newidx + _get_offset(newit))
    push!(cooᵥ, i == j ? v : v * newit.invscaling)
    return
end
#endregion

@inline function _primal_add_constr!(state::PrimalMomentSolver{I,K,V,C,B}, indval::Indvals{K,V}, max_con::I,
    ::Val{checkfree}=Val(false)) where {I,K,V,checkfree,C,B}
    it = psd_indextype(state.parent)
    @inbounds for (k, v) in indval
        if checkfree && !haskey(state.mon_eq, FastKey(k))
            # We did not see this monomial before in a PSD/nonnegative constraint, but now we need it in a fixed constraint.
            # This can happen due to basis reduction and still lead to a solvable problem if there are multiple well-balanced
            # fixed constraints. We would need to create free variables, which we solvers in this form usually don't support;
            # so create two nonnegatives instead.
            refₘ = -2
            refᵢ = state.max_nonneg += I(2)
            refⱼ = zero(B) # just to have it defined
            state.mon_eq[FastKey(k)] = (refₘ, refᵢ, refⱼ)
        else
            refₘ, refᵢ, refⱼ = state.mon_eq[FastKey(k)]
        end
        if refₘ > 0 # likely
            _primal_push_psd!(state.psds[refₘ], it, max_con, refᵢ, refⱼ, v)
        else # quite unlikely (would require size-1 moment matrix or the case above)
            if checkfree && refₘ == -2 # -1 for nonnegative variables, -2 for difference of two nonnegative variables
                push!(state.nonnegs[1], max_con, max_con)
                push!(state.nonnegs[2], refᵢ - one(I), refᵢ)
                push!(state.nonnegs[3], v, -v)
            else
                push!(state.nonnegs[1], max_con)
                push!(state.nonnegs[2], refᵢ)
                push!(state.nonnegs[3], v)
            end
        end
    end
    return
end

function add_constr_nonnegative!(state::PrimalMomentSolver{I,K,V,<:Any,B}, indvals::IndvalsIterator{K,V}) where {I,K,V,B}
    cooᵢ, cooⱼ, cooᵥ = state.nonnegs
    mon_eq = state.mon_eq
    req_elems = length(indvals)
    prepare_push!(cooᵢ, req_elems)
    prepare_push!(cooⱼ, req_elems)
    prepare_push!(cooᵥ, req_elems)
    max_con, max_nonneg = state.max_con, state.max_nonneg
    @inbounds for indval in indvals
        max_nonneg += one(I)
        if isone(length(indval))
            k, v = first(indval)
            if !haskey(mon_eq, FastKey(k)) # unlikely
                @assert(isone(v))
                state.mon_eq[FastKey(k)] = (-1, max_nonneg, zero(B))
                continue
            end
        end
        max_con += one(I)
        req_elems -= 1

        # first this element
        unsafe_push!(cooᵢ, max_con)
        unsafe_push!(cooⱼ, max_nonneg)
        unsafe_push!(cooᵥ, -one(V))
        # then the others
        _primal_add_constr!(state, indval, max_con)
        # check if the unlikely case happened: then, the previous command will have led to pushing into our coos.
        # As a consequence, our buffer might not be large enough any more.
        prepare_push!(cooᵢ, req_elems)
        prepare_push!(cooⱼ, req_elems)
        prepare_push!(cooᵥ, req_elems)
    end
    state.max_con = max_con
    state.max_nonneg = max_nonneg
    return
end

function add_constr_psd!(state::PrimalMomentSolver{I,K,V,C}, dim::Int, data::IndvalsIterator{K,V}) where {I,K,V,C}
    it = psd_indextype(state.parent)
    full = it isa PSDIndextypePrimal{:F}
    push!(state.psd_dim, dim)
    psd_index = Int(length(state.psd_dim))
    mon_eq = state.mon_eq
    cooᵢ = FastVec{I}(buffer=full ? dim^2 + 2length(rowvals(data)) : trisize(dim) + length(rowvals(data)))
    cooⱼ = C <: Tuple ? (similar(cooᵢ), similar(cooᵢ)) : similar(cooᵢ)
    cooᵥ = similar(cooᵢ, V)
    push!(state.psds, (cooᵢ, cooⱼ, cooᵥ))
    offset = C <: Tuple ? _get_offset(it) : one(I)
    col = offset
    row = offset - one(offset)
    max_con = state.max_con
    scaling = full || it isa PSDIndextypeMatrixCartesian ? -inv(V(2)) : -it.invscaling
    if !(C <: Tuple)
        j = _get_offset(it) - one(I)
        if full
            jtrans = _get_offset(it) - dim
        end
    end
    @inbounds for indval in data
        if it isa PSDIndextypePrimal{:U}
            if (row += one(row)) > col
                col += one(col)
                row = offset
            end
        else
            if (row += one(row)) ≥ dim + offset
                if full && !(C <: Tuple)
                    j += col # + 1 - offset, but if we are here, offset = 1
                    jtrans = j - dim + one(I)
                end
                col += one(col)
                row = col
            end
        end
        if C <: Tuple
            _i, _j = row, col
        else
            _i = (j += one(I))
            if full
                _j = (jtrans += dim)
            else
                _j = row == col
            end
        end

        if isone(length(indval))
            k, v = first(indval)
            if !haskey(mon_eq, FastKey(k))
                @assert(isone(v))
                state.mon_eq[FastKey(k)] = (psd_index, _i, _j)
                continue
            end
        end
        max_con += one(I)

        # first this element
        _primal_push_psd!((cooᵢ, cooⱼ, cooᵥ), it, max_con, _i, _j, -one(V))
        # then the others
        _primal_add_constr!(state, indval, max_con)
    end
    state.max_con = max_con
    return
end

function add_constr_fix!(state::PrimalMomentSolver{I,K,V}, ::Nothing, indvals::Indvals{K,V}, rhs::V) where {I,K,V}
    max_con = (state.max_con += one(I))
    _primal_add_constr!(state, indvals, max_con, Val(true))
    @inbounds if !iszero(rhs)
        push!(state.b[1], max_con)
        push!(state.b[2], rhs)
    end
    return
end

function fix_objective!(state::PrimalMomentSolver{<:Integer,K,V}, indvals::Indvals{K,V}) where {K,V}
    state.c = indvals
    return
end

_get_SI(::AbstractSparseMatrixSolver{I}) where {I<:Integer} = I
_get_SI(::AbstractAPISolver{<:Integer,T}) where {T} = T
_get_SI(::AnySolver) = Int

"""
    primal_moment_setup!(state::AbstractSolver, relaxation::AbstractRelaxation,
        groupings::RelaxationGroupings; verbose=false)

Sets up all the necessary moment matrices, constraints, and objective of a polynomial optimization problem in primal form,
i.e., for solvers which allow to declare monolithic semidefinite and nonnegative variables and linear constraints. While
usually, the SOS form would be more suitable for such solvers, forcing the moment form in primal representation ensures that
low-rank assumptions about the primal variable hold true, which can be exploited by some solvers.
This function returns a `Vector{<:Vector{<:Tuple{Symbol,Any}}}` that contains internal information on the problem. This
information is required to obtain dual variables and re-optimize the problem and should be stored in the `state`. It also
returns an internal state that is important for the reconstruction of the moment matrices using
[`MomentVector`](@ref MomentVector(::AbstractRelaxation, ::PrimalMomentSolver, ::Any, ::Any)).
The internal state has three properties that are part of the public interface:
- `num_con::Int` is the number of constraints
- `num_nonneg::Int` is the number of nonnegative variables
- `psd_dim::FastVec{I}` holds the side dimensions of the PSD variables (where `I` is the index type of `state`)

The following methods must be implemented by a solver to make this function work:
- [`mindex`](@ref)
- [`add_var_nonnegative!`](@ref add_var_nonnegative!(::AbstractSolver{<:Integer,V}, ::Int, ::Int, ::SparseMatrixCOO{I,I,V}, ::Tuple{FastVec{I},FastVec{V}}) where {I,V}),
  which is called no more than once
- [`add_var_psd!`](@ref add_var_psd!(::AbstractSolver{<:Integer,V}, ::Int, ::I, ::SparseMatrixCOO{I,I,V}, ::Union{Nothing,Tuple{FastVec{I},FastVec{V}}}) where {I,V})
- [`psd_indextype`](@ref)
- [`objective_indextype`](@ref)
- [`fix_constraints!`](@ref fix_constraints!(::AbstractSolver{<:Integer,V}, ::Int, ::Indvals{<:Integer,V}) where {V})

!!! warning
    During the reformulation, the function is able to detect a certain class of unbounded or infeasible problem formulations.
    If this is the case, it will return `missing` without invoking the solver.
"""
function primal_moment_setup!(state::AbstractSolver{K,V}, relaxation::AbstractRelaxation, groupings::RelaxationGroupings;
    verbose::Bool=false) where {K,V}
    objective_it = objective_indextype(state)
    I = _get_SI(state)
    it = psd_indextype(state)
    obj_vect = objective_it isa PSDIndextypeCOOVectorized
    offset = _get_offset(it)
    obj_offset = _get_offset(objective_it)

    conversion_time = @elapsed begin
        pstate = PrimalMomentSolver{I,K,V}(state)
        info = moment_setup!(pstate, relaxation, groupings; representation=RepresentationPSD())

        # Now we must call our solver with the data; but we also need the objective for this. So lets first convert it.
        C = Vector{Union{Nothing,Tuple{obj_vect ? FastVec{I} : Tuple{FastVec{I},FastVec{I}},FastVec{V}}}}(
            undef, length(pstate.psds)
        )
        c_lin_data = (FastVec{I}(), FastVec{V}())
        try
            @inbounds for (k, v) in pstate.c
                refₘ, refᵢ, refⱼ = pstate.mon_eq[FastKey(k)]
                if refₘ > 0
                    dim = pstate.psd_dim[refₘ]
                    if isassigned(C, refₘ)
                        spₘ = C[refₘ]
                    else
                        C[refₘ] = spₘ = (obj_vect ? FastVec{I}() : (FastVec{I}(), FastVec{I}()), FastVec{V}())
                    end
                    # In principle, even the triangle could change for the objective...
                    _primal_push_psd_obj!(spₘ, it, objective_it, dim, refᵢ, refⱼ, v)
                elseif refₘ == -1
                    push!(c_lin_data[1], refᵢ + (obj_offset - offset))
                    push!(c_lin_data[2], v)
                else
                    @assert(refₘ == -2)
                    refᵢ += obj_offset - offset
                    push!(c_lin_data[1], refᵢ - one(I), refᵢ)
                    push!(c_lin_data[2], v, -v)
                end
            end
        catch e
            if e isa KeyError
                # This can happen if the objective contains monomials that were not present before due to a smaller basis. In
                # principle, we should add a free variable for each of those monomials; however, these variables will never
                # occur anywhere else. Therefore, the problem is naturally unbounded (unless it is infeasible, which we simply
                # disregard here) and we don't need to solve it.
                @verbose_info("Detected unbounded objective during problem construction; skipping solver")
                return missing
            else
                rethrow(e)
            end
        end
    end

    @verbose_info("Conversion into primal format complete in ", conversion_time, " seconds, now setting up solver")

    num_con = pstate.num_con
    if pstate.max_nonneg ≥ offset
        nonneg_arg = (num_con, Int(pstate.max_nonneg - offset) +1, SparseMatrixCOO(pstate.nonnegs..., offset))
        add_var_nonnegative!(state, nonneg_arg..., c_lin_data)
    end

    @inbounds for (i, (dim, psd)) in enumerate(zip(pstate.psd_dim, pstate.psds))
        psd_arg = it isa PSDIndextypeCOOVectorized ? SparseMatrixCOO(psd[1], psd[2], psd[3], offset) : psd
        add_var_psd!(state, num_con, dim, psd_arg, isassigned(C, i) ? C[i] : nothing)
    end

    fix_constraints!(state, num_con, Indvals(pstate.b...))

    return info, pstate
end

"""
    MomentVector(relaxation::AbstractRelaxation, primal_moment_state, mm, mm_lin)

Given the moment data obtained from an optimization using [`primal_moment_setup!`](@ref), convert it to a
[`MomentVector`](@ref). Here, `primal_moment_state` must be the internal state that is returned as the second output of
[`primal_moment_setup!`](@ref). `mm` must be a vector of indexable objects (according to the indexing scheme specified by
[`psd_indextype`](@ref)), and `mm_lin` must be indexable; they should return the value of the primal PSD and nonnegative
variables. Note that this function is not fully type-stable, as the result may be based either on a dense or sparse vector
depending on the relaxation.
"""
function MomentVector(relaxation::AbstractRelaxation{<:Problem{<:SimplePolynomial{<:Any,Nr,Nc}}},
    state::PrimalMomentSolver{I,K,V}, mm, mm_lin) where {Nr,Nc,I<:Integer,K<:Integer,V<:Real}
    mon_eq = state.mon_eq
    y = Vector{V}(undef, length(mon_eq))
    mon_pos = convert.(K, keys(mon_eq))
    it = psd_indextype(state.parent)
    if it isa PSDIndextypeMatrixCartesian
        @inbounds for (j, (psdᵢ, row, col)) in enumerate(values(mon_eq))
            if psdᵢ > 0
                y[j] = mm[psdᵢ][row, col]
            elseif psdᵢ == -1
                y[i] = mm_lin[row]
            else
                y[i] = mm_lin[row-one(I)] - mm_lin[row]
            end
        end
    else
        @inbounds for (j, (psdᵢ, idx, _)) in enumerate(values(mon_eq))
            if psdᵢ > 0
                y[j] = mm[psdᵢ][idx]
            elseif psdᵢ == -1
                y[i] = mm_lin[idx]
            else
                y[i] = mm_lin[idx-one(I)] - mm_lin[idx]
            end
        end
    end
    sort_along!(mon_pos, y)
    max_mons = mon_pos[end]
    if length(y) == max_mons
        solution = y
    elseif 3length(y) < max_mons
        solution = SparseVector(max_mons, mon_pos, y)
    else
        solution = fill(NaN, max_mons)
        copyto!(@view(solution[mon_pos]), y)
    end
    return MomentVector(relaxation, ExponentsAll{Nr+2Nc,K}(), solution)
end

"""
    coo_to_csc!(nCols, s::SparseMatrixCOO)

Converts a COO matrix into a CSC matrix, where the three vectors (colptr, rowvals, nzvals) are returned. The offset of `s` is
respected. Note that some of the vectors in `s` are modified by this function.
"""
function coo_to_csc!(nCols::Integer, s::SparseMatrixCOO{I,I,V,offset}) where {I,V,offset}
    i, j, v = s.rows, s.cols, s.vals
    sort_along!(j, i, v, relevant=2) # let's put it into normal form with ordered rowvals also
    nCols > j[end] - offset || throw(ArgumentError("Too few columns"))
    j[begin] ≥ offset || throw(ArgumentError("Invalid offset"))

    colptr = Vector{I}(undef, nCols +1)
    colptridx = 1
    jidx = iterate(j)
    pos = offset
    curcol = offset
    @inbounds while colptridx ≤ length(colptr)
        if isnothing(jidx)
            # no column follow, fill up with whatever remains
            @label done
            fill!(@view(colptr[colptridx:end]), pos)
            return colptr, finish!(i), finish!(v)
        else
            fill!(@view(colptr[colptridx:colptridx+(jidx[1]-curcol)]), pos)
            colptridx += jidx[1] - curcol +1
            curcol = jidx[1]
            while jidx[1] == curcol
                pos += one(pos)
                jidx = iterate(j, jidx[2])
                isnothing(jidx) && @goto done
            end
            curcol += one(curcol)
        end
    end
    @assert(false)
end

"""
    coo_to_csr!(nRows, s::SparseMatrixCOO)

Converts a COO matrix into a CSR matrix, where the three vectors (rowptr, colvals, nzvals) are returned. The offset of `s` is
respected. Note that some of the vectors in `s` are modified by this function.
"""
coo_to_csr!(nRows::Int, s::SparseMatrixCOO{I,I,<:Real,offset}) where {I,offset} =
    coo_to_csc!(nRows, SparseMatrixCOO(s.cols, s.rows, s.vals, offset))