mutable struct StateSOS{K<:Integer,V<:Real} <: AbstractSparseMatrixSolver{Int,K,V}
    Aᵀcoo::SparseMatrixCOO{Int,K,V,1} # during buildup: all but the free ones; afterwards: all
    const Aᵀcoo_free::SparseMatrixCOO{Int,K,V,1}
    const c::Tuple{FastVec{Int},FastVec{V}}
    num_frees::Int # we don't need it, as size(Aᵀcoo_free, 1) contains the same information, but after merging it is useful
    slack::K
    const psds::FastVec{Int}
    b::Tuple{Vector{K},Vector{V}}
    info::Vector{<:Vector{<:Tuple{Symbol,Any}}}

    StateSOS{K,V}() where {K<:Integer,V<:Real} = new{K,V}(
        SparseMatrixCOO{Int,K,V,1}(), SparseMatrixCOO{Int,K,V,1}(),
        (FastVec{Int}(), FastVec{V}()),
        0, K <: Signed ? -one(K) : typemax(K), FastVec{Int}()
    )
end

Solver.issuccess(::Val{:SpecBMSOS}, status::Symbol) = status === :Optimal

Solver.psd_indextype(::StateSOS{<:Integer,V}) where {V} = PSDIndextypeVector(:L, sqrt(V(2)))

@counter_alias(StateSOS, :nonnegative, :psd)

function Solver.add_var_nonnegative!(state::StateSOS{K,V}, indvals::IndvalsIterator{K,V}) where {K<:Integer,V<:Real}
    append!(state.Aᵀcoo, indvals)
    append!(state.psds, Iterators.repeated(1, length(indvals)))
    return
end

function Solver.add_var_psd!(state::StateSOS{K,V}, dim::Int, data::IndvalsIterator{K,V}) where {K<:Integer,V<:Real}
    append!(state.Aᵀcoo, data)
    push!(state.psds, dim)
    return
end

function Solver.add_var_free_prepare!(state::StateSOS, num::Int)
    # Those is just a lower bound, as the number of constraints is multiplied by the individual index count. But better than
    # nothing.
    prepare_push!(state.Aᵀcoo_free, num)
    return
end

function Solver.add_var_free!(state::StateSOS{K,V}, ::Nothing, indvals::Indvals{K,V}, obj::V) where {K<:Integer,V<:Real}
    v = append!(state.Aᵀcoo_free, indvals)
    if !iszero(obj)
        push!(state.c[1], v)
        push!(state.c[2], -obj)
    end
    return
end

function Solver.fix_constraints!(state::StateSOS{K,V}, indvals::Indvals{K,V}) where {K<:Integer,V<:Real}
    state.b = (indvals.indices, indvals.values)
    return
end

function Solver.poly_optimize(::Val{:SpecBMSOS}, relaxation::AbstractRelaxation, groupings::RelaxationGroupings;
    representation, verbose::Bool=false, customize=_ -> nothing, parameters...)
    setup_time = @elapsed begin
        K = _get_I(eltype(monomials(poly_problem(relaxation).objective)))
        V = real(coefficient_type(poly_problem(relaxation).objective))
        state = StateSOS{K,V}()

        state.info = sos_setup!(state, relaxation, groupings; representation)
        customize(state)

        state.num_frees = num_frees = size(state.Aᵀcoo_free, 1)
        if isempty(state.Aᵀcoo_free.rowinds)
            n = size(state.Aᵀcoo, 1)
        elseif length(state.Aᵀcoo_free) ≤ length(state.Aᵀcoo)
            state.Aᵀcoo.rowinds .+= num_frees
            n = size(state.Aᵀcoo, 1)
            # will be sorted according to column anyway, so we don't care
            prepare_push!(state.Aᵀcoo, length(state.Aᵀcoo_free))
            unsafe_append!(state.Aᵀcoo.rowinds, state.Aᵀcoo_free.rowinds)
            unsafe_append!(state.Aᵀcoo.moninds, state.Aᵀcoo_free.moninds)
            unsafe_append!(state.Aᵀcoo.nzvals, state.Aᵀcoo_free.nzvals)
        else
            prepare_push!(state.Aᵀcoo_free, length(state.Aᵀcoo))
            @inbounds for x in state.Aᵀcoo.rowinds
                unsafe_push!(state.Aᵀcoo_free.rowinds, x + num_frees)
            end
            unsafe_append!(state.Aᵀcoo_free.moninds, state.Aᵀcoo.moninds)
            unsafe_append!(state.Aᵀcoo_free.nzvals, state.Aᵀcoo.nzvals)
            state.Aᵀcoo = state.Aᵀcoo_free
            n = size(state.Aᵀcoo, 1)
        end

        c = zeros(V, n)
        copyto!(@view(c[state.c[1]]), state.c[2])

        moncount, (Aᵀcolptr, Aᵀrowind, Aᵀnzval), b = coo_to_csc!(state.Aᵀcoo, state.b)
        A = transpose(SparseMatrixCSC{V,Int}(length(c), moncount, Aᵀcolptr, Aᵀrowind, Aᵀnzval))

        # Figure out some proper default arguments. Let us rely on the adaptive strategies in the solver to do this instead of
        # some sophisticated analysis that might give us better values provided the problem features certain constraints.
        args = merge!(Dict{Symbol,Any}(:ρ => V(10), :adaptiveρ => !haskey(parameters, :ρ), :r_past => 0, :r_current => 3,
            :ϵ => V(1//100_000)), parameters)
    end
    @verbose_info("Setup complete in ", setup_time, " seconds")

    result = specbm_primal(A, b, c; num_frees, psds=finish!(state.psds), verbose, args...)
    status = result.status
    value = -result.objective
    @verbose_info("Optimization complete, retrieving moments")

    rmul!(result.y, -one(V))
    return (state, result), status, value
end

Solver.extract_moments(relaxation::AbstractRelaxation, (state, result)::Tuple{StateSOS,Any}) =
    MomentVector(relaxation, result.y, state.slack, state.Aᵀcoo)

Solver.extract_sos(relaxation::AbstractRelaxation, (state, result)::Tuple{StateSOS,Any}, type::Val,
    index::AbstractUnitRange, ::Nothing) = @view(result.x[index])

Solver.psd_indextype(::Tuple{StateSOS{<:Integer,V},Any}) where {V} = PSDIndextypeVector(:L, sqrt(V(2)))