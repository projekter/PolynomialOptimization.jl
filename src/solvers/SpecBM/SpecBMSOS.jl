mutable struct StateSOS{K<:Integer,V<:Real} <: AbstractSparseMatrixSolver{Int,K,V}
    const Aᵀcoo::SparseMatrixCOO{Int,K,V,1}
    const c::Tuple{FastVec{Int},FastVec{V}}
    num_frees::Int
    slack::K
    const psds::FastVec{Int}
    b::Tuple{Vector{K},Vector{V}}
    info::Vector{<:Vector{<:Tuple{Symbol,Any}}}

    StateSOS{K,V}() where {K<:Integer,V<:Real} = new{K,V}(
        SparseMatrixCOO{Int,K,V,1}(),
        (FastVec{Int}(), FastVec{V}()),
        0, K <: Signed ? -one(K) : typemax(K), FastVec{Int}()
    )
end

Solver.issuccess(::Val{:SpecBMSOS}, status::Symbol) = status === :Optimal

Solver.psd_indextype(::StateSOS) = PSDIndextypeVector(:L)

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
    prepare_push!(state.Aᵀcoo, num)
    return
end

function Solver.add_var_free!(state::StateSOS{K,V}, ::Nothing, indvals::Indvals{K,V}, obj::V) where {K<:Integer,V<:Real}
    v = append!(state.Aᵀcoo, indvals)
    if !iszero(obj)
        push!(state.c[1], v)
        push!(state.c[2], -obj)
    end
    state.num_frees += 1
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

        # Now we have all the transposed data in COO form. The reason for this choice is that we were able to assign arbitrary
        # row indices - i.e., we could just use the monomial index. However, now we have to modify the row indices to make them
        # consecutive, removing all monomials that do not occur. We already know that no entry will ever occur twice, so we can
        # make our own optimized COO -> CSC function.
        c = zeros(V, size(state.Aᵀcoo, 1))
        copyto!(@view(c[state.c[1]]), state.c[2])

        moncount, (Aᵀcolptr, Aᵀrowind, Aᵀnzval), b = coo_to_csc!(state.Aᵀcoo, state.b)
        A = transpose(SparseMatrixCSC{V,Int}(length(c), moncount, Aᵀcolptr, Aᵀrowind, Aᵀnzval))

        # Figure out some proper default arguments. Let us rely on the adaptive strategies in the solver to do this instead of
        # some sophisticated analysis that might give us better values provided the problem features certain constraints.
        args = merge!(Dict{Symbol,Any}(:ρ => V(10), :adaptiveρ => !haskey(parameters, :ρ), :r_past => 0, :r_current => 3,
            :ϵ => V(1//100_000)), parameters)
    end
    @verbose_info("Setup complete in ", setup_time, " seconds")

    result = specbm_primal(A, b, c; state.num_frees, psds=finish!(state.psds), verbose, args...)
    status = result.status
    value = -result.objective
    @verbose_info("Optimization complete, retrieving moments")

    rmul!(result.y, -one(V))
    return (state, result), status, value
end

Solver.extract_moments(relaxation::AbstractRelaxation, (state, result)::Tuple{StateSOS,Any}) =
    MomentVector(relaxation, result.y, state.Aᵀcoo)

Solver.extract_sos(relaxation::AbstractRelaxation, (state, result)::Tuple{StateSOS,Any}, type::Val,
    index::AbstractUnitRange, ::Nothing) = @view(result.x[index])

Solver.psd_indextype(::Tuple{StateSOS,Any}) = PSDIndextypeVector(:L)