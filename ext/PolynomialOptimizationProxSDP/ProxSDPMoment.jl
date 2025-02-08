mutable struct StateMoment{K<:Integer} <: AbstractSparseMatrixSolver{Int64,K,Float64}
    const Acsc::Tuple{FastVec{Int64},FastVec{Int64},FastVec{Float64}}
    const cones::FastVec{ProxSDP.SDPSet}
    const c::Tuple{FastVec{Int64},FastVec{Float64}}
    b::Vector{Float64}
    info::Vector{<:Vector{<:Tuple{Symbol,Any}}}
    data

    StateMoment{K}() where {K<:Integer} = new{K}(
        (FastVec{Int64}(), FastVec{Int64}(), FastVec{Float64}()),
        FastVec{ProxSDP.SDPSet}(),
        (FastVec{Int64}(), FastVec{Float64}())
    )
end

Solver.issuccess(::Val{:ProxSDPMoment}, status::Int) = status == 1

Solver.psd_indextype(::StateMoment) = PSDIndextypeCOOVectorized(:U, true, 1)

@counter_atomic(StateMoment, :psd)

function Solver.add_var_nonnegative!(state::StateMoment, m::Int, n::Int, data::SparseMatrixCOO{Int64,Int64,Float64,one(Int64)},
    obj::Tuple{FastVec{Int64},FastVec{Float64}})
    lastprevvar = length(state.Acsc[1]) # index of last preceding variable before this cone starts
                                        # (should be 0, as nonnegs are the first)
    csc = coo_to_csc!(n, data)
    if !isempty(state.Acsc[2])
        csc[1] .+= length(state.Acsc[2])
    end
    append!.(state.Acsc, csc)
    Base._deleteend!(state.Acsc[1], 1) # the last item will only be added at the very end
    if !isnothing(obj)
        append!(state.c[1], Iterators.map(x -> x + lastprevvar, obj[1]))
        append!(state.c[2], obj[2])
    end
    return
end

function Solver.add_var_psd!(state::StateMoment, m::Int, dim::Int, data::SparseMatrixCOO{Int64,Int64,Float64,one(Int64)},
    obj::Union{Nothing,Tuple{FastVec{Int64},FastVec{Float64}}})
    ts = trisize(dim)
    mat_i = Vector{Int}(undef, ts) # will hold the position in the _lower_ triangle when matrix is seen as full vector
                                   # (even though data is generated in upper triangular form)
    i = 1
    mat_idx = 1
    @inbounds for col in 1:dim
        @simd for row in 0:dim-col
            mat_i[i+row] = mat_idx + row
        end
        i += dim - col +1
        mat_idx += dim +1
    end
    lastprevvar = length(state.Acsc[1]) # index of last preceding variable before this cone starts
    push!(state.cones, ProxSDP.SDPSet(collect(lastprevvar+1:lastprevvar+ts), mat_i, ts, dim))

    csc = coo_to_csc!(ts, data)
    if !isempty(state.Acsc[2])
        csc[1] .+= length(state.Acsc[2])
    end
    append!.(state.Acsc, csc)
    Base._deleteend!(state.Acsc[1], 1) # the last item will only be added at the very end

    if !isnothing(obj)
        append!(state.c[1], Iterators.map(x -> x + lastprevvar, obj[1]))
        append!(state.c[2], obj[2])
    end
    return
end

function Solver.fix_constraints!(state::StateMoment, m::Int, indvals::Indvals{Int64,Float64})
    state.b = b = zeros(Float64, m)
    for (i, v) in indvals
        b[i] = v
    end
    return
end

function Solver.poly_optimize(::Val{:ProxSDPMoment}, relaxation::AbstractRelaxation, groupings::RelaxationGroupings;
    representation, verbose::Bool=false, customize=(state) -> nothing, parameters...)
    setup_time = @elapsed @inbounds begin
        K = _get_I(eltype(monomials(poly_problem(relaxation).objective)))

        if haskey(parameters, :log_verbose)
            options = ProxSDP.Options(; parameters...)
        else
            options = ProxSDP.Options(; log_verbose=verbose, parameters...)
        end
        state = StateMoment{K}()
        primal_data = primal_moment_setup!(state, relaxation, groupings; verbose)
        ismissing(primal_data) && return missing, 4, typemin(V)
        push!(state.Acsc[1], length(state.Acsc[2]) +1) # complete the CSC form
        state.info, state.data = primal_data
        customize(state)

        num_con, num_var, num_nonneg = primal_data[2].num_con, length(state.Acsc[1]) -1, primal_data[2].num_nonneg
        A = SparseMatrixCSC{Float64,Int64}(num_con, num_var, finish!.(state.Acsc)...)
        c = zeros(Float64, num_var)
        @inbounds for (i, v) in zip(state.c...)
            c[i] = v
        end
        if iszero(num_nonneg)
            G = spzeros(Float64, Int64, (0, num_var))
        else
            G = let colptr=Vector{Int64}(undef, num_var +1)
                copyto!(@view(colptr[1:num_nonneg]), one(Int64):Int64(num_nonneg))
                fill!(@view(colptr[num_nonneg+1:end]), Int64(num_nonneg +1))
                SparseMatrixCSC{Float64,Int64}(num_nonneg, num_var, colptr, collect(one(Int64):Int64(num_nonneg)),
                    fill(-1., num_nonneg))
            end
        end
        h = zeros(Float64, size(G, 1))

        aff = ProxSDP.AffineSets(num_var, num_con, length(h), 0, A, G, state.b, h, c)
        con = ProxSDP.ConicSets(finish!(state.cones), ProxSDP.SOCSet[])
    end
    @verbose_info("Setup complete in ", setup_time, " seconds")
    solution = ProxSDP.chambolle_pock(aff, con, options)
    status = solution.status
    value = solution.objval
    @verbose_info("Optimization complete")

    return (state, solution), status, value
end

Solver.extract_moments(relaxation::AbstractRelaxation, (state, solution)::Tuple{StateMoment,ProxSDP.Result}) =
    MomentVector(relaxation, state.data, [@view(solution.primal[i.vec_i[begin]:i.vec_i[end]]) for i in state.cones],
        @view(solution.primal[1:state.data.num_nonneg]))

Solver.extract_sos(::AbstractRelaxation, (_, solution)::Tuple{StateMoment,ProxSDP.Result}, ::Val{:nonnegative},
    index::AbstractUnitRange, ::Nothing) = @view(solution.dual_in[index])

Solver.extract_sos(::AbstractRelaxation, (state, solution)::Tuple{StateMoment,ProxSDP.Result}, ::Val{:psd},
    index::Integer, ::Nothing) = SPMatrix(state.cones[index].sq_side,
        @view(solution.dual_cone[state.cones[index].vec_i[begin]:state.cones[index].vec_i[end]]), :U)