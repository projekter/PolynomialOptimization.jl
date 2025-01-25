mutable struct StateMoment{K<:Integer} <: AbstractSparseMatrixSolver{Cint,K,Cdouble}
    # while this is in the end an AbstractAPISolver, the add_ functions are matrix-based
    const A::Vector{Tuple{Vector{Cint},Vector{Cint},Vector{Cdouble}}}
    A_lin::Tuple{Vector{Cint},Vector{Cint},Vector{Cdouble}}
    b::Vector{Cdouble}
    info::Vector{<:Vector{<:Tuple{Symbol,Any}}}
    data

    StateMoment{K}() where {K<:Integer} = new{K}(
        Tuple{Vector{Cint},Vector{Cint},Vector{Cdouble}}[]
    )
end

Solver.issuccess(::Val{:LoRADSMoment}, status::LoRADS.Status) =
    status ∈ (LoRADS.ASDP_PRIMAL_DUAL_OPTIMAL, LoRADS.ASDP_PRIMAL_OPTIMAL)

Solver.psd_indextype(::StateMoment) = PSDIndextypeCOOVectorized(:L, 2., zero(Cint))

@counter_atomic(StateMoment, :psd)

function Solver.add_var_nonnegative!(state::StateMoment, m::Int, n::Int, data::SparseMatrixCOO{Cint,Cint,Cdouble,zero(Cint)},
    obj::Tuple{FastVec{Cint},FastVec{Cdouble}})
    @assert(!isdefined(state, :A_lin))
    data.rows .+= one(Cint)
    append!(data.rows, Iterators.repeated(zero(Cint), length(obj[1])))
    append!(data.cols, obj[1])
    append!(data.vals, obj[2])
    state.A_lin = coo_to_csr!(m +1, data)
    return
end

function Solver.add_var_psd!(state::StateMoment, m::Int, dim::Cint, data::SparseMatrixCOO{Cint,Cint,Cdouble,zero(Cint)},
    obj::Union{Nothing,Tuple{FastVec{Cint},FastVec{Cdouble}}})
    data.rows .+= one(Cint)
    if !isnothing(obj)
        append!(data.rows, Iterators.repeated(zero(Cint), length(obj[1])))
        append!(data.cols, obj[1])
        append!(data.vals, obj[2])
    end
    push!(state.A, coo_to_csr!(m +1, data))
    return
end

function Solver.fix_constraints!(state::StateMoment, m::Int, indvals::Indvals{Cint,Cdouble})
    state.b = b = zeros(Cdouble, m)
    for (i, v) in indvals
        b[i+1] = v
    end
    return
end

function Solver.poly_optimize(::Val{:LoRADSMoment}, relaxation::AbstractRelaxation, groupings::RelaxationGroupings;
    representation, verbose::Bool=false, customize=(state) -> nothing,
    timesLogRank::Real=2., phase1Tol::Real=1e-3, timeLimit::Real=10_000.,
    rhoMax::Real=5000., rhoFreq::Integer=5, rhoFactor::Real=1.2, rho::Real=0., maxIter::Integer=10_000,
    strategy::LoRADS.Strategy=LoRADS.STRATEGY_DEFAULT, admmStrategy::LoRADS.Strategy=LoRADS.STRATEGY_MIN_BISECTION,
    tau::Real=0., gamma::Real=0., heuristicFactor::Real=1., rankFactor::Real=1.5)
    setup_time = @elapsed @inbounds begin
        K = _get_I(eltype(monomials(poly_problem(relaxation).objective)))

        state = StateMoment{K}()
        primal_data = primal_moment_setup!(state, relaxation, groupings; verbose)
        ismissing(primal_data) && return missing, LoRADS.ASDP_INFEAS_OR_UNBOUNDED, -Inf
        isempty(primal_data[2].psd_dim) && error("LoRADS requires at least one semidefinite variable")
        state.info, state.data = primal_data
        customize(state)

        @verbose_info("Initializing solver")
        solver = LoRADS.ASDP()
        num_con, num_nonneg, psd_dim = primal_data[2].num_con, primal_data[2].num_nonneg, primal_data[2].psd_dim
        LoRADS.init_solver(solver, num_con, psd_dim, num_nonneg; rho, rhoMax, maxIter, strategy)
        LoRADS.set_dual_objective(solver, state.b)
        if isdefined(state, :A_lin)
            LoRADS.set_lp_cone(solver, num_con, num_nonneg, state.A_lin...) # this will copy all the A data
            empty!.(state.A_lin)
        end
        for (i, dim, data) in zip(Iterators.countfrom(zero(Cint)), psd_dim, state.A)
            LoRADS.set_cone(solver, i, LoRADS.conedata_to_userdata(LoRADS.ASDP_CONETYPE_DENSE_SDP, num_con, dim, data...))
        end
        LoRADS.preprocess(solver, psd_dim)
        empty!(state.A) # the data was transferred to an internal buffer in preprocess
        LoRADS.determine_rank(solver, psd_dim, timesLogRank, 0)
        LoRADS.detect_sparsity_sdp_coeff(solver)
        LoRADS.init_bm_vars(solver, psd_dim, num_nonneg)
        LoRADS.init_admm_vars(solver, psd_dim, num_nonneg)
        LoRADS.scale(solver)
        pre_mainiter = Ref(zero(Cint))
        pre_miniter = Ref(zero(Cint))
    end
    @verbose_info("Setup complete in ", setup_time, " seconds")
    solver_time = @elapsed begin
        time = LoRADS.get_timestamp()
        is_rank_max = 1e-10 > timesLogRank
        local bm_ret
        while true
            bm_ret = LoRADS.bm_optimize(solver, phase1Tol, -.001, 1e-16, 1e-10, time, is_rank_max, pre_mainiter, pre_miniter,
                timeLimit)
            bm_ret == LoRADS.RETCODE_RANK || break
            if !LoRADS.check_all_rank_max(solver, rankFactor)
                if (is_rank_max = LoRADS.aug_rank(solver, state.psd_dim, rankFactor))
                    @verbose_info("Restarting BM with maximum rank")
                else
                    timesLogRank *= rankFactor
                    @verbose_info("Restarting BM with updated rank (now ", timesLogRank, " logm)")
                end
            end
        end
        if bm_ret != LoRADS.RETCODE_EXIT
            bta_time = @elapsed LoRADS.bm_to_admm(solver, rhoMax)
            @verbose_info("Converted BM to ADMM in ", bta_time, " seconds")
            LoRADS.optimize(solver, rhoFreq, rhoFactor, admmStrategy, tau, gamma, 0., time, timeLimit)
        end
        LoRADS.end_program(solver) # this will print a lot of stuff, but it will also be responsible for calculating the errors
                                   # properly, so we do this even if verbose is not set
    end

    @verbose_info("Optimization complete in ", solver_time, " seconds")

    return (state, solver), solver.AStatus, solver.pObjVal
end

struct GetX
    solver::LoRADS.ASDP
    mm::Vector{Matrix{Cdouble}}
    trinumbers::Vector{Vector{Int}}

    function GetX(solver::LoRADS.ASDP, dims::AbstractVector{I}) where {I<:Integer}
        mm = Vector{Matrix{Cdouble}}(undef, length(dims))
        trinumbers = similar(mm, Vector{Int})
        trinumbers_lookup = Dict{I,Vector{Int}}()
        @inbounds for (i, dim) in enumerate(dims)
            trinumbers[i] = get!(trinumbers_lookup, dim) do
                tn = Vector{Int}(undef, dim)
                tn[1] = 0
                rest = dim
                for i in 2:dim
                    tn[i] = tn[i-1] + rest
                    rest -= 1
                end
                tn
            end
        end
        new(solver, mm, trinumbers)
    end
end

@inline function Base.getindex(x::GetX, i::Integer)
    @inbounds begin
        if !isassigned(x.mm, i)
            x.mm[i] = LoRADS.get_X(x.solver, i)
        end
        return GetXI(x.solver, x.mm[i], x.trinumbers[i])
    end
end

struct GetXI
    solver::LoRADS.ASDP
    mm::Matrix{Cdouble}
    trinumbers::Vector{Int}
end

@inline function Base.getindex(x::GetXI, idx::Integer)
    @inbounds begin
        col = searchsortedlast(x.trinumbers, idx)
        row = idx - x.trinumbers[col] + col
        return x.mm[row, col]
    end
end

Solver.extract_moments(relaxation::AbstractRelaxation, (state, solver)::Tuple{StateMoment,LoRADS.ASDP}) =
    # Dilemma: Do we want to calculate the full matrix product U Uᵀ (or one triangle) and extract the relevant entries, which
    # would be favorable if we have a large moment matrix that contains everything? Or do we want to calculate the entries
    # on-the-fly doing vector products, which in total is not cache-optimal, but would be better if we only need few entries
    # from the matrix, perhaps because it represents a moment sub-matrix where a lot of entries have already been obtained by
    # other matrices? Let's do the first: We assume that this solver is used for large moment matrices, and we assume that
    # there won't be too much overlap between sparse groupings.
    MomentVector(relaxation, state.data, GetX(solver, state.data.psd_dim), LoRADS.get_Xlin(solver))

Solver.extract_sos(::AbstractRelaxation, ::Tuple{StateMoment,LoRADS.ASDP}, ::Val, ::Any, ::Nothing) =
    error("LoRADS does not provide useful information about the dual variables; no SOS data can be extracted.")

# While the fields are there, they are not populated by LoRADS (and how could they be?). And the slack variable calculation,
# which is also used for the dual infeasibility calculation, is just way off.
#=
function Solver.extract_sos(::AbstractRelaxation, (_, solver)::Tuple{StateMoment,LoRADS.ASDP}, ::Val,
    index::AbstractUnitRange, ::Nothing)
    vlagLp = unsafe_load(solver.vlagLp)
    x = Vector{Cdouble}(undef, length(index))
    unsafe_copyto!(pointer(x), vlagLp.matElem + (first(index) -1) * sizeof(Cdouble), length(index) * sizeof(Cdouble))
    x .^= 2
    return x
end

function Solver.extract_sos(::AbstractRelaxation, (_, solver)::Tuple{StateMoment,LoRADS.ASDP}, ::Val{:psd}, index::Integer,
    ::Nothing)
    return LoRADS.get_S(solver, index)
    # vlag = unsafe_load(unsafe_load(solver.Vlag, index))
    # s = Matrix{Cdouble}(undef, vlag.nRows, vlag.nRows)
    # BLAS.syrk!('L', 'N', true, unsafe_wrap(Array, vlag.matElem, (vlag.nRows, vlag.rank)), false, s)
    # return Symmetric(s, :L)
end
=#