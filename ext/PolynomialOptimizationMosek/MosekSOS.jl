mutable struct StateSOS{K<:Integer} <: AbstractAPISolver{K,Int32,Float64}
    const task::Mosek.Task
    num_vars::Int32
    num_bar_vars::Int32
    num_solver_cons::Int32 # total number of constraints available in the solver
    num_used_cons::Int32 # number of constraints already used for something (might include scratch constraints)
    const mon_to_solver::Dict{FastKey{K},Int32}
    info::Vector{<:Vector{<:Tuple{Symbol,Any}}}

    StateSOS{K}(task::Mosek.Task) where {K<:Integer} = new{K}(
        task, zero(Int32), zero(Int32), zero(Int32), zero(Int32), Dict{FastKey{K},Int32}()
    )
end

Solver.issuccess(::Val{:MosekSOS}, status::Mosek.Solsta) = status === Mosek.MSK_SOL_STA_OPTIMAL

function Base.append!(state::StateSOS, key)
    if state.num_used_cons == state.num_solver_cons
        newcon = overallocation(state.num_solver_cons + one(state.num_solver_cons))
        Mosek.@MSK_appendcons(state.task.task, newcon - state.num_solver_cons)
        state.num_solver_cons = newcon
    end
    return state.mon_to_solver[FastKey(key)] = let uc=state.num_used_cons
        state.num_used_cons += one(state.num_used_cons)
        uc
    end
end

Solver.supports_rotated_quadratic(::StateSOS) = true

Solver.psd_indextype(::StateSOS) = PSDIndextypeMatrixCartesian(:L, zero(Int32))

@counter_alias(StateSOS, (:nonnegative, :quadratic, :rotated_quadratic), :free)
@counter_atomic(StateSOS, :psd)

function Solver.add_constr_slack!(state::StateSOS, num::Int)
    if state.num_used_cons + num > state.num_solver_cons
        newnum = overallocation(state.num_used_cons + Int32(num))
        Mosek.@MSK_appendcons(state.task.task, newnum - state.num_solver_cons)
        state.num_solver_cons = newnum
    end
    result = state.num_used_cons:state.num_used_cons+Int32(num -1)
    state.num_used_cons += num
    return result
end

function Solver.add_var_nonnegative!(state::StateSOS, indvals::Indvals{Int32,Float64})
    task = state.task.task
    Mosek.@MSK_appendvars(task, 1)
    Mosek.@MSK_putacol(task, state.num_vars, length(indvals), indvals.indices, indvals.values)
    Mosek.@MSK_putvarbound(task, state.num_vars, MSK_BK_LO.value, 0., Inf)
    state.num_vars += one(state.num_vars)
    return
end

function Solver.add_var_nonnegative!(state::StateSOS, indvals::IndvalsIterator{Int32,Float64})
    task = state.task.task
    Mosek.@MSK_appendvars(task, length(indvals))
    id = state.num_vars
    curvar = id
    for indval in indvals
        @capture(Mosek.@MSK_putacol(task, $curvar, length(indval), indval.indices, indval.values))
        curvar += one(curvar)
    end
    Mosek.@MSK_putvarboundsliceconst(task, id, curvar, MSK_BK_LO.value, 0., Inf)
    state.num_vars = curvar
    return
end

function Solver.add_var_quadratic!(state::StateSOS, indvals::IndvalsIterator{Int32,Float64}, ::Val{rotated}=Val(false)) where {rotated}
    task = state.task.task
    conedim = length(indvals)
    rhsdim = conedim -2
    Mosek.@MSK_appendvars(task, conedim)
    id = state.num_vars
    if rotated
        Mosek.@MSK_putvarboundsliceconst(task, id, id + Int32(2), MSK_BK_LO.value, 0., Inf)
    else
        Mosek.@MSK_putvarbound(task, id, MSK_BK_LO.value, 0., Inf)
    end
    Mosek.@MSK_putvarboundsliceconst(task, id + Int32(rotated ? 2 : 1), id + Int32(conedim), MSK_BK_FR.value, -Inf, Inf)
    Mosek.@MSK_appendconeseq(task, (rotated ? MSK_CT_RQUAD : MSK_CT_QUAD).value, 0., conedim, id)
    curvar = id
    for indval in indvals
        @capture(Mosek.@MSK_putacol(task, $curvar, length(indval), indval.indices, indval.values))
        curvar += one(curvar)
    end
    state.num_vars = curvar
    return
end

Solver.add_var_rotated_quadratic!(state::StateSOS, indvals::IndvalsIterator{Int32,Float64}) =
    add_var_quadratic!(state, indvals, Val(true))

function Solver.add_var_psd!(state::StateSOS, dim::Int, data::PSDMatrixCartesian{Int32,Float64})
    task = state.task.task
    Mosek.@MSK_appendbarvars(task, 1, Ref(Int32(dim)))
    outidx = Ref{Int64}()
    oneref = Ref{Float64}(1.)
    @inbounds for (midx, (symrows, symcols, vals)) in data
        Mosek.@MSK_appendsparsesymmat(task, dim, length(symrows), symrows, symcols, vals, outidx)
        Mosek.@MSK_putbaraij(task, midx, state.num_bar_vars, 1, outidx, oneref)
    end
    state.num_bar_vars += one(state.num_bar_vars)
    return
end

function Solver.add_var_free_prepare!(state::StateSOS, num::Int)
    task = state.task.task
    Mosek.@MSK_appendvars(task, num)
    Mosek.@MSK_putvarboundsliceconst(task, state.num_vars, state.num_vars + num, MSK_BK_FR.value, -Inf, Inf)
    return
end

function Solver.add_var_free!(state::StateSOS, ::Nothing, indvals::Indvals{Int32,Float64}, obj::Float64)
    task = state.task.task
    Mosek.@MSK_putacol(task, state.num_vars, length(indvals), indvals.indices, indvals.values)
    iszero(obj) || Mosek.@MSK_putcj(task, state.num_vars, obj)
    state.num_vars += 1
    return
end

function Solver.fix_constraints!(state::StateSOS, indvals::Indvals{Int32,Float64})
    task = state.task.task
    len = length(indvals)
    Mosek.@MSK_putconboundsliceconst(task, 0, state.num_used_cons, MSK_BK_FX, 0., 0.)
    Mosek.@MSK_putconboundlist(task, len, indvals.indices, fill(MSK_BK_FX.value, len), indvals.values, indvals.values)
    return
end

function Solver.poly_optimize(::Val{:MosekSOS}, relaxation::AbstractRelaxation, groupings::RelaxationGroupings;
    representation, verbose::Bool=false, customize=(state) -> nothing, parameters...)
    task = Mosek.Task(msk_global_env::Env)
    setup_time = @elapsed begin
        K = _get_I(eltype(monomials(poly_problem(relaxation).objective)))

        verbose && putstreamfunc(task, MSK_STREAM_LOG, printstream)
        for (k, v) in parameters
            putparam(task, string(k), string(v))
        end
        Mosek.@MSK_putobjsense(task.task, MSK_OBJECTIVE_SENSE_MAXIMIZE.value)

        state = StateSOS{K}(task)
        state.info = sos_setup!(state, relaxation, groupings; representation)
        customize(state)
    end
    @verbose_info("Setup complete in ", setup_time, " seconds")

    optimize(task)
    @verbose_info("Optimization complete")

    return state, getsolsta(task, MSK_SOL_ITR), getprimalobj(task, MSK_SOL_ITR)
end

function Solver.extract_moments(relaxation::AbstractRelaxation, state::StateSOS)
    y = Vector{Float64}(undef, length(state.mon_to_solver))
    Mosek.@MSK_getyslice(state.task.task, MSK_SOL_ITR.value, 0, length(y), y)
    return MomentVector(relaxation, y, state)
end

function Solver.extract_sos(::AbstractRelaxation, state::StateSOS, ::Val, index::AbstractUnitRange, ::Nothing)
    x = Vector{Float64}(undef, length(index))
    Mosek.@MSK_getxxslice(state.task.task, MSK_SOL_ITR.value, first(index) -1, last(index), x)
    return x
end

function Solver.extract_sos(::AbstractRelaxation, state::StateSOS, ::Val{:psd}, index::Integer, ::Nothing)
    dim = Ref{Int32}()
    Mosek.@MSK_getdimbarvarj(state.task.task, index -1, dim)
    s = Vector{Float64}(undef, trisize(dim[]))
    Mosek.@MSK_getbarxj(state.task.task, MSK_SOL_ITR.value, index -1, s)
    return SPMatrix(dim[], s, :L)
end