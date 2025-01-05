mutable struct StateMoment{K<:Integer} <: AbstractAPISolver{K,Int32,Float64}
    const task::Mosek.Task
    num_solver_vars::Int32 # total number of variables available in the solver
    num_used_vars::Int32 # number of variables already used for something (might include scratch variables)
    num_cons::Int32
    num_afes::Int64
    const mon_to_solver::Dict{FastKey{K},Int32}
    info::Vector{<:Vector{<:Tuple{Symbol,Any}}}

    StateMoment{K}(task::Mosek.Task) where {K<:Integer} = new{K}(
        task, zero(Int32), zero(Int32), zero(Int32), zero(Int64), Dict{FastKey{K},Int32}()
    )
end

Solver.issuccess(::Val{:MosekMoment}, status::Mosek.Solsta) = status === Mosek.MSK_SOL_STA_OPTIMAL

function Base.append!(state::StateMoment, key)
    if state.num_used_vars == state.num_solver_vars
        newnum = overallocation(state.num_solver_vars + one(state.num_solver_vars))
        Mosek.@MSK_appendvars(state.task.task, newnum - state.num_solver_vars)
        state.num_solver_vars = newnum
    end
    return state.mon_to_solver[FastKey(key)] = let uv=state.num_used_vars
        state.num_used_vars += one(state.num_used_vars)
        uv
    end
end

Solver.supports_rotated_quadratic(::StateMoment) = true

Solver.psd_indextype(::StateMoment) = PSDIndextypeVector(:L)

@counter_alias(StateMoment, :nonnegative, :fix)
@counter_atomic(StateMoment, :quadratic)
@counter_atomic(StateMoment, (:rotated_quadratic, :psd), :quadratic)

function Solver.add_var_slack!(state::StateMoment, num::Int)
    if state.num_used_vars + num > state.num_solver_vars
        newnum = overallocation(state.num_used_vars + Int32(num))
        Mosek.@MSK_appendvars(state.task.task, newnum - state.num_solver_vars)
        state.num_solver_vars = newnum
    end
    result = state.num_used_vars:state.num_used_vars+Int32(num -1)
    state.num_used_vars += num
    return result
end

function Solver.add_constr_nonnegative!(state::StateMoment, indvals::Indvals{Int32,Float64})
    task = state.task.task
    Mosek.@MSK_appendcons(task, 1)
    Mosek.@MSK_putarow(task, state.num_cons, length(indvals), indvals.indices, indvals.values)
    Mosek.@MSK_putconbound(task, state.num_cons, MSK_BK_LO.value, 0., Inf)
    state.num_cons += one(state.num_cons)
    return
end

function Solver.add_constr_nonnegative!(state::StateMoment, indvals::IndvalsIterator{Int32,Float64})
    task = state.task.task
    Mosek.@MSK_appendcons(task, length(indvals))
    id = state.num_cons
    curcon = id
    for indval in indvals
        @capture(Mosek.@MSK_putarow(task, $curcon, length(indval), indval.indices, indval.values))
        curcon += one(curcon)
    end
    Mosek.@MSK_putconboundsliceconst(task, id, curcon, MSK_BK_LO.value, 0., Inf)
    state.num_cons = curcon
    return
end

function Solver.add_constr_quadratic!(state::StateMoment, indvals::IndvalsIterator{Int32,Float64}, ::Val{rotated}=Val(false)) where {rotated}
    task = state.task.task
    conedim = length(indvals)
    Mosek.@MSK_appendafes(task, conedim)
    id = state.num_afes
    curafe = id
    for indval in indvals
        @capture(Mosek.@MSK_putafefrow(task, $curafe, length(indval), indval.indices, indval.values))
        curafe += one(curafe)
    end
    state.num_afes = curafe
    # TODO: benchmark which approach is better, subsequence putafefrow as here or one putafefrowlist with preprocessing
    if !rotated && conedim == 3
        Mosek.@MSK_appendaccseq(task, 1, 3, id, C_NULL)
    elseif rotated && conedim == 3
        Mosek.@MSK_appendaccseq(task, 2, 3, id, C_NULL)
    elseif rotated && conedim == 4
        Mosek.@MSK_appendaccseq(task, 3, 4, id, C_NULL)
    else
        dom = Ref{Int64}()
        if rotated
            Mosek.@MSK_appendrquadraticconedomain(task, conedim, dom)
        else
            Mosek.@MSK_appendquadraticconedomain(task, conedim, dom)
        end
        Mosek.@MSK_appendaccseq(task, dom[], conedim, id, C_NULL)
    end
    return
end

Solver.add_constr_rotated_quadratic!(state::StateMoment, indvals::IndvalsIterator{Int32,Float64}) =
    add_constr_quadratic!(state, indvals, Val(true))

function Solver.add_constr_psd!(state::StateMoment, dim::Int, data::IndvalsIterator{Int32,Float64})
    task = state.task.task
    n = trisize(dim)
    Mosek.@MSK_appendafes(task, n)
    curafe = state.num_afes
    # putafefrowlist would require us to create collect(state.num_afes:state.num_afes+length(data)) for the afeidx parameter
    # and accumulate Base.index_lengths(data) for ptrrow
    for indvals in data
        @capture(Mosek.@MSK_putafefrow(task, $curafe, length(indvals), indvals.indices, indvals.values))
        curafe += 1
    end
    # for sure, the number of PSD cones will not be the bottleneck, so let's save ourselves the trouble of caching already
    # existing domains
    cone = Ref{Int64}()
    Mosek.@MSK_appendsvecpsdconedomain(task, n, cone)
    Mosek.@MSK_appendaccseq(task, cone[], n, state.num_afes, C_NULL)
    state.num_afes = curafe
    return
end

function Solver.add_constr_fix_prepare!(state::StateMoment, num::Int)
    task = state.task.task
    Mosek.@MSK_appendcons(task, num)
    # most of the constraints will be zero, so we can save some calls if we all set them to zero and only change the few
    # nonzero ones
    Mosek.@MSK_putconboundsliceconst(task, state.num_cons, state.num_cons + num, MSK_BK_FX.value, 0., 0.)
    return
end

function Solver.add_constr_fix!(state::StateMoment, ::Nothing, indvals::Indvals{Int32,Float64}, rhs::Float64)
    task = state.task.task
    Mosek.@MSK_putarow(task, state.num_cons, length(indvals), indvals.indices, indvals.values)
    iszero(rhs) || Mosek.@MSK_putconbound(task, state.num_cons, MSK_BK_FX.value, rhs, rhs)
    state.num_cons += one(state.num_cons)
    return
end

function Solver.fix_objective!(state::StateMoment, indvals::Indvals{Int32,Float64})
    Mosek.@MSK_putclist(state.task.task, length(indvals), indvals.indices, indvals.values)
    return
end

function Solver.poly_optimize(::Val{:MosekMoment}, relaxation::AbstractRelaxation, groupings::RelaxationGroupings;
    representation, verbose::Bool=false, customize=_ -> nothing, parameters...)
    task = Mosek.Task(msk_global_env::Env)
    setup_time = @elapsed begin
        K = _get_I(eltype(monomials(poly_problem(relaxation).objective)))

        taskptr = task.task
        verbose && putstreamfunc(task, MSK_STREAM_LOG, printstream)
        for (k, v) in parameters
            putparam(task, string(k), string(v))
        end
        Mosek.@MSK_putobjsense(taskptr, MSK_OBJECTIVE_SENSE_MINIMIZE.value)
        # just set up some domains, this is cheap and we don't have to query them afterwards
        domidx_ = Ref{Int64}()
        Mosek.@MSK_appendrplusdomain(taskptr, 1, domidx_)
        Mosek.@MSK_appendquadraticconedomain(taskptr, 3, domidx_)
        Mosek.@MSK_appendrquadraticconedomain(taskptr, 3, domidx_)
        Mosek.@MSK_appendrquadraticconedomain(taskptr, 4, domidx_)

        state = StateMoment{K}(task)
        state.info = moment_setup!(state, relaxation, groupings; representation)
        Mosek.@MSK_putvarboundsliceconst(taskptr, 0, state.num_used_vars, MSK_BK_FR.value, -Inf, Inf)
        customize(state)
    end
    @verbose_info("Setup complete in ", setup_time, " seconds")

    optimize(task)
    @verbose_info("Optimization complete")

    return state, getsolsta(task, MSK_SOL_ITR), getprimalobj(task, MSK_SOL_ITR)
end

function Solver.extract_moments(relaxation::AbstractRelaxation, state::StateMoment)
    x = Vector{Float64}(undef, state.num_used_vars)
    Mosek.@MSK_getxxslice(state.task.task, MSK_SOL_ITR.value, 0, length(x), x)
    return MomentVector(relaxation, x, state)
end

function Solver.extract_sos(::AbstractRelaxation, state::StateMoment, ::Union{Val{:fix},Val{:nonnegative}},
    index::AbstractUnitRange, ::Nothing)
    y = Vector{Float64}(undef, length(index))
    Mosek.@MSK_getyslice(state.task.task, MSK_SOL_ITR.value, first(index) -1, last(index), y)
    return y
end

Solver.extract_sos(::AbstractRelaxation, state::StateMoment,
    ::Union{Val{:quadratic},Val{:rotated_quadratic},Val{:psd}}, index::Integer, ::Nothing) =
    Mosek.getaccdoty(state.task, MSK_SOL_ITR, index)