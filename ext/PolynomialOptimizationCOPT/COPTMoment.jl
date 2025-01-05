mutable struct StateMoment{K<:Integer} <: AbstractAPISolver{K,Cint,Float64}
    problem::COPTProb
    num_solver_vars::Cint # total number of variables available in the solver
    num_used_vars::Cint # number of variables already used for something (might include scratch variables)
    num_symmat::Cint
    const mon_to_solver::Dict{FastKey{K},Cint}
    info::Vector{<:Vector{<:Tuple{Symbol,Any}}}

    StateMoment{K}(task::COPTProb) where {K<:Integer} = new{K}(
        task, zero(Cint), zero(Cint), zero(Cint), Dict{FastKey{K},Cint}()
    )
end

Solver.issuccess(::Val{:COPTMoment}, status::Cint) = status == COPT_LPSTATUS_OPTIMAL # maybe also COPT_LPSTATUS_IMPRECISE

function Base.append!(state::StateMoment, key)
    if state.num_used_vars == state.num_solver_vars
        newnum = overallocation(state.num_solver_vars + one(state.num_solver_vars))
        Δ = newnum - state.num_solver_vars
        _check_ret(copt_env, COPT_AddCols(state.problem, Δ, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL,
            fill(-COPT_INFINITY, Δ), C_NULL, C_NULL))
        state.num_solver_vars = newnum
    end
    return state.mon_to_solver[FastKey(key)] = let uv=state.num_used_vars
        state.num_used_vars += one(state.num_used_vars)
        uv
    end
end

Solver.supports_rotated_quadratic(::StateMoment) = true

Solver.psd_indextype(::StateMoment) = PSDIndextypeMatrixCartesian(:L, zero(Cint))

@counter_alias(StateMoment, (:nonnegative, :quadratic, :rotated_quadratic), :fix)
@counter_atomic(StateMoment, :psd)

function Solver.add_var_slack!(state::StateMoment, num::Int)
    if state.num_used_vars + num > state.num_solver_vars
        newnum = overallocation(state.num_used_vars + Cint(num))
        Δ = newnum - state.num_solver_vars
        _check_ret(copt_env, COPT_AddCols(state.problem, Δ, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL,
            fill(-COPT_INFINITY, Δ), C_NULL, C_NULL))
        state.num_solver_vars = newnum
    end
    result = state.num_used_vars:state.num_used_vars+Cint(num -1)
    state.num_used_vars += num
    return result
end

function Solver.add_constr_nonnegative!(state::StateMoment, indvals::Indvals{Cint,Float64})
    _check_ret(copt_env, COPT_AddRow(state.problem, length(indvals), indvals.indices, indvals.values, 0, 0., COPT_INFINITY,
        C_NULL))
    return
end

function Solver.add_constr_nonnegative!(state::StateMoment, indvals::IndvalsIterator{Cint,Float64})
    N = length(indvals)
    rowMatBeg = Vector{Cint}(undef, N +1)
    @inbounds rowMatBeg[1] = 0
    @inbounds for (i, len) in zip(Iterators.countfrom(2), Base.index_lengths(indvals))
        rowMatBeg[i] = rowMatBeg[i-1] + len
    end
    _check_ret(copt_env, COPT_AddRows(state.problem, N, rowMatBeg, C_NULL, rowvals(indvals), nonzeros(indvals), C_NULL,
        fill(0., N), fill(COPT_INFINITY, N), C_NULL))
    return
end

function Solver.add_constr_quadratic!(state::StateMoment, indvals::IndvalsIterator{Cint,Float64}, rotated::Bool=false)
    N = length(indvals)
    # COPT does not support an arbitrary-content quadratic cone. Either we put the cone into the form <x, Qx> ≤ b or we
    # need to create more variables that wrap this functionality.
    if state.num_solver_vars < state.num_used_vars + N
        newnum = overallocation(state.num_solver_vars + Cint(N))
        Δ = newnum - state.num_solver_vars
        lb = fill(-COPT_INFINITY, Δ)
        @inbounds lb[1] = 0.
        rotated && @inbounds lb[2] = 0.
        _check_ret(copt_env, COPT_AddCols(state.problem, Δ, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL,
            lb, C_NULL, C_NULL))
        state.num_solver_vars += Δ
    end
    rowMatBeg = Vector{Cint}(undef, N +1)
    @inbounds rowMatBeg[1] = 0
    @inbounds for (i, len) in zip(Iterators.countfrom(2), Base.index_lengths(indvals))
        rowMatBeg[i] = rowMatBeg[i-1] + len +1
    end
    rowMatIdx = Vector{Cint}(undef, length(rowvals(indvals)) + N)
    rowMatElem = similar(rowMatIdx, Cdouble)
    i = 1
    @inbounds for (j, indval) in zip(Iterators.countfrom(Cint(0)), indvals)
        n = length(indval)
        copyto!(rowMatIdx, i, indval.indices, 1, n)
        copyto!(rowMatElem, i, indval.values, 1, n)
        i += n
        rowMatIdx[i] = state.num_used_vars + j
        rowMatElem[i] = -1.
        i += 1
    end
    zv = fill(0., N)
    _check_ret(copt_env, COPT_AddRows(state.problem, N, rowMatBeg, C_NULL, rowMatIdx, rowMatElem, C_NULL, zv, zv, C_NULL))
    _check_ret(copt_env, COPT_AddCones(state.problem, 1, Ref(Cint(rotated ? COPT_CONE_RQUAD : COPT_CONE_QUAD)), Ref(Cint(0)),
        Ref(Cint(N)), [state.num_used_vars + i for i in Cint(0):Cint(N -1)]))
    state.num_used_vars += N
    return
end

Solver.add_constr_rotated_quadratic!(state::StateMoment, indvals::IndvalsIterator{Cint,Float64}) =
    add_constr_quadratic!(state, indvals, true)

function Solver.add_constr_psd!(state::StateMoment, dim::Int, data::PSDMatrixCartesian{Cint,Float64})
    colIdx = data.rowind    # We can do variable piracy: those vectors have more than sufficient length, and we will only write
    symMatIdx = data.colind # to positions that were already traveled.
    start = state.num_symmat
    j = 1
    @inbounds for (midx, (symrows, symcols, vals)) in data
        _check_ret(copt_env, COPT_AddSymMat(state.problem, dim, length(symrows), symrows, symcols, vals))
        colIdx[j] = midx
        symMatIdx[j] = state.num_symmat
        state.num_symmat += 1
        j += 1
    end
    # It appears (undocumented) that we can use -1 to specify that there is no constant term present. Other ways of enforcing
    # this don't work:
    # - creation of a zero-element symmat fails
    # - creation of a one-element symmat with just a zero entry works; but as soon as a quadratic constraint is present, the
    #   solver fails with Invalid Data.
    _check_ret(copt_env, COPT_AddLMIConstr(state.problem, dim, j -1, colIdx, symMatIdx, Cint(-1), C_NULL))
    return
end

function Solver.add_constr_fix!(state::StateMoment, ::Nothing, indvals::Indvals{Cint,Float64}, rhs::Float64)
    _check_ret(copt_env, COPT_AddRow(state.problem, length(indvals), indvals.indices, indvals.values, 0, rhs, rhs, C_NULL))
    return
end

function Solver.fix_objective!(state::StateMoment, indvals::Indvals{Cint,Float64})
    _check_ret(copt_env, COPT_ReplaceColObj(state.problem, length(indvals), indvals.indices, indvals.values))
    return
end

function Solver.poly_optimize(::Val{:COPTMoment}, relaxation::AbstractRelaxation{<:Problem{P}}, groupings::RelaxationGroupings;
    representation, verbose::Bool=false, customize::Function=_ -> nothing, parameters=()) where {P}
    setup_time = @elapsed begin
        K = _get_I(eltype(monomials(poly_problem(relaxation).objective)))

        task = COPTProb(copt_env)
        _check_ret(copt_env, COPT_SetIntParam(task, COPT_INTPARAM_LOGTOCONSOLE, Cint(verbose)))
        for (k, v) in parameters
            if v isa Integer
                _check_ret(copt_env, COPT_SetIntParam(task, k, Cint(v)))
            elseif v isa AbstractFloat
                _check_ret(copt_env, COPT_SetDblParam(task, k, Cdouble(v)))
            else
                throw(ArgumentError("Parameter $k is not of type Integer or AbstractFloat"))
            end
        end

        _check_ret(copt_env, COPT_SetObjSense(task, COPT_MINIMIZE))

        state = StateMoment{K}(task)
        state.info = moment_setup!(state, relaxation, groupings; representation)
        customize(state)
    end
    @verbose_info("Setup complete in ", setup_time, " seconds")

    _check_ret(copt_env, COPT_Solve(task))
    status = Ref{Cint}()
    _check_ret(copt_env, COPT_GetIntAttr(task, COPT_INTATTR_LPSTATUS, status))
    value = Ref{Cdouble}()
    _check_ret(copt_env, COPT_GetDblAttr(task, COPT_DBLATTR_LPOBJVAL, value))
    @verbose_info("Optimization complete")

    if status[] ∈ (COPT_LPSTATUS_OPTIMAL, COPT_LPSTATUS_IMPRECISE)
        return state, status[], value[]
    else
        return missing, status[], value[]
    end
end

function Solver.extract_moments(relaxation::AbstractRelaxation, state::StateMoment)
    x = Vector{Cdouble}(undef, state.num_solver_vars)
    _check_ret(copt_env, COPT_GetLpSolution(state.problem, x, C_NULL, C_NULL, C_NULL))
    return MomentVector(relaxation, resize!(x, state.num_used_vars), state)
end

function Solver.extract_sos_prepare(relaxation::AbstractRelaxation, state::StateMoment)
    num = Ref{Cint}()
    _check_ret(copt_env, COPT_GetIntAttr(state.problem, COPT_INTATTR_ROWS, num))
    need_slack = state.num_used_vars > length(state.mon_to_solver)
    x = Vector{Cdouble}(undef, need_slack ? state.num_solver_vars : 0)
    z = Vector{Cdouble}(undef, num[])
    _check_ret(copt_env, COPT_GetLpSolution(state.problem, need_slack ? x : C_NULL, C_NULL, z, C_NULL))
    return x, z
end

Solver.extract_sos(relaxation::AbstractRelaxation, state::StateMoment, ::Val, index::AbstractUnitRange, (_, z)) =
    @view(z[index])

function Solver.extract_sos(relaxation::AbstractRelaxation, state::StateMoment, ::Val{:psd}, index::Integer, _)
    dim = Ref{Cint}()
    _check_ret(copt_env, COPT_GetLMIConstr(state.problem, index -1, dim, C_NULL, C_NULL, C_NULL, C_NULL, 0, C_NULL))
    z = Matrix{Cdouble}(undef, dim[], dim[])
    _check_ret(copt_env, COPT_GetLMIConstrInfo(state.problem, COPT_DBLINFO_DUAL, index -1, z))
    return z
end