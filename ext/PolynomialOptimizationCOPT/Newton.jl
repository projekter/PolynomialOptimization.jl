let
    status = quote
        if verbose
            nextinfo = time_ns()
            if nextinfo - lastinfo > 1_000_000_000
                if verbose
                    print("\33[2KStatus update: ", progress, " of ", nc, " (removed ", 100removed ÷ progress,
                        "% so far)\r")
                    flush(stdout)
                    lastinfo = nextinfo
                end
            end
        end
    end
    for checkall in (false, true)
        @eval function Newton.preproc(::Val{:COPT}, mons, vertexindices::$(checkall ? :(Val{:all}) : :(Any)), verbose,
            singlethread; parameters...)
            nv = nvariables(mons)
            nc = length(mons)
            nvertices = $(checkall ? :nc : :(length(vertexindices)))
            required_exps = fill!(BitVector(undef, nc), true)
            lastinfo = time_ns()
            task = COPTProb(copt_env)
            @inbounds begin
                singlethread && _check_ret(copt_env, COPT_SetIntParam(task, COPT_INTPARAM_THREADS, one(Cint)))
                _check_ret(copt_env, COPT_SetIntParam(task, COPT_INTPARAM_LOGGING, zero(Cint)))
                for (k, v) in parameters
                    if v isa Integer
                        _check_ret(copt_env, COPT_SetIntParam(task, k, Cint(v)))
                    elseif v isa AbstractFloat
                        _check_ret(copt_env, COPT_SetDblParam(task, k, Cdouble(v)))
                    else
                        throw(ArgumentError("Parameter $k is not of type Integer or AbstractFloat"))
                    end
                end
                # basic initialization: every point will get a variable
                ri = Ref{Cint}()
                idxs = collect(Cint(0):Cint(max(nv, nvertices) -1))
                tmp = Vector{Float64}(undef, max(nv, nvertices))
                let
                    $checkall && fill!(@view(tmp[1:nv]), 0.)
                    # ^ since we always fix the point in question to be -1, the sum of all points must be zero (condition nv +1)
                    # we must make sure that we have all the rows available, else AddCol silently won't fill them
                    _check_ret(copt_env, COPT_LoadProb(task, 0, nv, COPT_MAXIMIZE, 0., C_NULL, C_NULL, C_NULL, C_NULL, C_NULL,
                        C_NULL, C_NULL, C_NULL, C_NULL, tmp, tmp, C_NULL, C_NULL))
                    for exps in veciter($(checkall ? :mons : :(@view(mons[vertexindices]))))
                        copyto!(tmp, exps)
                        _check_ret(copt_env, COPT_AddCol(task, 0., nv, idxs, tmp, COPT_CONTINUOUS, 0., COPT_INFINITY, C_NULL))
                    end
                    @inbounds fill!(@view(tmp[1:nvertices]), 1.)
                    _check_ret(copt_env, COPT_AddRow(task, nvertices, idxs, tmp, 0, $((checkall ? (0., 0.) : (1., 1.))...),
                        C_NULL))
                end

                removed = 0
                $(checkall ?
                    quote
                        varnum = nvertices
                        dropvars = FastVec{Cint}(buffer=20)
                        progress = 1
                        vm = Ref{Cdouble}(-1.)
                        vz = Ref{Cdouble}(0.)
                        vi = Ref{Cdouble}(COPT_INFINITY)
                        for i in nvertices-1:-1:0
                            # first enforce this variable to be fixed: all others must add up to this point
                            ri[] = i
                            _check_ret(copt_env, COPT_SetColLower(task, 1, ri, vm))
                            _check_ret(copt_env, COPT_SetColUpper(task, 1, ri, vm))
                            # then try to find a solution
                            _check_ret(copt_env, COPT_SolveLp(task))
                            _check_ret(copt_env, COPT_GetIntAttr(task, COPT_INTATTR_LPSTATUS, ri))
                            if ri[] == COPT_LPSTATUS_OPTIMAL
                                # this was indeed possible, so our point is redundant; remove it!
                                ri[] = i
                                _check_ret(copt_env, COPT_SetColLower(task, 1, ri, vz))
                                _check_ret(copt_env, COPT_SetColUpper(task, 1, ri, vz))
                                required_exps[i+1] = false
                                removed += 1
                                if varnum ≥ 200
                                    unsafe_push!(dropvars, i)
                                    # Deleting Mosek variables is expensive, but every once in a while, it may be worth the
                                    # effort
                                    if length(dropvars) == 20
                                        _check_ret(copt_env, COPT_DelCols(task, length(dropvars), dropvars))
                                        varnum -= 20
                                        empty!(dropvars)
                                    end
                                end
                            else
                                # it was not possible, we must keep this point
                                ri[] = i
                                _check_ret(copt_env, COPT_SetColLower(task, 1, ri, vz))
                                _check_ret(copt_env, COPT_SetColUpper(task, 1, ri, vi))
                            end
                            progress += 1
                            $status
                        end
                    end :
                    quote
                        vertexpos = 1
                        vertexindex = vertexindices[begin]
                        for (progress, exps) in enumerate(veciter(mons))
                            if progress == vertexindex
                                vertexpos += 1
                                if vertexpos ≤ nvertices
                                    vertexindex = vertexindices[vertexpos]
                                end
                                continue
                            end
                            copyto!(tmp, exps)
                            _check_ret(copt_env, COPT_SetRowLower(task, nv, idxs, tmp))
                            _check_ret(copt_env, COPT_SetRowUpper(task, nv, idxs, tmp))
                            _check_ret(copt_env, COPT_SolveLp(task))
                            _check_ret(copt_env, COPT_GetIntAttr(task, COPT_INTATTR_LPSTATUS, ri))
                            if ri[] == COPT_LPSTATUS_OPTIMAL
                                required_exps[progress] = false
                                removed += 1
                            end
                            $status
                        end
                    end)
            end
            return required_exps
        end
    end
end

function Newton.prepare(::Val{:COPT}, mons, num, verbose; parameters...)
    nv = nvariables(mons)
    nc = length(mons)
    # now we build a task with the minimal number of extremal points and try to find every possible part of the Newton polytope
    # for SOS polynomials
    task = COPTProb(copt_env)
    _check_ret(copt_env, COPT_SetIntParam(task, COPT_INTPARAM_LOGGING, zero(Cint)))
    for (k, v) in parameters
        if v isa Integer
            _check_ret(copt_env, COPT_SetIntParam(task, k, Cint(v)))
        elseif v isa AbstractFloat
            _check_ret(copt_env, COPT_SetDblParam(task, k, Cdouble(v)))
        else
            throw(ArgumentError("Parameter $k is not of type Integer or AbstractFloat"))
        end
    end
    let
        idxs = collect(Int32(0):Int32(max(nv, nc) -1))
        tmp = Vector{Float64}(undef, max(nv, nc))
        # we must make sure that we have all the rows available, else AddCol silently won't fill them
        _check_ret(copt_env, COPT_LoadProb(task, 0, nv, COPT_MAXIMIZE, 0., C_NULL, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL,
            C_NULL, C_NULL, C_NULL, tmp, tmp, C_NULL, C_NULL))
        for cf in veciter(mons)
            @inbounds @view(tmp[1:nv]) .= 0.5 .* cf
            _check_ret(copt_env, COPT_AddCol(task, 0., nv, idxs, tmp, COPT_CONTINUOUS, 0., COPT_INFINITY, C_NULL))
        end
        @inbounds fill!(@view(tmp[1:nc]), 1.)
        _check_ret(copt_env, COPT_AddRow(task, nc, idxs, tmp, 0, 1., 1., C_NULL))
    end
    if num < 10_000 || isone(nv)
        nthreads = 1
        secondtask = nothing
    else
        nthreads = Threads.nthreads()
        if nthreads > 1
            # we need to figure out if we have enough memory for all the threads. Let's assume that at least a second task can
            # be created.
            mem = @allocdiff begin
                secondtask = COPTProb(copt_env)
                _check_ret(copt_env, COPT_CreateCopy(task, secondtask))
                COPT_SetIntParam(secondtask, COPT_INTPARAM_THREADS, one(Cint))
                COPT_Solve(task)
            end
            if mem ≤ 0
                @verbose_info("Memory requirements of a single thread could not be determined, using all available threads")
            else
                @verbose_info("Memory requirements of a single thread: ", div(mem, 1024*1024, RoundUp), " MiB")
                nthreads = min(nthreads, Int(Sys.free_memory() ÷ mem +2))
            end
            # Note that this is potentially still an underestimation, as our candidates list will also grow. But this is
            # something that can potentially be swapped, so if swap space is available beyond the free_memory limit, then we
            # are still fine.
        else
            secondtask = nothing
        end
    end
    isone(nthreads) || COPT_SetIntParam(task, COPT_INTPARAM_THREADS, one(Cint)) # single-threaded for COPT itself
    return nthreads, task, secondtask
end

Newton.alloc_global(::Val{:COPT}, nv) = collect(Int32(0):Int32(nv -1))
Newton.alloc_local(::Val{:COPT}, nv) = Vector{Float64}(undef, nv)
function Newton.clonetask(t::COPTProb)
    secondtask = COPTProb(copt_env)
    _check_ret(copt_env, COPT_CreateCopy(task, secondtask))
    return secondtask
end

@inline function Newton.work(::Val{:COPT}, task, idxs, tmp, expiter, Δprogress, Δacceptance, add_callback, iteration_callback)
    for (idx, exponents) in expiter
        # check the previous exponent in the linear program and add it if possible
        copyto!(tmp, exponents)
        _check_ret(copt_env, COPT_SetRowLower(task, length(idxs), idxs, tmp))
        _check_ret(copt_env, COPT_SetRowUpper(task, length(idxs), idxs, tmp))
        status = Ref{Cint}()
        _check_ret(copt_env, COPT_SolveLp(task))
        _check_ret(copt_env, COPT_GetIntAttr(task, COPT_INTATTR_LPSTATUS, status))
        if status[] == COPT_LPSTATUS_OPTIMAL
            # this candidate is part of the Newton polytope
            @inline add_callback(idx)
            Δacceptance[] += 1
        end
        Δprogress[] += 1
        isnothing(iteration_callback) || @inline iteration_callback()
    end
    return
end