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
        @eval function Newton.preproc(::Val{:Mosek}, mons, vertexindices::$(checkall ? :(Val{:all}) : :(Any)), verbose,
            singlethread; parameters...)
            nv = nvariables(mons)
            nc = length(mons)
            nvertices = $(checkall ? :nc : :(length(vertexindices)))
            required_exps = fill!(BitVector(undef, nc), true)
            lastinfo = time_ns()
            task = Mosek.Task(msk_global_env::Env)
            @inbounds try
                singlethread && putintparam(task, MSK_IPAR_NUM_THREADS, 1)
                # verbose && putstreamfunc(task, MSK_STREAM_LOG, printstream)
                for (k, v) in parameters
                    putparam(task, string(k), v)
                end
                # basic initialization: every point will get a variable
                appendvars(task, nvertices)
                putvarboundsliceconst(task, 1, nvertices +1, MSK_BK_LO, 0., Inf)
                appendcons(task, nv +1)
                $(checkall ? :(putconboundsliceconst(task, 1, nv +2, MSK_BK_FX, 0., 0.)) :
                             :(tmp = Vector{Float64}(undef, max(nv, nvertices))))
                # ^ since we always fix the point in question to be -1, the sum of all points must be zero (condition nv +1)
                let
                    idxs = collect(Int32(0):Int32(max(nv, nvertices) -1))
                    $(checkall ? :(tmp = Vector{Float64}(undef, max(nv, nvertices))) : :(nothing))
                    for (i, exps) in zip(Iterators.countfrom(zero(Int32)),
                                         veciter($(checkall ? :mons : :(@view(mons[vertexindices])))))
                        copyto!(tmp, exps)
                        let task=task.task, i=i, nv=nv, idxs=idxs, tmp=tmp
                            Mosek.@MSK_putacol(task, i, nv, idxs, tmp)
                        end
                    end
                    @inbounds fill!(@view(tmp[1:nvertices]), 1.)
                    let task=task.task, nv=nv, nvertices=nvertices, idxs=idxs, tmp=tmp
                        Mosek.@MSK_putarow(task, nv, nvertices, idxs, tmp)
                    end
                    $checkall || putconbound(task, nv +1, MSK_BK_FX, 1., 1.)
                end

                removed = 0
                $(checkall ?
                    quote
                        varnum = nvertices
                        dropvars = FastVec{Int32}(buffer=20)
                        progress = 1
                        for i in nvertices:-1:1
                            # first enforce this variable to be fixed: all others must add up to this point
                            putvarbound(task, i, MSK_BK_FX, -1., -1.)
                            # then try to find a solution
                            optimize(task)
                            if getsolsta(task, MSK_SOL_BAS) == MSK_SOL_STA_OPTIMAL
                                # this was indeed possible, so our point is redundant; remove it!
                                putvarbound(task, i, MSK_BK_FX, 0., 0.)
                                required_exps[i] = false
                                removed += 1
                                if varnum ≥ 200
                                    unsafe_push!(dropvars, i -1)
                                    # Deleting Mosek variables is expensive, but every once in a while, it may be worth the
                                    # effort
                                    if length(dropvars) == 20
                                        let task=task.task, dropvars=dropvars
                                            Mosek.@MSK_removevars(task, length(dropvars), dropvars)
                                        end
                                        varnum -= 20
                                        empty!(dropvars)
                                    end
                                end
                            else
                                # it was not possible, we must keep this point
                                putvarbound(task, i, MSK_BK_LO, 0., Inf)
                            end
                            progress += 1
                            $status
                        end
                    end :
                    quote
                        fx = fill(MSK_BK_FX.value, nv)
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
                            Mosek.@MSK_putconboundslice(task.task, 0, nv, fx, tmp, tmp)
                            optimize(task)
                            if getsolsta(task, MSK_SOL_BAS) == MSK_SOL_STA_OPTIMAL
                                required_exps[progress] = false
                                removed += 1
                            end
                            $status
                        end
                    end)
            finally
                deletetask(task)
            end
            return required_exps
        end
    end
end

function Newton.prepare(::Val{:Mosek}, mons, num, verbose; parameters...)
    nv = nvariables(mons)
    nc = length(mons)
    # now we build a task with the minimal number of extremal points and try to find every possible part of the Newton polytope
    # for SOS polynomials
    task = Mosek.Task(msk_global_env::Env)
    # verbose && putstreamfunc(task, MSK_STREAM_LOG, printstream)
    for (k, v) in parameters
        putparam(task, string(k), v)
    end
    appendvars(task, nc)
    putvarboundsliceconst(task, 1, nc +1, MSK_BK_LO, 0., Inf)
    appendcons(task, nv +1)
    let
        idxs = collect(Int32(0):Int32(max(nv, nc) -1))
        tmp = Vector{Float64}(undef, max(nv, nc))
        for (i, cf) in zip(Iterators.countfrom(zero(Int32)), veciter(mons))
            @inbounds @view(tmp[1:nv]) .= 0.5 .* cf
            let task=task.task, i=i, nv=nv, idxs=idxs, tmp=tmp
                Mosek.@MSK_putacol(task, i, nv, idxs, tmp)
            end
        end
        @inbounds fill!(@view(tmp[1:nc]), 1.)
        let task=task.task, nv=nv, nc=nc, idxs=idxs, tmp=tmp
            Mosek.@MSK_putarow(task, nv, nc, idxs, tmp)
        end
    end
    putconbound(task, nv +1, MSK_BK_FX, 1., 1.)
    if num < 10_000 || isone(nv)
        nthreads = 1
        secondtask = nothing
    else
        nthreads = Threads.nthreads()
        if nthreads > 1
            # we need to figure out if we have enough memory for all the threads. Unfortunately, getmemusage() seems to
            # always return 1081, "No available information about the space usage." - so we need to fetch the information
            # ourselves. Let's assume that at least a second task can be created.
            mem = @allocdiff begin
                secondtask = Mosek.Task(task)
                putintparam(secondtask, MSK_IPAR_NUM_THREADS, 1)
                optimize(secondtask)
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
    isone(nthreads) || putintparam(task, MSK_IPAR_NUM_THREADS, 1) # single-threaded for Mosek itself
    return nthreads, task, secondtask
end

Newton.alloc_global(::Val{:Mosek}, nv) = fill(MSK_BK_FX.value, nv)
Newton.alloc_local(::Val{:Mosek}, nv) = Vector{Float64}(undef, nv)
Newton.clonetask(t::Mosek.Task) = Mosek.Task(t)

@inline function Newton.work(::Val{:Mosek}, task, bk, tmp, expiter, Δprogress, Δacceptance, add_callback, iteration_callback)
    for (idx, exponents) in expiter
        # check the previous exponent in the linear program and add it if possible
        copyto!(tmp, exponents)
        let task=task.task, bk=bk, tmp=tmp
            Mosek.@MSK_putconboundslice(task, 0, length(bk), bk, tmp, tmp)
        end
        optimize(task)
        if getsolsta(task, MSK_SOL_BAS) == MSK_SOL_STA_OPTIMAL
            # this candidate is part of the Newton polytope
            @inline add_callback(idx)
            Δacceptance[] += 1
        end
        Δprogress[] += 1
        isnothing(iteration_callback) || @inline iteration_callback()
    end
    return
end