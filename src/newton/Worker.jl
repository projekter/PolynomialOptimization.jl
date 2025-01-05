function verbose_worker(num, init_progress, init_time, threadprogress, threadacceptance)
    self_running = Ref(true)
    Threads.@spawn(:interactive,
        while $self_running[]
            allprogress = sum(threadprogress, init=0)
            allacceptance = sum(threadacceptance, init=0)
            Δt = allprogress == $init_progress ? 1 : allprogress - init_progress
            # ^ if a finished job is started, this might happen
            rem_sec = round(Int, ((time_ns() - init_time) / 1_000_000_000Δt) * (num - allprogress))
            @printf("\33[2KStatus update: %.2f%%, acceptance: %.2f%%, remaining time: %02d:%02dmin\r",
                100allprogress / num, 100allacceptance / allprogress, rem_sec ÷ 60, rem_sec % 60)
            flush(stdout)
            sleep(1)
        end)
    return self_running
end

function step_callback(verbose, ::Nothing, lastinfo, candidates, progress, acceptance)
    if verbose
        return () -> let nextinfo=time_ns()
            if nextinfo - lastinfo[] > 1_000_000_000
                yield()
                lastinfo[] = nextinfo
            end
        end
    else
        return nothing
    end
end

step_callback(verbose, (fileprogress, fileout), lastinfo, candidates, progress, acceptance) = () -> let nextinfo=time_ns()
    if nextinfo - lastinfo[] > 1_000_000_000
        write(fileout, candidates)
        flush(fileout)
        seekstart(fileprogress)
        write(fileprogress, progress[], acceptance[])
        flush(fileprogress)
        empty!(candidates)
        verbose && yield()
        lastinfo[] = nextinfo
    end
end

function execute_taskfun(V, tid, task, iter, nv, data_global, cond, progresses, acceptances, allcandidates, verbose, filestuff)
    lastinfo = Ref(time_ns())
    data_local = alloc_local(V, nv)
    candidates = FastVec{UInt}()
    progress = Ref(progresses, tid)
    acceptance = Ref(acceptances, tid)
    cb = step_callback(verbose, filestuff, lastinfo, candidates, progress, acceptance)
    try
        work(V, task, data_global, data_local, iter, progress, acceptance, @capture(idx -> push!($candidates, idx)), cb)
        if isnothing(filestuff)
            lock(cond)
            try
                append!(allcandidates, candidates)
            finally
                unlock(cond)
            end
        else
            write(filestuff[2], candidates)
            seekstart(filestuff[1])
            write(filestuff[1], progress[], acceptance[])
        end
    finally
        finalize(task)
        if !isnothing(filestuff)
            close.(filestuff)
        end
    end
end

function execute(V, verbose, mons, task, filepath)
    @verbose_info("Preparing to determine Newton polytope (single-threaded)")
    nv = nvariables(mons)
    num = length(mons)

    # While we precalculate the size of the list exactly, we don't pre-allocate the output candidates - we hope to eliminate a
    # lot of exponents by the polytope containment, so we might overallocate so much memory that we hit a resource constraint
    # here. If instead, we grow the list dynamically, we pay the price in speed, but the impossible might become feasible.
    candidates = FastVec{UInt}()

    progress = Ref(0)
    acceptance = Ref(0)
    if isnothing(filepath)
        filestuff = nothing
        iter = veciter(mons, Val(true))
    else
        fileprogress = open("$filepath.prog", read=true, write=true, create=true, lock=false)
        fs = filesize(fileprogress)
        if fs == 2sizeof(Int)
            progress[] = read(fileprogress, Int)
            acceptance[] = read(fileprogress, Int)
        elseif !iszero(fs)
            error("Unknown progress file format - please delete existing files.")
        end
        iter = veciter(@view(mons[progress[]+1:end]), Val(true))
        if isempty(iter)
            close(fileprogress)
            return progress[], acceptance[]
        end
        fileout = open("$filepath.out", append=true, lock=false)
        filestuff = (fileprogress, fileout)
    end
    if verbose
        self_running = verbose_worker(num, progress[], time_ns(), progress, acceptance)
    end
    lastinfo = Ref(time_ns())
    cb = step_callback(verbose, filestuff, lastinfo, candidates, progress, acceptance)
    try
        work(V, task, alloc_global(V, nv), alloc_local(V, nv), iter, progress, acceptance,
            @capture(idx -> push!($candidates, idx)), cb)
        verbose && print("\33[2K")
        if !isnothing(filepath)
            write(fileout, candidates)
            seekstart(fileprogress)
            write(fileprogress, progress[], acceptance[])
        end
    finally
        if verbose
            self_running[] = false
        end
        finalize(task)
        if !isnothing(filepath)
            close(fileout)
            close(fileprogress)
        end
    end
    if isnothing(filepath)
        return finish!(candidates)
    else
        return progress[], acceptance[]
    end
end

function execute(V, verbose, mons, nthreads::Integer, task, secondtask, filepath)
    if isone(nthreads)
        @assert(isnothing(secondtask))
        return execute(V, verbose, mons, task, filepath)
    end

    nv = nvariables(mons)
    num = length(mons)
    threadsize = div(num, nthreads, RoundUp)
    @verbose_info("Preparing to determine Newton polytope using ", nthreads, " threads, each checking about ", threadsize,
        " candidates")
    data_global = alloc_global(V, nv)
    # While we precalculate the size of the list exactly, we don't pre-allocate the output candidates - we hope to eliminate a
    # lot of exponents by the polytope containment, so we might overallocate so much memory that we hit a resource constraint
    # here. If instead, we grow the list dynamically, we pay the price in speed, but the impossible might become feasible.
    candidates = FastVec{UInt}()
    iterators = Vector{Base.promote_op(veciter, typeof(@view(mons[begin:end])), Val{true})}(undef, nthreads)
    let start=1
        for i in 1:nthreads
            @inbounds iterators[i] = veciter(@view(mons[start:min(num, start+threadsize-1)]), Val(true))
            start += threadsize
        end
    end

    threadprogress = zeros(Int, nthreads)
    threadacceptance = zeros(Int, nthreads)
    if !isnothing(filepath)
        fileprogresses = Vector{IOStream}(undef, nthreads)
        fileouts = Vector{IOStream}(undef, nthreads)
        @inbounds try
            for i in 1:nthreads
                fileprogresses[i] = fileprogress = open("$filepath-$i.prog", read=true, write=true, create=true, lock=false)
                fileouts[i] = open("$filepath-$i.out", append=true, lock=false)
                fs = filesize(fileprogress)
                if fs == 2sizeof(Int)
                    seekstart(fileprogress)
                    threadprogress[i] = prog = read(fileprogress, Int)
                    threadacceptance[i] = read(fileprogress, Int)
                    iterators[i] = veciter(@view(parent(iterators[i])[prog+1:end]), Val(true))
                elseif !iszero(fs)
                    error("Unknown progress file format - please delete existing files.")
                end
            end
        catch
            for i in 1:nthreads
                isassigned(fileprogresses, i) || break
                close(fileprogresses[i])
                isassigned(fileouts, i) && close(fileouts[i])
            end
            rethrow()
        end
    end
    cond = Threads.SpinLock()

    if verbose
        self_running = verbose_worker(num, sum(threadprogress, init=0), time_ns(), threadprogress, threadacceptance)
    end
    ccall(:jl_enter_threaded_region, Cvoid, ())
    try
        # To avoid naming confusion with Mosek's Task, we call the parallel Julia tasks threads.
        threads = Vector{Union{Task,Nothing}}(undef, nthreads)
        # We can already start all the tasks; this main task that must still feed the data will continue running until we yield
        # to the scheduler.
        @inbounds for tid in nthreads:-1:1
            if tid ≥ 3
                taskₜ = clonetask(task)
            elseif tid == 2
                taskₜ = secondtask
            else
                taskₜ = task
            end
            iter = iterators[tid]
            if isnothing(filepath)
                filestuff = nothing
            elseif isempty(iter)
                threads[tid] = nothing
                close(fileouts[tid])
                close(fileprogresses[tid])
                continue
            else
                filestuff = (fileprogresses[tid], fileouts[tid])
            end
            # secondtask has a solution, so we just use task (better than deletesolution).
            # We must create the copy in the main thread; Mosek will crash occasionally if the copies are created in parallel,
            # even if we make sure not to modify the base task until all copies are done.
            threads[tid] = Threads.@spawn execute_taskfun($V, $tid, $taskₜ, $iter, $nv, $data_global, $cond, $threadprogress,
                $threadacceptance, $candidates, $verbose, $filestuff)
        end
        for thread in threads
            isnothing(thread) || wait(thread)
        end
    finally
        if verbose
            self_running[] = false
        end
        ccall(:jl_exit_threaded_region, Cvoid, ())
    end
    @verbose_info("\33[2KAll tasks have finished, sorting the output")
    # We need to return the deglex monomial order, but due to the partitioning and multithreading, our output is unordered (not
    # completely, we have ordered slices of varying length, but does this help?).
    # TODO: Maybe defer appending to the main thread, then we can already do it in the correct order...
    if isnothing(filepath)
        return finish!(candidates)
    else
        return sum(threadprogress, init=0), sum(threadacceptance, init=0)
    end
end