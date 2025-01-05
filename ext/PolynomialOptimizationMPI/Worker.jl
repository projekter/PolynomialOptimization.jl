function verbose_worker(num, (otherprogress, otheracceptance), init_time, comm, ::RootRank, threadprogress, threadacceptance)
    init_progress = otherprogress + sum(threadprogress, init=0)
    self_running = Ref(true)
    # We take care of the gathering of information in a separate task that must run on the main thread (due to MPI funneling).
    # This requires the other tasks (and the main one) to yield regularly.
    workers = Ref(MPI.Comm_size(comm) -1)
    progress_acceptances = [Vector{Int}(undef, 2) for _ in 1:workers[]]
    buffers = MPI.Buffer.(progress_acceptances)
    reqs = MPI.MultiRequest(workers[])
    for (worker, (req, Δbuffer)) in enumerate(zip(reqs, buffers))
        MPI.Irecv!(Δbuffer, comm, req, tag=1, source=worker)
    end
    curotherprogress = Ref(otherprogress)
    curotheracceptance = Ref(otheracceptance)
    check = @capture(print -> let
        δprogress = 0
        δacceptance = 0
        # We need to give MPI some time to collect all the data; but we don't want to completely block using Waitall for
        # performance reasons. So we wait at most 10ms (= 1% of the calucation time) before we interrupt and resume the
        # calculations.
        # Note that the strange scheme where we receive on all nodes seems to be the only way not to completely loose some
        # packages. If we just receive without any binding to a source, this surprisingly is very unreliable.
        init_wait = time_ns()
        @inbounds while time_ns() - init_wait < 10_000_000
            _, idx = MPI.Testany($reqs)
            if !isnothing(idx)
                progress_acceptance = $progress_acceptances[idx]
                if progress_acceptance[1] < 0
                    $workers[] -= 1
                else
                    δprogress += progress_acceptance[1]
                    δacceptance += progress_acceptance[2]
                    MPI.Irecv!($buffers[idx], $comm, reqs[idx], tag=1, source=idx)
                end
            end
        end
        $curotherprogress[] += δprogress
        $curotheracceptance[] += δacceptance
        print || return
        allprogress = sum(threadprogress, init=0) + curotherprogress[]
        allacceptance = sum(threadacceptance, init=0) + curotheracceptance[]
        #println(string("Root: others ", curotherprogress[], ", this ", threadprogress, ", total ", allprogress))
        # no need to lock, we only modify this here
        iszero(allprogress) && (allprogress = 1) # we divide, so if this task runs too early, we might be in trouble
        Δt = allprogress == $init_progress ? 1 : allprogress - $init_progress # if a finished job is started, this might happen
        rem_sec = round(Int, ((time_ns() - $init_time) / 1_000_000_000Δt) * ($num - allprogress))
        @printf("\33[2KStatus update: %.2f%%, acceptance: %.2f%%, remaining time: %02d:%02dmin\r",
            100allprogress / num, 100allacceptance / allprogress, rem_sec ÷ 60, rem_sec % 60)
        flush(stdout)
    end)
    Threads.@spawn(:interactive, begin
        while !iszero($workers[]) || $self_running[]
            $check(true)
            sleep(1)
        end
        $check(false) # make sure to empty all messages, but no need to print again, we are done anyway
    end)
    return self_running
end

function verbose_worker(num, (initprogress, initacceptance), init_time, comm, rank::OtherRank, threadprogress, threadacceptance)
    self_running = Ref(true)
    progress_acceptance = Array{Int}(undef, 2)
    lastprogress = Ref(initprogress)
    lastacceptance = Ref(initacceptance)
    buffer = MPI.Buffer(progress_acceptance)
    # We take care of the gathering of information in a separate task that must run on the main thread (due to MPI funneling).
    # This requires the other tasks (and the main one) to yield regularly.
    check = @capture(() -> begin
        newprogress = sum($threadprogress, init=0)
        @inbounds $progress_acceptance[1] = newprogress - $lastprogress[]
        lastprogress[] = newprogress
        newacceptance = sum($threadacceptance, init=0)
        @inbounds progress_acceptance[2] = newacceptance - $lastacceptance[]
        lastacceptance[] = newacceptance
        #println(string("Worker ", convert(Int, rank), ": ", threadprogress, " delta ", progress_acceptance[1]))
        MPI.Send($buffer, $comm, dest=root, tag=1)
    end)
    Threads.@spawn(:interactive, begin
        while $self_running[]
            $check()
            sleep(1)
        end
        $check() # make sure to have a final report
        # and then signal that we are done
        $progress_acceptance[1] = -1
        MPI.Send($buffer, $comm, dest=root, tag=1)
    end)
    return self_running
end

function execute(V, verbose, mons, task, filepath, comm, rank::MPIRank)
    nv = nvariables(mons)
    num = length(mons)
    nworkers = MPI.Comm_size(comm)
    workersize = div(num, nworkers, RoundUp)
    isroot(rank) && @verbose_info("Preparing to determine Newton polytope using ", nworkers,
        " single-threaded workers, each checking about ", workersize, " candidates")

    # While we precalculate the size of the list exactly, we don't pre-allocate the output candidates - we hope to eliminate a
    # lot of exponents by the polytope containment, so we might overallocate so much memory that we hit a resource constraint
    # here. If instead, we grow the list dynamically, we pay the price in speed, but the impossible might become feasible.
    candidates = FastVec{UInt}()

    workerprogress = Ref(0)
    workeracceptance = Ref(0)
    if isroot(rank)
        otherprogress = otheracceptance = 0
    end
    if isnothing(filepath)
        filestuff = nothing
        iter = let start=convert(Int, rank) * workersize +1
            veciter(@view(mons[start:min(num, start+workersize-1)]), Val(true))
        end
    else
        prefix = "$filepath-n$(convert(Int, rank))"
        fileprogress = open("$prefix.prog", read=true, write=true, create=true, lock=false)
        success = fill(false, nworkers)
        try
            fs = filesize(fileprogress)
            if fs == 2sizeof(Int)
                workerprogress[] = read(fileprogress, Int)
                workeracceptance[] = read(fileprogress, Int)
            elseif !iszero(fs)
                error()
            end
        catch
            close(fileprogress)
            MPI.Allgather!(MPI.UBuffer(success, 1), comm)
            isroot(rank) && error("Unknown progress file format - please delete existing files.")
            return
        end
        @inbounds success[convert(Int, rank)+1] = true
        MPI.Allgather!(MPI.UBuffer(success, 1), comm)
        if !all(success)
            close(fileprogress)
            isroot(rank) && error("Unknown progress file format - please delete existing files.")
            return
        end
        iter = let start=convert(Int, rank) * workersize +1
            veciter(@view(mons[start+workerprogress[]:min(num, start+workersize-1)]), Val(true))
        end
        if isroot(rank)
            worker_pas = Vector{Int}(undef, 2nworkers)
            MPI.Gather!(MPI.IN_PLACE, MPI.UBuffer(worker_pas, 2), comm; root)
            otherprogress = sum(@view(worker_pas[3:2:end]), init=0)
            otheracceptance = sum(@view(worker_pas[4:2:end]), init=0)
        else
            MPI.Gather!([workerprogress[], workeracceptance[]], nothing, comm; root)
        end
        fileout = open("$prefix.out", append=true, lock=false)
        filestuff = (fileprogress, fileout)
    end
    if verbose
        self_running = verbose_worker(num,
            isroot(rank) ? (otherprogress, otheracceptance) : (workerprogress[], workeracceptance[]), time_ns(), comm, rank,
            workerprogress, workeracceptance)
    end
    lastinfo = Ref(time_ns())
    cb = step_callback(verbose, filestuff, lastinfo, candidates, workerprogress, workeracceptance)
    try
        work(V, task, alloc_global(V, nv), alloc_local(V, nv), iter, workerprogress, workeracceptance,
            @capture(idx -> push!($candidates, idx)), cb)
        if !isnothing(filepath)
            write(fileout, candidates)
            seekstart(fileprogress)
            write(fileprogress, workerprogress[], workeracceptance[])
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
    return isnothing(filepath) ? finish!(candidates) : nothing
end

function execute(V, verbose, mons, nthreads::Integer, task, secondtask, filepath, comm, rank::MPIRank)
    if isone(nthreads)
        @assert(isnothing(secondtask))
        return execute(V, verbose, mons, task, filepath, comm, rank)
    end

    nv = nvariables(mons)
    num = length(mons)
    nworkers = MPI.Comm_size(comm)
    workersize = div(num, nworkers, RoundUp)
    threadsize = div(workersize, nthreads, RoundUp)
    if threadsize == workersize
        finalize(secondtask)
        return execute(V, verbose, mons, task, filepath, comm, rank)
    end
    isroot(rank) && @verbose_info("Preparing to determine Newton polytope using ", nworkers, " workers, each with ", nthreads,
        " threads checking about ", threadsize, " candidates")
    data_global = alloc_global(V, nv)
    # While we precalculate the size of the list exactly, we don't pre-allocate the output candidates - we hope to eliminate a
    # lot of exponents by the polytope containment, so we might overallocate so much memory that we hit a resource constraint
    # here. If instead, we grow the list dynamically, we pay the price in speed, but the impossible might become feasible.
    candidates = FastVec{UInt}()
    iterators = Vector{Base.promote_op(veciter, typeof(@view(mons[begin:end])), Val{true})}(undef, nthreads)
    let start=convert(Int, rank) * workersize +1, stop=min(num, start+workersize-1)
        for i in 1:nthreads
            @inbounds iterators[i] = veciter(@view(mons[start:min(stop, start+threadsize-1)]), Val(true))
            start += threadsize
        end
    end

    threadprogress = zeros(Int, nthreads)
    threadacceptance = zeros(Int, nthreads)
    otherprogress = otheracceptance = 0
    if !isnothing(filepath)
        prefix = "$filepath-n$(convert(Int, rank))"
        fileprogresses = Vector{IOStream}(undef, nthreads)
        fileouts = Vector{IOStream}(undef, nthreads)
        success = fill(false, nworkers)
        @inbounds try
            for i in 1:nthreads
                fileprogresses[i] = fileprogress = open("$prefix-$i.prog", read=true, write=true, create=true, lock=false)
                fileouts[i] = open("$prefix-$i.out", append=true, lock=false)
                fs = filesize(fileprogress)
                if fs == 2sizeof(Int)
                    threadprogress[i] = prog = read(fileprogress, Int)
                    threadacceptance[i] = read(fileprogress, Int)
                    iterators[i] = veciter(@view(parent(iterators[i])[prog+1:end]), Val(true))
                elseif !iszero(fs)
                    error()
                end
            end
        catch
            for i in 1:nthreads
                isassigned(fileprogresses, i) || break
                close(fileprogresses[i])
                isassigned(fileouts, i) && close(fileouts[i])
            end
            MPI.Allgather!(MPI.UBuffer(success, 1), comm)
            isroot(rank) && error("Unknown progress file format - please delete existing files.")
            return
        end
        @inbounds success[convert(Int, rank)+1] = true
        MPI.Allgather!(MPI.UBuffer(success, 1), comm)
        if !all(success)
            for (fp, fo) in zip(fileprogresses, fileouts)
                close(fp)
                close(fo)
            end
            isroot(rank) && error("Unknown progress file format - please delete existing files.")
            return
        end
        if isroot(rank)
            worker_pas = Vector{Int}(undef, 2nworkers)
            MPI.Gather!(MPI.IN_PLACE, MPI.UBuffer(worker_pas, 2), comm; root)
            otherprogress = sum(@view(worker_pas[3:2:end]), init=0)
            otheracceptance = sum(@view(worker_pas[4:2:end]), init=0)
        else
            otherprogress = sum(threadprogress, init=0)
            otheracceptance = sum(threadacceptance, init=0)
            MPI.Gather!([otherprogress, otheracceptance], nothing, comm; root)
        end
    end
    cond = Threads.SpinLock()

    if verbose
        self_running = verbose_worker(num, (otherprogress, otheracceptance), time_ns(), comm, rank, threadprogress,
            threadacceptance)
    end
    ccall(:jl_enter_threaded_region, Cvoid, ())
    try
        # To avoid naming confusion with Mosek's Task, we call the parallel Julia tasks threads.
        threads = Vector{Union{Task,Nothing}}(undef, nthreads)
        # We can already start all the tasks; this main task that must still feed the data will continue running until we
        # yield to the scheduler.
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
            # We must create the copy in the main thread; Mosek will crash occasionally if the copies are created in
            # parallel, even if we make sure not to modify the base task until all copies are done.
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
    return isnothing(filepath) ? finish!(candidates) : nothing
end