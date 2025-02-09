function preproc_prequick(V, mons::SimpleMonomialVector{Nr,0}, verbose; parameters...) where {Nr}
    @inbounds begin
        vertexindices = fill(1, 2Nr)
        lowestidx = @view(vertexindices[1:Nr])
        highestidx = @view(vertexindices[Nr+1:2Nr])
        lowestvals = collect(exponents(mons[begin]))
        highestvals = copy(lowestvals)
        # we might also add the points with the smallest/largest sum of all coordinates, or differences (but there are 2^nv
        # ways to combine, so let's skip it)
        for (j, exps) in enumerate(veciter(mons))
            for (i, expᵢ) in enumerate(exps)
                if lowestvals[i] > expᵢ
                    lowestidx[i] = j
                    lowestvals[i] = expᵢ
                end
                if highestvals[i] < expᵢ
                    highestidx[i] = j
                    highestvals[i] = expᵢ
                end
            end
        end
        Base._groupedunique!(sort!(vertexindices))
        return preproc(V, mons, vertexindices, false, verbose; parameters...)
    end
end

function preproc_randomized_taskfun(V, mons, required_exps, subset, done, event, stop; parameters...)
    @inbounds while !stop[]
        dropped = 0
        for (cfi, rem) in zip(subset, preproc(V, @view(mons[subset]), Val(:all), false, true; parameters...))
            required_exps[cfi] = rem
            if !rem
                dropped += 1
            end
        end
        put!(done, dropped)
        stop[] && return
        wait(event)
    end
end

function preproc_randomized(V, mons, verbose; parameters...)
    nc = length(mons)
    nthreads = Threads.nthreads()
    subset_size = min(1000, nc ÷ 20) # samples too small -> no chance of success; samples too large -> takes too long
    @assert(subset_size ≥ 1)
    if nthreads * subset_size > nc
        nthreads = nc ÷ subset_size # we could squeeze in another incomplete thread, but this would make the code more
        # complicated than necessary without a real benefit
    end
    if isone(nthreads)
        required_exps = fill!(BitVector(undef, nc), true)
        let
            _subset = sample(required_exps, subset_size, nc)
            @inbounds while true
                totaldrop = 0
                for (cfi, rem) in zip(_subset, preproc(V, @view(mons[_subset]), Val(:all), false, false); parameters...)
                    required_exps[cfi] = rem
                    if !rem
                        totaldrop += 1
                    end
                end
                nc -= totaldrop
                if verbose
                    print("\33[2KStatus update: ", nc, " remaining extremal points, last drop was ", totaldrop, "\r")
                    flush(stdout)
                end
                (totaldrop < 20 || 10totaldrop < subset_size) && break
                subset_size > nc && break
                sample!(required_exps, _subset, nc)
            end
        end
    else
        required_exps = fill(true, nc) # BitVector is not thread-safe; and we cannot guarantee that the samples are chosen
                                       # with nonoverlapping regions of size 64
        ccall(:jl_enter_threaded_region, Cvoid, ())
        try
            # We divide our variables in randomized, disjoint (in the case of multi-threading) subsets and try to find linear
            # dependencies in them. Depending on our success rate, we'll try to increase the size of the subsets.
            allsubs = sample(required_exps, nthreads * subset_size, nc)
            done = Base.Channel{Int}(1)
            events = [Base.Event(true) for _ in 1:nthreads]
            stop = Ref(false)
            # Initialize all threads; as everything is already set up, they will directly start
            threads = Vector{Task}(undef, nthreads)
            for (tid, start) in enumerate(1:subset_size:subset_size*nthreads)
                @inbounds threads[tid] = Threads.@spawn preproc_randomized_taskfun($V, $mons, $required_exps,
                    $(@view(allsubs[start:start+subset_size-1])), $done, $(events[tid]), $stop; $parameters...)
            end
            while true
                # We wait for all threads to finish - if we have a sufficiently good coverage of the whole set with all
                # available threads, it would make no sense to directly re-start a thread, it would just get the same set
                # again. Instead, when all threads are finished (which happens approximately at the same time), we can rehash
                # the full set.
                totaldrop = 0
                for _ in 1:nthreads
                    totaldrop += take!(done)
                end
                nc -= totaldrop
                if verbose
                    print("\33[2KStatus update: ", nc, " remaining extremal points, last drop was ", totaldrop, "\r")
                    flush(stdout)
                end
                totalsize = nthreads * subset_size
                (totaldrop < 20 || 10totaldrop < totalsize) && break
                if totalsize > nc
                    nthreads = nc ÷ subset_size
                    iszero(nthreads) && break
                    stop[] = true
                    for event in events
                        notify(event)
                    end
                    totalsize = nthreads * subset_size
                end
                subs = @view(allsubs[1:totalsize])
                sample!(required_exps, subs, nc)
                for event in @view(events[1:nthreads])
                    notify(event) # don't start those beyond nthreads
                end
            end
            stop[] = true
            for event in events
                notify(event) # to properly finish all tasks
            end
        finally
            ccall(:jl_exit_threaded_region, Cvoid, ())
        end
    end
    return required_exps
end

function preproc(V, objective::SimplePolynomial{<:Any,Nr,0}; verbose::Bool=false,
    zero::AbstractVector{<:SimplePolynomial{<:Any,Nr,0}}, nonneg::AbstractVector{<:SimplePolynomial{<:Any,Nr,0}},
    psd::AbstractVector{<:AbstractMatrix{<:SimplePolynomial{<:Any,Nr,0}}}, groupings::RelaxationGroupings{Nr,0},
    preprocess::Union{Nothing,Bool}=nothing, preprocess_quick::Bool=true, preprocess_randomized::Bool=false,
    preprocess_fine::Bool=false, warn_disable_randomization::Bool=true, parameters...) where {Nr}
    if !isnothing(preprocess)
        preprocess_quick = preprocess_randomized = preprocess_fine = preprocess
    end
    @verbose_info("Determining Newton polytope (quick preprocessing: ", preprocess_quick, ", randomized preprocessing: ",
        preprocess_randomized, ", fine preprocessing: ", preprocess_fine, ")")
    mons = merge_constraints(objective, zero, nonneg, psd, groupings, verbose,
        preprocess_quick | preprocess_randomized | preprocess_fine)
    if preprocess_quick
        @verbose_info("Removing redundancies from the convex hull - quick heuristic, ", length(mons), " initial candidates")
        preproc_time = @elapsed begin
            mons = SimplePolynomials.keepat!!(mons, preproc_prequick(V, mons, verbose; parameters...))
        end
        @verbose_info("Found ", length(mons), " potential extremal points of the convex hull in ", preproc_time, " seconds")
    end
    if preprocess_randomized
        if length(mons) ≥ 100
            @verbose_info("Removing redundancies from the convex hull - randomized, ", length(mons), " initial candidates")
            preproc_time = @elapsed begin
                mons = SimplePolynomials.keepat!!(mons, preproc_randomized(V, mons, verbose; parameters...))
            end
            @verbose_info("Found ", length(mons), " extremal points of the convex hull via randomization in ", preproc_time, " seconds")
        else
            warn_disable_randomization &&
                @info("Removing redundancies from the convex hull via randomization was requested, but skipped due to the small size of the problem")
        end
    end
    if preprocess_fine
        # eliminate all the coefficients that are redundant themselves to make the linear system smaller
        @verbose_info("Removing redundancies from the convex hull - fine, ", length(mons), " initial candidates")
        preproc_time = @elapsed begin
            mons = SimplePolynomials.keepat!!(mons, preproc(V, mons, Val(:all), verbose, false; parameters...))
        end
        @verbose_info("Found ", length(mons), " extremal points of the convex hull in ", preproc_time, " seconds")
    end
    return parameters, mons
end

function analyze(mons)
    # This is some quick preprocessing to further restrict the potential degrees (actually, this is what SumOfSquares.jl calls
    # Newton polytope)
    nv = nvariables(mons)

    maxdeg, mindeg = zero(Int), typemax(Int)
    maxmultideg, minmultideg = fill(zero(Int), nv), fill(typemax(Int), nv)
    # We expect mons to be an indexed monomial vector (as it comes from a user-defined polynomial), so optimizations that check
    # whether a whole exponent range is covered are most likely useless.
    for exps in veciter(mons)
        deg = 0
        @inbounds for (i, expᵢ) in enumerate(exps)
            deg += expᵢ
            if expᵢ > maxmultideg[i]
                maxmultideg[i] = expᵢ
            end
            if expᵢ < minmultideg[i]
                minmultideg[i] = expᵢ
            end
        end
        if deg > maxdeg
            maxdeg = deg
        end
        if deg < mindeg
            mindeg = deg
        end
    end
    maxmultideg .= div.(maxmultideg, 2, RoundDown)
    minmultideg .= div.(minmultideg, 2, RoundUp)
    maxdeg = div(maxdeg, 2, RoundDown)
    mindeg = div(mindeg, 2, RoundUp)

    # let's just always go for UInt as the index type. We need to feed this to a linear solver after all, and we cannot expect
    # a linear solver to cope with more than typemax(UInt) constraints.
    return mindeg, maxdeg, minmultideg, maxmultideg
end