# This is an adaptation of StatsBase.sample! with replace=false, however specifically adapted to the case where we are
# sampling indices from `a` whose value is `true`, we will accumulate them into `x` and set them to `false` in `a`.
# `total` must be the number of `true`s occuring in `a`, else the function output will be undefined.
function sample!(a::AbstractVector{Bool}, x::AbstractVector{<:Integer}, total::Integer=sum(a, init=0))
    k = length(x)
    iszero(k) && return x
    0 < total ≤ length(a) || error("Invalid total number")
    k ≤ total || error("Cannot draw more samples without replacement.")
    if k == 1
        @inbounds x[1] = sample!(a, total)
    elseif k == 2
        @inbounds (x[1], x[2]) = samplepair!(a, total)
    elseif total < 24k
        fisher_yates_sample!(a, x, total)
    else
        self_avoid_sample!(a, x, total)
    end
    return x
end

sample(a::AbstractVector{Bool}, n::Integer, total::Integer=sum(a, init=0)) = sample!(a, Vector{Int}(undef, n), total)

function sample!(a::AbstractVector{Bool}, total::Integer)
    idx = rand(1:total)
    for (i, aᵢ) in enumerate(a)
        if aᵢ
            if !isone(idx)
                idx -= 1
            else
                @inbounds a[i] = false
                return i
            end
        end
    end
    error("No valid sample vector")
end

function samplepair!(a::AbstractVector{Bool}, total::Integer)
    idx1 = rand(1:total)
    idx2 = rand(1:total-1)
    if idx1 == idx2
        idx2 = total
    elseif idx1 > idx2
        idx1, idx2 = idx2, idx1
    end
    idx2 -= idx1
    i₁ = 0
    for (i, aᵢ) in enumerate(a)
        if aᵢ
            if !isone(idx1)
                idx1 -= 1
            elseif !iszero(i₁)
                if !isone(idx2)
                    idx2 -= 1
                else
                    @inbounds a[i] = false
                    return i₁, i
                end
            else
                @inbounds a[i] = false
                i₁ = i
            end
        end
    end
    error("No valid sample vector")
end

function fisher_yates_sample!(a::AbstractVector{Bool}, x::AbstractVector{<:Integer}, total::Integer)
    0 < total ≤ length(a) || error("Invalid total number")
    k = length(x)
    k ≤ total || error("length(x) should not exceed total")
    inds = let inds=FastVec{Int}(buffer=total)
        for (i, aᵢ) in enumerate(a)
            aᵢ && unsafe_push!(inds, i)
        end
        finish!(inds)
    end
    @inbounds for i in 1:k
        j = rand(i:total)
        t = inds[j]
        inds[j] = inds[i]
        inds[i] = t
        x[i] = t
        a[t] = false
    end
    return x
end

function self_avoid_sample!(a::AbstractVector{Bool}, x::AbstractVector{<:Integer}, total::Integer)
    0 < total ≤ length(a) || error("Invalid total number")
    k = length(x)
    k ≤ total || error("length(x) should not exceed total")
    curidx = findlast(a)::eltype(eachindex(a))
    currelidx = total
    @inbounds for i in 1:k
        newidx = rand(1:total)
        # newidx is the relative index; search from the current position until we find it
        while newidx < currelidx
            curidx -= 1
            if a[curidx]
                currelidx -= 1
            end
        end
        while newidx > currelidx
            curidx += 1
            if a[curidx]
                currelidx += 1
            end
        end
        # then set the output
        x[i] = curidx
        a[curidx] = false
        total -= 1
        iszero(total) && break
        # and make sure that currelidx indeed corresponds to curidx: for this, we must increase curidx unless we are at the
        # last item.
        if currelidx != total
            while !a[curidx]
                curidx += 1
            end
        else
            while !a[curidx]
                curidx -= 1
            end
            currelidx -= 1
        end
    end
    return x
end