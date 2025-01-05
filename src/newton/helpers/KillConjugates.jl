struct KillConjugates{Nr,P}
    parent::P
end

Base.IteratorSize(::Type{<:KillConjugates{<:Any,P}}) where {P} = Base.IteratorSize(P)
Base.IteratorEltype(::Type{<:KillConjugates{<:Any,P}}) where {P} = Base.HasEltype()
Base.length(k::KillConjugates) = length(k.parent)
Base.eltype(::Type{<:KillConjugates{<:Any,P}}) where {P} = eltype(P)
function Base.iterate(k::KillConjugates{Nr}) where {Nr}
    result = iterate(k.parent)
    isnothing(result) && return nothing
    return result[1], (Nr, result[2])
end
function Base.iterate(k::KillConjugates, (take, state))
    result = iterate(k.parent, state)
    isnothing(result) && return nothing
    if iszero(take)
        return zero(result[1]), (one(take), state)
    else
        return result[1], (take - one(take), state)
    end
end