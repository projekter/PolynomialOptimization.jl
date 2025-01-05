export Counters, addtocounter!, @counter_alias, @counter_atomic

"""
    Counters

Mutable struct that keeps track of the number of conic constraints/variables currently in use. Every possible symbol listed in
the documentation for [`extract_sos`](@ref) has an `Int` field of the same name.
"""
mutable struct Counters
    fix::Int
    nonnegative::Int
    quadratic::Int
    rotated_quadratic::Int
    psd::Int
    psd_complex::Int
    dd::Int
    dd_complex::Int
    lnorm::Int
    lnorm_complex::Int
    sdd::Int
    sdd_complex::Int

    Counters() = new(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
end

"""
    addtocounter!(state, counters::Counters, ::Val{type}, dim::Integer) where {type} -> Int
    addtocounter!(state, counters::Counters, ::Val{type}, num::Integer, dim::Integer) where {type} -> UnitRange{Int}

Increments the internal information about how many constraints/variables of a certain type were already added. Usually, a
solver implementation does not have to overwrite the default implementation; but it might be useful to do so if some types are
for example counted in a single counter or bunch of added elements counts as a single entry.
`dim` is the length of the conic constraint or variable, while `num` indicates that multiple such constraints or variables of
the same type are added.
For a list of possible symbols `type`, see the documentation for [`extract_sos`](@ref).

See also [`Counters`](@ref), [`@counter_alias`](@ref), [`@counter_atomic`](@ref).
"""
function addtocounter! end

@inline function addtocounter!(state::AbstractSolver, counters::Counters, ::Val{type}, dim::Integer) where {type}
    v = getproperty(counters, type)
    return v+1:setproperty!(counters, type, v + dim)
end

@inline function addtocounter!(state::AbstractSolver, counters::Counters, ::Val{type}, num::Integer, dim::Integer) where {type}
    v = getproperty(counters, type)
    return v+1:setproperty!(counters, type, v + num * dim)
end

function _parse_counter(counter)
    if counter isa QuoteNode
        return :(Val{$counter})
    elseif counter isa Symbol
        counter === :Any && return :(Val)
    else
        counter.head === :tuple && !isempty(counter.args) && all(x isa QuoteNode for x in counter.args) &&
            return :(Union{$((:(Val{$x}) for x in counter.args)...),})
    end
    throw(ArgumentError("@counter_atomic expects Symbol, tuple of Symbols, or `Any` as a second parameter, got $counter"))
end

"""
    @counter_atomic(::Type{<:AbstractSolver}, counter[, alias])

Defines the [`addtocounter!`](@ref) function in such a way that a multi-dimensional constraint of type `counter` is always
counted as a single entry.
This macro may be called at most once for each counter when defining a solver.
As a consequence, [`extract_sos`](@ref) will may have an `Int` or a `UnitRange{Int}` (if multiple constraints are required,
which will never be the case for matrix cones) as `index` parameter for the `counter` type.

If the `alias` parameter is present, whenever the type `counter` is encountered, the counter of type `alias` is incremented
instead. Do not use [`@counter_alias`](@ref) together with this macro on the same `counter`.

`counter` may either be a Symbol, a tuple of Symbols, or the value `Any`. Note that `Any` has weakest precedence, irrespective
of when the macro is called.
"""
macro counter_atomic(S, counter, alias::QuoteNode=counter)
    m = Base.parentmodule(addtocounter!)
    f = alias.value
    c = _parse_counter(counter)
    esc(quote
        @inline $m.addtocounter!(::$S, counters::$Counters, ::$c, dim::Integer) =
            counters.$f += 1
        @inline $m.addtocounter!(::$S, counters::$Counters, ::$c, num::Integer, dim::Integer) =
            counters.$f+1:(counters.$f += num)
    end)
end

"""
    @counter_alias(::Type{<:AbstractSolver}, counter, alias)

Defines the [`addtocounter!`](@ref) function in such a way that contraints of type `counter` instead only affect the counter
`alias`. Do not use [`@counter_atomic`](@ref) together with this macro on the same `counter`.

`counter` may either be a Symbol, a tuple of Symbols, or the value `Any`. Note that `Any` has weakest precendence, irrespective
of when the macro was called.

!!! warning
    Regardless of whether `counter` or `alias` where made atomic before, after this macro, `counter` will not be so (although
    it shares the same counter as `alias`). This may or may not be desirable (most likely not), so always make atomic counters
    explicit using [`@counter_atomic`](@ref), which allows the definition of aliases.
"""
macro counter_alias(S, counter, alias::QuoteNode)
    m = Base.parentmodule(addtocounter!)
    f = alias.value
    c = _parse_counter(counter)
    esc(quote
        @inline $m.addtocounter!(::$S, counters::$Counters, ::$c, dim::Integer) =
            counters.$f+1:(counters.$f += dim)
        @inline $m.addtocounter!(::$S, counters::$Counters, ::$c, num::Integer, dim::Integer) =
            counters.$f+1:(counters.$f += num * dim)
    end)
end