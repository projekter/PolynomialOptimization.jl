```@meta
CurrentModule = PolynomialOptimization
```

# Reference of auxilliaries
This reference page lists some additional functions and type that are required for the package to work properly. They are not
part of the public API and are not exported; however, they may be useful even when separated from the purpose of polynomial
optimization.

```@docs
@allocdiff
keepcol!
issubset_sorted
```

## Chordal graphs
The following functions provide a thin wrapper to [`CliqueTrees.jl`](https://github.com/AlgebraicJulia/CliqueTrees.jl).
```@docs
Relaxation.chordal_cliques!
Relaxation.chordal_cliques
```

## Sorting of multiple vectors
To simplify the task of sorting one vector and at the same time sorting multiple other vectors according to the same order - a
common task that usually requires first computing a `sortperm` and then indexing all vectors with this permutation -
`PolynomialOptimization` provides a helper function.
```@docs
sort_along!
```

## FastVector
To improve the speed in some implementation details, `PolynomialOptimization` provides a "fast" vector type. This is basically
just a wrapper around the stdlib `Vector`; however, it actually takes proper advantage of sizehints. The fact that Julia does
this badly has been known for quite some time ([#24909](https://github.com/JuliaLang/julia/issues/24909)), but the default
implementation has not changed (this will be different starting from Julia 1.11). Our own `FastVec` is a bit more specific than
the [PushVector](https://github.com/tpapp/PushVectors.jl), but also allows for more aggressive optimization.
`FastVec` can directly be used in `ccall` with a pointer element type.
```@meta
CurrentModule = PolynomialOptimization.FastVector
```
```@docs
FastVec
Base.sizehint!(::FastVec, ::Integer)
Base.resize!(::FastVec, ::Integer)
Base.empty!(::FastVec)
prepare_push!
Base.push!(::FastVec{V}, ::Any) where {V}
unsafe_push!
Base.insert!(::FastVec{V}, ::Integer, ::Any) where {V}
unsafe_insert!
Base.append!(::FastVec, ::Any)
unsafe_append!
Base.prepend!(::FastVec, ::Any)
unsafe_prepend!
Base.splice!(::FastVec, ::Integer, ::Any)
Base.similar(::FastVec{V}) where {V}
Base.copyto!(::FastVec, ::Integer, ::FastVec, ::Integer, ::Integer)
Base.deleteat!(::FastVec, ::Integer)
finish!
```
