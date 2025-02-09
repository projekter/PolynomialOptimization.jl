"""
    preproc(V, mons::SimpleMonomialVector, vertexindices, verbose::Bool, singlethread::Bool;
        parameters...) -> AbstractVector{Bool}

Checks for convex dependencies of exponents from the list of monomials.
`vertexindices` might either be a `Vector{Int}` that indexes `mons` or is it `Val(:all)`.
In the former case, the convex polytope is spanned by the exponents of all monomials in `mons[vertexindices]` (and it can be
assumed that those are independent, and that `vertexindices` is fairly short compared to `mons`); in the latter case, the
convex polytope is potentially spanned by the exponents of all monomials in `mons`.
The implementation has to return an `AbstractVector{Bool}` that for every monomial in `mons` indicates whether the convex hull
spanned by the elements just described (excluding the monomial in question) already contains the monomial. An entry must be
`false` if and only if this is possible, i.e., it is redundant.
`singlethread` will be `true` if this function will be called in parallel by multiple threads, so that the linear solver itself
should be single-threaded.
"""
function preproc end

"""
    prepare(V, vertices::SimpleMonomialVector, num::Int, verbose::Bool; parameters...) ->
        (nthreads::Int, userdata, userdata_clone)

This function is responsible for creating an optimization task that can be used to check membership in the Newton halfpolytope.
The vertices of the polytope are given by the exponents of `mons`. The total number of monomials that have to be checked is
given by `num`.
The function must return the number of threads that will be used to carry out the optimization (which is not the number of
threads that the optimizer uses internally, but determines how `PolynomialOptimization` will distribute the jobs), an internal
optimization task that is passed on to [`work`](@ref), and a copy of this task for use in a second thread if the number of
threads is greater than one, else `nothing`. More copies will be created as required by [`clonetask`](@ref); however, assuming
that setting up the task will potentially require more resources than cloning a task allows the function to estimate the
required memory (and therefore a sensible number of threads) better by already performing one clone.

See also [`@allocdiff`](@ref var"@allocdiff").
"""
function prepare end

"""
    alloc_global(V, nv)

This function is called once in the main thread before [`work`](@ref) is executed (which, due to multithreading, might occur
more than once). It can create some shared data that is used in a read-only manner by all workers at the same time.
The default implementation of this function does nothing.
"""
alloc_global(_, _) = nothing

"""
    alloc_local(V, nv)

This function is called once in every computation thread before [`work`](@ref) is executed (which, due to task splitting,
might occur more than once). It can create some shared data that is available for reading and writing by every worker.
The default implementation of this function does nothing.
"""
alloc_local(_, _) = nothing

"""
    clonetask(t)

This function must create a copy of the optimization task `t` that can run in parallel to `t` in a different thread.
"""
function clonetask end

"""
    work(V, task, data_global, data_local, expiter::SimpleMonomialVectorIterator{true},
        Δprogress::Ref{Int}, Δacceptance::Ref{Int}, add_callback::Function,
        iteration_callback::Union{Nothing,Function})

Iterates through `expiter` (which gives a tuple consisting of the index and the exponents) and for every exponent checks
whether this it can be obtained by a convex combination of the coefficients as they are set up in `task`. If yes,
`add_callback` must be called with the exponent index as a parameter, and `Δacceptance` should be incremented. In any case,
`Δprogress` should be incremented. Additionally, after every check, `iteration_callback` should be called with no parameters,
if it is a function.
The `data` parameters contain the custom data that was previously generated using [`alloc_global`](@ref) and
[`alloc_local`](@ref). Only `data_local` may be mutated.
"""
function work end
