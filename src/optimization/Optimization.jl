export poly_optimize, optimality_certificate, RepresentationPSD, RepresentationDD, RepresentationSDD, RepresentationIAs,
    IterateRepresentation

include("./Result.jl")
include("./MomentMatrix.jl")
include("./OptimalityCertificate.jl")
include("./solver/Solver.jl")
include("./SOSCertificate.jl")
using .Solver: default_solver_method, monomial_count, RepresentationMethod, RepresentationPSD, RepresentationDD,
    RepresentationSDD, RepresentationIAs
import .Solver: poly_optimize
using StandardPacked: SPMatrix, SPMatrixUpper, SPMatrixLower
import PositiveFactorizations
import LinearAlgebra: UpperOrUnitUpperTriangular, LowerOrUnitLowerTriangular

_val_of(::Val{S}) where {S} = S

"""
    poly_optimize([method, ]relaxation::AbstractRelaxation; verbose=false,
        representation=RepresentationPSD(), kwargs...)

Optimize a relaxed polynomial optimization problem that was construced via [`poly_problem`](@ref) and then wrapped into an
[`AbstractRelaxation`](@ref). Returns a [`Result`](@ref) object.

Instead of modeling the moment/SOS matrices as positive semidefinite, other representations such as the (scaled) diagonally
dominant description are also possible. The `representation` parameter can be used to define a representation that is employed
for the individual groupings. This may either be an instance of a [`RepresentationMethod`](@ref) - which requires the method to
be independent of the dimension of the grouping - or a callable. In the latter case, it will be passed as a first parameter an
identifier[^1] of the current conic variable, and as a second parameter the side dimension of its matrix. The method must then
return a [`RepresentationMethod`](@ref) instance.

`verbose=true` will enable logging; this will print basic information about the relaxation itself as well as instruct the
solver to output a detailed log. The PSD block sizes reported accurately represent the side dimensions of semidefinite
variables and how many of these variables appear. The free block sizes are only very loose upper bounds on the maximal number
of equality constraints that will be constructed by multiplying two elements from a block, as duplicates will be ignored.
Any additional keyword argument is passed on to the solver.

For a list of supported methods, see [the solver reference](@ref solvers_poly_optimize). If `method` is omitted, the default
solver is used. Note that this depends on the loaded solver packages, and possibly also their loading order if no preferred
solver has been loaded.

See also [`RepresentationIAs`](@ref).

[^1]: This identifier will be a tuple, where the first element is a symbol - either `:objective`, `:nonneg`, or `:psd` - to
      indicate the general reason why the variable is there. The second element is an `Int` denoting the index of the
      constraint (and will be undefined for the objective, but still present to avoid extra compilation). The last element
      is an `Int` denoting the index of the grouping within the constraint/objective.
"""
function poly_optimize(@nospecialize(v::Val), relaxation::AbstractRelaxation; verbose::Bool=false,
    representation=RepresentationPSD(), kwargs...)
    otime = @elapsed begin
        @verbose_info("Beginning optimization...")
        groups = groupings(relaxation) # This is instantaneous, as the groupings were already calculated when the relaxation
                                       # was constructed.
        if verbose
            bs = StatsBase.countmap(length.(groups.obj))
            for constr in groups.nonnegs
                mergewith!(+, bs, StatsBase.countmap(length.(constr)))
            end
            for (constr, constr_mat) in zip(groups.psds, poly_problem(relaxation).constr_psd)
                mergewith!(+, bs, StatsBase.countmap(length.(constr) .* size(constr_mat, 1)))
            end
            print("PSD block sizes:\n  ", sort!(collect(bs), rev=true))
            if !isempty(groups.zeros)
                empty!(bs)
                for constr in groups.zeros
                    mergewith!(+, bs, StatsBase.countmap(length.(constr)))
                end
                print("\nFree block sizes:\n  ", sort!(collect(bs), rev=true))
            end
            println("\nStarting solver...")
        end
        result = poly_optimize(v, relaxation, groups; verbose, representation, kwargs...)
    end
    return Result(relaxation, _val_of(v), otime, result...)
end

"""
    poly_optimize([method, ]problem::Problem[, degree::Int]; kwargs...)

Construct a [`Relaxation.Dense`](@ref) by default.
"""
poly_optimize(v::Val, problem::Problem, rest...; kwargs...) =
    poly_optimize(v, Relaxation.Dense(problem, rest...); kwargs...)

poly_optimize(s::Symbol, rest...; kwrest...) = poly_optimize(Val(s), rest...; kwrest...)

function poly_optimize(args...; kwargs...)
    if !isempty(args) && args[1] isa Val
        error("Unknown solver method specified. Are the required solver packages loaded?")
    end
    method = default_solver_method()
    @info("No solver method specified: choosing $method")
    poly_optimize(Val(method), args...; kwargs...)
end

struct IterateRepresentation
    keep_structure::Bool

    @doc """
        IterateRepresentation(; keep_structure=false)

    Default iteration method for DD and SDD representations. This is will perform a Cholesky decomposition of the old SOS
    matrix and use it as the new rotation, ensuring that results never get worse (at least in theory; since a positive definite
    SOS matrix is only guaranteed up to a certain tolerance, bad things could still happen).

    Note that the resulting rotation matrix will be upper triangular, which may break a previous structure. By setting
    `keep_structure` to `true`, the structure will be preserved (if it was diagonal, this would mean keeping only the diagonal
    of the Cholesky factor, with no theoretical guarantees, not even about convergence; if it was lower triangular the adjoint
    will be taken, which will _not_ give any convergence guarantees, as the rotated DD/SDD cone is implemented with respect to
    the upper triangular factorization).

    See also [`poly_optimize`](@ref poly_optimize(::Result))
    """
    IterateRepresentation(; keep_structure::Bool=false) = new(keep_structure)
end

function (i::IterateRepresentation)((type, index, grouping), _, oldrep::Type{<:RepresentationMethod}, oldsos)
    oldrep <: RepresentationPSD && return RepresentationPSD()
    # We need to perfom a Cholesky decomposition in the form oldsos = Uᵀ * U, where U is upper triangular. However, oldsos
    # will have tiny negative eigenvalues as numerical artifacts of the solver. Our regular BLAS cholesky with disabled checks
    # may yield results that are good for improvement nevertheless; but it might instead also stall or even be detrimental.
    # So we'll use are more stable version that is still quite efficient, which is the implementation in
    # PositiveFactorizations. While success is not guaranteed any more, as Uᵀ * U might not hold and therefore the convergence
    # proof fails, it is still quite likely. And in fact, nothing that we could do (even truncating the negative eigenvalues
    # using an eigendecomposition, then Choleskying, very expensive...) would give this guarantee, so why not.
    # To have it actually efficient, we unwrap the symmetric structure, although it would work on any AbstractMatrix.
    if oldsos isa (Symmetric{T,Matrix{T}} where {T<:Real}) || oldsos isa (Hermitian{T,Matrix{T}} where {T})
        fullm = LinearAlgebra.copytri!(parent(oldsos), oldsos.uplo, true)
        inplace = false
    elseif oldsos isa Matrix
        fullm = oldsos
        inplace = false
    elseif oldsos isa SPMatrixUpper
        fullm = LinearAlgebra.copytri!(convert(Matrix{eltype(oldsos)}, oldsos), 'U', true)
        inplace = true
    elseif oldsos isa SPMatrixLower
        fullm = LinearAlgebra.copytri!(convert(Matrix{eltype(oldsos)}, oldsos), 'L', true)
        inplace = true
    else
        fullm = Matrix(oldsos)::Matrix{eltype(oldsos)}
        inplace = true
    end
    newrot = (inplace ? cholesky! : cholesky)(PositiveFactorizations.Positive, fullm).U
    (!i.keep_structure || oldrep <: RepresentationMethod{<:UpperOrUnitUpperTriangular}) && return oldrep(newrot)
    oldrep <: RepresentationMethod{<:Union{<:UniformScaling,<:Diagonal}} && return oldrep(Diagonal(newrot)) # bad
    oldrep <: RepresentationMethod{<:LowerOrUnitLowerTriangular} && return oldrep(newrot') # also bad!
    dense = parent(newrot)
    @inbounds for j in axes(dense, 2)
        fill!(@view(dense[j+1:size(dense, 1)]), zero(eltype(dense)))
    end
    return oldrep(dense)
end

"""
    poly_optimize(result::Result; [representation=IterateRepresentation(), ]kwargs...)

Re-optimizes a previously optimized polynomial optimization problem. This is usually pointless, as the employed optimizers will
find globally optimal solutions. However, this method allows to change the representation used for the constraints (or
objective). If `representation` is a callable, it will now receive as a third parameter the type of the
[`RepresentationMethod`](@ref) used before for this constraint[^2], and as a fourth parameter the associated SOS matrix from
the previous optimization. For efficiency reasons, this should only be used for changes that preserve the structure of the
representation (i.e., whether it was PSD/DD/SDD and if its rotation was diagonal, triangular, or dense). If a
structure non-preserving change is made, the problem needs to be constructed from scratch. For non-diagonal rotations, consider
using [`RepresentationIAs`](@ref) in the first optimization.

!!! warning
    The internal state of the previous solver run will be re-used whenever possible. Therefore, no further data may be queried
    from the previous result afterwards, unless a re-optimization from scratch was necessary. While `result` will still be able
    to offer information about the relaxation, method, time, status, and objective value, moment matrices can only be accessed
    if they were already cached (i.e., accessed) before. Existing SOS certificates of the previous result will still be
    available, but new ones may not be constructed.

See also [`IterateRepresentation`](@ref).

[^2]: Roughly, as the exact type is not known. For sure, it will be possible to distinguish between
      [`RepresentationPSD`](@ref), [`RepresentationDD`](@ref), and [`RepresentationSDD`](@ref). The matrix type will not be
      concrete, but either `Union{<:UniformScaling,<:Diagonal}` if a diagonal representation was used before,
      `UpperOrUnitUpperTriangular`, `LowerOrUnitLowerTriangular`, or `Matrix` else. The complex identification will be `true`
      if a complex-valued cone was used and `false` else (where during specification, it could also have been `true` for
      real-valued data, which would simply be ignored). In any case, the third parameter can be used as a constructor accepting
      (unless it is for [`RepresentationPSD`](@ref)) the new rotation matrix as parameter. This is recommended, as in this way
      the `complex` value cannot change back to `true` for real-valued data, which would be interpreted as a change in
      structure, even if it is not.
"""
function poly_optimize(result::Result; representation=IterateRepresentation(), verbose::Bool=false, kwargs...)
    ismissing(result.state) && throw(ArgumentError("The given result does not contain any data."))
    @verbose_info("Beginning re-optimization...")
    oldstate = result.state
    if !(representation isa RepresentationMethod)
        representation = Solver.Rerepresent([[s for (s, _) in x] for x in Solver.extract_info(result.state)],
            SOSCertificate(result), representation, true)
    end
    relaxation = result.relaxation
    local otime
    try
        otime = @elapsed begin
            new_result = poly_optimize(Val(result.method), result.state, relaxation, Relaxation.groupings(relaxation);
                verbose, representation, kwargs...)
        end
    catch e
        if e isa Solver.RepresentationChangedError && representation isa Solver.Rerepresent
            @warn("The representation of at least one grouping changed, either in type or diagonally. The problem has to be \
                   set-up from the beginning. If this is due to a change in diagonality, consider using \
                   RepresentationIAs as initial value.")
            representation = Solver.Rerepresent(representation, false)
            otime = @elapsed begin
                new_result = Solver._Copied(poly_optimize(Val(result.method), relaxation, Relaxation.groupings(relaxation),
                    verbose, representation, kwargs...))
            end
        else
            rethrow(e)
        end
    end
    if new_result isa Solver._Copied
        new_result = new_result.data
    else
        result.state = missing
    end
    return Result(relaxation, result.method, otime, new_result...)
end