"""
    RepresentationPSD <: RepresentationMethod

Model the constraint "σ ⪰ 0" as a positive semidefinite cone membership, `σ ∈ PSD`. This is the strongest possible model, but
the most resource-intensive.
"""
struct RepresentationPSD end

"""
    RepresentationSDD([u]; complex=true) <: RepresentationMethod

Model the constraint "σ ⪰ 0" as a membership in the scaled diagonally dominant cone, ``\\sigma = u^\\dagger Q u`` for some
`Q ∈ SDD`. The matrix `u` is by default an identity of any dimension; however, usually, care must be taken to have a matrix
of suitable dimension.
The membership `Q ∈ SDD` is achieved using the scaled diagonally dominant (dual) cone directly or (rotated) quadratic cones.

If σ is a Hermitian matrix, a complex-valued scaled diagonally dominant (dual) cone will be used, if supported. If not,
fallbacks to the (rotated) quadratic cones are used; however, if `complex=false` and the ordinary scaled diagonally dominant
(dual) cone is supported, rewrite the matrix as a real one and then use the real-valued cone. This is usually never advisable,
as the rotated quadratic cone always works on the complex data.
Note that if rewritten, `u` must be real-valued and have twice the side dimension of the complex-valued matrix.
"""
struct RepresentationSDD{M,Complex}
    u::M

    RepresentationSDD(u::M=I; complex::Bool=true) where {M} = new{M,complex}(u)
end

(::Type{<:RepresentationSDD{<:Any,complex}})(u) where {complex} = RepresentationSDD(u; complex)

@doc raw"""
    RepresentationDD([u]; complex=true) <: RepresentationMethod

Model the constraint "σ ⪰ 0" as a membership in the diagonally dominant cone, ``\sigma = u^\dagger Q u`` for some `Q ∈ DD`.
The matrix `u` is by default an identity of any dimension; however, usually, care must be taken to have a matrix of suitable
dimension.
The membership `Q ∈ DD` is achieved using the diagonally dominant (dual) cone directly, ``\ell_1``- or ``\ell_\infty``-norm
cones or linear inequalities; slack variables will be added as necessary.

If σ is a Hermitian matrix, a complex-valued diagonally dominant (dual) cone will be used, if supported. If not, fallbacks will
first try quadratic cones on the complex-valued data, and if this is also not supported, rewrite the matrix as a real one and
then apply the real-valued DD constraint.
By setting the `complex` parameter to `false`, the rewriting to a real matrix will always be used, regardless of complex-valued
solver support.
Note that if rewritten, `u` must be real-valued and have twice the side dimension of the complex-valued matrix.
"""
struct RepresentationDD{M,Complex}
    u::M

    RepresentationDD(u::M=I; complex::Bool=true) where {M} = new{M,complex}(u)
end

(::Type{<:RepresentationDD{<:Any,complex}})(u) where {complex} = RepresentationDD(u; complex)

"""
    RepresentationMethod{M,Complex}

Union type that defines how the optimizer constraint σ ⪰ 0 is interpreted. Usually, "⪰" means positive semidefinite;
however, there are various other possibilities giving rise to weaker results, but scale more favorably. The following methods
are supported:
- [`RepresentationPSD`](@ref)
- [`RepresentationSDD`](@ref)
- [`RepresentationDD`](@ref)

See also [`RepresentationIAs`](@ref).
"""
const RepresentationMethod{M,Complex} = Union{RepresentationPSD,RepresentationSDD{M,Complex},RepresentationDD{M,Complex}}

"""
    RepresentationIAs(r::Type{RepresentationDD,RepresentationSDD},
        m::Type{<:AbstractMatrix}=UpperTriangular; complex=true)

Default callable that will instantiate a correctly-sized representation of type `r` with an identity rotation that is, however,
not recognized as a diagonal rotation but as type `m` instead. Use this type in the first call to [`poly_optimize`](@ref) if
you want to re-optimize the problem afterwards with rotations of type `m`. The default for `m`, `UpperTriangular`, is suitable
for the automatic Cholesky-based reoptimization.
"""
struct RepresentationIAs{R<:Union{RepresentationDD,RepresentationSDD},M,Complex}
    RepresentationIAs(r::Type{<:Union{RepresentationDD,RepresentationSDD}}, m::Type{<:AbstractMatrix}; complex::Bool=true) =
        new{r,m,complex}()
end

(::RepresentationIAs{R,M,complex})(_, dim) where {R,M<:Matrix,complex} = R(Matrix(I, dim, dim); complex)
(::RepresentationIAs{R,M,complex})(_, dim) where {R,M,complex} = R(M(Diagonal(I, dim)); complex)

struct RepresentationChangedError <: Exception end

struct Rerepresent{C,F}
    info::Vector{Vector{Symbol}}
    cert::C
    fn::F
    fix_structure::Bool
end

Rerepresent(r::Rerepresent, fix_structure::Bool) = Rerepresent(r.info, r.cert, r.fn, fix_structure)

function (r::Rerepresent)((type, index, grouping), dim)
    idx = PolynomialOptimization.id_to_index(r.cert.relaxation, (type, index, grouping))
    oldtype = r.info[idx][grouping]
    if oldtype in INFO_PSD
        oldrep = RepresentationPSD
    else
        complex = oldtype in INFO_COMPLEX
        if oldtype in INFO_DIAG
            oldmat = Union{<:UniformScaling,<:Diagonal}
        elseif oldtype in INFO_TRIU
            oldmat = UpperOrUnitUpperTriangular
        elseif oldtype in INFO_TRIL
            oldmat = LowerOrUnitLowerTriangular
        else
            oldmat = Matrix
        end
        if oldtype in INFO_SDD
            oldrep = RepresentationSDD{<:oldmat,complex}
        else
            @assert(oldtype in INFO_DD)
            oldrep = RepresentationDD{<:oldmat,complex}
        end
    end
    newrep = r.fn((type, index, grouping), dim, oldrep, r.cert.data[idx][grouping])::RepresentationMethod
    r.fix_structure && !(newrep isa oldrep) && throw(RepresentationChangedError())
    return newrep
end