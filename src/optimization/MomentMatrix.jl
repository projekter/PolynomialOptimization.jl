export moment_matrix

moment_matrix(moments::MomentVector{<:Any,V,Nr,Nc}, ::Val{:symmetric}, rows, cols,
    prefix::Union{<:SimpleMonomialOrConj{Nr,Nc},<:SimpleVariable{Nr,Nc}}...) where {Nr,Nc,V} =
    V[i ≤ j ? get(moments, let n=V(NaN); () -> n end, rows[i], cols[j], prefix...) : V(NaN)
      for i in 1:length(rows), j in 1:length(cols)]

moment_matrix(moments::MomentVector{<:Any,V,Nr,Nc}, ::Val{:throw}, rows, cols,
    prefix::Union{<:SimpleMonomialOrConj{Nr,Nc},<:SimpleVariable{Nr,Nc}}...) where {Nr,Nc,V} =
    V[get(moments, () -> throw(MonomialMissing()), rows[i], cols[j], prefix...)
      for i in 1:length(rows), j in 1:length(cols)]

"""
    moment_matrix(problem::Result; max_deg=Inf, prefix=1)

After a problem has been optimized, this function assembles the associated moment matrix (possibly by imposing a degree bound
`max_deg`, and possibly multiplying each monomial by the monomial or variable `prefix`, which does not add to `max_deg`).
Note that `prefix` has to be a valid `SimpleMonomial` or `SimpleVariable` of appropriate type.

See also [`poly_optimize`](@ref), [`poly_optimize`](@ref).
"""
function moment_matrix(result::Result{R}; max_deg=Inf,
    prefix::Union{<:SimpleMonomialOrConj{Nr,Nc},<:SimpleVariable{Nr,Nc},Nothing}=nothing) where {Nr,Nc,R<:AbstractRelaxation{<:Problem{<:SimplePolynomial{<:Any,Nr,Nc}}}}
    relaxation = result.relaxation
    b = basis(relaxation)
    if max_deg < Inf
        b = Relaxation.truncate_basis(b, max_deg)
    end
    return (isreal(relaxation) ? Symmetric : Hermitian)(
        isnothing(prefix) ? moment_matrix(result.moments, Val(:symmetric), b, conj(b)) :
                            moment_matrix(result.moments, Val(:symmetric), b, conj(b), prefix))
end