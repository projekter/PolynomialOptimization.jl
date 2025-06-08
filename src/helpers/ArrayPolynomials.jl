for f in (:variables, :effective_variables)
    @eval function MultivariatePolynomials.$f(a::AbstractArray{<:AbstractPolynomialLike})
        isempty(a) && return variable_union_type(eltype(a))[]
        result = Set{variable_union_type(eltype(a))}()
        for x in a
            union_!(result, $f(x))
        end
        return sort!(collect(result), rev=true)
    end
end

function MultivariatePolynomials.effective_variables(a::AbstractArray{<:IntPolynomial{<:Any,Nr,Nc}};
    rettype::Type{V}=Vector, by::F=identity) where {Nr,Nc,V<:Union{Vector,Set},F}
    ET = Base.promote_op(by, IntVariable{Nr,Nc,IntPolynomials.smallest_unsigned(Nr+2Nc)})
    isempty(a) && return V <: Vector ? ET[] : Set{ET}()
    a, rest = Iterators.peel(a)
    @inbounds result = effective_variables(a; rettype=Set, by)
    for x in rest
        union_!(result, effective_variables(x; rettype=Set, by))
    end
    V <: Vector && return sort!(collect(result), rev=true)
    return result
end

MultivariatePolynomials.monomials(a::AbstractArray{<:AbstractPolynomialLike}) = merge_monomial_vectors(monomials.(a))
# For the type-stable implementation let's assume that if all of the polynomials share a common exponent type, then they'll
# actually share the same exponents object.
MultivariatePolynomials.monomials(a::AbstractArray{<:IntPolynomial{<:Any,Nr,Nc,<:IntMonomialVector{Nr,Nc,E}}}) where {Nr,Nc,E} =
    merge_monomial_vectors(Val(Nr), Val(Nc), monomials(first(a)).e, monomials.(Iterators.flatten(a)))

MultivariatePolynomials.coefficients(a::AbstractArray{<:AbstractPolynomialLike}) = coefficients.(a, (monomials(a),))
MultivariatePolynomials.coefficients(a::AbstractArray{<:AbstractPolynomialLike}, X::AbstractVector) = coefficients.(a, (X,))

for f in (:mindegree, :mindegree_complex, :minhalfdegree, :maxdegree, :maxdegree_complex, :maxhalfdegree)
    mini = startswith(string(f), "min")
    @eval MultivariatePolynomials.$f(a::AbstractArray{<:AbstractPolynomialLike}, args...) =
        $(mini ? :minimum : :maximum)(($f(p, args...) for p in a), init=$(mini ? typemax(Int) : 0))::Int
end
for f in (:extdegree, :extdegree_complex, :exthalfdegree)
    @eval function MultivariatePolynomials.$f(a::AbstractArray{<:AbstractPolynomialLike}, args...)
        l = typemax(Int)
        u = 0
        for p in a
            (newl, newu) = $f(p, args...)
            newl < l && (l = newl)
            newu > u && (u = newu)
        end
        return l, u
    end
end

IntPolynomials.effective_variables_in(a::AbstractArray{<:AbstractPolynomialLike}, in) =
    all(Base.Fix2(effective_variables_in, in), a)

Base.isreal(a::AbstractArray{<:AbstractPolynomialLike}) = all(isreal, a)