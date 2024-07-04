for f in (:variables, :effective_variables)
    @eval begin
        function MultivariatePolynomials.$f(m::AbstractMatrix{<:AbstractPolynomialLike})
            isempty(m) && return variable_union_type(eltype(m))[]
            result = Set{variable_union_type(eltype(m))}()
            for x in m
                union!(result, $f(x))
            end
            return sort!(collect(result))
        end
    end
end

function MultivariatePolynomials.effective_variables(m::AbstractMatrix{<:SimplePolynomial{<:Any,Nr,Nc}};
    rettype::Type{V}=Vector, by::F=identity) where {Nr,Nc,V<:Union{Vector,Set},F}
    ET = Base.promote_op(by, SimpleVariable{Nr,Nc,SimplePolynomials.smallest_unsigned(Nr+2Nc)})
    isempty(m) && return V <: Vector ? ET[] : Set{ET}()
    a, rest = Iterators.peel(m)
    @inbounds result = effective_variables(a; rettype=Set, by)
    for x in rest
        union!(result, effective_variables(x; rettype=Set, by))
    end
    V <: Vector && return sort!(collect(result))
    return result
end

MultivariatePolynomials.monomials(m::AbstractMatrix{<:AbstractPolynomialLike}) = merge_monomial_vectors(monomials.(m))
# For the type-stable implementation let's assume that if all of the polynomials share a common exponent type, then they'll
# actually share the same exponents object.
MultivariatePolynomials.monomials(m::AbstractMatrix{<:SimplePolynomial{<:Any,Nr,Nc,<:SimpleMonomialVector{Nr,Nc,E}}}) where {Nr,Nc,E} =
    merge_monomial_vectors(Val(Nr), Val(Nc), monomials(first(m)).e, monomials.(Iterators.flatten(m)))

MultivariatePolynomials.coefficients(m::AbstractMatrix{<:AbstractPolynomialLike}) = coefficients.(m, (monomials(m),))
MultivariatePolynomials.coefficients(m::AbstractMatrix{<:AbstractPolynomialLike}, X::AbstractVector) = coefficients.(m, (X,))

MultivariatePolynomials.mindegree(m::AbstractMatrix{<:AbstractPolynomialLike}, args...) =
    minimum((mindegree(p, args...) for p in m), init=0)::Int
MultivariatePolynomials.mindegree_complex(m::AbstractMatrix{<:AbstractPolynomialLike}, args...) =
    minimum((mindegree_complex(p, args...) for p in m), init=0)::Int
MultivariatePolynomials.minhalfdegree(m::AbstractMatrix{<:AbstractPolynomialLike}, args...) =
    minimum((minhalfdegree(p, args...) for p in m), init=0)::Int
MultivariatePolynomials.maxdegree(m::AbstractMatrix{<:AbstractPolynomialLike}, args...) =
    maximum((maxdegree(p, args...) for p in m), init=0)::Int
MultivariatePolynomials.maxdegree_complex(m::AbstractMatrix{<:AbstractPolynomialLike}, args...) =
    maximum((maxdegree_complex(p, args...) for p in m), init=0)::Int
MultivariatePolynomials.maxhalfdegree(m::AbstractMatrix{<:AbstractPolynomialLike}, args...) =
    maximum((maxhalfdegree(p, args...) for p in m), init=0)::Int
function MultivariatePolynomials.extdegree(m::AbstractMatrix{<:AbstractPolynomialLike}, args...)
    l = typemax(Int)
    u = 0
    for p in m
        (newl, newu) = extdegree(p, args...)
        newl < l && (l = newl)
        newu > u && (u = newu)
    end
    return l, u
end
function MultivariatePolynomials.extdegree_complex(m::AbstractMatrix{<:AbstractPolynomialLike}, args...)
    l = typemax(Int)
    u = 0
    for p in m
        (newl, newu) = extdegree_complex(p, args...)
        newl < l && (l = newl)
        newu > u && (u = newu)
    end
    return l, u
end
function MultivariatePolynomials.exthalfdegree(m::AbstractMatrix{<:AbstractPolynomialLike}, args...)
    l = typemax(Int)
    u = 0
    for p in m
        (newl, newu) = exthalfdegree(p, args...)
        newl < l && (l = newl)
        newu > u && (u = newu)
    end
    return l, u
end

SimplePolynomials.effective_variables_in(m::AbstractMatrix{<:AbstractPolynomialLike}, in) =
    all(Base.Fix2(effective_variables_in, in), m)

Base.isreal(m::AbstractMatrix{<:AbstractPolynomialLike}) = all(isreal, m)