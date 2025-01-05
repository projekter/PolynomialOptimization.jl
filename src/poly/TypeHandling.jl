# These are simplified versions, relying on default fill-ins. A much more intricate way of dealing with incomplete type
# specifications was present previously (see commit 04ab29d91c42045d3d10bf54344ec34eeb1c0c68).
MultivariatePolynomials.variable_union_type(::XorTX{<:Union{<:SimpleVariable{Nr,Nc},
                                                            <:SimpleMonomialOrConj{Nr,Nc},
                                                            <:Term{<:Any,<:SimpleMonomial{Nr,Nc}},
                                                            <:SimplePolynomial{<:Any,Nr,Nc}}}) where {Nr,Nc} =
    SimpleVariable{Nr,Nc,smallest_unsigned(Nr + 2Nc)}

MultivariatePolynomials.monomial_type(::XorTX{<:SimpleVariable{Nr,Nc}}) where {Nr,Nc} =
    SimpleMonomial{Nr,Nc,UInt,ExponentsAll{Nr+2Nc,UInt}}

MultivariatePolynomials.term_type(::XorTX{<:Union{SimpleMonomial{Nr,Nc,I,E},<:SimpleMonomial{Nr,Nc,I},<:SimpleMonomial{Nr,Nc,<:Integer,E}}},
    ::Type{C}) where {Nr,Nc,I<:Integer,E<:AbstractExponents{<:Any,I},C} =
    Term{C,SimpleMonomial{Nr,Nc,I,@isdefined(E) ? E : ExponentsAll{Nr+2Nc,I}}}
MultivariatePolynomials.term_type(::XorTX{<:SimplePolynomial{Cold,Nr,Nc,<:SimpleMonomialVector{Nr,Nc,I,<:Any,M}}},
    ::Type{C}=Cold) where {Cold,C,Nr,Nc,I<:Integer,M<:SimpleMonomial{Nr,Nc,I}} =
    Term{C,monomial_type(M)}

MultivariatePolynomials.polynomial_type(::XorTX{<:Term{C,<:SimpleMonomial{Nr,Nc,I}}}) where {C,Nr,Nc,I<:Integer} =
    SimplePolynomial{C,Nr,Nc,SimpleMonomialVector{Nr,Nc,I,Tuple{ExponentsAll{Nr+2Nc,I},Vector{I}},SimpleMonomial{Nr,Nc,I,ExponentsAll{Nr+2Nc,I}}}}
MultivariatePolynomials.polynomial_type(::XorTX{<:SimplePolynomial{<:Any,Nr,Nc,M}}, ::Type{C}) where {C,Nr,Nc,M<:SimpleMonomialVector{Nr,Nc}} =
    SimplePolynomial{C,Nr,Nc,M}
MultivariatePolynomials.polynomial_type(::XorTX{M}, ::Type{C}) where {C,Nr,Nc,M<:SimpleMonomialVector{Nr,Nc}} =
    SimplePolynomial{C,Nr,Nc,M}

MultivariatePolynomials.coefficient_type(::XorTX{SimplePolynomial{C}}) where {C} = C

Base.convert(::Type{T}, x::Q) where {T<:Union{<:SimpleVariable,<:SimpleMonomial,<:SimpleMonomialVector,<:SimplePolynomial},Q<:T} = x

MultivariatePolynomials.promote_rule_constant(::Type{C}, M::Type{<:SimpleMonomial}) where {C} =
    term_type(M, promote_type(C, Int))
MultivariatePolynomials.promote_rule_constant(::Type{C}, T::Type{<:Term{Cold,<:SimpleMonomial}}) where {C,Cold} =
    term_type(T, promote_type(C, Cold))
MultivariatePolynomials.promote_rule_constant(::Type{C}, P::Type{<:SimplePolynomial{Cold}}) where {C,Cold} =
    polynomial_type(P, promote_type(T, Cold))
MultivariatePolynomials.promote_rule_constant(::Type, ::Type{<:Union{<:Term{<:Any,<:SimpleMonomial},<:SimplePolynomial}}) = Any

function Base.promote_rule(::Type{<:Term{C1,<:SimpleMonomial{Nr1,Nc1,I1,E1}}}, ::Type{<:Term{C2,<:SimpleMonomial{Nr2,Nc2,I2,E2}}}) where
    {C1,Nr1,Nc1,I1<:Integer,E1<:AbstractExponents,C2,Nr2,Nc2,I2<:Integer,E2<:AbstractExponents}
    # make sure this will always be an error and not recurse infinitely
    (Nr1 === Nr2 && Nc1 === Nc2) || throw(MethodError(promote, (SimpleMonomial{Nr1,Nc1}, SimpleMonomial{Nr2,Nc2})))
    I = promote_type(I1, I2)
    return Term{promote_type(C1, C2),
                <:(E1 isa ExponentsAll && E2 isa ExponentsAll ? SimpleMonomial{Nr1,Nc1,I,ExponentsAll{Nr1+2Nc1,I}} :
                                                                SimpleMonomial{Nr1,Nc1,I,<:AbstractExponents{Nr1+2Nc1,I}})}
end
function Base.promote_rule(::Type{<:Term{C,<:SimpleMonomial{Nr1,Nc1,I1,E1}}}, ::Type{SimpleMonomial{Nr2,Nc2,I2,E2}}) where
    {C,Nr1,Nc1,I1<:Integer,E1<:AbstractExponents,Nr2,Nc2,I2<:Integer,E2<:AbstractExponents}
    # make sure this will always be an error and not recurse infinitely
    (Nr1 === Nr2 && Nc1 === Nc2) || throw(MethodError(promote, (SimpleMonomial{Nr1,Nc1}, SimpleMonomial{Nr2,Nc2})))
    I = promote_type(I1, I2)
    return Term{C,<:(E1 isa ExponentsAll && E2 isa ExponentsAll ? SimpleMonomial{Nr1,Nc1,I,ExponentsAll{Nr1+2Nc1,I}} :
                                                                  SimpleMonomial{Nr1,Nc1,I,<:AbstractExponents{Nr1+2Nc1,I}})}
end
function Base.promote_rule(::Type{<:Term{C,<:SimpleMonomial{Nr1,Nc1,I,E}}}, ::Type{<:SimpleVariable{Nr2,Nc2,<:Unsigned}}) where
    {C,Nr1,Nc1,I<:Integer,E<:AbstractExponents,Nr2,Nc2}
    # make sure this will always be an error and not recurse infinitely
    (Nr1 === Nr2 && Nc1 === Nc2) || throw(MethodError(promote, (SimpleMonomial{Nr1,Nc1}, SimpleVariable{Nr2,Nc2})))
    return Term{C,SimpleMonomial{Nr1,Nc1,I,E}}
end
function Base.promote_rule(::Type{<:SimpleMonomial{Nr1,Nc1,I,E}}, ::Type{<:SimpleVariable{Nr2,Nc2,<:Unsigned}}) where
    {Nr1,Nc1,I<:Integer,E<:AbstractExponents,Nr2,Nc2}
    # make sure this will always be an error and not recurse infinitely
    (Nr1 === Nr2 && Nc1 === Nc2) || throw(MethodError(promote, (SimpleMonomial{Nr1,Nc1}, SimpleVariable{Nr2,Nc2})))
    return SimpleMonomial{Nr1,Nc1,I,E}
end
function Base.promote_rule(::Type{SimpleVariable{Nr1,Nc1,I1}}, ::Type{SimpleVariable{Nr2,Nc2,I2}}) where
    {Nr1,Nc1,I1<:Unsigned,Nr2,Nc2,I2<:Unsigned}
    # make sure this will always be an error and not recurse infinitely
    (Nr1 === Nr2 && Nc1 === Nc2) || throw(MethodError(promote, (SimpleVariable{Nr1,Nc1}, SimpleVariable{Nr2,Nc2})))
    @assert I1 === I2
    return SimpleVariable{Nr1,Nc1,I1}
end