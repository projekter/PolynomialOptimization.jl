struct Custom{P<:Problem,MV<:IntMonomialVector,G<:RelaxationGroupings} <: AbstractRelaxationBasis{P}
    problem::P
    degree::Int
    basis::MV
    groupings::G

    @doc """
        Custom(problem::Problem, basis)

    Constructs a relaxation out of a polynomial optimization problem for the case in which a suitable basis is already known.
    """
    function Custom(problem::P,
        basis::Union{<:IntMonomialVector{Nr,Nc},<:AbstractVector{<:AbstractMonomialLike}}) where {Nr,Nc,I<:Integer,Poly<:IntPolynomial{<:Any,Nr,Nc,<:IntMonomialVector{Nr,Nc,I}},P<:Problem{Poly}}
        if !(basis isa IntMonomialVector)
            basis = IntMonomialVector{I}(basis; vars=problem.original_variables)
            if !iszero(Nc)
                for v in effective_variables(basis)
                    isconj(v) && throw(ArgumentError("Any custom basis must not contain explicit conjugates"))
                end
            end
        end
        degree = maxdegree(basis) # ≡ maxdegree_complex, but maxdegree cannot break branch prediction
        gr = groupings(problem, basis, degree)
        new{P,typeof(basis),typeof(gr)}(problem, Int(degree), basis, gr)
    end
end

default_solution_method(::Custom) = :heuristic