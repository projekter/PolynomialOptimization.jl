struct Dense{P<:Problem,MV<:SimpleMonomialVector,G<:RelaxationGroupings} <: AbstractRelaxationBasis{P}
    problem::P
    degree::Int
    basis::MV
    groupings::G

    @doc """
        Dense(problem::Problem[, degree])

    Constructs a full dense relaxation out of a polynomial optimization problem. This is the largest possible representation
    for a given degree bound, giving the best bounds. It is wasteful at the same time, as a Newton relaxation gives equally
    good bounds; but contrary to the Newton one, solution reconstruction works much better with a dense basis.
    `degree` is the degree of the Lasserre relaxation, which must be larger or equal to the halfdegree of all polynomials that
    are involved. If `degree` is omitted, the minimum required degree will be used.
    Specifying a degree larger than the minimal only makes sense if there are inequality or PSD constraints present, else it
    needlessly complicates calculations without any benefit.
    """
    function Dense(problem::P,
        degree::Integer=(@info("Automatically selecting minimal degree cutoff $(problem.mindegree)"); problem.mindegree)) where
        {Nr,Nc,I<:Integer,Poly<:SimplePolynomial{<:Any,Nr,Nc,<:SimpleMonomialVector{Nr,Nc,I}},P<:Problem{Poly}}
        degree < problem.mindegree && throw(ArgumentError("The minimally required degree is $(problem.mindegree)"))
        if iszero(Nc)
            basis = monomials(Val(Nr), Val(Nc), 0:degree; I)
        else
            maxmultideg = Vector{Int}(undef, Nr + 2Nc)
            @inbounds fill!(@view(maxmultideg[1:Nr]), degree)
            idx = Nr +1
            @inbounds for i in 1:Nc
                maxmultideg[idx] = degree
                maxmultideg[idx+1] = 0
                idx += 2
            end
            basis = monomials(Val(Nr), Val(Nc), 0:degree; maxmultideg, I)
        end
        gr = groupings(problem, basis, degree, nothing)
        new{P,typeof(basis),typeof(gr)}(problem, Int(degree), basis, gr)
    end
end

default_solution_method(::Dense) = :mvhankel