struct FacialReduction{P<:Problem,G<:RelaxationGroupings} <: AbstractRelaxationSparse{P}
    parent::AbstractRelaxation{P} # no specialization
    groupings::G

    @doc """
        FacialReduction(relaxation::AbstractRelaxation; [method,] parameters...)

    Constructs a relaxation based on facial reduction applied to another relaxation of a polynomial optimization problem. This
    will be a superset of the largest possible representation for a given degree bound, with no negative consequences for
    finding the optimum. It can be much smaller than a dense basis, but solution reconstruction may be harder. It is the most
    effective solution-neutral method, surpassing the Newton polytope; however, its construction can take significantly longer.

    When constraints are present, the polytope is determined based on the Putinar reformulation, where all constraints and the
    objective are moved to one side (comprising a new virtual objective). The prefactors for the constraints are determined by
    the previous relaxation method. The fact that nonnegative constraints must yield SOS prefactors is taken into account by
    this method. PSD constraints are not supported.

    The `method` determines which solver to use for performing the facial reduction. If omitted, this will be the default
    solver.
    The `parameters` are passed on to
    [`FacialReduction.facial_reduction!!`](@ref PolynomialOptimization.FacialReduction.facial_reduction!!).
    """
    function FacialReduction(relaxation::AbstractRelaxation{P};
        method::Symbol=PolynomialOptimization.FacialReduction.default_reduction_method(), verbose::Bool=false, parameters...) where
        {P<:Problem}
        @verbose_info("Performing facial reduction")
        gentime = @elapsed(gr = PolynomialOptimization.FacialReduction.facial_reduction(method, relaxation; verbose,
            parameters...))
        @verbose_info("Facial reduction done in ", gentime, " seconds")

        new{P,typeof(gr)}(relaxation, gr)
    end
end

# All the parameters are only for choosing a different method, which might impact runtime behavior; but all methods valid for
# a certain problem must always give the same result. So there's no use in doing the work all over, even if the parameters
# have changed.
FacialReduction(relaxation::FacialReduction; kwargs...) = relaxation

default_solution_method(::FacialReduction) = :heuristic