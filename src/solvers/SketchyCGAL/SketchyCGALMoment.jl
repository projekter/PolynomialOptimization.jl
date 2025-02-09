mutable struct StateMoment{I<:Integer,K<:Integer,V<:Real} <: AbstractSparseMatrixSolver{I,K,V}
    const n::FastVec{Int}
    const A::FastVec{Transpose{V,SparseMatrixCSC{V,I}}}
    const C::FastVec{SparseMatrixCSC{V,I}}
    b::Vector{V}
    info::Vector{<:Vector{<:Tuple{Symbol,Any}}}
    data

    StateMoment{I,K,V}() where {I<:Integer,K<:Integer,V<:Real} = new{I,K,V}(
        FastVec{Int}(),
        FastVec{Transpose{V,SparseMatrixCSC{V,I}}}(),
        FastVec{SparseMatrixCSC{V,I}}()
    )
end

Solver.issuccess(::Val{:SketchyCGALMoment}, status::Symbol) = status === :optimal

Solver.psd_indextype(::StateMoment) = PSDIndextypeCOOVectorized(:F, 1)

Solver.objective_indextype(::StateMoment) = PSDIndextypeMatrixCartesian(:F, 1)

@counter_atomic(StateMoment, :psd)

Solver.add_var_nonnegative!(state::StateMoment{I,<:Integer,V}, m::Int, n::Int, data::SparseMatrixCOO{I,I,V,Offset},
    obj::Tuple{FastVec{I},FastVec{V}}) where {I<:Integer,V<:Real,Offset} =
    error("Size-1 constraints are not supported by SketchyCGAL")

function Solver.add_var_psd!(state::StateMoment{I,<:Integer,V}, m::Int, dim::Int, data::SparseMatrixCOO{I,I,V,Offset},
    obj::Union{Nothing,Tuple{Tuple{FastVec{I},FastVec{I}},FastVec{V}}}) where {I<:Integer,V<:Real,Offset}
    isone(Offset) || error("Internal error")
    dsq = dim^2
    # we often need row slices, so A is constructed as a CSR matrix
    push!(state.A, transpose(SparseMatrixCSC(dsq, m, coo_to_csr!(m, data)...)))
    push!(state.C, isnothing(obj) ? spzeros(V, I, dim, dim) :
                                    SparseMatrixCSC(dim, dim, coo_to_csc!(dim, SparseMatrixCOO(obj[1]..., obj[2], Offset))...))
    return
end

function Solver.fix_constraints!(state::StateMoment{I,<:Integer,V}, m::Int, indvals::Indvals{I,V}) where {I<:Integer,V<:Real}
    state.b = b = zeros(V, m)
    for (i, v) in indvals
        b[i] = v
    end
    return
end

function Solver.poly_optimize(::Val{:SketchyCGALMoment}, relaxation::AbstractRelaxation, groupings::RelaxationGroupings;
    representation, verbose::Bool=false, customize=(state) -> nothing, α::Union{Tuple{Real,Real},Real}, parameters...)
    if α isa Real
        α = (zero(α), α)
    end
    setup_time = @elapsed @inbounds begin
        K = _get_I(eltype(monomials(poly_problem(relaxation).objective)))
        V = real(coefficient_type(poly_problem(relaxation).objective))

        state = StateMoment{Int,K,V}()
        primal_data = primal_moment_setup!(state, relaxation, groupings; verbose)
        ismissing(primal_data) && return missing, :unknown, typemin(V)
        state.info, state.data = primal_data
        customize(state)
    end
    @verbose_info("Setup complete in ", setup_time, " seconds")
    solver_time = @elapsed begin
        status, value, Xs = SketchyCGAL.sketchy_cgal(finish!(state.A), state.b, finish!(state.C); verbose, α, parameters...)
    end

    @verbose_info("Optimization complete in ", solver_time, " seconds")

    return (state, Xs), status, value
end

Solver.extract_moments(relaxation::AbstractRelaxation, (state, solution)::Tuple{StateMoment,Any}) =
    MomentVector(relaxation, state.data, collect(map(Matrix, solution)), nothing)