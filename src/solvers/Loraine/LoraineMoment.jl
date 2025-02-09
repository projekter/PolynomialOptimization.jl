mutable struct StateMoment{I<:Integer,K<:Integer,V<:Real} <: AbstractSparseMatrixSolver{I,K,V}
    const A::Vector{SparseMatrixCSC{V,I}}
    const C::Vector{SparseMatrixCSC{V,I}}
    A_lin::SparseMatrixCSC{V,I}
    c_lin::SparseVector{V,I}
    b::Vector{V}
    info::Vector{<:Vector{<:Tuple{Symbol,Any}}}
    data

    StateMoment{I,K,V}() where {I<:Integer,K<:Integer,V<:Real} = new{I,K,V}(
        SparseMatrixCSC{V,I}[],
        SparseMatrixCSC{V,I}[]
    )
end

Solver.issuccess(::Val{:LoraineMoment}, status::Int) = status == 1

Solver.psd_indextype(::StateMoment) = PSDIndextypeCOOVectorized(:F, 1)

Solver.objective_indextype(::StateMoment) = PSDIndextypeMatrixCartesian(:F, 1)

@counter_atomic(StateMoment, :psd)

function Solver.add_var_nonnegative!(state::StateMoment{I,<:Integer,V}, m::Int, n::Int, data::SparseMatrixCOO{I,I,V,Offset},
    obj::Tuple{FastVec{I},FastVec{V}}) where {I<:Integer,V<:Real,Offset}
    isone(Offset) || error("Internal error")
    @assert(!isdefined(state, :A_lin) && !isdefined(state, :c_lin))
    state.A_lin = SparseMatrixCSC(m, n, coo_to_csc!(n, data)...)
    state.c_lin = SparseVector(n, finish!.(obj)...)
    return
end

function Solver.add_var_psd!(state::StateMoment{I,<:Integer,V}, m::Int, dim::Int, data::SparseMatrixCOO{I,I,V,Offset},
    obj::Union{Nothing,Tuple{Tuple{FastVec{I},FastVec{I}},FastVec{V}}}) where {I<:Integer,V<:Real,Offset}
    isone(Offset) || error("Internal error")
    dsq = dim^2
    push!(state.A, SparseMatrixCSC(m, dsq, coo_to_csc!(dsq, data)...))
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

function Solver.poly_optimize(::Val{:LoraineMoment}, relaxation::AbstractRelaxation, groupings::RelaxationGroupings;
    representation, verbose::Bool=false, customize=(state) -> nothing, parameters...)
    setup_time = @elapsed @inbounds begin
        K = _get_I(eltype(monomials(poly_problem(relaxation).objective)))
        V = real(coefficient_type(poly_problem(relaxation).objective))

        state = StateMoment{Int,K,V}()
        primal_data = primal_moment_setup!(state, relaxation, groupings; verbose)
        ismissing(primal_data) && return missing, Loraine.STATUS_INFEASIBLE_OR_UNBOUNDED, typemin(V)
        state.info, state.data = primal_data
        customize(state)

        @verbose_info("Initializing solver")
        solver = Loraine.Solver(Loraine.Model(; state.A, state.C, state.b,
            A_lin=isdefined(state, :A_lin) ? state.A_lin : nothing,
            c_lin=isdefined(state, :c_lin) ? state.c_lin : nothing,
            coneDims=finish!(primal_data[2].psd_dim), check=false);
            verb=verbose ? Loraine.VERBOSITY_FULL : Loraine.VERBOSITY_NONE, parameters...)
    end
    @verbose_info("Setup complete in ", setup_time, " seconds")
    solver_time = @elapsed Loraine.solve!(solver)

    @verbose_info("Optimization complete in ", solver_time, " seconds")

    return (state, solver), solver.status, dot(state.b, solver.y)
end

Solver.extract_moments(relaxation::AbstractRelaxation, (state, solver)::Tuple{StateMoment,Loraine.Solver}) =
    MomentVector(relaxation, state.data, solver.X, solver.X_lin)

Solver.extract_sos(::AbstractRelaxation, (_, solver)::Tuple{StateMoment,Loraine.Solver}, ::Val{:nonnegative},
    index::AbstractUnitRange, ::Nothing) = @view(solver.S_lin[index])

Solver.extract_sos(::AbstractRelaxation, (_, solver)::Tuple{StateMoment,Loraine.Solver}, ::Val{:psd},
    index::Integer, ::Nothing) = solver.S[index]