mutable struct StateMoment{K<:Integer} <: AbstractSparseMatrixSolver{Cint,K,Float64}
    const solver::ConicBundle.CBMatrixCBSolver
    const opAt::ConicBundle.CBSparseCoeffmatMatrix
    const C::ConicBundle.CBSparseCoeffmatMatrix
    const cmats::FastVec{ConicBundle.CBCMsymsparse}
    const Cmats::FastVec{ConicBundle.CBSparsesym}
    psds::ConicBundle.CBPSCAffineFunction
    info::Vector{<:Vector{<:Tuple{Symbol,Any}}}
    data

    StateMoment{K}(solver::ConicBundle.CBMatrixCBSolver) where {K<:Integer} = new{K}(
        solver,
        ConicBundle.CBSparseCoeffmatMatrix(),
        ConicBundle.CBSparseCoeffmatMatrix(),
        FastVec{ConicBundle.CBCMsymsparse}(),
        FastVec{ConicBundle.CBSparsesym}(),
    )
end

Solver.issuccess(::Val{:ConicBundleMoment}, status::Int) = status == 1

Solver.psd_indextype(::StateMoment) = PSDIndextypeMatrixCartesian(:U, zero(Cint))

@counter_atomic(StateMoment, :psd)

Solver.add_var_nonnegative!(state::StateMoment, m::Int, n::Int, data::SparseMatrixCOO{Cint,Cint,Float64,zero(Cint)},
    obj::Tuple{FastVec{Cint},FastVec{Float64}}) = error("Size-1 constraints are not supported by the implementation")
    # perhaps use BoxOracle?

function Solver.add_var_psd!(state::StateMoment, m::Int, dim::Cint,
    data::Tuple{FastVec{Cint},Tuple{FastVec{Cint},FastVec{Cint}},FastVec{Float64}},
    obj::Union{Nothing,Tuple{Tuple{FastVec{Cint},FastVec{Cint}},FastVec{Float64}}})
    i = 1
    indm = ConicBundle.CBIndexmatrix(1, 1, dim)
    tmp = ConicBundle.CBSparseCoeffmatMatrix(indm, m)
    prepare_push!(state.cmats, m + !isnothing(obj))
    len = length(data[1])
    rmul!(data[3], -1.)
    @inbounds while i ≤ len
        j = i
        while j < len && data[1][j+1] == data[1][i]
            j += 1
        end
        mat = @views ConicBundle.CBSparsesym(dim, j - i +1, data[2][1][i:j], data[2][2][i:j], data[3][i:j])
        cmat = ConicBundle.CBCMsymsparse(mat)
        ConicBundle.cb_destroy!(mat)
        unsafe_push!(state.cmats, cmat)
        ConicBundle.cb_set!(tmp, 0, data[1][i], cmat)
        i = j +1
    end
    ConicBundle.cb_append_blocks!(state.opAt, tmp)
    tmp = ConicBundle.CBSparseCoeffmatMatrix(indm, 1)
    if !isnothing(obj)
        rmul!(obj[2], -1.)
        mat = ConicBundle.CBSparsesym(dim, length(obj[1][1]), obj[1]..., obj[2])
        push!(state.Cmats, mat)
        cmat = ConicBundle.CBCMsymsparse(mat)
        unsafe_push!(state.cmats, cmat)
        ConicBundle.cb_set!(tmp, 0, 0, cmat)
    end
    ConicBundle.cb_append_blocks!(state.C, tmp)
    ConicBundle.cb_destroy!(tmp)
    ConicBundle.cb_destroy!(indm)
    return
end

function Solver.fix_constraints!(state::StateMoment, m::Int, indvals::Indvals{Cint,Float64})
    mat = ConicBundle.CBMatrix(m, 1, 0.)
    for (i, v) in indvals
        mat[i, 0] = v
    end
    ConicBundle.cb_init_problem!(state.solver, m, nothing, nothing, nothing, mat)
    ConicBundle.cb_destroy!(mat)
    return
end

function Solver.poly_optimize(::Val{:ConicBundleMoment}, relaxation::AbstractRelaxation, groupings::RelaxationGroupings;
    representation, verbose::Bool=false, customize=(state) -> nothing, bound::Real=-1.)
    setup_time = @elapsed @inbounds begin
        K = _get_I(eltype(monomials(poly_problem(relaxation).objective)))

        solver = ConicBundle.CBMatrixCBSolver(verbose ? 1 : 0)
        ConicBundle.cb_set_augvalfailslimit!(ConicBundle.cb_get_terminator(ConicBundle.cb_get_solver(solver)), 100)
        state = StateMoment{K}(solver)
        primal_data = primal_moment_setup!(state, relaxation, groupings; verbose)
        ismissing(primal_data) && return missing, -1, typemin(V)
        state.info, state.data = primal_data
        customize(state)

        length(state.Cmats) > 1 && @warn("More than a single PSD constraint present: not well-supported")
        # In this case, we should use a BlockPSCPrimal of GramSparsePSCPrimal, but this is currently not implemented in the
        # Julia interface, as the data is passed as a std::map.
        state.psds = psds = ConicBundle.CBPSCAffineFunction(state.C, state.opAt,
            isone(length(state.Cmats)) ? ConicBundle.CBGramSparsePSCPrimal(state.Cmats[1]) : nothing)
        ConicBundle.cb_add_function!(solver, psds, abs(bound), bound > 0 ? cbft_constant_penalty_function :
                                                                           cbft_adaptive_penalty_function)
        ConicBundle.cb_set_term_relprec!(solver, 1e-8)
        ConicBundle.cb_set_active_bounds_fixing!(solver, true)
        finalizer(state) do _
            ConicBundle.cb_destroy!(solver, psds)
            for cm in state.cmats
                ConicBundle.cb_destroy!(cm)
            end
            for Cm in state.Cmats
                ConicBundle.cb_destroy!(Cm)
            end
        end
    end
    @verbose_info("Setup complete in ", setup_time, " seconds")
    ret = ConicBundle.cb_solve!(solver)
    status = ConicBundle.cb_termination_code(solver)
    value = -ConicBundle.cb_get_objval(solver)
    @verbose_info("Optimization complete")

    return state, status, value
end

function Solver.extract_moments(relaxation::AbstractRelaxation, state::StateMoment)
    isone(length(state.Cmats)) || error("Moment matrix extraction not implemented for more than a single PSD cone")
    # Better: implement our own type (and eliminate the OffsetArrays dependency)
    data = ConicBundle.CBGramSparsePSCPrimal(ConicBundle.cb_get_approximate_primal(state.solver, state.psds))
    gram = ConicBundle.cb_get_grammatrix(data)
    gramm = unsafe_wrap(Array, ConicBundle.cb_get_store(gram), (ConicBundle.cb_rowdim(gram), ConicBundle.cb_coldim(gram)))
    mom = OffsetArrays.Origin(0)(Symmetric(BLAS.syrk('U', 'N', 1., gramm), :U))
    mv = MomentVector(relaxation, state.data, (mom,), nothing)
end