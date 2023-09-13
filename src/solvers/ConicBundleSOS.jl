struct CBStateSOS{M}
    solver::ConicBundle.CBMatrixCBSolver
    var_dim::ConicBundle.CBIndexmatrix
    constraint_ops::Dict{M,Dict{Int,Tuple{Vector{Int32},Vector{Int32},Vector{Float64}}}}
    obj_ops::Dict{Int,Tuple{Vector{Int32},Vector{Int32},Vector{Float64}}}
end

get_constraint!(cb::CBStateSOS{M}, mon::M) where {M} =
    get!(cb.constraint_ops, mon) do
        return valtype(cb.constraint_ops)()
    end

const unitonedata = ConicBundle.CBSymmatrix(1, 1.)
const unitone = ConicBundle.CBCMsymdense(unitonedata)

function sos_matrix!(cb::CBStateSOS{M}, gröbner_basis, grouping::AbstractVector{M}, constraint::P) where {P,M}
    lg = length(grouping)
    ConicBundle.cb_enlarge_below!(cb.var_dim, 1, lg)
    mat_idx = ConicBundle.cb_rowdim(cb.var_dim) -1
    for exp2 in 1:lg
        for exp1 in exp2:lg
            sqr = grouping[exp1] * grouping[exp2]
            for mon_constr in constraint
                @inbounds for term in rem(sqr * mon_constr, gröbner_basis)
                    if isconstant(term)
                        constr_rs, constr_cs, constr_vs = get!(() -> (Int32[], Int32[], Float64[]), cb.obj_ops, mat_idx)
                        push!(constr_vs, -Float64(coefficient(term)))
                    else
                        constr = get_constraint!(cb, monomial(term))
                        constr_rs, constr_cs, constr_vs = get!(() -> (Int32[], Int32[], Float64[]), constr, mat_idx)
                        push!(constr_vs, Float64(coefficient(term)))
                    end
                    push!(constr_rs, exp2 -1)
                    push!(constr_cs, exp1 -1)
                end
            end
        end
    end
    return
end

function sos_matrix!(cb::CBStateSOS{M}, gröbner_basis, grouping::AbstractVector{M}, constraint::AbstractMatrix{P}) where {P,M}
    @assert size(constraint, 1) == size(constraint, 2)
    block_size = size(constraint, 1)
    if block_size == 1
        @inbounds return sos_matrix!(cb, gröbner_basis, grouping, constraint[1, 1])
    end
    lg = length(grouping)
    dim = lg * block_size
    ConicBundle.cb_enlarge_below!(cb.var_dim, 1, dim)
    mat_idx = ConicBundle.cb_rowdim(cb.var_dim) -1

    for exp2 in 1:lg
        exp2idx = (exp2 -1) * block_size
        for block_j in 1:block_size
            for exp1 in exp2:lg
                sqr = grouping[exp1] * grouping[exp2]
                exp1idx = (exp1 -1) * block_size
                for block_i in (exp1 == exp2 ? block_j : 1):block_size
                    @inbounds for mon_constr in constraint[block_i, block_j]
                        @inbounds for term in rem(sqr * mon_constr, gröbner_basis)
                            if isconstant(term)
                                constr_rs, constr_cs, constr_vs = get!(() -> (Int32[], Int32[], Float64[]), cb.obj_ops, mat_idx)
                                push!(constr_vs, -Float64(coefficient(term)))
                            else
                                constr = get_constraint!(cb, monomial(term))
                                constr_rs, constr_cs, constr_vs = get!(() -> (Int32[], Int32[], Float64[]), constr, mat_idx)
                                push!(constr_vs, Float64(coefficient(term)))
                            end
                            push!(constr_rs, exp2idx + block_j -1)
                            push!(constr_cs, exp1idx + block_i -1)
                        end
                    end
                end
            end
        end
    end
    return
end

function sparse_optimize(::Union{Val{:ConicBundleSOS},Val{:CBSOS}}, problem::PolyOptProblem{P,M,V},
    groupings::Vector{<:Vector{<:AbstractVector{M}}}; verbose::Bool=false, customize::Function=(cb, so) -> nothing,
    parameters...) where {P,M,V}
    @assert(!problem.complex)
    @assert isone(problem.prefactor)

    setup_time = @elapsed begin
        cb = CBStateSOS{M}(
            ConicBundle.CBMatrixCBSolver(verbose ? 1 : 0),
            ConicBundle.CBIndexmatrix(0, 1),
            Dict{M,Dict{Int,Tuple{Vector{Int32},Vector{Int32},Vector{Float64}}}}(),
            Dict{Int,Tuple{Vector{Int32},Vector{Int32},Vector{Float64}}}())
        solver = cb.solver
        term = ConicBundle.cb_get_terminator(ConicBundle.cb_get_solver(solver))
        ConicBundle.cb_set_augvalfailslimit!(term, 100)

        # SOS term for objective
        for grouping in groupings[1]
            sos_matrix!(cb, problem.gröbner_basis, sort(grouping, by=degree),
                polynomial(constant_monomial(problem.objective)))
        end
        # localizing matrices
        for (groupings, constr) in zip(Iterators.drop(groupings, 1), problem.constraints)
            if constr.type == pctNonneg || constr.type == pctPSD
                for grouping in groupings
                    sos_matrix!(cb, problem.gröbner_basis, sort(grouping, by=degree), constr.constraint)
                end
            elseif constr.type == pctEqualityNonneg
                for grouping in groupings
                    let sg = sort(grouping, by=degree)
                        sos_matrix!(cb, problem.gröbner_basis, sg, constr.constraint)
                        sos_matrix!(cb, problem.gröbner_basis, sg, -constr.constraint)
                    end
                end
            elseif constr.type == pctEqualityGröbner
                for grouping in groupings
                    sos_matrix_eq!(cb, EmptyGröbnerBasis{P}(), sort(grouping, by=degree), constr.constraint)
                end
            elseif constr.type == pctEqualitySimple
                for grouping in groupings
                    sos_matrix_eq!(cb, problem.gröbner_basis, sort(grouping, by=degree), constr.constraint)
                end
            else
                @assert(false)
            end
        end

        ConicBundle.cb_enlarge_below!(cb.var_dim, 1, 1) # enlarge with one innocent scalar for trace purposes
        # turn dict constraints into sparse matrix objects (need to keep them all in a vector for GC purposes)
        sortops = sort(cb.constraint_ops, lt=(a, b) -> begin
            if degree(a) == degree(b)
                for v in @view(problem.variables[end:-1:begin])
                    dav, dbv = degree(a, v), degree(b, v)
                    dav == dbv || return dav < dbv
                end
                return false
            else
                return degree(a) < degree(b)
            end
        end)
        opAt = ConicBundle.CBSparseCoeffmatMatrix(cb.var_dim, length(sortops))
        cmats = FastVec{ConicBundle.CBCMsymsparse}(buffer=sum(length ∘ last, sortops, init=0) + length(cb.obj_ops))
        for (i, (_, constr)) in zip(Iterators.countfrom(0), sortops)
            for (j, (constr_rs, constr_cs, constr_vs)) in constr
                @assert(length(constr_rs) == length(constr_cs) == length(constr_vs))
                mat = ConicBundle.CBSparsesym(cb.var_dim[j], length(constr_rs), constr_rs, constr_cs, constr_vs)
                cmat = ConicBundle.CBCMsymsparse(mat)
                ConicBundle.cb_destroy!(mat)
                unsafe_push!(cmats, cmat)
                ConicBundle.cb_set!(opAt, j, i, cmat)
            end
        end
        C = ConicBundle.CBSparseCoeffmatMatrix(cb.var_dim, 1)
        for (j, (constr_rs, constr_cs, constr_vs)) in cb.obj_ops
            @assert(length(constr_rs) == length(constr_cs) == length(constr_vs))
            mat = ConicBundle.CBSparsesym(cb.var_dim[j], length(constr_rs), constr_rs, constr_cs, constr_vs)
            cmat = ConicBundle.CBCMsymsparse(mat)
            ConicBundle.cb_destroy!(mat)
            unsafe_push!(cmats, cmat)
            ConicBundle.cb_set!(C, j, 0, cmat)
        end

        # add objective (which is already mod gröbner_basis)
        bmat_data = ∘(-, Float64).(coefficients(problem.objective, collect(keys(sortops))))
        bmat = ConicBundle.CBMatrix(length(sortops), 1, bmat_data)
        ConicBundle.cb_init_problem!(solver, length(sortops), nothing, nothing, bmat)

        # construct function
        psds = ConicBundle.CBPSCAffineFunction(C, opAt)
        ConicBundle.cb_add_function!(solver, psds, 100., ConicBundle.cbft_adaptive_penalty_function)
        ConicBundle.cb_set_term_relprec!(solver, 1e-8)
        ConicBundle.cb_set_active_bounds_fixing!(solver, true)
    end
    @verbose_info("Setup complete in ", setup_time, " seconds")

    customize(cb, sortops)

    ret = ConicBundle.cb_solve!(solver)
    status = ConicBundle.cb_termination_code(solver)
    value = coefficient(problem.objective, constant_monomial(problem.objective)) + ConicBundle.cb_get_objval(solver)
    @verbose_info("Optimization complete, code ", ret)

    empty!(problem.last_moments)
    sizehint!(problem.last_moments, length(cb.constraint_ops))
    mon_vals = ConicBundle.CBMatrix()
    ConicBundle.cb_get_center(solver, mon_vals)
    for (i, mon) in enumerate(keys(cb.constraint_ops))
        push!(problem.last_moments, mon => mon_vals[i])
    end

    ConicBundle.cb_destroy!(mon_vals, solver, psds, bmat, cmats..., cb.var_dim)
    return status, value
end