include("SpecBM/SpecBM.jl")

struct SpecBMState{M,R}
    psds::FastVec{Int}
    A::Tuple{FastVec{Int},FastVec{Int},FastVec{R}}
    constraint_mappings::Dict{M,Int}
end

get_constraint!(sbm::SpecBMState{M}, mon::M) where {M} =
    get!(sbm.constraint_mappings, mon, length(sbm.constraint_mappings) +1)

function sos_matrix!(sbm::SpecBMState{M,R}, gröbner_basis, grouping::AbstractVector{M}, constraint::P, psd::Union{Val{true},Val{false}}) where {P,M,R}
    lg = length(grouping)
    psd isa Val{true} && push!(sbm.psds, lg)
    items = (lg * (lg + 1)) >> 1
    i = isempty(sbm.A[2]) ? 1 : last(sbm.A[2]) +1
    prepare_push!.(sbm.A, items * nterms(constraint))
    for exp2 in 1:lg
        for exp1 in exp2:lg
            sqr = grouping[exp1] * grouping[exp2]
            for mon_constr in constraint
                @inbounds for term in rem(sqr * mon_constr, gröbner_basis)
                    if gröbner_basis isa EmptyGröbnerBasis
                        unsafe_push!(sbm.A[1], get_constraint!(sbm, monomial(term)))
                        unsafe_push!(sbm.A[2], i)
                        unsafe_push!(sbm.A[3], psd isa Val{false} || exp2 == exp1 ? R(coefficient(term)) : sqrt(R(2)) * coefficient(term))
                    else
                        push!(sbm.A[1], get_constraint!(sbm, monomial(term)))
                        push!(sbm.A[2], i)
                        push!(sbm.A[3], psd isa Val{false} || exp2 == exp1 ? R(coefficient(term)) : sqrt(R(2)) * coefficient(term))
                    end
                end
            end
            i += 1
        end
    end
    return
end

function sos_matrix!(sbm::SpecBMState{M,R}, gröbner_basis, grouping::AbstractVector{M}, constraint::AbstractMatrix{P}, ::Val{true}) where {P,M,R}
    block_size = LinearAlgebra.checksquare(constraint)
    isone(block_size) && @inbounds return sos_matrix!(sbm, gröbner_basis, grouping, constraint[1, 1])
    lg = length(grouping)
    dim = lg * block_size
    push!(sbm.psds, dim)
    items = (dim * (dim +1)) >> 1
    @assert(!isempty(sbm.A[2])) # the moment matrix is already there
    i = last(sbm.A[2]) +1
    prepare_push!.(sbm.A, items * sum(nterms, constraint, init=0))
    for exp2 in 1:lg
        for block_j in 1:block_size
            for exp1 in exp2:lg
                sqr = grouping[exp1] * grouping[exp2]
                for block_i in (exp1 == exp2 ? block_j : 1):block_size
                    for mon_constr in constraint[block_i, block_j]
                        for term in rem(sqr * mon_constr, gröbner_basis)
                            if gröbner_basis isa EmptyGröbnerBasis
                                unsafe_push!(sbm.A[1], get_constraint!(sbm, monomial(term)))
                                unsafe_push!(sbm.A[2], i)
                                unsafe_push!(sbm.A[3], exp2 == exp1 && block_j == block_i ? R(coefficient(term)) :
                                    sqrt(R(2)) * coefficient(term))
                            else
                                push!(sbm.A[1], get_constraint!(sbm, monomial(term)))
                                push!(sbm.A[2], i)
                                push!(sbm.A[3], exp2 == exp1 && block_j == block_i ? R(coefficient(term)) :
                                    sqrt(R(2)) * coefficient(term))
                            end
                        end
                    end
                    i += 1
                end
            end
        end
    end
    return
end

function sparse_optimize(::Val{:SpecBMSOS}, problem::PolyOptProblem{P,M,V}, groupings::Vector{<:Vector{<:AbstractVector{M}}};
    verbose::Bool=false, kwargs...) where {P,M,V}
    @assert(!problem.complex)
    R = coefficient_type(problem.objective)
    R <: AbstractFloat || (R = Float64)
    setup_time = @elapsed begin
        sbm = SpecBMState(
            FastVec{Int}(buffer=length(groupings)),
            (FastVec{Int}(), FastVec{Int}(), FastVec{R}()),
            Dict{M,Int}()
        )
        # add lower bound. This must be done in the beginning, so that the free variable for this is created first.
        if isone(problem.prefactor)
            push!(sbm.A[1], get_constraint!(sbm, constant_monomial(problem.objective)))
            push!(sbm.A[2], 1)
            push!(sbm.A[3], one(R))
        else
            append!(sbm.A[1], get_constraint!(sbm, monomial(x)) for x in problem.prefactor)
            append!(sbm.A[2], Iterators.repeated(1, length(problem.prefactor)))
            append!(sbm.A[3], R.(coefficients(problem.prefactor)))
        end

        # We need to put equality constraints (-> free variables) to the front, so we just take care of them beforehand,
        # avoiding moving data. This is also why we don't need a sos_matrix_eq! function - the only difference is the
        # characterization of the free variables
        for (groupings, constr) in zip(Iterators.drop(groupings, 1), problem.constraints)
            if constr.type == pctEqualityGröbner
                for grouping in groupings
                    @assert(issorted(grouping, by=degree))
                    sos_matrix!(sbm, EmptyGröbnerBasis{P}(), grouping, constr.constraint, Val(false))
                end
            elseif constr.type == pctEqualitySimple
                for grouping in groupings
                    @assert(issorted(grouping, by=degree))
                    sos_matrix!(sbm, problem.gröbner_basis, grouping, constr.constraint, Val(false))
                end
            end
        end
        num_frees = last(sbm.A[2])
        # SOS term for objective
        for grouping in groupings[1]
            @assert(issorted(grouping, by=degree))
            sos_matrix!(sbm, problem.gröbner_basis, grouping, polynomial(constant_monomial(problem.objective)), Val(true))
        end
        # localizing matrices
        for (groupings, constr) in zip(Iterators.drop(groupings, 1), problem.constraints)
            if constr.type == pctNonneg || constr.type == pctPSD
                for grouping in groupings
                    sos_matrix!(sbm, problem.gröbner_basis, grouping, constr.constraint, Val(true))
                end
            elseif constr.type == pctEqualityNonneg
                for grouping in groupings
                    @assert(issorted(grouping, by=degree))
                    sos_matrix!(sbm, problem.gröbner_basis, grouping, constr.constraint, Val(true))
                    sos_matrix!(sbm, problem.gröbner_basis, grouping, -constr.constraint, Val(true))
                end
            end
        end
        nvars = last(sbm.A[2])
        nconstrs = length(sbm.constraint_mappings)

        A, At = let A_coo=finish!.(sbm.A)
            @assert(nconstrs == maximum(A_coo[1]))
            coolen = length(A_coo[1])
            csrrowptr = Vector{Int}(undef, nconstrs +1)
            csrcolval = Vector{Int}(undef, coolen)
            csrnzval = Vector{R}(undef, coolen)
            A = SparseArrays.sparse!(A_coo..., nconstrs, nvars, +, Vector{Int}(undef, nvars), csrrowptr, csrcolval, csrnzval,
                A_coo...)
            # Now the csr data already contain the CSR representation of A - which, when seen as CSC, corresponds to the
            # transpose. So we get Aᵀ almost for free - however, the columns are still unsorted, so we have to do the
            # sorting.
            for (from, toplus1) in zip(csrrowptr, Iterators.drop(csrrowptr, 1))
                to = toplus1 -1
                sort_along!(csrcolval, from, to, Base.Forward, csrnzval)
            end
            A, SparseMatrixCSC(nvars, nconstrs, csrrowptr, csrcolval, csrnzval)
        end
        c = SparseVector{R,Int}(nvars, [1], [-one(R)])
        # Will we work with sparse vectors at all?
        if 5length(problem.objective) < nconstrs
            b_idx = FastVec{Int}(buffer=length(problem.objective))
            b_val = FastVec{R}(buffer=length(problem.objective))
            for x in problem.objective
                unsafe_push!(b_idx, get_constraint!(sbm, monomial(x)))
                unsafe_push!(b_val, R(coefficient(x)))
            end
            b = sparsevec(finish!(b_idx), finish!(b_val), nconstrs, +)
        else
            b = zeros(nconstrs)
            for x in problem.objective
                @inbounds b[get_constraint!(sbm, monomial(x))] = coefficient(x)
            end
        end
    end
    @verbose_info("Setup complete in ", setup_time, " seconds")

    # Now we need to figure out some proper default arguments. Let us rely on the adaptive strategies in the solver to do this
    # instead of some sophisticated analysis that might give us better values provided the problem features certain
    # constraints.
    args = merge!(Dict{Symbol,Any}(:ρ => R(10), :adaptiveρ => !haskey(kwargs, :ρ), :r_past => 0, :r_current => 3,
        :ϵ => R(1e-5)), kwargs)

    solve_time = @elapsed begin
        result = specbm_primal(A, b, c; num_frees, psds=finish!(sbm.psds), verbose, At, args...)
    end
    @verbose_info("Optimization complete in ", solve_time, " seconds")

    mon_vals = result.y
    empty!(problem.last_moments)
    sizehint!(problem.last_moments, nconstrs +1)
    problem.last_moments[constant_monomial(problem.objective)] = 1.
    for (mon, i) in sbm.constraint_mappings
        push!(problem.last_moments, mon => -Float64(mon_vals[i]))
    end
    return result.status, -result.objective
end

for sp in SparsityAny
    @eval default_solution_method(::$sp, ::Val{:SpecBMSOS}) = :heuristic
end