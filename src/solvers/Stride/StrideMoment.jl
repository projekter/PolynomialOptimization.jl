mutable struct StateMoment{K<:Integer,V<:Real} <: AbstractSparseMatrixSolver{Int,K,V}
    const A::FastVec{SparseMatrixCSC{V,Int}}
    const C::FastVec{SPMatrix{V,SparseVector{V,Int},:LS}}
    b::SparseVector{V,Int}
    info::Vector{<:Vector{<:Tuple{Symbol,Any}}}
    data

    StateMoment{K,V}() where {K<:Integer,V<:Real} = new{K,V}(
        FastVec{SparseMatrixCSC{V,Int}}(),
        FastVec{SPMatrix{V,SparseVector{V,Int},:LS}}(),
    )
end

Solver.issuccess(::Val{:StrideMoment}, status::Symbol) = status === :ok

Solver.psd_indextype(::StateMoment{<:Integer,V}) where {V<:Real} = PSDIndextypeCOOVectorized(:L, sqrt(V(2)), 1)

@counter_atomic(StateMoment, :psd)

Solver.add_var_nonnegative!(state::StateMoment{<:Integer,V}, m::Int, dim::Int, data::SparseMatrixCOO{Int,Int,V,1},
    obj::Union{Nothing,Tuple{FastVec{Int},FastVec{V}}}) where {V<:Real} =
    error("Size-1 constraints are not supported by Stride")

function Solver.add_var_psd!(state::StateMoment{<:Integer,V}, m::Int, dim::Int, data::SparseMatrixCOO{Int,Int,V,1},
    obj::Union{Nothing,Tuple{FastVec{Int},FastVec{V}}}) where {V<:Real}
    dtr = trisize(dim)
    push!(state.A, SparseMatrixCSC(m, dtr, coo_to_csc!(dtr, data)...))
    cvec = isnothing(obj) ? spzeros(V, Int, dtr) : SparseVector(dtr, finish!.(obj)...)
    push!(state.C, SPMatrix(dim, cvec, :LS))
    return
end

function Solver.fix_constraints!(state::StateMoment{<:Integer,V}, m::Int, indvals::Indvals{Int,V}) where {V<:Real}
    state.b = SparseVector(m, finish!(indvals.indices), finish!(indvals.values))
    return
end

function Solver.poly_optimize(::Val{:StrideMoment}, relaxation::AbstractRelaxation, groupings::RelaxationGroupings;
    representation, verbose::Bool=false, customize=_ -> nothing, parameters...)
    setup_time = @elapsed begin
        K = _get_I(eltype(monomials(poly_problem(relaxation).objective)))
        V = real(coefficient_type(poly_problem(relaxation).objective))

        state = StateMoment{K,V}()
        primal_data = primal_moment_setup!(state, relaxation, groupings; verbose)
        ismissing(primal_data) && return missing, :infeasible_or_unbounded, typemin(V)
        state.info, state.data = primal_data
        customize(state)

        num_con = primal_data[2].num_con
        A = finish!(state.A)
        Astar = similar(A)
        AAstars = similar(A)
        @inbounds for (i, Aᵢ) in enumerate(A)
            Astar[i] = transpose(Aᵢ)
            AAstars[i] = Aᵢ * Astar[i]
        end
        # We must add together all the AAstar in one matrix, which is then factorized. However, these are CSC matrices whose
        # addition is highly inefficient. But they may still be large, so we don't want to densify them, either. So we follow a
        # multi-step path, where we first calculate the sparsity pattern of the output before we do the actual addition.
        AAstarinv = @inbounds let
            AAstarcolsize = 0
            for AAstarᵢ in AAstars
                AAstarcolsize += AAstarᵢ.colptr[end] -1
            end
            # now we have certain upper bounds on how large the individual columns can become.
            AAstarrowval = Vector{Int}(undef, AAstarcolsize)
            AAstarnzval = similar(AAstarrowval, Float64)
            AAstarcolptr = Vector{Int}(undef, num_con +1)
            j = 1
            for col in 1:num_con
                AAstarcolptr[col] = j
                # First, we add the elements all consecutively, don't care about the ordering or possible duplicates.
                firstpos = j
                for AAstarᵢ in AAstars
                    colptr = AAstarᵢ.colptr
                    srcstart = colptr[col]
                    len = colptr[col+1] - srcstart
                    copyto!(AAstarrowval, j, AAstarᵢ.rowval, srcstart, len)
                    copyto!(AAstarnzval, j, AAstarᵢ.nzval, srcstart, len)
                    j += len
                end
                j ≤ firstpos +1 && continue # j == firstpos must be caught, and j == firstpos +1 is a short-circuit
                # Then, we sort all the elements
                sort_along!(AAstarrowval, AAstarnzval, lo=firstpos, hi=j -1)
                # Finally, we take care of duplicates
                for p in firstpos+1:j-1
                    if AAstarrowval[firstpos] == AAstarrowval[p]
                        AAstarnzval[firstpos] += AAstarnzval[p]
                    else
                        firstpos += 1
                        AAstarrowval[firstpos] = AAstarrowval[p]
                        AAstarnzval[firstpos] = AAstarnzval[p]
                    end
                end
                j = firstpos +1
            end
            AAstarcolptr[num_con+1] = j
            resize!(AAstarrowval, j -1)
            resize!(AAstarnzval, j -1)
            # We cannot simply use cholesky, as A A* will probably be singular.
            # ldlt might work, but can also lead to zero pivot values, so this will sometimes also fail.
            AAstarinv = try
                EfficientCholmod(ldlt(SparseMatrixCSC(num_con, num_con, AAstarcolptr, AAstarrowval, AAstarnzval)))
            catch
                qr(SparseMatrixCSC(m, m, AAstarcolptr, AAstarrowval, AAstarnzval))
            end
        end

        psds = primal_data[2].psd_dim
        ssd = Stride.Data(
            num_con, A, Astar, AAstarinv, finish!(state.C), (I, SparseVector{V,Int}(num_con, Int[], V[])),
            state.b, zero(V),
            [SPMatrix{V}(undef, nᵢ, :LS) for nᵢ in psds],
            [SPMatrix{V}(undef, nᵢ, :LS) for nᵢ in psds],
            [SPMatrix{V}(undef, nᵢ, :LS) for nᵢ in psds],
            [SPMatrix{V}(undef, nᵢ, :LS) for nᵢ in psds],
            Vector{V}(undef, num_con),
            Vector{V}(undef, num_con),
            Vector{V}(undef, max(trisize(maximum(psds, init=0)), num_con)),
        )
    end
    @verbose_info("Setup complete in ", setup_time, " seconds")

    status, value = stride_solve(ssd, relaxation, primal_data[2]; verbose, parameters...)

    return (state, ssd), status, value
end

Solver.extract_moments(relaxation::AbstractRelaxation, (state, ssd)::Tuple{StateMoment,Stride.Data}) =
    MomentVector(relaxation, state.data, ssd.X, nothing)

Solver.extract_sos(::AbstractRelaxation, (_, ssd)::Tuple{StateMoment,Stride.Data}, ::Val{:psd},
    index::Integer, ::Nothing) = ssd.S[index]