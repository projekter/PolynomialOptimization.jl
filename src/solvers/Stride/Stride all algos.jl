# This is an implementation of the STRIDE solver, doi:10.1007/s10107-022-01912-6, tightly integrated with the
# PolynomialOptimization framework
# This implementation still contains remnants (all? commented out) of constructing the local problem with Optim, which is not
# a good idea.
using Printf
import SuiteSparse: CHOLMOD
import StaticPolynomials
import Manifolds
import Manopt

include("./PackedMatrices.jl")

const bfgs_optim = true

bfgs_optim && import Optim

mutable struct StrideState{M,R}
    const monomial_mapping::Dict{M,Tuple{Int,Int}}
    const psds::Vector{Int}
    constraints::Int
    const A::Vector{Tuple{Vector{Int},Vector{Int},Vector{R}}}
end

struct StrideData{R<:Real,F<:Factorization}
    d::Int
    A::Vector{SparseMatrixCSC{R,Int}}
    Astar::Vector{SparseMatrixCSC{R,Int}}
    AAstarinv::F
    C::Vector{PackedMatrix{R,SparseVector{R,Int}}}
    X::Vector{PackedMatrix{R,Vector{R}}} # each dimensions nᵢ × nᵢ
    S::Vector{PackedMatrix{R,Vector{R}}}
    tmpm1::Vector{PackedMatrix{R,Vector{R}}}
    tmpm2::Vector{PackedMatrix{R,Vector{R}}}
    y::Vector{R} # length d
    tmpv::Vector{R}
    tmpvlong::Vector{R} # length max(nᵢ)*(max(nᵢ)+1)÷2
end

# Begin enhancement of SuiteSparse: use the factorization functions that don't allocate
struct EfficientCholmod{T,F<:Factorization{T}} <: Factorization{T}
    cholmod::F
    Y::Ref{Ptr{CHOLMOD.cholmod_dense}}
    E::Ref{Ptr{CHOLMOD.cholmod_dense}}

    function EfficientCholmod(f::CHOLMOD.Factor)
        result = new{eltype(f),typeof(f)}(f, Ref{Ptr{CHOLMOD.cholmod_dense}}(C_NULL), Ref{Ptr{CHOLMOD.cholmod_dense}}(C_NULL))
        finalizer(f) do _
            let result = result
                common = CHOLMOD.COMMONS[Threads.threadid()]
                CHOLMOD.cholmod_free_dense(result.Y, common)
                CHOLMOD.cholmod_free_dense(result.E, common)
            end
        end
        return result
    end
end

function cholmod_dense_wrap(A::StridedVecOrMat{T}) where {T}
    out = CHOLMOD.cholmod_dense()
    out.nrow = size(A, 1)
    out.ncol = size(A, 2)
    out.d = Base.stride(A, 2)
    out.nzmax = out.d * out.ncol
    out.x = pointer(A)
    out.z = C_NULL
    out.xtype = CHOLMOD.xtyp(T)
    out.dtype = CHOLMOD.COMMONS[Threads.threadid()][].dtype
    return out
end

function LinearAlgebra.ldiv!(out::StridedVecOrMat{T}, lhs::EfficientCholmod{T,<:Factorization{T}},
    rhs::StridedVecOrMat{T}) where {T<:CHOLMOD.VTypes}
    F = lhs.cholmod
    if size(F,1) != size(rhs,1)
        throw(DimensionMismatch("LHS and RHS should have the same number of rows. " *
            "LHS has $(size(F,1)) rows, but RHS has $(size(rhs,1)) rows."))
    end
    if size(F,2) != size(out,1) # F is square, but anyway...
        throw(DimensionMismatch("Output must have the same number of rows as LHS has columns. " *
            "LHS has $(size(F,2)) columns, but output has $(size(out,1)) rows."))
    end
    if size(out,2) != size(rhs,2)
        throw(DimensionMismatch("Output and RHS should have the same number of column." *
            "Output has $(size(out,2)) columns, but RHS has $(size(rhs,2)) columns."))
    end
    if !CHOLMOD.issuccess(F)
        s = unsafe_load(pointer(F))
        if s.is_ll == 1
            throw(LinearAlgebra.PosDefException(s.minor))
        else
            throw(LinearAlgebra.ZeroPivotException(s.minor))
        end
    end
    # Casting via CHOLMOD.Dense{T}(rhs) will actually copy around all the data.
    B = cholmod_dense_wrap(rhs)
    X = cholmod_dense_wrap(out)
    pX = Ref(X)
    GC.@preserve pX CHOLMOD.cholmod_l_solve2(
        CHOLMOD.CHOLMOD_A, # system to solve
        F,                 # factorization to use
        Ref(B), C_NULL,    # right-hand-side, dense and sparse
        Ref(Base.unsafe_convert(Ptr{CHOLMOD.cholmod_dense}, pX)),
        C_NULL,            # solution, dense and sparse
        lhs.Y, lhs.E,      # workspace
        CHOLMOD.COMMONS[Threads.threadid()]
    ) == 0 && error("Error in cholmod")
    return out
end

LinearAlgebra.ldiv!(lhs::EfficientCholmod{T,<:Factorization{T}}, rhs::StridedVecOrMat{T}) where {T<:CHOLMOD.VTypes} =
    ldiv!(rhs, lhs, rhs) # this appears to work, although it is not documented
# End enhancement of SuiteSparse

# Begin enhancement of StaticPolynomials: provide a Hessian for polynomials systems by directly taking products with a vector
#=@generated function StaticPolynomials._hessian!(U::AbstractMatrix{T}, λ::AbstractVector{T}, F::StaticPolynomials.PolynomialSystem{N, NVars}, x...) where {T, N, NVars}
    quote
        @boundscheck checkbounds(U, 1:$NVars, 1:$NVars)
        @boundscheck checkbounds(λ, 1:$N)
        @myinbounds begin
            $(map(1:N) do i
                quote
                    $(Symbol("H", i)) = StaticPolynomials.hessian(F.polys[$i], x...)
                    @myinbounds @simd for k=1:$(NVars*NVars)
                        U[k] += λ[$i] * $(Symbol("H", i))[k]
                    end
                end
            end...)
        end
        U
    end
end

Base.@propagate_inbounds StaticPolynomials.hessian!(U, λ::AbstractVector, F::StaticPolynomials.PolynomialSystem, x::AbstractVector) =
    StaticPolynomials._hessian!(U, λ, F, x)=#
# End enhancement of StaticPolynomials

macro myinbounds(expr)
    esc(expr)
end

function get_varrepr!(sst::StrideState{M}, gröbner_basis, grouping::AbstractVector{M}, row, col, constraint::P) where {P,M}
    @myinbounds rowcol = grouping[row] * grouping[col]
    constr_terms = mergewith(+, (monomial(term) => coefficient(term)
                                 for mon_constr in constraint
                                 for term in rem(rowcol * mon_constr, gröbner_basis)), M, Float64)
    cur_index = rowcol_to_vec(row, col)
    cur_A = last(sst.A)
    # We will first call this on the moment matrix; and since everything is degree-ordered, we are guaranteed to already have
    # defined all the necessary single-monomial entries. However, note that the moment matrix may be split in various parts due
    # to our sparsity pattern. TODO: This is sad. Sparsity may invalidate the approach, as not every monomial may occur in one
    # of the moment matrices, but they might still occur in localizing matrices, probably as terms in a polynomial - which will
    # lead to an error.
    # Note that the value of A(X) can be calculated using only the upper triangle, which is what we do and therefore do not
    # have factors of 1/2. However, the effective coefficient matrices Aᵢ are supposed to be symmetric, which means that we
    # have to put the 1/2 factors in A*.
    @myinbounds if length(constr_terms) == 1
        mon, val = first(constr_terms)
        if haskey(sst.monomial_mapping, mon)
            sst.constraints += 1
            push!(cur_A[1], sst.constraints)
            push!(cur_A[2], cur_index)
            push!(cur_A[3], -1.)
            ref_psd, ref_index = sst.monomial_mapping[mon]
            ref_A = sst.A[ref_psd]
            push!(ref_A[1], sst.constraints)
            push!(ref_A[2], ref_index)
            push!(ref_A[3], val)
        else
            @assert(isone(val))
            sst.monomial_mapping[mon] = (length(sst.psds), cur_index)
        end
    else
        sst.constraints += 1
        push!(cur_A[1], sst.constraints)
        push!(cur_A[2], cur_index)
        push!(cur_A[3], -1.)
        for (mon, val) in constr_terms
            ref_psd, ref_index = sst.monomial_mapping[mon]
            ref_A = sst.A[ref_psd]
            push!(ref_A[1], sst.constraints)
            push!(ref_A[2], ref_index)
            push!(ref_A[3], val)
        end
    end
    return
end

function moment_matrix!(sst::StrideState{M}, gröbner_basis, grouping::AbstractVector{M}, constraint::P) where {P,M}
    lg = length(grouping)
    items = (lg * (lg +1)) >> 1
    push!(sst.psds, UInt(lg))
    # These sizehints are just rough estimates, we cannot (with reasonable effort) precompute the actual sizes.
    push!(sst.A, (sizehint!(Int[], items), sizehint!(Int[], items), sizehint!(Float64[], items)))
    if length(sst.psds) == 1
        cur_A = last(sst.A)
        push!(cur_A[1], 1)
        push!(cur_A[2], 1)
        push!(cur_A[3], 1.) # constant monomial
    end
    for exp2 in 1:lg
        for exp1 in 1:exp2
            get_varrepr!(sst, gröbner_basis, grouping, exp1, exp2, constraint)
        end
    end
    return
end

function moment_matrix!(sst::StrideState{M}, gröbner_basis, grouping::AbstractVector{M}, constraint::AbstractMatrix{P}) where {P,M}
    block_size = LinearAlgebra.checksquare(constraint)
    if block_size == 1
        @inbounds return moment_matrix!(sst, gröbner_basis, grouping, constraint[1, 1])
    end
    lg = length(grouping)
    dim = lg * block_size
    items = (dim * (dim +1)) >> 1
    push!(sst.psds, UInt(lg))
    push!(sst.A, (sizehint!(Int[], items), sizehint!(Int[], items), sizehint!(Float64[], items)))
    for exp2 in 1:lg, block_j in 1:block_size, exp1 in 1:exp2, block_i in 1:(exp1 == exp2 ? block_j : block_size)
        @inbounds get_varrepr!(sst, gröbner_basis, grouping, exp1, exp2, constraint[block_i, block_j])
    end
    return
end

function moment_matrix_eq!(sst::StrideState{M}, gröbner_basis, grouping::AbstractVector{M}, constraint::P) where {P,M}
    lg = length(grouping)
    i = 1
    @myinbounds for exp2 in 1:lg
        for exp1 in 1:exp2
            rowcol = grouping[exp1] * grouping[exp2]
            constr_terms = mergewith(+, (monomial(term) => coefficient(term)
                                         for mon_constr in constraint
                                         for term in rem(rowcol * mon_constr, gröbner_basis)),
                M, Float64)
            sst.constraints += 1
            for (mon, val) in constr_terms
                ref_psd, ref_index = sst.monomial_mapping[mon]
                ref_A = sst.A[ref_psd]
                push!(ref_A[1], sst.constraints)
                push!(ref_A[2], ref_index)
                push!(ref_A[3], val)
            end
            i += 1
        end
    end
    return
end

function sort_along!(v::AbstractVector, lo::Integer, hi::Integer, o::Base.Ordering, along::AbstractVector)
    @myinbounds while lo < hi
        if hi - lo <= Base.SMALL_THRESHOLD
            # function sort!(v::AbstractVector, lo::Integer, hi::Integer, ::InsertionSortAlg, o::Ordering)
            lo_plus_1 = (lo +1)::Integer
            for i in lo_plus_1:hi
                j = i
                x, x_along = v[i], along[i]
                while j > lo
                    y = v[j-1]
                    if !(Base.lt(o, x, y)::Bool)
                        break
                    end
                    v[j] = y
                    along[j] = along[j-1]
                    j -= 1
                end
                v[j] = x
                along[j] = x_along
            end
            return v, along
            # end
        end
        # function partition!(v::AbstractVector, lo::Integer, hi::Integer, o::ordering)
            # function selectpivot!(v::AbstractVector, lo::Integer, hi::Integer, o::Ordering)
            mi = midpoint(lo, hi)
            if Base.lt(o, v[lo], v[mi])
                v[mi], v[lo] = v[lo], v[mi]
                along[mi], along[lo] = along[lo], along[mi]
            end
            if Base.lt(o, v[hi], v[lo])
                if lt(o, v[hi], v[mi])
                    v[hi], v[lo], v[mi] = v[lo], v[mi], v[hi]
                    along[hi], along[lo], along[mi] = along[lo], along[mi], along[hi]
                else
                    v[hi], v[lo] = v[lo], v[hi]
                    along[hi], along[lo] = along[lo], along[hi]
                end
            end
            pivot = v[lo]
            pivot_along = along[lo]
            # end
            i, j = lo, hi
            while true
                i += 1; j -= 1
                while Base.lt(o, v[i], pivot); i += 1; end;
                while Base.lt(o, pivot, v[j]); j -= 1; end;
                i >= j && break
                v[i], v[j] = v[j], v[i]
                along[i], along[j] = along[j], along[i]
            end
            v[j], v[lo] = pivot, v[j]
            along[j], along[lo] = pivot_along, along[j]
        # end
        if j - lo < hi - j
            lo < (j -1) && sort_along!(v, lo, j -1, o, along)
            lo = j +1
        else
            j +1 < hi && sort_along!(v, j +1, hi, o, along)
            hi = j -1
        end
    end
    return v, along
end

function prepare_relaxation(problem::PolyOptProblem{P,M,V}, groupings::Vector{<:Vector{<:AbstractVector{M}}},
    verbose::Bool) where {P,M,V}
    # How may constraints will we need? Every monomial that arises somewhere in the relaxation must be assigned the same value.
    # This means that the first occurrence is free, and all others give one constraint each.
    sst = StrideState{M,Float64}(Dict{M,Tuple{Int,Int,Int}}(), Int[], 1, Tuple{Vector{Int},Vector{Int},Vector{Float64}}[])
    # moment matrix
    @verbose_info("Assembling primal linear constraint function")
    for grouping in groupings[1]
        moment_matrix!(sst, problem.gröbner_basis, sort(grouping, by=degree), polynomial(constant_monomial(problem.objective)))
    end
    for (groupings, constr) in zip(Iterators.drop(groupings, 1), problem.constraints)
        if constr.type == pctNonneg || constr.type == pctPSD
            for grouping in groupings
                moment_matrix!(sst, problem.gröbner_basis, sort(grouping, by=degree), constr.constraint)
            end
        elseif constr.type == pctEqualityNonneg
            for grouping in groupings
                let sg = sort(grouping, by=degree)
                    moment_matrix!(sst, problem.gröbner_basis, sg, constr.constraint)
                    moment_matrix!(sst, problem.gröbner_basis, sg, -constr.constraint)
                end
            end
        elseif constr.type == pctEqualityGröbner
            for grouping in groupings
                moment_matrix_eq!(sst, EmptyGröbnerBasis{P}(), sort(grouping, by=degree), constr.constraint)
            end
        elseif constr.type == pctEqualitySimple
            for grouping in groupings
                moment_matrix_eq!(sst, problem.gröbner_basis, sort(grouping, by=degree), constr.constraint)
            end
        else
            @assert(false)
        end
    end
    npsd = length(sst.psds)
    maxdim = maximum(sst.psds) # should be first(sst.psds), but with sparsity, who knows...
    A = Vector{SparseMatrixCSC{Float64,Int}}(undef, npsd)
    Astar = Vector{SparseMatrixCSC{Float64,Int}}(undef, npsd)
    AAstars = Vector{SparseMatrixCSC{Float64,Int}}(undef, npsd)
    @verbose_info("Constructing sparse representations, duals, and factorizations")
    @myinbounds let tmp = Vector{Int}(undef, maxdim * (maxdim +1) ÷ 2), m = sst.constraints
        for (i, (Adata, n)) in enumerate(zip(sst.A, sst.psds))
            n = n * (n +1) ÷ 2
            coolen = length(Adata[1])
            csrrowptr = Vector{Int}(undef, m +1)
            csrcolval = Vector{Int}(undef, coolen)
            csrnzval = Vector{Float64}(undef, coolen)
            A[i] = SparseArrays.sparse!(Adata..., m, n, +, tmp, csrrowptr, csrcolval, csrnzval, Adata...)
            Avalscaled = copy(Adata[3])
            # Now the csr data already contain the CSR representation of A[i] - which, when seen as CSC, corresponds to the
            # transpose. So we get A* almost for free - however, the columns are still unsorted, so we have to do the sorting.
            for (from, toplus1) in zip(csrrowptr, Iterators.drop(csrrowptr, 1))
                to = toplus1 -1
                # We implement QuickSort directly, for then we can sort csrnzval along the way without requiring any
                # temporaries.
                sort_along!(csrcolval, from, to, Base.Forward, csrnzval)
                # We must also put factors of 1/2 to the off-diagonal elements - we didn't do this in A, as the 1/2 would be
                # compensated by an addition of the symmetric counterpart, and all that we care for in A is the sum. However,
                # for the dual, we are interested in the actual _matrix elements_.
                curcol = 1
                nextdiag = 1
                for j in from:to
                    el = csrcolval[j]
                    while el > nextdiag
                        curcol += 1
                        nextdiag += curcol
                    end
                    el == nextdiag || (csrnzval[j] *= .5)
                end
            end
            Astar[i] = SparseMatrixCSC(n, m, csrrowptr, csrcolval, csrnzval)
            AAstars[i] = A[i] * Astar[i]
        end
    end
    # We must add together all the AAstar in one matrix, which is then factorized. However, those are CSC matrices whose
    # addition is highly inefficient. But they may still be large, so we don't want to densify them, either. So we follow a
    # multi-step path, where we first calculate the sparsity pattern of the output before we do the actual addition.
    local AAstarinv
    @myinbounds let m = sst.constraints, AAstarcolsizes = zeros(Int, m +1) # one more which we don't need now, but use later
        for AAstarᵢ in AAstars
            colptr = AAstarᵢ.colptr
            for col in 1:m
                AAstarcolsizes[col] += colptr[col+1] - colptr[col]
            end
        end
        # now we have certain upper bounds on how large the individual columns can become.
        AAstarrowval = Vector{Int}(undef, sum(AAstarcolsizes))
        AAstarnzval = similar(AAstarrowval, Float64)
        AAstarcolptr = fill!(AAstarcolsizes, 0)
        j = 1
        for col in 1:m
            AAstarcolptr[col] = j
            # First, we add the elements all consecutively, don't care about the ordering or possible duplicates.
            firstpos = j
            for AAstarᵢ in AAstars
                colptr = AAstarᵢ.colptr
                srcstart = colptr[col]
                len = colptr[col+1] - colptr[col]
                copyto!(AAstarrowval, j, AAstarᵢ.rowval, srcstart, len)
                copyto!(AAstarnzval, j, AAstarᵢ.nzval, srcstart, len)
                j += len
            end
            # Then, we sort all the elements
            sort_along!(AAstarrowval, firstpos, j -1, Base.Forward, AAstarnzval)
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
        AAstarcolptr[m+1] = j
        resize!(AAstarrowval, j -1)
        resize!(AAstarnzval, j -1)
        # We cannot simply use cholesky, was A A* will probably be singular.
        AAstarinv = EfficientCholmod(ldlt(SparseMatrixCSC(m, m, AAstarcolptr, AAstarrowval, AAstarnzval)))
    end
    @verbose_info("Assembling objective constraint matrix")
    # Finally, we need our primal objective matrix. We take a simple approach in which we just use one entry in the moment
    # matrix instead of averaging over all the terms that the constraints enforced to be the same.
    local C
    let
        I = [Int[] for _ in 1:npsd]
        Vs = [Float64[] for _ in 1:npsd]
        lastidx = 0
        @myinbounds for term in problem.objective
            ref_psd, ref_idx = sst.monomial_mapping[monomial(term)]
            push!(I[ref_psd], ref_idx)
            # We have to be aware of whether this index occurs on the diagonal or not. We have two different scenarios: In the
            # scalar product ⟨C, X⟩, we could just do dot(C.data, X.data) on the upper triangle and would be fine. However, we
            # also have the operations M .+= C, where we have to be careful. Since the latter occurs far more often, we do
            # scale our matrix appropriately.
            push!(Vs[ref_psd], isdiag_triu(ref_idx) ? coefficient(term) : .5coefficient(term))
            lastidx < ref_psd && (lastidx = ref_psd)
        end
        C = Vector{PackedMatrix{Float64,SparseVector{Float64,Int}}}(undef, lastidx)
        for (i, dim, Iᵢ, Vᵢ) in zip(1:lastidx, sst.psds, I, Vs)
            @myinbounds C[i] = PackedMatrix(dim, sparsevec(Iᵢ, Vᵢ, dim * (dim +1) ÷ 2))
        end
    end
    maxn = maximum(sst.psds)
    return StrideData(
        sst.constraints, A, Astar, AAstarinv, C,
        [PackedMatrix{Float64}(undef, nᵢ) for nᵢ in sst.psds],
        [PackedMatrix{Float64}(undef, nᵢ) for nᵢ in sst.psds],
        [PackedMatrix{Float64}(undef, nᵢ) for nᵢ in sst.psds],
        [PackedMatrix{Float64}(undef, nᵢ) for nᵢ in sst.psds],
        Vector{Float64}(undef, sst.constraints),
        Vector{Float64}(undef, sst.constraints),
        Vector{Float64}(undef, max(maxn*(maxn+1)÷2, sst.constraints))
    ), sst.monomial_mapping
end

function sparse_optimize(::Val{:Stride}, problem::PolyOptProblem{P,M,V}, groupings::Vector{<:Vector{<:AbstractVector{M}}};
    verbose::Bool=false, kwargs...) where {P,M,V}
    # σ: step size for projected gradient descent
    @assert(!problem.complex && isone(problem.prefactor))
    # The STRIDE solver is an SDP solver specifically designed for problems from polynomial optimization. Hence, we need to
    # be able to quickly evaluate the constraints and objectives as well as construct the semidefinite relaxation.
    @verbose_info("Constructing relaxation data")
    setup_relaxation = @elapsed begin
        ssd, monomial_mapping = prepare_relaxation(problem, groupings, verbose)
    end
    @verbose_info("Got relaxation data in ", setup_relaxation, " seconds. Generating initial point")
    # now we got compiled functions A, Astar, AAstarinv that do the desired operations; but they live in a newer world age than
    # our function call. We therefore cannot call them at all! Just use one slow invokelatest to avoid this problem.
    #return Base.invokelatest(stride, ssd, problem, monomial_mapping, groupings; verbose, kwargs...)
    return stride(ssd, problem, monomial_mapping, groupings; verbose, kwargs...)
end

if bfgs_optim
    mutable struct PolyLBFGSState{R,S<:Optim.LBFGSState} <: Optim.AbstractOptimizerState
        const state::S
        const ssd::StrideData{R}
        const tol::R
        success::Bool
    end

    Optim.reset!(method, state::PolyLBFGSState, obj, x) = Optim.reset!(method, state.state, obj, x)
    Optim.update_state!(d, state::PolyLBFGSState, method::Optim.LBFGS) = Optim.update_state!(d, state.state, method)
    Optim.update_h!(d, state::PolyLBFGSState, method::Optim.LBFGS) = Optim.update_h!(d, state.state, method)
    Optim.trace!(tr, d, state::PolyLBFGSState, iteration, method::Optim.LBFGS, options, curr_time=time()) =
        Optim.trace!(tr, d, state.state, iteration, method, options, curr_time)
    function Base.getproperty(state::PolyLBFGSState, prop::Symbol)
        if prop === :state || prop === :ssd || prop === :tol || prop === :success
            return getfield(state, prop)
        else
            return getproperty(getfield(state, :state), prop)
        end
    end
    function Base.setproperty!(state::PolyLBFGSState, prop::Symbol, v)
        if prop === :success
            return setfield!(state, prop, v)
        else
            return setfield!(getfield(state, :state), prop, v)
        end
    end

    function Optim.assess_convergence(state::PolyLBFGSState{R}, d, ::Optim.Options) where {R}
        ssd = state.ssd
        Astar, Z, W, ξ = ssd.Astar, ssd.X, ssd.S, state.state.x
        for (Wᵢ, Zᵢ, Astarᵢ) in zip(W, Z, Astar)
            copyto!(Wᵢ, Zᵢ)
            mul!(Wᵢ, Astarᵢ, ξ, -one(R), -one(R))
            psdproject!(Wᵢ)
        end
        # x_converged, f_converged, g_converged, f_increased - we disregard all this and use the criterion from Algorithm 6
        state.success = converged = ηproj(ssd, W, ξ) ≤ state.tol
        return converged, converged, converged, false
    end
end

function LinearAlgebra.isposdef(M::PackedMatrix, tol)
    M = copy(M)
    return isposdef(cholesky!(M, shift=tol, check=false))
end

function psdproject!(m::PackedMatrix)
    eigs = eigen!(m, 0., Inf)
    fill!(m, zero(eltype(m)))
    for (eval, evec) in zip(eigs.values, eachcol(eigs.vectors))
        spr!(eval, evec, m)
    end
    return m
end

function stride(ssd::StrideData{R}, problem::PolyOptProblem{P,M,V}, monomial_mapping::Dict{M,Tuple{Int,Int}}, groupings;
    opti_round::Union{Function,Nothing}=nothing, opti_local::Union{Function,Nothing}=nothing, verbose::Bool=false,
    tol::R=1e-8, maxiter::Integer=10,
    tol_init::R=1e-4, maxiter_init::Integer=1000, verbose_init::Bool=false,
    tol_sgsapg::R=1e-12, maxiter_sgsapg::Integer=1000, verbose_sgsapg::Bool=false,
    tol_lbfgs::R=1e-12, mem_lbfgs::Integer=10, maxiter_lbfgs::Integer=1000, verbose_lbfgs::Bool=false,
    σ::R=10., options_nlp::Optim.Options=Optim.Options()) where {R<:Real,P,M,V}
    @assert(tol_init > 0 && maxiter_init > 0 && tol > 0 && σ > 0)
    # A: output, [X₁, ...], α, β -> α output + β A(X₁, ...)
    # Astar: output, y, i, α, β -> α output + β Xᵢ(y)
    # AAstarinv: output, y -> ̃y
    d, A, Astar, C, X, S, y, tmpv = ssd.d, ssd.A, ssd.Astar, ssd.C, ssd.X, ssd.S, ssd.y, ssd.tmpv
    tmpm = ssd.tmpm1
    # Input: Generate
    # - an initial point (X₀, y₀, S₀) via spadmm (i.e., Algorithm 4). Not necessarily feasible.
    spadmm!(ssd, σ, (1 + sqrt(5)) / 2, tol_init, maxiter_init, verbose_init)
    @assert(all(M -> isposdef(M, tol), S))
    # - a nondecreasing positive sequence {σₖ}; choose σₖ ≥ O(√k) to get convergence according to Theorem 2.
    #   > choose σₖ = sqrt(k) ... according to the IEEE paper, choose 10 (constant)
    # - a nonnegative sequence {ϵₖ} such that {kϵₖ} is summable,
    #   > looks like it is not required
    # - a stopping criterion η: 𝕊₊ⁿ × ℝᵐ × 𝕊₊ⁿ → ℝ₊ with a tolerance Tol > 0
    # - V = {X₀} if X₀ is feasible, V = ∅ otherwise. We don't store the Xs, but just their objectives
    #   > make it feasible by projecting onto the positive part as a first step.
    @verbose_info("Post-validating the initial point")
    Vopt = zero(R)
    for (Xᵢ, Cᵢ) in zip(X, C)
        Vopt += dot(Cᵢ, psdproject!(Xᵢ))
    end
    #   > now we have positive X, but feasibility also requires the linear constraints to be satisfied. Check it:
    for (i, (Aᵢ, Xᵢ)) in enumerate(zip(A, X))
        mul!(tmpv, Aᵢ, Xᵢ, true, i != 1)
    end
    @myinbounds tmpv[1] -= 1.
    ηp = norm(tmpv) / R(2)
    if ηp > tol
        # no, we are infeasible. So Vopt is not an appropriate value for the last good result.
        @verbose_info("Initial candidate generated by ADMM+ is infeasible with error ", ηp)
        Vopt = R(Inf)
    else
        @verbose_info("Initial candidate generated by ADMM+ is feasible")
    end
    CX = Vopt
    # - a positive integer r ∈ [1, n] - this is the rank
    # - a positive constant ϵ > 0. According to the IEEE paper, choose 1e-12

    if isnothing(opti_local)
        @verbose_info("No local optimization specified; automatically constructing one. ",
           "Compiling polynomials for fast evaluation")
        setup_static = @elapsed let status
            nvars = length(problem.variables)
            #=df = Optim.TwiceDifferentiable(
                static_objective,
                (storage, x) -> StaticPolynomials.gradient!(storage, static_objective, x),
                (storage, x) -> StaticPolynomials.evaluate_and_gradient!(storage, static_objective, x),
                (storage, x) -> StaticPolynomials.hessian!(storage, static_objective, x),
                Vector{R}(undef, nvars)
            )
            dfc = Optim.TwiceDifferentiableConstraints(
                (storage, x) -> StaticPolynomials.evaluate!(storage, static_constraints, x),
                (storage, x) -> StaticPolynomials.jacobian!(storage, static_constraints, x),
                (storage, x, λ) -> StaticPolynomials.hessian!(storage, λ, static_constraints, x),
                fill(R(-Inf), nvars), fill(R(Inf), nvars), fill(zero(R), length(problem.constraints)),
                [c.type == pctNonneg ? R(Inf) : zero(R) for c in problem.constraints]
            )
            ip = Optim.IPNewton()=#
            # First, we need to construct the manifold. Usually, this is just ℝ^n, but in the case of PSD constraints, we
            # additionally create the necessary variables for the matrices which are then in 𝕊^d and use equality constraints
            # to force them to their values.
            equality_constraints = 0
            inequality_constraints = 0
            psd_constraints = 0
            npsd = 0
            for c in problem.constraints
                if c.constraint isa AbstractMatrix
                    s = size(c.constraint, 1)
                    psd_constraints += (s * (s +1)) >> 1
                    npsd += 1
                elseif c.type == pctNonneg
                    inequality_constraints += 1
                else
                    equality_constraints += 1
                end
            end
            psd_variables = variable_union_type(problem)[similar_variable(problem.objective, gensym("psd"))
                                                         for _ in 1:psd_constraints]
            constrs_eq = Vector{polynomial_type(problem)}(undef, equality_constraints + psd_constraints)
            constrs_ineq = Vector{polynomial_type(problem)}(undef, inequality_constraints)
            manis = Manifolds.AbstractManifold{Manifolds.ℝ}[Manifolds.Euclidean(nvars)]
            sizehint!(manis, 1 + npsd)
            i_eq, i_ineq, i_psd = 1, 1, 1
            @myinbounds for c in problem.constraints
                if c.constraint isa AbstractMatrix
                    for mj in 1:size(c.constraint, 1), mi in 1:mj
                        constrs_eq[i] = c.constraint[mi, mj] - psd_variables[i_psd]
                        i_eq += 1
                        i_psd += 1
                    end
                    push!(manis, Manifolds.SymmetricPositiveDefinite(size(c.constraints, 1)))
                elseif c.type == pctNonneg
                    constrs_ineq[i_ineq] = c.constraint
                    i_ineq += 1
                else
                    constrs_eq[i_eq] = c.constraint
                    i_eq += 1
                end
            end
            all_variables = [problem.variables; psd_variables]
            mani = Manifolds.ProductManifold(manis...)
            point = rand(mani)
            point_assembled = Vector{R}(undef, length(all_variables))
            point_dissembled = similar(point_assembled)
            assemble! = p::Manifolds.ProductRepr -> let i = 1
                for part in Manifolds.submanifold_components(p)
                    if part isa AbstractVector
                        copyto!(@view(point_assembled[i:i+length(part)-1]), part)
                        i += length(part)
                    elseif part isa AbstractMatrix # symmetric
                        s = LinearAlgebra.checksquare(part)
                        s = (s * (s +1)) >> 1
                        trttp!('U', part, @view(point_assembled[i:i+s-1]))
                        i += s
                    else
                        @assert(false)
                    end
                end
                return point_assembled
            end
            dissemble! = to::Manifolds.ProductRepr -> let i = 1
                for part in Manifolds.submanifold_components(to)
                    if part isa AbstractVector
                        copyto!(part, @view(point_dissembled[i:i+length(part)-1]))
                        i += length(part)
                    elseif part isa AbstractMatrix # symmetric
                        s = LinearAlgebra.checksquare(part)
                        s = (s * (s +1)) >> 1
                        tpttr!('U', @view(point_dissembled[i:i+s-1]), part)
                        tpttr!('L', @view(point_dissembled[i:i+s-1]), part)
                        i += S
                    else
                        @assert(false)
                    end
                end
                return to
            end
            static_objective = StaticPolynomials.Polynomial(problem.objective, variables=all_variables)
            # We might not have all constraints, and an empty PolynomialSystem is not supported directly.
            static_constrs_eq = StaticPolynomials.Polynomial[StaticPolynomials.Polynomial(c, variables=all_variables)
                                                             for c in constrs_eq]
            static_constrs_ineq = StaticPolynomials.Polynomial[StaticPolynomials.Polynomial(-c, variables=all_variables)
                                                               for c in constrs_ineq]
            local alms, epms, mp
            # We reproduce all the defaults from the method! instantiation, as we want to keep the state alive.
            # There are two possibilities when constraints are present: augmented_Lagrangian_method! and exact_penalty_method!.
            # The latter is much better with infeasible starts (only with the standard smoothing), but less precise in total.
            # For this reason, we run the penalty method before switching to the lagrangian.
            let μ = ones(inequality_constraints), λ = ones(equality_constraints + psd_constraints), ρ = 1.0, u = 1e-1,
                evaluation = Manopt.InplaceEvaluation(), smoothing = Manopt.LogarithmicSumOfExponentials(),
                _objective = Manopt.Manopt.ConstrainedManifoldObjective(
                    (M, p) -> static_objective(assemble!(p)),
                    (M, x, p) -> (StaticPolynomials.gradient!(point_dissembled, static_objective, assemble!(p)); dissemble!(x)),
                    isempty(static_constrs_ineq) ? nothing : Function[(M, p) -> c(assemble!(p)) for c in static_constrs_ineq],
                    isempty(static_constrs_ineq) ? nothing : Function[
                        (M, x, p) -> (StaticPolynomials.gradient!(point_dissembled, c, assemble!(p)); dissemble!(x))
                        for c in static_constrs_ineq
                    ],
                    isempty(static_constrs_eq) ? nothing : Function[(M, p) -> c(assemble!(p)) for c in static_constrs_eq],
                    isempty(static_constrs_eq) ? nothing : Function[
                        (M, x, p) -> (StaticPolynomials.gradient!(point_dissembled, c, assemble!(p)); dissemble!(x))
                        for c in static_constrs_eq
                    ];
                    evaluation
                ),
                sub_state = Manopt.QuasiNewtonState(
                    mani,
                    copy(mani, point);
                    initial_vector = Manifolds.zero_vector(mani, point),
                    direction_update = Manopt.QuasiNewtonLimitedMemoryDirectionUpdate(
                        mani, copy(mani, point), Manopt.InverseBFGS(), 30
                    ),
                    stopping_criterion = Manopt.StopAfterIteration(300) |
                                            Manopt.StopWhenGradientNormLess(1e-3) |
                                            Manopt.StopWhenStepsizeLess(1e-8),
                    stepsize = Manopt.default_stepsize(mani, Manopt.QuasiNewtonState),
                )
                alms = Manopt.AugmentedLagrangianMethodState(
                    mani, _objective, point, Manopt.DefaultManoptProblem(
                        mani, Manopt.ManifoldGradientObjective(
                            Manopt.AugmentedLagrangianCost(_objective, ρ, μ, λ),
                            Manopt.AugmentedLagrangianGrad(_objective, ρ, μ, λ);
                            evaluation
                        )
                    ), sub_state; μ, λ, ρ
                )
                epms = Manopt.ExactPenaltyMethodState(
                    mani, point, Manopt.DefaultManoptProblem(
                        mani, Manopt.ManifoldGradientObjective(
                            Manopt.ExactPenaltyCost(_objective, ρ, u; smoothing),
                            Manopt.ExactPenaltyGrad(_objective, ρ, u; smoothing);
                            evaluation
                        )
                    ), sub_state; u, ρ
                )
                obj = Manopt.decorate_objective!(mani, _objective)
                mp = Manopt.DefaultManoptProblem(mani, obj)
            end
            opti_local = (x0) -> begin
                # Optim uses an interior point algorithm. We might need a phase I algorithm if we are not at the interior. So
                # first check the constraints at the initial point.
                #=boundary = false
                for (c, p) in zip(problem.constraints, static_constraints.polys)
                    val = p(x0)
                    if c.type == pctNonneg
                        if abs(val) < tol
                            boundary = true
                        elseif val < -tol
                            @verbose_info("Initial point is infeasible")
                            return R(Inf), x0
                        end
                    elseif abs(val) > tol
                        @verbose_info("Initial point is infeasible")
                        return R(Inf), x0
                    end
                end
                if boundary
                    @verbose_info("Initial point is on the boundary")
                    return static_objective(x0), x0
                end=#
                # we just reset our objective instead of reconstructing things all over
                # x0 contains the data for problem.variables
                copyto!(Manifolds.submanifold_components(point)[1], x0)
                copyto!(@view(point_assembled[1:nvars]), x0) # copy the points
                fill!(@view(point_assembled[nvars+1:end]), 0.)
                i_man = 2
                i_eq = 1
                for c in problem.constraints
                    if c.constraint isa AbstractMatrix
                        point_man = Manifolds.submanifold_components(point)[i_man]
                        for mj in 1:size(c.constraint, 1), mi in 1:mj
                            point_man[mj, mi] = point_man[mi, mj] = static_constrs_eq.polys[i_eq](point_assembled)
                            i_eq += 1
                        end
                        i_man += 1
                    elseif c.type != pctNonneg
                        i_eq += 1
                    end
                end
                copyto!(point, Manopt.retract(mani, point, point))
                # ^ the direction is irrelevant for our manifolds, so just give some vector of the appropriate size
                # copyto!(point, Manopt.get_solver_return(Manopt.solve!(mp, epms))) # exact penalty method
                # fill!(alms.λ, 1.)
                # fill!(alms.μ, 1.)
                # alms.penalty = Inf
                # minim = assemble!(Manopt.get_solver_return(Manopt.solve!(mp, alms))) # augmented Lagrangian method
                minim = assemble!(Manopt.get_solver_return(Manopt.solve!(mp, epms)))
                obj = static_objective(minim)
                #=Optim.clear!(df)
                result = Optim.optimize(df, dfc, x0, ip, options_nlp)
                violation = zero(R)
                minim = Optim.minimizer(result)=#
                violation = sum(p -> abs(p(minim)), static_constrs_eq, init=zero(R)) +
                    sum(p -> let v = p(minim); v < zero(R) ? -v : zero(R) end, static_constrs_ineq, init=zero(R))
                @verbose_info(@sprintf("NLP: constraint violation: %.8g, cost: %.8g\n", violation, obj))
                return violation < tol ? obj : R(Inf), minim
                #=if verbose
                    if Optim.converged(result)
                        status = "converged"
                    elseif Optim.iteration_limit_reached(result)
                        status = "iteration limit"
                    elseif Optim.f_increased(result)
                        status = "objective increased"
                    elseif isa(result.ls_success, Bool) && !result.ls_success
                        status = "line search failed"
                    elseif Optim.time_run(result) > Optim.time_limit(result)
                        status = "time limit"
                    else
                        status = "unknown failure"
                    end
                    @printf("NLP: %s, constraint violation: %.8g, cost: %.8g\n", status, violation, Optim.minimum(result))
                end
                return violation < tol ? Optim.minimum(result) : R(Inf), Optim.minimizer(result)=#
            end
        end
        @verbose_info("Compilation finished in ", setup_static, " seconds.")
    end

    if bfgs_optim
        local lbfgs_fg!
        let Z = X
            # ϕ(ξ) := ‖Π(A* ξ + Z)‖^2/2 - ⟨b, ξ⟩
            # ∇ϕ(ξ) = AΠ(A* ξ + Z) - b
            lbfgs_fg! = Optim.only_fg!((F, G, ξ) -> begin
                ϕ = zero(R)
                for (i, (tmpmᵢ, Zᵢ, Aᵢ, Astarᵢ)) in enumerate(zip(tmpm, Z, A, Astar))
                    copyto!(tmpmᵢ, Zᵢ)
                    mul!(tmpmᵢ, Astarᵢ, ξ, true, true)
                    psdproject!(tmpmᵢ)
                    isnothing(F) || (ϕ += norm(tmpmᵢ)^2)
                    isnothing(G) || mul!(G, Aᵢ, tmpmᵢ, true, i != 1)
                end
                isnothing(G) || (G[1] -= one(eltype(G)))
                isnothing(F) || return R(1/2)*ϕ - ξ[1]
                nothing
            end)
        end
        lbfgs = Optim.LBFGS(m = mem_lbfgs, linesearch=Optim.LineSearches.MoreThuente())
        lbfgs_options = Optim.Options(g_tol=tol_lbfgs, show_trace=verbose_lbfgs, show_every=50, iterations=maxiter_lbfgs)
        lbfgs_d = Optim.OnceDifferentiable(lbfgs_fg!, y, zero(R))
        lbfgs_state = PolyLBFGSState(Optim.LBFGSState(
            y, # x
            similar(y), # x_previous
            similar(y), # g_previous
            Vector{R}(undef, mem_lbfgs), # rho
            [similar(y) for _ in 1:mem_lbfgs], # dx_history
            [similar(y) for _ in 1:mem_lbfgs], # dg_history
            R(NaN)*y, # dx
            R(NaN)*y, # # dg
            R(NaN)*y, # u
            R(NaN), # f_x_previous
            similar(y), # twoloop_q
            Vector{R}(undef, mem_lbfgs), # twoloop_alpha
            0, # pseudo_iteration
            similar(y), # s
            similar(y), # x_ls
            one(R)
        ), ssd, tol_lbfgs, false)
    else
        bfgsmemα = Vector{Float64}(undef, mem_lbfgs)
        bfgsmemρ = Vector{Float64}(undef, mem_lbfgs)
        bfgsmemS = Matrix{Float64}(undef, d, mem_lbfgs)
        bfgsmemY = Matrix{Float64}(undef, d, mem_lbfgs)
    end
    mons = problem.last_moments
    empty!(mons)
    sizehint!(mons, length(monomial_mapping))
    status = :max_iter
    σₖ = zero(R)
    # Iterate the following steps for k = 1, ...:
    @myinbounds for k in 1:maxiter
        σₖ += σ
        # Step 1 (Projection). Compute (̄Xₖ, yₖ, Sₖ) satisfying (5) and (6), where ̄Xₖ is an inexact projection of Xₖ₋₁ - σₖ C onto
        # the primal feasible set of the SDP, by running algorithm sgsapg and algorithm mlbfgs, warmstarted by (yₖ₋₁, Sₖ₋₁).
        # → Section 4
        @myinbounds for (Xᵢ, Cᵢ) in zip(X, C)
            for (idx, val) in zip(rowvals(Cᵢ.data), nonzeros(Cᵢ.data))
                @myinbounds Xᵢ[idx] -= σₖ * val
            end
        end
        sgsapg!(ssd, tol_sgsapg, maxiter_sgsapg, verbose_sgsapg)
        if bfgs_optim
            verbose_lbfgs && @verbose_info("Entering limited-memory BFGS method")
            # we have to duplicate the behavior of initial_state here
            Optim.reset!(lbfgs, lbfgs_state, lbfgs_d, y)
            copyto!(lbfgs_state.state.x_previous, y) # not strictly necessary (else it would be part of reset!)
            lbfgs_state.success = false
            Optim.optimize(lbfgs_d, y, lbfgs, lbfgs_options, lbfgs_state)
            if !lbfgs_state.success
                for (Wᵢ, Zᵢ, Astarᵢ) in zip(S, X, Astar)
                    copyto!(Wᵢ, Zᵢ)
                    mul!(Wᵢ, Astarᵢ, y, -one(R), -one(R))
                    psdproject!(Wᵢ)
                end
                ηproj(ssd, S, y) # to calculate X (which is stored in tmpm2)
            end
            # else checking the termination criteria already constructed the state
        else
            # the supposed parameters are completely unspecified.
            mlbfgs!(ssd, bfgsmemα, bfgsmemρ, bfgsmemS, bfgsmemY, Inf, tol_lbfgs, 1/10, 1/2, one(R), one(R), maxiter_lbfgs,
                verbose_lbfgs)
        end
        # From the solution ξ of L-BFGS, we output
        # W = Π(-A*ξ - Z) [already done in mlbfgs!]
        # X(W, ξ) = A*ξ + Z + W [the ηproj step at the end of mlbfgs! stores this value in tmpm2]
        # (̄Xₖ, yₖ, Sₖ) = (X(W, ξ), 1/σₖ ξ, 1/σₖ W) [note that ξ and y are aliased, as well as W and S]
        copyto!.(X, ssd.tmpm2)
        σinvₖ = inv(σₖ)
        lmul!(σinvₖ, y)
        lmul!.(σinvₖ, S)

        # Step 2 (Certification). If η(̄Xₖ, yₖ, Sₖ) < Tol, output (̄Xₖ, yₖ, Sₖ) and stop.
        CX, by, ηp, ηd, ηg = η(ssd, Val{true})
        (verbose_lbfgs || verbose_sgsapg || !isnothing(opti_local)) &&
            @verbose_info("Iteration | Primal objective | Dual objective | Primal infeasibility | Dual infeasibility | Relative duality gap\n")
        @verbose_info(@sprintf("%9d | %16g | %14g | %20g | %18g | %20g", k, CX, by, ηp, ηd, ηg))
        if max(ηp, ηg, ηd) ≤ tol
            @verbose_info("Successfully converged.")
            Vopt = CX
            status = :ok
            break
        end

        # Step 3 (Acceleration): Compute ̂Xₖ via rounding, local NLP optimization and lifting, starting from ̄Xₖ:
        if !isnothing(opti_local)
            # ̄Xₖ = ∑ᵢ λᵢ vᵢ vᵢᵀ (truncate after r)
            # ̄xₖᵢ = rounding(vᵢ). Problem-dependent, but can first normalize vᵢ such that the constant component is 1, then
            #                     extract the order-1 entries. If easily possible, project to the feasible set.
            # Problem: we now have multiple X, contrary to the Yang paper. If those just came from constraints, we would simply
            # disregard all but the moment matrix. However, we allow for generic sparsity structure, i.e., there may be
            # multiple sub-moment matrices with possibly overlapping variables. To account for this, we will use our solution
            # extraction heuristic. This is different from the eigenvector decomposition method!
            bound = min(Vopt, CX) # This seems to be a bit unfair. CX may well be generated by something infeasible, but
                                  # depending on the NLP algorithm, these optimization results might always be feasible.
                                  # However, we cannot guarantee this - and how to then decide whether a slightly smaller value
                                  # is better than a slightly more feasible point?
            @verbose_info("Entering local optimization")
            for (mon, (i, idx)) in monomial_mapping
                mons[mon] = X[i][idx]
            end
            solutions = isnothing(opti_round) ? poly_solutions_heuristic(problem; verbose) : opti_round(problem; verbose)
            # ̂xₖᵢ = nlp(̄xₖᵢ)
            # ̂xₖ = argmin p(̂xₖᵢ) over all { ̂xₖᵢ : i = 1, ..., r }
            best = R(Inf)
            local best_solution
            for solution in solutions
                val, pos = opti_local(solution) # note: we REQUIRE solutions returned by opti_local to be feasible (should
                                                # return (Inf, arbitrary) if not)!
                if val < best
                    best_solution = pos
                    best = val
                end
            end
            # ̂Xₖ = proj(monomial lifting of ̂xₖ)
            # Step 4 (Policy). Choose the better candidate in {̄Xₖ, ̂Xₖ} to update Xₖ:
            # Xₖ = ̂Xₖ if ⟨C, ̂Xₖ⟩ < min( ⟨C, ̄Xₖ⟩, min_{X ∈ V} { ⟨C, X⟩ : X ∈ V } ) - ϵ ∧ ̂Xₖ feasible
            if best < bound - 1e-12
                @verbose_info("Found better optimization candidate with objective ", best)
                varmap = problem.var_map
                for (grouping, Xᵢ) in zip(Iterators.flatten(groupings), X)
                    idx = 1
                    for (col, mon_col) in enumerate(grouping)
                        val_col = one(R)
                        for (var, exp) in powers(mon_col)
                            val_col *= best_solution[varmap[var]]^exp
                        end
                        for mon_row in Iterators.take(grouping, col)
                            val = val_col
                            for (var, exp) in powers(mon_row)
                                val *= best_solution[varmap[var]]^exp
                            end
                            Xᵢ[idx] = val
                            idx += 1
                        end
                    end
                end
                # If Xₖ = ̂Xₖ, set V ← V ∪ { ̂Xₖ }
                Vopt = best
            else
                @verbose_info("Optimization did not improve")
            end
        end
        #    = ̄Xₖ otherwise (nop, we already did this)
    end
    @verbose_info("Optimization completed")
    for (mon, (i, idx)) in monomial_mapping
        mons[mon] = X[i][idx]
    end
    λmins = Vector{R}(undef, length(groupings))
    for (i, (tmpmᵢ, Cᵢ, Astarᵢ)) in enumerate(zip(tmpm, C, Astar))
        copyto!(tmpmᵢ, Cᵢ)
        mul!(tmpmᵢ, Astarᵢ, y, -one(R), true)
        @myinbounds λmins[i] = eigvals!(tmpmᵢ, 1:1)[1]
    end
    return status, isinf(Vopt) ? CX : Vopt, x -> let by = y[1], groupings = groupings, problem = problem,
        varmap = problem.var_map, λmins = λmins
        # This is a function that returns the relative suboptimality gap when passed a solution that is extracted in some way
        # from the moment matrix.
        # ηs = |p(̂x) - (⟨b, y⟩ + Mb λₘᵢₙ(C - A*y))|/(1 + |p(̂x)| + |⟨b, y⟩ + Mb λₘᵢₙ(C - A*y)|)
        px = problem.objective(problem.variables => x)
        # ̂x is feasible approximate solution that is rounded from leading eigenvector of X
        # Mb ≥ tr(X)
        Mbλₘᵢₙ = zero(R)
        for (grouping, λmin) in zip(Iterators.flatten(groupings), λmins)
            Mbᵢ = zero(R)
            for mon in grouping
                val_col = one(R)
                for (var, exp) in powers(mon)
                    val_col *= x[varmap[var]]^(2exp)
                end
                Mbᵢ += val_col
            end
            Mbλₘᵢₙ += Mbᵢ * λmin
        end
        return abs(px - (by + Mbλₘᵢₙ)) / (1 + abs(px) + abs(by + Mbλₘᵢₙ))
    end
end

function spadmm!(ssd::StrideData{R}, σ::R, γ::R, tol::R, maxiter::Integer, verbose::Bool) where {R<:Real}
    @assert(0 ≤ γ ≤ 2 && σ > 0 && tol ≥ 0 && maxiter > 0)
    A, Astar, AAstarinv, C, X, S, y = ssd.A, ssd.Astar, ssd.AAstarinv, ssd.C, ssd.X, ssd.S, ssd.y
    tmpm = ssd.tmpm1
    @verbose_info("Entering ADMM+ with penalty σ = ", σ, ", step length γ = ", γ, ", and tolerance ", tol,
        "\nIteration | Primal objective | Dual objective | Primal infeasibility | Dual infeasibility | Relative duality gap")
    lC = length(C)
    # Initialization: Initial points X₀ = S₀ = 0 ∈ 𝕊ⁿ
    fill!.(X, zero(R))
    fill!.(S, zero(R))
    σ⁻¹ = inv(σ)
    local y
    # Iterate the following steps for k = 1, ...
    @myinbounds for k in 1:maxiter
        # Step 1: Compute
        #    ̂yₖ₊₁ = (A A*)⁻¹(b/σ - A(Xₖ/σ + Sₖ - C))
        for (i, (Xᵢ, Sᵢ, Aᵢ)) in enumerate(zip(X, S, A))
            # We don't need Sᵢ any more, so let's overwrite it
            if i ≤ lC
                Sᵢ .+= Xᵢ .* σ⁻¹ .- C[i]
            else
                Sᵢ .+= Xᵢ .* σ⁻¹
            end
            mul!(y, Aᵢ, Sᵢ, -one(R), i != 1)
        end
        y[1] += σ⁻¹
        ldiv!(AAstarinv, y)
        # Step 2: Compute
        #    Sₖ₊₁ = ( Π(Xₖ + σ(A* ̂yₖ₊₁ - C)) - (Xₖ + σ(A* ̂yₖ₊₁ - C)) ) / σ
        #             ^ projection onto PSD, use LAPACK.syevx
        for (i, (Xᵢ, Sᵢ, tmpmᵢ, Astarᵢ)) in enumerate(zip(X, S, tmpm, Astar))
            mul!(tmpmᵢ, Astarᵢ, y, σ, false)
            if i ≤ lC
                tmpmᵢ .+= Xᵢ .- σ .* C[i]
            else
                tmpmᵢ .+= Xᵢ
            end
            copyto!(Sᵢ, tmpmᵢ)
            # We effectively project onto the (absolute of the) negative part, which is done by a positive projection, since
            # the positive part is expected to have low rank
            psdproject!(Sᵢ)
            Sᵢ .= (Sᵢ .- tmpmᵢ) .* σ⁻¹
        end
        # Step 3: Compute
        #    yₖ₊₁ = (A A*)⁻¹(b/σ - A(Xₖ/σ + Sₖ₊₁ - C))
        for (i, (Xᵢ, Sᵢ, tmpmᵢ, Aᵢ)) in enumerate(zip(X, S, tmpm, A))
            # This time, we need the Sᵢ again (as well as the Xs), so we need the temporary storage
            if i ≤ lC
                tmpmᵢ .= Xᵢ .* σ⁻¹ .+ Sᵢ .- C[i]
            else
                tmpmᵢ .= Xᵢ .* σ⁻¹ .+ Sᵢ
            end
            mul!(y, Aᵢ, tmpmᵢ, -one(R), i != 1)
        end
        y[1] += σ⁻¹
        ldiv!(AAstarinv, y)
        # Step 4: Compute
        #    Xₖ₊₁ = Xₖ + γ σ (Sₖ₊₁ + A* yₖ₊₁ - C)
        for (i, (Xᵢ, Sᵢ, Astarᵢ)) in enumerate(zip(X, S, Astar))
            if i ≤ lC
                Xᵢ .+= (γ * σ) .* (Sᵢ .- C[i])
            else
                Xᵢ .+= (γ * σ) .* Sᵢ
            end
            mul!(Xᵢ, Astarᵢ, y, γ * σ, true)
        end
        # Until η(Xₖ₊₁, yₖ₊₁, Sₖ₊₁) ≤ tol (eq. 28)
        CX, by, ηp, ηd, ηg = η(ssd, Val{true})
        success = max(ηp, ηg, ηd) ≤ tol
        if verbose && (success || k == maxiter || isone(k % 20))
            @verbose_info(@sprintf("%9d | %16g | %14g | %20g | %18g | %20g", k, CX, by, ηp, ηd, ηg))
        end
        # Output Xₖ₊₁, yₖ₊₁, Sₖ₊₁
        if success
            @verbose_info("Found an initial point up to the given tolerance")
            return
        end
    end
    @verbose_info("Terminated initial point finding due to maxiter")
end

function sgsapg!(ssd::StrideData{R}, tol::R, maxiter::Integer, verbose::Bool) where {R<:Real}
    @assert(tol ≥ 0) # also all(W .⪰ 0), but this is too expensive to check
    A, Astar, AAstarinv, Z, W, ξ = ssd.A, ssd.Astar, ssd.AAstarinv, ssd.X, ssd.S, ssd.y
    tmpm1, tmpm2, tmpv = ssd.tmpm1, ssd.tmpm2, ssd.tmpv
    @verbose_info("Entering sGS-based accelerated proximal gradient method\nIteration |     residue")
    # Initialization: Set ̃W₁ = W₀ and t₁ = 1
    copyto!.(tmpm1, W)
    t = 1.
    # Iterate the following steps for k = 1, ...
    k = 1
    @myinbounds while true
        # Step 1 (sGS update): Compute
        #    ̃ξₖ = (A A*)⁻¹(b - A(Z) - A(̃Wₖ))
        for (i, (Aᵢ, Zᵢ, tmpm1ᵢ)) in enumerate(zip(A, Z, tmpm1))
            tmpm1ᵢ .+= Zᵢ
            mul!(tmpv, Aᵢ, tmpm1ᵢ, -one(R), i != 1)
        end
        tmpv[1] += one(R)
        copyto!.(tmpm1, W) # we need to back up the previous W
        ldiv!(ξ, AAstarinv, tmpv)
        #    Wₖ = Π(-A*̃ξₖ - Z) = Π(A*ξₖ + Z) - (A*ξₖ + Z)
        for (i, (Zᵢ, Wᵢ, tmpm2ᵢ, Astarᵢ)) in enumerate(zip(Z, W, tmpm2, Astar))
            copyto!(Wᵢ, Zᵢ)
            mul!(Wᵢ, Astarᵢ, ξ, true, true)
            copyto!(tmpm2ᵢ, Wᵢ)
            psdproject!(Wᵢ)
            Wᵢ .-= tmpm2ᵢ
        end
        #    ξₖ = (A A*)⁻¹(b - A(Z) - A(Wₖ))
        for (i, (Aᵢ, Zᵢ, Wᵢ)) in enumerate(zip(A, Z, W))
            mul!(tmpv, Aᵢ, Zᵢ, -one(R), i != 1)
            mul!(tmpv, Aᵢ, Wᵢ, -one(R), true)
        end
        tmpv[1] += one(R)
        ldiv!(ξ, AAstarinv, tmpv)

        # Until ηproj(Wₖ, ξₖ) ≤ tol
        ηp = ηproj(ssd, W, ξ) # ηproj will use tmpm2, tmpvlong, and tmpv as temporaries
        success = ηp ≤ tol
        if verbose && (success || k == maxiter || isone(k % 20))
            @verbose_info(@sprintf("%9d | %11g", k, ηp))
        end
        if success
            return # Output Wₖ, ξₖ
        elseif k == maxiter
            @verbose_info("Maximum iteration count reached")
            return
        end

        # Step 2: Compute
        #    tₖ₊₁ = (1 + sqrt(1 + 4tₖ^2))/2
        newt = (one(R) + sqrt(one(R) + R(4)*t^2)) / R(2)
        #    ̃Wₖ₊₁ = Wₖ + (tₖ-1)/tₖ₊₁ (Wₖ - Wₖ₋₁)
        frac = (t - one(R)) / newt
        t = newt
        # now Wₖ is W, Wₖ₋₁ is tmpm1, and we need to fill ̃Wₖ₊₁, which will be tmpm1
        for (Wᵢ, tmpm1ᵢ) in zip(W, tmpm1)
            tmpm1ᵢ .= (one(R) + frac) .* Wᵢ .- frac .* tmpm1ᵢ
        end

        k += 1
    end
end

function mlbfgs!(ssd::StrideData{R}, α::Vector{R}, ρs::Vector{R}, S::Matrix{R}, Y::Matrix{R}, K::R, tol::R, μ::R, ρ::R, τ₁::R,
    τ₂::R, maxiter::Integer, verbose::Bool) where {R<:Real}
    @assert(K > 0 && tol > 0 && 0 < μ < 1/2 && μ < 1 && 0 < ρ < 1 && τ₁ > 0 && τ₂ > 0)
    @verbose_info("Entering modified limited-memory BFGS method\nIteration |     Residue")
    mem = length(α)
    A, Astar, Z, W, ξ = ssd.A, ssd.Astar, ssd.X, ssd.S, ssd.y
    @myinbounds tmpm, ∇ϕ, tmpv = ssd.tmpm1, ssd.tmpv, view(ssd.tmpvlong, 1:ssd.d)
    # Iterate the following steps for k = 1, ...
    k = 1
    for (tmpmᵢ, Zᵢ, Astarᵢ) in zip(tmpm, Z, Astar)
        copyto!(tmpmᵢ, Zᵢ)
        mul!(tmpmᵢ, Astarᵢ, ξ, true, true)
        psdproject!(tmpmᵢ)
    end
    local idx
    @myinbounds while true
        # ϕ(ξ) := ‖Π(A* ξ + Z)‖^2/2 - ⟨b, ξ⟩
        # ∇ϕ(ξ) = AΠ(A* ξ + Z) - b
        # tmpm will already contain the up-to-date Π(A* ξ + Z)
        ϕ = sum(normsq, tmpm, init=zero(R)) / R(2) - ξ[1]
        for (i, (Aᵢ, tmpmᵢ)) in enumerate(zip(A, tmpm))
            mul!(∇ϕ, Aᵢ, vec(tmpmᵢ), true, i != 1)
        end
        ∇ϕ[1] -= one(R)
        @views if k != 1
            # leftover from the last iteration (Y[:, idx] already contains last ∇ϕ)
            # Step 3 (Update memory):
            #     Compute and save uₖ = ξₖ₊₁ - ξₖ and wₖ = ∇ϕ(ξₖ₊₁) - ∇ϕ(ξₖ).
            Y[:, idx] .= ∇ϕ .- Y[:, idx]
            ρs[idx] = one(R) / dot(S[:, idx], Y[:, idx]) # required for the two-loop recursion
        end
        # Step 1 (Search direction):
        #     Choose Qₖ₀ ≻ 0 (we always choose Iₘ - actually, due to Nocedal, Wright, we choose a proportionality)
        #            βₖ := τ₁ ‖∇ϕ(ξₖ)‖^(τ₂)
        β = τ₁*norm(∇ϕ)^τ₂
        #     and compute dₖ = -βₖ ∇ϕ(ξₖ) - gₖ
        #     where gₖ := Qₖ ∇ϕ(ξₖ) with Qₖ ⪰ 0 is obtained via the two-loop recursion as in [60, Algorithm 7.4]
        copyto!(tmpv, ∇ϕ)
        @views for i in k-1:-1:max(1, k-mem)
            idx = mod1(i, mem)
            α[idx] = ρs[idx] * dot(S[:, idx], tmpv)
            tmpv .-= α[idx] .* Y[:, idx]
        end
        idx = mod1(k -1, mem)
        @views k > 1 && lmul!(dot(S[:, idx], Y[:, idx]) / norm(Y[:, idx])^2, tmpv) # this is the proportionality factor
        @views for i = max(1, k-mem):k-1
            idx = mod1(i, mem)
            tmpv .+= S[:, idx] .* (α[idx] - ρs[idx] * dot(Y[:, idx], tmpv))
        end
        tmpv .+= β .* ∇ϕ # tmpv now corresponds to -dₖ
        #     If ‖dₖ‖ ≥ K, then choose dₖ = -βₖ ∇ϕ(ξₖ) (i.e., set Qₖ = 0)
        if norm(tmpv) ≥ K
            tmpv .= β .* ∇ϕ
        end
        # Step 2 (Line search):
        #     Set αₖ = ρ^(mₖ) where mₖ is the smallest nonnegative integer m such that
        #     ϕ(ξₖ + ρ^k dₖ) ≤ ϕ(ξₖ) + μ ρ^m ⟨∇ϕ(ξₖ), dₖ⟩.
        ϕmiddle = 0.
        mρᵏ = -ρ^k
        for (tmpmᵢ, Zᵢ, Astarᵢ) in zip(tmpm, Z, Astar)
            copyto!(tmpmᵢ, Zᵢ)
            mul!(tmpmᵢ, Astarᵢ, ξ, true, true)
            mul!(tmpmᵢ, Astarᵢ, tmpv, mρᵏ, true)
            psdproject!(tmpmᵢ)
            ϕmiddle += norm(tmpmᵢ)^2
        end
        ϕmiddle = ϕmiddle/R(2) - (ξ[1] + mρᵏ * tmpv[1])
        m = max(1, ceil(Int, log(ρ, max(ρ, (ϕ - ϕmiddle) / (μ * dot(∇ϕ, tmpv))))))
        idx = mod1(k, mem)
        α[idx] = ρ^m
        # Step 3 (Update memory):
        #     Compute and save uₖ = ξₖ₊₁ - ξₖ = αₖ dₖ and wₖ = ∇ϕ(ξₖ₊₁) - ∇ϕ(ξₖ).
        S[:, idx] .= .-α[idx] .* tmpv
        copyto!(@view(Y[:, idx]), ∇ϕ) # just store, the difference is calculated in the next iteration
        # [Step2]    Compute ξₖ₊₁ = ξₖ + αₖ dₖ
        ξ .+= @view(S[:, idx])
        #     If k > mem, discard the vector {uₖ₋ₘₑₘ, wₖ₋ₘₑₘ} from storage.
        # this is a no-op, as we cycle through our storage
        # Until: ηproj(Wₖ₊₁, ξₖ₊₁) ≤ tol with Wₖ₊₁ = Π(-A* ξₖ₊₁ - Z)  (but the positive part is expected be of low rank)
        for (i, (Wᵢ, tmpmᵢ, Zᵢ, Astarᵢ)) in enumerate(zip(W, tmpm, Z, Astar))
            copyto!(Wᵢ, Zᵢ)
            mul!(Wᵢ, Astarᵢ, ξ, true, true)
            copyto!(tmpmᵢ, Wᵢ)
            psdproject!(tmpmᵢ)
            Wᵢ .= tmpmᵢ .- Wᵢ
        end
        ηp = ηproj(ssd, W, ξ, false) # ηproj will use tmpm2, tmpvlong, and tmpv as temporaries
        success = ηp ≤ tol
        (success || k == maxiter || isone(k % 50)) && @verbose_info(@sprintf("%9d | %11g", k, ηp))
        success && return
        if !iszero(maxiter) && k == maxiter
            @verbose_info("Maximum iteration count reached")
            return
        end
        # Output Wₖ₊₁, ξₖ₊₁
        k += 1
    end
end

function η(ssd::StrideData{R}, details::Union{Type{Val{false}},Type{Val{true}}}=Val{false}) where {R}
    A, Astar, C, X, S, y, tmpm, tmpv = ssd.A, ssd.Astar, ssd.C, ssd.X, ssd.S, ssd.y, ssd.tmpm1, ssd.tmpv
    # here, note that b = (1, 0, ..., 0)
    lC = length(C)
    # ηp(X) = ‖A(X) - b‖/(1 + ‖b‖)
    for (i, (Aᵢ, Xᵢ)) in enumerate(zip(A, X))
        mul!(tmpv, Aᵢ, Xᵢ, true, i != 1)
    end
    @myinbounds tmpv[1] -= 1.
    ηp = norm(tmpv) / R(2)
    # ηd(y, S) = ‖A*(y) + S - C‖/(1 + ‖C‖)
    # ηg(X, y) = |⟨C, X⟩ - ⟨b, y⟩|/(1 + |⟨C, X⟩| + |⟨b, y⟩|)
    ηdnum, ηdden, CX, by = zero(R), zero(R), zero(R), y[1]
    @myinbounds for (i, (Sᵢ, Xᵢ, tmpmᵢ, Astarᵢ)) in enumerate(zip(S, X, tmpm, Astar))
        mul!(tmpmᵢ, Astarᵢ, y, true, false)
        if i ≤ lC
            Cᵢ = C[i]
            tmpmᵢ .+= Sᵢ .- Cᵢ
            CX += dot(Cᵢ, Xᵢ)
        else
            tmpmᵢ .+= Sᵢ
        end
        ηdnum += LinearAlgebra.norm2(tmpmᵢ)^2
    end
    ηd = sqrt(ηdnum) / (1 + sqrt(ηdden))
    ηg = abs(CX - by) / (1 + abs(CX) + abs(by))
    return details == Val{true} ? (CX, by, ηp, ηd, ηg) : max(ηp, ηd, ηg)
end

function ηproj(ssd::StrideData{R}, W, ξ, recompute_AstarXiplusZ::Bool=true) where {R}
    A, Astar, Z, X = ssd.A, ssd.Astar, ssd.X, ssd.tmpm2
    tmpmvec, tmpv = ssd.tmpvlong, ssd.tmpv
    # X(W, ξ) := A*ξ + W + Z
    # ηproj := max( ‖A(X(W, ξ)) - b‖, ‖X(W, ξ) - Π(A*ξ + Z)‖ )
    norm2 = zero(R)
    if recompute_AstarXiplusZ
        @myinbounds for (Xᵢ, Wᵢ, Zᵢ, Astarᵢ) in zip(X, W, Z, Astar)
            tmpmᵢ = let s = size(Xᵢ, 1)
                PackedMatrix(s, view(tmpmvec, 1:s*(s+1)÷2))
            end
            copyto!(tmpmᵢ, Zᵢ)
            mul!(tmpmᵢ, Astarᵢ, ξ, true, true)
            Xᵢ .= tmpmᵢ .+ Wᵢ
            psdproject!(tmpmᵢ)
            j = 1
            @simd for q in 1:size(tmpmᵢ, 2)
                for _ in 1:q-1
                    norm2 += 2(Xᵢ[j] - tmpmᵢ[j])^2
                    j += 1
                end
                norm2 += (Xᵢ[j] - tmpmᵢ[j])^2
                j += 1
            end
        end
    else
        # if this is called from the BFGS function, the algorithm guarantees that norm2 is indeed exactly zero. But we still
        # need X.
        @myinbounds for (Xᵢ, Wᵢ, Zᵢ, Astarᵢ) in zip(X, W, Z, Astar)
            copyto!(Xᵢ, Zᵢ)
            mul!(Xᵢ, Astarᵢ, ξ, true, true)
            Xᵢ .+= Wᵢ
        end
    end
    for (i, (Aᵢ, Xᵢ)) in enumerate(zip(A, X))
        mul!(tmpv, Aᵢ, Xᵢ, true, i != 1)
    end
    tmpv[1] -= one(R)
    norm1 = norm(tmpv)
    return max(norm1, sqrt(norm2))
end