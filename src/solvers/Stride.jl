struct StrideState{V,R}
    psds::FastVec{Int}
    Astar::FastVec{Tuple{Vector{Int},Vector{Int},Vector{R}}}
    Astar_eq::Tuple{FastVec{Int},FastVec{Int},FastVec{R}}
    b::Tuple{Vector{Int},Vector{R}}
    Ccoeffs::FastVec{R}
    c_eq::FastVec{R}
    variables::Vector{V}
end

struct StrideData{R<:Real,F<:Factorization,M<:AbstractMonomial}
    d::Int
    A::Vector{SparseMatrixCSC{R,Int}}
    Astar::Vector{SparseMatrixCSC{R,Int}}
    AAstarinv::F
    C::Vector{PackedMatrix{R,SparseVector{R,Int}}}
    equality_subspace::Tuple{SparseMatrixCSC{R,Int},SparseVector{R,Int}}
    b::SparseVector{R,Int}
    offset::R
    X::Vector{PackedMatrix{R,Vector{R}}} # each dimensions nᵢ × nᵢ
    S::Vector{PackedMatrix{R,Vector{R}}}
    tmpm1::Vector{PackedMatrix{R,Vector{R}}}
    tmpm2::Vector{PackedMatrix{R,Vector{R}}}
    y::Vector{R} # length d
    tmpv::Vector{R}
    tmpvlong::Vector{R} # length packedsize(max(nᵢ))
    monomial_map::Vector{M}
end

# calculates the index of a given monomial in a deglex ordering (lex according to vars)
function monomial_index(mon, vars)
    nvars = length(vars)
    mostvars = @view(vars[end:-1:2])
    mondeg = degree(mon)
    mindex = binomial(mondeg + nvars -1, nvars) # how many monomials were there with a lower total degree?
    for var in mostvars
        vardeg = degree(mon, var)
        #= for i in vardeg+1:mondeg # check for all possible higher degrees that the current variable may have had
            mindex += binomial((n -1) + (mondeg - i) -1, (n -1) -1) # and add the number of possible monomials
                                                                    # that remain in this subspace
        end =#
        # The above loop is exactly equivalent to the following
        mondeg -= vardeg
        nvars -= 1
        @assert(nvars ≥ 1)
        mindex += ((nvars + mondeg -1) * binomial(nvars + mondeg -2, nvars -1)) ÷ nvars
    end
    return mindex # no +1: we never have the constant monomial
end

function moment_matrix!(sst::StrideState{V,R}, gröbner_basis, grouping, constraint) where {V,R}
    lg = length(grouping)
    unsafe_push!(sst.psds, lg)
    # we construct the dual map A*, which in Lasserre's notation corresponds to the B matrices (though we work with the vec
    # representation)
    # let us first guess how many items we need in our Astar. This is exact in the case of no Gröbner basis; but if we do have
    # one, we cannot efficiently make this guess - we'd have to carry out the products and remainders twice.
    items = let len = lg * nterms(constraint)
        (len * (len +1)) >> 1
    end
    rows = FastVec{Int}(buffer=items)
    cols = FastVec{Int}(buffer=items)
    vals = FastVec{R}(buffer=items)
    cval = zero(R)
    cur_index = 1
    vars = sst.variables
    # now we populate Astar
    for exp2 in 1:lg
        for exp1 in 1:exp2
            for term in rem(grouping[exp1] * grouping[exp2] * constraint, gröbner_basis)
                # every y corresponds to one monomial. We first assume that every monomial occurs (dense basis). This allows us
                # to directly calculate the index of a monomial without having to store them in a dictionary whose size we
                # cannot predict (unless we take the upper bound, of course, which defeats the purpose of sparsity). We count
                # all the monomials that were actually used and in the end remove empty ranges
                # Note that it is still not too easy to calculate the index in a basis with degree bound (instead of maxdegree
                # bound - but this would lead to much bigger numbers, potentially overflowing, so we spend the extra effort).
                # Further note that we don't have to care about the actual monomial order. We just impose deglex here and are
                # self-consistent; the indices are arbitrary anyways.
                if isconstant(term)
                    cval += coefficient(term)
                elseif gröbner_basis isa EmptyGröbnerBasis
                    unsafe_push!(rows, cur_index)
                    unsafe_push!(cols, monomial_index(monomial(term), vars))
                    unsafe_push!(vals, coefficient(term))
                else
                    # same, but safeguarded: the vector might need to grow
                    push!(rows, cur_index)
                    push!(cols, monomial_index(monomial(term), vars))
                    push!(vals, coefficient(term))
                end
            end
            cur_index += 1
        end
    end
    unsafe_push!(sst.Astar, finish!.((rows, cols, vals)))
    unsafe_push!(sst.Ccoeffs, cval)
    return
end

function moment_matrix!(sst::StrideState{V,R}, gröbner_basis, grouping, constraint::AbstractMatrix) where {V,R}
    block_size = LinearAlgebra.checksquare(constraint)
    if block_size == 1
        @inbounds return moment_matrix!(sst, gröbner_basis, grouping, constraint[1, 1])
    end
    lg = length(grouping)
    # to precompute the size, we need to take into account the number of terms in the constraint matrix. However, for the
    # "outer" diagonal blocks, we are interested only in on triangle, where for the "outer" off-diagonal blocks, we need the
    # whole matrix.
    mat_size_full, mat_size_tri = 0
    for j in 1:block_size
        δ = 0
        for i in 1:j-1
            @inbounds δ += nterms(constraints[i, j])
        end
        @inbounds δ2 = nterms(constraints[j, j])
        mat_size_full += 2δ + δ2
        mat_size_tri += δ + δ2
    end
    items = ((lg * (lg -1)) >> 1) * mat_size_full + lg * mat_size_tri
    unsafe_push!(sst.psds, lg)
    rows = FastVec{Int}(buffer=items)
    cols = FastVec{Int}(buffer=items)
    vals = FastVec{R}(buffer=items)
    cval = zero(R)
    cur_index = 1
    for exp2 in 1:lg, block_j in 1:block_size, exp1 in 1:exp2, block_i in 1:(exp1 == exp2 ? block_j : block_size)
        for term in rem(grouping[exp1] * grouping[exp2] * constraint[block_i, block_j], gröbner_basis)
            if isconstant(term)
                cval += coefficient(term)
            elseif gröbner_basis isa EmptyGröbnerBasis
                unsafe_push!(rows, cur_index)
                unsafe_push!(cols, monomial_index(monomial(term), vars))
                unsafe_push!(vals, coefficient(term))
            else
                push!(rows, cur_index)
                push!(cols, monomial_index(monomial(term), vars))
                push!(vals, coefficient(term))
            end
        end
        cur_index += 1
    end
    unsafe_push!(sst.Astar, finish!.((rows, cols, vals)))
    unsafe_push!(sst.Ccoeffs, cval)
    return
end

function moment_matrix_eq!(sst::StrideState{V,R}, gröbner_basis, grouping, constraint) where {V,R}
    lg = length(grouping)
    items = let len = lg * nterms(constraint)
        (len * (len +1)) >> 1
    end
    c = sst.c_eq
    cur_index = length(c) +1
    rows, cols, vals = sst.Astar_eq
    prepare_push!(rows, items)
    prepare_push!(cols, items)
    prepare_push!(vals, items)
    prepare_push!(c, (lg * (lg +1)) >> 1)
    vars = sst.variables
    for exp2 in 1:lg
        for exp1 in 1:exp2
            cval = zero(R)
            for term in rem(grouping[exp1] * grouping[exp2] * constraint, gröbner_basis)
                if isconstant(term)
                    cval += coefficient(term)
                elseif gröbner_basis isa EmptyGröbnerBasis
                    unsafe_push!(rows, cur_index)
                    unsafe_push!(cols, monomial_index(monomial(term), vars))
                    unsafe_push!(vals, coefficient(term))
                else
                    push!(rows, cur_index)
                    push!(cols, monomial_index(monomial(term), vars))
                    push!(vals, coefficient(term))
                end
            end
            unsafe_push!(c, cval)
            cur_index += 1
        end
    end
    return
end

# All our equality constraints basically just introduce additional linear equality constraints into the problem. We can instead
# solve this underdetermined system of linear equations, thereby reducing the number of variables that enter the solver
# firsthand. However, we have to make sure that we preserve sparsity of the matrices for fast operations.
# The following algorithm follows a suggestion by Federico Polloni (https://mathoverflow.net/a/253997).
function reduce_variables(A::SparseMatrixCSC, b, check::Bool)
    # We need to restrict this function to sparse matrices as we need to be able to access the permutation done by the QR
    # procedure; however, this is different for different matrix types.
    @inbounds begin
        n = size(A, 2)
        qrA = qr(A) # A[qrA.prow,qrA.pcol] = qrA.Q * qrA.R
        r = rank(qrA)
        R1 = UpperTriangular(@view(qrA.R[1:r, 1:r]))
        b = @view(b[qrA.prow])
        if check
            QTb = transpose(Matrix(qrA.Q)) * b # without the explicit conversion to a matrix (which has memory drawbacks), this
                                               # can easily take forever, even for small problems.
            any(x -> abs(x) > 100eps(eltype(QTb)), @view(QTb[r+1:end])) && throw(LinearAlgebra.SingularException(-1))
            shifttmp = R1 \ @view(QTb[1:r])
        else
            shifttmp = R1 \ (transpose(@view(qrA.Q[:, 1:r])) * b)
        end
        R2 = @view(qrA.R[1:r, r+1:end])
        ns = [-R1 \ R2; I(n - r)]
        permute!(ns, invperm(qrA.pcol), 1:size(ns, 2))
        # Do we return a sparse or dense vector? If we make this dependent on the number of nonzeros, the function becomes
        # type-unstable. Better use the more likely case of a sparse vector.
        #=shift = zeros(eltype(shifttmp), n)
        for (i, v) in zip(qrA.pcol, shifttmp)
            shift[i] = v
        end=#
        pcol = resize!(qrA.pcol, r)
        sort_along!(pcol, shifttmp)
        shift = SparseVector(n, pcol, shifttmp)
        droptol!(shift, 100eps(eltype(shift)))
        return ns, shift
    end
end

function monomial_from_powers(::Type{<:AbstractMonomial}, vars, powers)
    @assert(!isempty(vars))
    @assert(length(vars) == length(powers))
    result = vars[1] ^ powers[1]
    for (v, p) in zip(Iterators.drop(vars, 1), Iterators.drop(powers, 1))
        map_exponents!(+, result, v^p)
    end
    return result
end

monomial_from_powers(m::Type{<:DynamicPolynomials.Monomial}, vars, powers) = m(vars, copy(powers))

function prepare_relaxation(problem::PolyOptProblem{P,M,V}, groupings::Vector{<:Vector{<:AbstractVector{M}}},
    verbose::Bool, check_eq_consistency::Bool, warn_check_eq_consistency::Bool) where {P,M,V}
    # How many constraints will we need? Every monomial that arises somewhere in the relaxation must be assigned the same
    # value. This means that the first occurrence is free, and all others give one constraint each.
    @verbose_info("Initializing problem relaxation data")
    npsd = length(groupings[1]) + sum(begin
        if constr.type == pctNonneg || constr.type == pctPSD
            length(grs)
        elseif constr.type == pctEqualityNonneg
            2length(grs)
        else
            0
        end
    end for (grs, constr) in zip(Iterators.drop(groupings, 1), problem.constraints); init=0)
    sst = StrideState{V,Float64}(
        FastVec{Int}(buffer=npsd),
        FastVec{Tuple{Vector{Int},Vector{Int},Vector{Float64}}}(buffer=npsd),
        (FastVec{Int}(), FastVec{Int}(), FastVec{Float64}()),
        (monomial_index.(monomials(problem.objective), (problem.variables,)), collect(coefficients(problem.objective))),
        FastVec{Float64}(buffer=npsd),
        FastVec{Float64}(),
        problem.variables
    )
    sort_along!(sst.b...)
    @inbounds if iszero(sst.b[1][1])
        popfirst!(sst.b[1])
        offset = popfirst!(sst.b[2])
    else
        offset = 0.
    end
    # note that STRIDE's dual formulation is
    # max {⟨b, y⟩ : A*y + S = C, S ⪰ 0}
    # whereas Lasserre's is
    # min {∑ᵢ pᵢ yᵢ : M(gⱼ y) ⪰ 0}
    # So while we have max vs. min, we also need to adjust A*y + S = C ⇔ C - A*y ⪰ 0 ⇔ C + A* (-y) ⪰ 0.
    # Therefore, we take A* to give the _actual_ moment matrices, but interpret y as the negative moments. Consequently, what
    # we solve is -min {⟨b, -y⟩ : C + A*(-y) ⪰ 0} = -min {⟨b, ̃y⟩ : C + A*̃y ⪰ 0}
    # moment matrix
    @verbose_info("Assembling moment matrix data")
    for grouping in groupings[1]
        @assert(issorted(grouping, by=degree))
        moment_matrix!(sst, problem.gröbner_basis, grouping, polynomial(constant_monomial(problem.objective)))
    end
    for (constr_groupings, constr) in zip(Iterators.drop(groupings, 1), problem.constraints)
        if constr.type == pctNonneg || constr.type == pctPSD
            for grouping in constr_groupings
                @assert(issorted(grouping, by=degree))
                moment_matrix!(sst, problem.gröbner_basis, grouping, constr.constraint)
            end
        elseif constr.type == pctEqualityNonneg
            for grouping in constr_groupings
                @assert(issorted(grouping, by=degree))
                moment_matrix!(sst, problem.gröbner_basis, grouping, constr.constraint)
                moment_matrix!(sst, problem.gröbner_basis, grouping, -constr.constraint)
            end
        elseif constr.type == pctEqualityGröbner
            for grouping in constr_groupings
                @assert(issorted(grouping, by=degree))
                moment_matrix_eq!(sst, EmptyGröbnerBasis{P}(), grouping, constr.constraint)
            end
        elseif constr.type == pctEqualitySimple
            for grouping in constr_groupings
                @assert(issorted(grouping, by=degree))
                moment_matrix_eq!(sst, problem.gröbner_basis, grouping, constr.constraint)
            end
        else
            @assert(false)
        end
    end
    @assert(length(sst.psds) == npsd)
    # We also need the (dual) objective vector b

    @verbose_info("Compressing moment matrix data")
    # now we have the Astar data - but we always assumed a dense basis when we assigned the y indices. Now, we want to remove
    # all unused regions from the columns. For this, we first start by sorting the data columnwise
    for (rows, cols, vals) in sst.Astar
        sort_along!(cols, rows, vals)
    end
    Astar_eq_coo = finish!.(sst.Astar_eq)
    sort_along!(Astar_eq_coo[2], Astar_eq_coo[1], Astar_eq_coo[3])
    # now we can simultaneously loop over all the matrices at the same time. Note that we must also take the objective vector
    # b into account, as it also references columns! And our equality constraints are stored separately.
    # During this loop, we remove all connection to the actual monomial that was represented. Therefore, we must keep track of
    # this re-indexing here. Unfortunately, we don't know now many monomials will actually be present before running through
    # the loop. So either we allocate the maximum - binom(nvars + 2problem.degree, nvars) - which might be a killer for large
    # sparse problems or we have to live with growing the buffer.
    vars = problem.variables
    nvars = length(vars)
    monomial_map = FastVec{M}(buffer=min(nvars + 2problem.degree, nvars,
        maximum(x -> trunc(Int, sum(length, x, init=0)^1.35), groupings))) # just a crude initial guess
    current_index = 0
    current_monomial = zeros(Int, nvars)
    current_monomial_degree = 0
    current_monomial_left = 1
    is = ones(Int, npsd +2)
    δ = 0
    remaining = length(is)
    if isempty(Astar_eq_coo[2])
        @inbounds is[2] = 0
        remaining -= 1
    end
    while true
        @assert(remaining > 0)
        new_index = typemax(Int)
        @inbounds for k in 1:npsd +2
            i = is[k]
            iszero(i) && continue # is this Astar already finished?
            if k > 2 # most likely first
                cols = sst.Astar[k-2][2]
            elseif isone(k)
                cols = sst.b[1]
            else
                cols = Astar_eq_coo[2]
            end
            maxi = length(cols)
            # a particular column may come multiple times in a row, we need to get the next entry
            while cols[i] - δ == current_index
                cols[i] -= δ # fix the δ to the entry
                i += 1 # and go on
                if i > maxi
                    i = 0
                    remaining -= 1
                    break
                end
            end
            is[k] = i
            iszero(i) && continue # maybe we finished in the previous loop
            this_index = cols[i] - δ
            @assert(this_index > current_index)
            if new_index > this_index
                new_index = this_index
            end
        end
        if new_index == typemax(Int)
            @assert(iszero(remaining))
            break
        else
            @assert(!iszero(remaining))
        end
        # now we know what the next index is - new_index. The last one was current_index. So unless
        # new_index = current_index +1, we can now actually start to subtract. But all our subtractions are virtual...
        current_index += 1
        thisδ = new_index - current_index
        δ += thisδ
        # go to the next monomial, skipping all those that were left out (ordering must be compatible with monomial_index)
        # TODO: is there a more direct way instead of iterating over them?
        for _ in 0:thisδ
            if current_monomial_left > 1
                @inbounds current_monomial[current_monomial_left] -= 1
                current_monomial_left -= 1
                @inbounds current_monomial[current_monomial_left] += 1
            else
                @inbounds current_monomial[1] = 0
                current_monomial_left = findfirst(x -> x > 0, current_monomial)
                if isnothing(current_monomial_left)
                    current_monomial_degree += 1
                    @inbounds current_monomial[end] = current_monomial_degree
                    current_monomial_left = nvars
                else
                    @inbounds current_monomial[current_monomial_left] -= 1
                    current_monomial_left -= 1
                    @inbounds current_monomial[current_monomial_left] = current_monomial_degree -
                        sum(@view(current_monomial[current_monomial_left+1:end]), init=0)
                end
            end
        end
        push!(monomial_map, monomial_from_powers(M, vars, current_monomial))
    end
    # now, current_index = max(max.(<all columns>)), so we now also know how large our matrices have to be (note that this
    # is the number of dual variables y!). The columns are, seen over all matrices, contiguous. We can convert all the Astar
    # vectors into proper sparse matrices.
    constraints = current_index

    @verbose_info("Eliminating equality constraints")
    # Instead of adding the equality constraints, we eliminate them by working in the smaller subspace that is dictated by
    # them. For this, we first need to convert the equality constraint data to a CSC matrix.
    local subspace_map, negsubspace_affine
    tmp = Vector{Int}(undef, constraints)
    let n = constraints
        # Finally, the equality constraints. These work on vectors, not matrices.
        coolen = length(Astar_eq_coo[1])
        if iszero(coolen)
            subspace_map = SparseMatrixCSC{Float64,Int}(I(n))
            negsubspace_affine = SparseVector{Float64,Int}(n, Int[], Float64[])
        else
            !check_eq_consistency && warn_check_eq_consistency &&
                @warn("Checking consistency of equality constraints is disabled.")
            m = length(sst.c_eq)
            Astar_eq_csc = SparseArrays.sparse!(Astar_eq_coo..., m, n, +, tmp, Vector{Int}(undef, m +1),
                Vector{Int}(undef, coolen), Vector{Float64}(undef, coolen), Astar_eq_coo...)
            # Now we need to reduce this linear system into a map between subspaces
            subspace_map, negsubspace_affine = reduce_variables(Astar_eq_csc, sst.c_eq, check_eq_consistency)
            rmul!(negsubspace_affine, -1.)
            constraints = size(subspace_map, 2)
        end
    end

    @verbose_info("Constructing sparse representations, duals, and factorizations")
    A = Vector{SparseMatrixCSC{Float64,Int}}(undef, npsd)
    Astar = Vector{SparseMatrixCSC{Float64,Int}}(undef, npsd)
    AAstars = Vector{SparseMatrixCSC{Float64,Int}}(undef, npsd)
    C = Vector{PackedMatrix{Float64,SparseVector{Float64,Int}}}(undef, npsd)
    @inbounds for (i, (Astardata, mside, Ccoeff)) in enumerate(zip(sst.Astar, sst.psds, sst.Ccoeffs))
        m = packedsize(mside)
        coolen = length(Astardata[1])
        csrrowptr = Vector{Int}(undef, m +1)
        csrcolval = Vector{Int}(undef, coolen)
        csrnzval = Vector{Float64}(undef, coolen)
        Astar[i] = SparseArrays.sparse!(Astardata..., m, size(subspace_map, 1), +, tmp, csrrowptr, csrcolval, csrnzval,
            Astardata...)
        # Due to equality constraints, we don't have the full y available, but only a subset:
        # y = (subspace_map * ̃y - negsubspace_affine)
        # This can also lead to us having more than just the (1,1) component in C.
        C[i] = PackedMatrix(mside, Astar[i] * negsubspace_affine)
        C[i][1] += Ccoeff
        # For this, we now already form the transpose, although the reinterpretation has not taken place yet. But in
        # principle, what would have been ideal is to carry out the multiplication on Astardata (which is inefficient/
        # impossible), i.e., before sparse! - so it would be reflected as-is in A[i].
        # In the case of no equality constraints, we already have the transpose available, so let's short-circuit:
        if isempty(Astar_eq_coo[1])
            # We could instead check whether the subspace_map is the identity (which would then also take the short-circuit
            # in the case of redundant equality constraints), but this unlikely case is probably not worth it.

            # Now the csr data already contain the CSR representation of Astar[i] - which, when seen as CSC, corresponds to the
            # transpose. So we get A almost for free - however, the columns are still unsorted, so we have to do the sorting.
            for (from, toplus1) in zip(csrrowptr, Iterators.drop(csrrowptr, 1))
                to = toplus1 -1
                sort_along!(csrcolval, from, to, Base.Forward, csrnzval)
            end
            A[i] = SparseMatrixCSC(constraints, m, csrrowptr, csrcolval, csrnzval)
        else
            Astar[i] = Astar[i] * subspace_map
            # We temporarily create a very invalid CSC matrix, but transpose! overwrites it anyway. As the number of
            # columns of A does not change, no adjustments are needed.
            A[i] = SparseMatrixCSC(size(subspace_map, 2), m, csrrowptr, csrcolval, csrnzval)
            transpose!(A[i], Astar[i]) # use the already allocated memory, we don't need it any more
        end
        # Finally, note that we only ever work with the vectorized triangle. This is fine for A*, which is supposed to
        # return the matrix (or its vec representation). However, the role of each A is to be multiplied with another vec
        # representation and give a scalar - for this to be equal to the trace of the matrix product, we need to double the
        # off-diagonal elements. Note that csrrowptr will correspond to the _columns_ (since we are feeding into a CSC
        # matrix), and a column in turn corresponds to a whole matrix by means of vec; so indeed, we need to scale whole
        # batches of columns, i.e., entries pointed to in csrrowptr.
        curcol = 2
        nextdiag = 1
        while curcol ≤ m
            @assert(curcol + nextdiag ≤ m)
            lmul!(2, @view(csrnzval[csrrowptr[curcol]:csrrowptr[curcol+nextdiag]-1]))
            nextdiag += 1
            curcol += nextdiag
        end
        AAstars[i] = A[i] * Astar[i]
    end
    # We must add together all the AAstar in one matrix, which is then factorized. However, these are CSC matrices whose
    # addition is highly inefficient. But they may still be large, so we don't want to densify them, either. So we follow a
    # multi-step path, where we first calculate the sparsity pattern of the output before we do the actual addition.
    local AAstarinv
    @inbounds let m = constraints
        AAstarcolsize = 0
        for AAstarᵢ in AAstars
            AAstarcolsize += AAstarᵢ.colptr[end] -1
        end
        # now we have certain upper bounds on how large the individual columns can become.
        AAstarrowval = Vector{Int}(undef, AAstarcolsize)
        AAstarnzval = similar(AAstarrowval, Float64)
        AAstarcolptr = Vector{Int}(undef, m +1)
        j = 1
        for col in 1:m
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
        # We cannot simply use cholesky, as A A* will probably be singular.
        # ldlt might work, but can also lead to zero pivot values, so this will sometimes also fail.
        try
            AAstarinv = EfficientCholmod(ldlt(SparseMatrixCSC(m, m, AAstarcolptr, AAstarrowval, AAstarnzval)))
        catch
            AAstarinv = qr(SparseMatrixCSC(m, m, AAstarcolptr, AAstarrowval, AAstarnzval))
        end
    end
    maxn = maximum(sst.psds)
    bvec = SparseVector(size(subspace_map, 1), sst.b...)
    return StrideData(
        constraints, A, Astar, AAstarinv, C, (subspace_map, negsubspace_affine),
        transpose(subspace_map) * bvec, offset + dot(bvec, negsubspace_affine),
        [PackedMatrix{Float64}(undef, nᵢ) for nᵢ in sst.psds],
        [PackedMatrix{Float64}(undef, nᵢ) for nᵢ in sst.psds],
        [PackedMatrix{Float64}(undef, nᵢ) for nᵢ in sst.psds],
        [PackedMatrix{Float64}(undef, nᵢ) for nᵢ in sst.psds],
        Vector{Float64}(undef, constraints),
        Vector{Float64}(undef, constraints),
        Vector{Float64}(undef, max(packedsize(maxn), constraints)),
        finish!(monomial_map)
    )
end

function sparse_optimize(::Val{:Stride}, problem::PolyOptProblem{P,M,V}, groupings::Vector{<:Vector{<:AbstractVector{M}}};
    verbose::Bool=false, check_eq_consistency=false, warn_eq_consistency=true, kwargs...) where {P,M,V}
    # σ: step size for projected gradient descent
    @assert(!problem.complex && isone(problem.prefactor))
    # The STRIDE solver is an SDP solver specifically designed for problems from polynomial optimization. Hence, we need to
    # be able to quickly evaluate the constraints and objectives as well as construct the semidefinite relaxation.
    @verbose_info("Constructing relaxation data")
    setup_relaxation = @elapsed begin
        ssd = prepare_relaxation(problem, groupings, verbose, check_eq_consistency, warn_eq_consistency)
    end
    @verbose_info("Got relaxation data in ", setup_relaxation, " seconds. Generating initial point")
    # The relaxation was now prepared in the way that:
    # - the bound on the objective is ssd.offset - stride(ssd)
    # - the duals (moments) are the negative of the actual moments; and note that y₀ = 1 is not part of the variable list
    return stride(ssd, problem, groupings; verbose, kwargs...)
end

function psdproject!(m::PackedMatrix)
    eigs = eigen!(m, 0., Inf)
    fill!(m, zero(eltype(m)))
    for (eval, evec) in zip(eigs.values, eachcol(eigs.vectors))
        spr!(eval, evec, m)
    end
    return m
end

function stride(ssd::StrideData{R}, problem::PolyOptProblem{P,M,V}, groupings; verbose::Bool=false,
    tol::R=1e-8, maxiter::Integer=10,
    tol_init::R=1e-4, maxiter_init::Integer=1000, verbose_init::Bool=false,
    tol_sgsapg::R=1e-12, maxiter_sgsapg::Integer=1000, verbose_sgsapg::Bool=false,
    σ::R=10.,
    verbose_local::Bool=false, kwargs_local::Dict=Dict(),
    opti_local::Union{Function,Nothing}=begin
        @verbose_info("No local optimization specified; automatically constructing one (disable by passing opti_local=nothing).")
        setup_static = @elapsed(ol = lancelot_solve(problem; verbose=verbose_local, kwargs_local...))
        @verbose_info("Local optimization constructed in ", setup_static, " seconds.")
        ol
    end,
    opti_round::Union{Function,Nothing}=nothing) where {R<:Real,P,M,V}
    @assert(tol_init > 0 && maxiter_init > 0 && tol > 0 && σ > 0)
    # A: output, [X₁, ...], α, β -> α output + β A(X₁, ...)
    # Astar: output, y, i, α, β -> α output + β Xᵢ(y)
    # AAstarinv: output, y -> ̃y
    Astar, b, C, X, S, y, tmpm, mmap = ssd.Astar, ssd.b, ssd.C, ssd.X, ssd.S, ssd.y, ssd.tmpm1, ssd.monomial_map
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
    for (Xᵢ, Sᵢ) in zip(X, S)
        psdproject!(Xᵢ)
        psdproject!(Sᵢ)
    end
    Vopt = dot(b, ssd.y)
    #   > now we have positive X and S, but feasibility also requires the linear constraints to be satisfied. Check it for the
    #     dual problem:
    ηd = zero(R)
    for (Astarᵢ, Sᵢ, Cᵢ, tmpmᵢ) in zip(Astar, S, C, tmpm)
        copyto!(tmpmᵢ, Sᵢ)
        mul!(tmpmᵢ, Astarᵢ, y, true, true)
        axpy!(-one(R), Cᵢ, tmpmᵢ)
        ηd += norm(tmpmᵢ)
    end
    if ηd > tol
        # no, we are infeasible. So Vopt is not an appropriate value for the last good result.
        @verbose_info("Initial candidate generated by ADMM+ is infeasible with error ", ηd)
        Vopt = R(Inf)
    else
        @verbose_info("Initial candidate generated by ADMM+ is feasible")
    end
    by = Vopt
    # - a positive integer r ∈ [1, n] - this is the rank
    # - a positive constant ϵ > 0. According to the IEEE paper, choose 1e-12

    mons = problem.last_moments
    empty!(mons)
    sizehint!(mons, length(mmap) +1)
    mons[constant_monomial(M)] = one(R)
    full_y = Vector{Float64}(undef, length(mmap))
    subspace_map_inv = qr(ssd.equality_subspace[1])
    yinv_large = similar(full_y) # To solve the overdetermined system, we need a buffer that is as large as full_y
    new_y = @view(yinv_large[1:length(y)])

    status = :max_iter
    σₖ = zero(R)
    # Iterate the following steps for k = 1, ...:
    @inbounds for k in 1:maxiter
        σₖ += σ
        # Step 1 (Projection). Compute (̄Xₖ, yₖ, Sₖ) satisfying (5) and (6), where ̄Xₖ is an inexact projection of Xₖ₋₁ - σₖ C onto
        # the primal feasible set of the SDP, by running algorithm sgsapg and algorithm mlbfgs, warmstarted by (yₖ₋₁, Sₖ₋₁).
        # → Section 4
        @inbounds for (Xᵢ, Cᵢ) in zip(X, C)
            axpy!(-σₖ, Cᵢ, Xᵢ)
        end
        sgsapg!(ssd, tol_sgsapg, maxiter_sgsapg, verbose_sgsapg)
        # We got an updated y and S, now we have to construct X = Π(A*ξ + Z), where Z = Xₖ₋₁ - σₖ C ≡ X

        for (Xᵢ, Astarᵢ) in zip(X, Astar)
            mul!(Xᵢ, Astarᵢ, y, true, true)
            psdproject!(Xᵢ)
        end

        # (̄Xₖ, yₖ, Sₖ) = (X(W, ξ), 1/σₖ ξ, 1/σₖ W) [note that ξ and y are aliased, as well as W and S]
        σinvₖ = inv(σₖ)
        rmul!(y, σinvₖ)
        rmul!.(S, σinvₖ)

        # Step 2 (Certification). If η(̄Xₖ, yₖ, Sₖ) < Tol, output (̄Xₖ, yₖ, Sₖ) and stop.
        CX, by, ηp, ηd, ηg = η(ssd, Val{true})
        (verbose_sgsapg || !isnothing(opti_local) || isone(k)) &&
            @verbose_info("Iteration | Primal objective | Dual objective | Primal infeasibility | Dual infeasibility | Relative duality gap")
        @verbose_info(@sprintf("%9d | %16g | %14g | %20g | %18g | %20g", k, CX, by, ηp, ηd, ηg))
        if max(ηp, ηg, ηd) ≤ tol
            @verbose_info("Successfully converged.")
            Vopt = by
            status = :ok
            break
        end

        # Step 3 (Acceleration): Compute ̂Xₖ via rounding, local NLP optimization and lifting, starting from ̄Xₖ:
        if !isnothing(opti_local)
            # Original description:
            # ̄Xₖ = ∑ᵢ λᵢ vᵢ vᵢᵀ (truncate after r)
            # ̄xₖᵢ = rounding(vᵢ). Problem-dependent, but can first normalize vᵢ such that the constant component is 1, then
            #                     extract the order-1 entries. If easily possible, project to the feasible set.
            # This would be problematic, to the claim in the Yang paper that "the algorithms...can be extended...in a
            # straightforward way to solve sparse relaxations. We now have multiple X, contrary to the Yang paper. If those
            # just came from constraints, we would simply disregard all but the moment matrix. However, we allow for generic
            # sparsity structure, i.e., there may be multiple sub-moment matrices with possibly overlapping variables. To
            # account for this, we will use our solution extraction heuristic. This is different from the eigenvector
            # decomposition method!
            # (Also note that the poor-man's option to add a dense first-order moment matrix for each clique suggested
            # previously by Wang et al. [CS-TSSOS] appears to be rather pointless.)
            # However, note that we are in a different situation than Yang. We created the whole problem based on the moments,
            # and we made the correspondence with their dual problem (whereas they somehow attribute the SOS for to the dual).
            # Consequently, we get all of our moments from the y (and S is the moment matrix), and it is also this y and S that
            # will then be overwritten with a potentially better candidate. While this is unbalanced - indeed, we do not modify
            # the primal X, it is also unbalanced in Yang's algorithm, where they do not modify the duals.
            @verbose_info("Entering local optimization")
            # Note that S = C - A* y; if we already take into account that
            # y_full = equality_subspace[1] * y - equality_subspace[2]
            # then S = Cᵣ - A*ᵣ y_full, where Cᵣ contains the (1, 1) -> 1 contribution only and A*ᵣ merely distributes the
            # moments in matrix form (the B matrices with Lasserre). However, sgsapg will not yield perfect results, therefore
            # this equality does not hold exactly, and while S ⪰ 0, we don't have C - A* y ⪰ 0. But then, S does not satisfy
            # the structural constraints of a moment matrix. So here, we decide to get the "moment matrix" from y, although it
            # is not PSD.
            copyto!(full_y, ssd.equality_subspace[2])
            mul!(full_y, ssd.equality_subspace[1], y, -one(R), true)
            for (i, mon) in enumerate(mmap)
                mons[mon] = full_y[i]
            end
            solutions = isnothing(opti_round) ? poly_solutions_heuristic(problem; verbose) : opti_round(problem; verbose)
            # ̂xₖᵢ = nlp(̄xₖᵢ)
            # ̂xₖ = argmin p(̂xₖᵢ) over all { ̂xₖᵢ : i = 1, ..., r }
            best = R(Inf)
            local best_solution
            for solution in solutions
                for i in eachindex(solution)
                    if isnan(solution[i])
                        solution[i] = rand()
                    end
                end
                val, pos = opti_local(solution) # note: we REQUIRE solutions returned by opti_local to be feasible (should
                                                # return (Inf, arbitrary) if not)!
                # opti_local is allowed to modify its parameter, e.g., by writing pos into solution.
                if val < best
                    best_solution = pos
                    best = val
                end
            end
            # ̂Xₖ = proj(monomial lifting of ̂xₖ)
            # Step 4 (Policy). Choose the better candidate in {̄Xₖ, ̂Xₖ} to update Xₖ:
            # Xₖ = ̂Xₖ if ⟨C, ̂Xₖ⟩ < min( ⟨C, ̄Xₖ⟩, min_{X ∈ V} { ⟨C, X⟩ : X ∈ V } ) - ϵ ∧ ̂Xₖ feasible
            # We decide to be a bit more careful here, which is necessary due to the absence of BFGS. Our SGSAPG result may be
            # very infeasible, so it is not fair to compare to this value unless both feasibilities are sufficienty well under
            # control.
            if ηp < 100tol && ηd < 100tol
                bound = min(Vopt, by)
            else
                bound = Vopt
            end
            if best < bound - 1e-12
                @verbose_info("Found better optimization candidate with objective ", best)
                substitution = problem.variables => best_solution
                for (i, mon) in enumerate(mmap)
                    full_y[i] = mon(substitution)
                end
                full_y .+= ssd.equality_subspace[2]
                ldiv!(yinv_large, subspace_map_inv, full_y)
                # now yinv_large aka new_y is a least squares solution, but we don't know whether the solution was actually
                # possible (it should be, if the local solver obeyed the constraints). Maybe checking
                # norm(equality_subspace[1] * y - equality_subspace[2]) < ϵ would be a good idea?
                # Let's update our data.
                copyto!(y, new_y)
                for (Sᵢ, Astarᵢ, Cᵢ) in zip(S, Astar, C)
                    mul!(Sᵢ, Astarᵢ, y, -one(R), false)
                    axpy!(one(R), Cᵢ, Sᵢ)
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

    copyto!(full_y, ssd.equality_subspace[2])
    mul!(full_y, ssd.equality_subspace[1], y, -one(R), true)
    for (i, mon) in enumerate(mmap)
        mons[mon] = full_y[i]
    end
    return status, Vopt
    #=
    λmins = Vector{R}(undef, length(groupings))
    for (i, (tmpmᵢ, Cᵢ, Astarᵢ)) in enumerate(zip(tmpm, C, Astar))
        copyto!(tmpmᵢ, Cᵢ)
        mul!(tmpmᵢ, Astarᵢ, y, -one(R), true)
        @inbounds λmins[i] = eigvals!(tmpmᵢ, 1:1)[1]
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
=#
end

function spadmm!(ssd::StrideData{R}, σ::R, γ::R, tol::R, maxiter::Integer, verbose::Bool) where {R<:Real}
    @assert(0 ≤ γ ≤ 2 && σ > 0 && tol ≥ 0 && maxiter > 0)
    A, Astar, AAstarinv, C, b, X, S, y = ssd.A, ssd.Astar, ssd.AAstarinv, ssd.C, ssd.b, ssd.X, ssd.S, ssd.y
    tmpm = ssd.tmpm1
    @verbose_info("Entering ADMM+ with penalty σ = ", σ, ", step length γ = ", γ, ", and tolerance ", tol,
        "\nIteration | Primal objective | Dual objective | Primal infeasibility | Dual infeasibility | Relative duality gap")
    # Initialization: Initial points X₀ = S₀ = 0 ∈ 𝕊ⁿ
    fill!.(X, zero(R))
    fill!.(S, zero(R))
    σ⁻¹ = inv(σ)
    γσ = γ * σ
    local y
    # Iterate the following steps for k = 1, ...
    @inbounds for k in 1:maxiter
        # Step 1: Compute
        #    ̂yₖ₊₁ = (A A*)⁻¹(b/σ - A(Xₖ/σ + Sₖ - C))
        for (i, (Xᵢ, Sᵢ, Cᵢ, Aᵢ)) in enumerate(zip(X, S, C, A))
            # We don't need Sᵢ any more, so let's overwrite it
            Sᵢ .+= Xᵢ .* σ⁻¹
            axpy!(-one(R), Cᵢ, Sᵢ) # we use axpy! for all sparse changes instead of broadcasting, as this can exploit the
                                   # sparsity pattern. In contrast, axpy! seems to be much worse than broadcasting for dense
                                   # vectors.
            mul!(y, Aᵢ, Sᵢ, -one(R), !isone(i))
        end
        axpy!(σ⁻¹, b, y)
        ldiv!(AAstarinv, y)
        # Step 2: Compute
        #    Sₖ₊₁ = ( Π(Xₖ + σ(A* ̂yₖ₊₁ - C)) - (Xₖ + σ(A* ̂yₖ₊₁ - C)) ) / σ
        #             ^ projection onto PSD, use LAPACK.syevx
        for (Xᵢ, Sᵢ, tmpmᵢ, Astarᵢ, Cᵢ) in zip(X, S, tmpm, Astar, C)
            mul!(tmpmᵢ, Astarᵢ, y, σ, false)
            axpy!(-σ, Cᵢ, tmpmᵢ)
            tmpmᵢ .+= Xᵢ
            copyto!(Sᵢ, tmpmᵢ)
            # We effectively project onto the (absolute of the) negative part, which is done by a positive projection, since
            # the positive part is expected to have low rank
            psdproject!(Sᵢ)
            Sᵢ .= (Sᵢ .- tmpmᵢ) .* σ⁻¹
        end
        # Step 3: Compute
        #    yₖ₊₁ = (A A*)⁻¹(b/σ - A(Xₖ/σ + Sₖ₊₁ - C))
        for (i, (Xᵢ, Sᵢ, Cᵢ, tmpmᵢ, Aᵢ)) in enumerate(zip(X, S, C, tmpm, A))
            # This time, we need the Sᵢ again (as well as the Xs), so we need the temporary storage
            tmpmᵢ .= Xᵢ .* σ⁻¹ .+ Sᵢ
            axpy!(-one(R), Cᵢ, tmpmᵢ)
            mul!(y, Aᵢ, tmpmᵢ, -one(R), !isone(i))
        end
        axpy!(σ⁻¹, b, y)
        ldiv!(AAstarinv, y)
        # Step 4: Compute
        #    Xₖ₊₁ = Xₖ + γ σ (Sₖ₊₁ + A* yₖ₊₁ - C)
        for (Xᵢ, Sᵢ, Cᵢ, Astarᵢ) in zip(X, S, C, Astar)
            Xᵢ .+= γσ .* Sᵢ
            axpy!(-γσ, Cᵢ, Xᵢ)
            mul!(Xᵢ, Astarᵢ, y, γσ, true)
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
    A, Astar, AAstarinv, b, Z, W, ξ = ssd.A, ssd.Astar, ssd.AAstarinv, ssd.b, ssd.X, ssd.S, ssd.y
    tmpm1, tmpm2, tmpv = ssd.tmpm1, ssd.tmpm2, ssd.tmpv
    @verbose_info("Entering sGS-based accelerated proximal gradient method\nIteration |     residue")
    # Initialization: Set ̃W₁ = W₀ and t₁ = 1
    copyto!.(tmpm1, W)
    t = 1.
    # Iterate the following steps for k = 1, ...
    k = 1
    @inbounds while true
        # Step 1 (sGS update): Compute
        #    ̃ξₖ = (A A*)⁻¹(b - A(Z) - A(̃Wₖ))
        for (i, (Aᵢ, Zᵢ, tmpm1ᵢ)) in enumerate(zip(A, Z, tmpm1))
            tmpm1ᵢ .+= Zᵢ
            mul!(tmpv, Aᵢ, tmpm1ᵢ, -one(R), !isone(i))
        end
        tmpv .+= b
        copyto!.(tmpm1, W) # we need to back up the previous W
        ldiv!(ξ, AAstarinv, tmpv)
        #    Wₖ = Π(-A*̃ξₖ - Z) = Π(A*ξₖ + Z) - (A*ξₖ + Z)
        for (Zᵢ, Wᵢ, tmpm2ᵢ, Astarᵢ) in zip(Z, W, tmpm2, Astar)
            copyto!(Wᵢ, Zᵢ)
            mul!(Wᵢ, Astarᵢ, ξ, true, true)
            copyto!(tmpm2ᵢ, Wᵢ)
            psdproject!(Wᵢ)
            Wᵢ .-= tmpm2ᵢ
        end
        #    ξₖ = (A A*)⁻¹(b - A(Z) - A(Wₖ))
        for (i, (Aᵢ, Zᵢ, Wᵢ)) in enumerate(zip(A, Z, W))
            mul!(tmpv, Aᵢ, Zᵢ, -one(R), !isone(i))
            mul!(tmpv, Aᵢ, Wᵢ, -one(R), true)
        end
        tmpv .+= b
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

function η(ssd::StrideData{R}, details::Union{Type{Val{false}},Type{Val{true}}}=Val{false}) where {R}
    A, Astar, C, b, X, S, y, tmpm, tmpv = ssd.A, ssd.Astar, ssd.C, ssd.b, ssd.X, ssd.S, ssd.y, ssd.tmpm1, ssd.tmpv
    # here, note that b = (1, 0, ..., 0)
    # ηp(X) = ‖A(X) - b‖/(1 + ‖b‖)
    for (i, (Aᵢ, Xᵢ)) in enumerate(zip(A, X))
        mul!(tmpv, Aᵢ, Xᵢ, true, !isone(i))
    end
    @inbounds tmpv .-= b
    ηp = norm(tmpv) / (1 + norm(b))
    # ηd(y, S) = ‖A*(y) + S - C‖/(1 + ‖C‖)
    # ηg(X, y) = |⟨C, X⟩ - ⟨b, y⟩|/(1 + |⟨C, X⟩| + |⟨b, y⟩|)
    ηdnumsq, normCsq, CX, by = zero(R), zero(R), zero(R), dot(b, y)
    @inbounds for (Sᵢ, Xᵢ, Cᵢ, tmpmᵢ, Astarᵢ) in zip(S, X, C, tmpm, Astar)
        mul!(tmpmᵢ, Astarᵢ, y, true, false)
        tmpmᵢ .+= Sᵢ
        axpy!(-one(R), Cᵢ, tmpmᵢ)
        CX += dot(Cᵢ, Xᵢ)
        ηdnumsq += LinearAlgebra.norm2(tmpmᵢ)^2
        normCsq += LinearAlgebra.norm2(Cᵢ)^2
    end
    ηd = sqrt(ηdnumsq) / (1 + sqrt(normCsq))
    ηg = abs(CX - by) / (1 + abs(CX) + abs(by))
    return details == Val{true} ? (ssd.offset - CX, ssd.offset - by, ηp, ηd, ηg) : max(ηp, ηd, ηg)
end

function ηproj(ssd::StrideData{R}, W, ξ, recompute_AstarξplusZ::Bool=true) where {R}
    A, Astar, b, Z, X = ssd.A, ssd.Astar, ssd.b, ssd.X, ssd.tmpm2
    tmpmvec, tmpv = ssd.tmpvlong, ssd.tmpv
    # X(W, ξ) := A*ξ + W + Z
    # ηproj := max( ‖A(X(W, ξ)) - b‖, ‖X(W, ξ) - Π(A*ξ + Z)‖ )
    norm2 = zero(R)
    if recompute_AstarξplusZ
        @inbounds for (Xᵢ, Wᵢ, Zᵢ, Astarᵢ) in zip(X, W, Z, Astar)
            tmpmᵢ = let s = size(Xᵢ, 1)
                PackedMatrix(s, view(tmpmvec, 1:packedsize(s)))
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
        @inbounds for (Xᵢ, Wᵢ, Zᵢ, Astarᵢ) in zip(X, W, Z, Astar)
            copyto!(Xᵢ, Zᵢ)
            mul!(Xᵢ, Astarᵢ, ξ, true, true)
            Xᵢ .+= Wᵢ
        end
    end
    for (i, (Aᵢ, Xᵢ)) in enumerate(zip(A, X))
        mul!(tmpv, Aᵢ, Xᵢ, true, !isone(i))
    end
    axpy!(-one(R), b, tmpv)
    norm1 = norm(tmpv)
    return max(norm1, sqrt(norm2))
end