function PolynomialOptimization.tighten_minimize_l1(::Val{:COPT}, spmat, rhs)
    task = COPTProb(copt_env)
    _check_ret(copt_env, COPT_SetIntParam(task, COPT_INTPARAM_LOGGING, zero(Cint)))
    lensp, len = size(spmat)
    sp_cols = Vector{Cint}(undef, 2len +1)
    sp_rows = Vector{Cint}(undef, nnz(spmat) + 4len)
    sp_nz = similar(sp_rows, Cdouble)
    oldcols = SparseArrays.getcolptr(spmat)
    oldrows = rowvals(spmat)
    oldnz = nonzeros(spmat)
    i = 1
    ci = lensp
    @inbounds for j in 1:len
        # Make sure the original variables satisfy spmat*x = rhs
        sp_cols[j] = i -1
        @simd for k in oldcols[j]:oldcols[j+1]-1
            sp_rows[i] = oldrows[k] -1
            sp_nz[i] = oldnz[k]
            i += 1
        end
        # For the ℓ₁ norm, we need the absolute values, so another len variables.
        # And the others are larger or equal to the absolute values.
        sp_rows[i] = ci
        sp_rows[i+1] = ci +1
        sp_nz[i] = -1.
        sp_nz[i+1] = 1.
        i += 2
        ci += 2
    end
    ci = lensp
    @inbounds for j in len+1:2len
        sp_cols[j] = i -1
        sp_rows[i] = ci
        sp_rows[i+1] = ci +1
        sp_nz[i] = 1.
        sp_nz[i+1] = 1.
        i += 2
        ci += 2
    end
    @inbounds sp_cols[2len+1] = i-1
    obj = Vector{Cdouble}(undef, 2len)
    fill!(@view(obj[1:len]), 0.)
    fill!(@view(obj[len+1:end]), 1.)
    bnd = Vector{Cdouble}(undef, 2len)
    fill!(@view(bnd[1:len]), -COPT_INFINITY)
    fill!(@view(bnd[len+1:end]), 0.)
    nrhsₗ = Vector{Cdouble}(undef, lensp + 2len)
    copyto!(nrhsₗ, rhs)
    fill!(@view(nrhsₗ[lensp+1:end]), 0.)
    nrhsᵤ = similar(nrhsₗ)
    copyto!(nrhsᵤ, rhs)
    fill!(@view(nrhsᵤ[lensp+1:end]), COPT_INFINITY)
    _check_ret(copt_env, COPT_LoadProb(
        task, 2len, lensp + 2len, COPT_MINIMIZE, 0., obj,
        sp_cols, C_NULL, sp_rows, sp_nz, C_NULL, bnd, C_NULL, C_NULL, nrhsₗ, nrhsᵤ,
        C_NULL, C_NULL
    ))
    status = Ref{Cint}()
    _check_ret(copt_env, COPT_SolveLp(task))
    _check_ret(copt_env, COPT_GetIntAttr(task, COPT_INTATTR_LPSTATUS, status))
    if status[] == COPT_LPSTATUS_OPTIMAL
        x = Vector{Cdouble}(undef, 2len)
        _check_ret(copt_env, COPT_GetLpSolution(task, x, C_NULL, C_NULL, C_NULL))
        return resize!(x, len)
    else
        throw(SingularException(0))
    end
end