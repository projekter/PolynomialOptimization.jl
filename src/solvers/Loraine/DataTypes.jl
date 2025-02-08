export Model, Solver

struct Model{T<:Real,I<:Integer}
    A::Vector{SparseMatrixCSC{T,I}}
    C::Vector{SparseMatrixCSC{T,I}}
    b::Vector{T}
    c_lin::SparseVector{T,I}
    A_lin::SparseMatrixCSC{T,I}
    coneDims::Vector{Int}
    n::Int
    nlin::Int
    nlmi::Int

    @doc raw"""
        Model(; A, C, b, [c_lin, A_lin,] coneDims, check=true)

    Constructs a model for `PolynomialOptimization`'s rewrite of the Loraine solver. This solves the problem
    ```math
        \min \vec c_{\mathrm{lin}} \cdot \vec x + \sum_j \langle C_j, X_j\rangle \\
        \text{such that} \\
        x_i \geq 0 \ \forall i \\
        X_j \succeq 0 \forall j \\
        A_{\mathrm{lin}, k} \cdot \vec x + \sum_j \langle \operatorname{mat}(A_{k, j}), X_k\rangle = b_k \ \forall k
    ```
    with the following representation in the variables:
    - `1 ≤ j ≤ length(coneDims) = length(A) = length(C)`
    - `A` is a vector of `SparseMatricCSC`s, where every matrix corresponds to a semidefinite optimization variable. Each row
      in the matrix corresponds to one constraint; when the row is reshaped into a matrix in a col-major order (which must then
      be symmetric), it defines the coefficient matrix ``\operatorname{mat}(A_{k, j})``. The side dimension of the matrix is
      stored in `coneDims`.
    - If present, `c_lin` is a `SparseVector` and `A_lin` a `SparseMatrixCSC` that define the objective and constraints
      coefficients for the nonnegative variables. They may be omitted only together.

    See also [`Solver`](@ref), [`solve!`](@ref).

    !!! info
        This function checks the validity of the variables; however, it is quite expensive to do the checks for symmetry. They
        can be disabled by setting `check` to `false`.
    """
    function Model(; A::Vector{SparseMatrixCSC{T,I}}, C::Vector{SparseMatrixCSC{T,I}}, b::Vector{T},
        c_lin::Union{SparseVector{T,I},Nothing}=nothing, A_lin::Union{Nothing,SparseMatrixCSC{T,I}}=nothing,
        coneDims::Vector{Int}, check::Bool=true) where {T<:Real,I<:Integer}
        isnothing(c_lin) == isnothing(A_lin) || throw(ArgumentError("c_lin and A_lin must both be present or both be omitted"))
        length(A) == length(C) == length(coneDims) || throw(ArgumentError("Inconsistent number of PSD matrices"))
        if isnothing(c_lin)
            c_lin = spzeros(T, I, 0)
            A_lin = spzeros(T, I, length(b), 0)
        else
            length(c_lin) == size(A_lin, 2) || throw(ArgumentError("Inconsistent number of linear variables"))
            length(b) == size(A_lin, 1) || throw(ArgumentError("Inconsistent number of constraints"))
        end
        for (Aᵢ, Cᵢ, s) in zip(A, C, coneDims)
            size(Aᵢ, 1) == length(b) || throw(ArgumentError("Inconsistent number of constraints"))
            (size(Aᵢ, 2) == LinearAlgebra.checksquare(Cᵢ)^2 && size(Cᵢ, 1) == s) || throw("Inconsistent PSD matrix dimensions")
        end
        if check
            all(issymmetric, C) || throw(ArgumentError("The objective matrix must be symmetric"))
            all(x -> all(r -> issymmetric(reshape(r, (x[2], x[2]))), eachrow(x[1])), zip(A, coneDims)) ||
                throw(ArgumentError("The constraint matrices must be symmetric"))
        end
        new{T,I}(A, C, b, c_lin, A_lin, coneDims, length(b), length(c_lin), length(coneDims))
    end
end

mutable struct Solver{T<:Real,I<:Integer}
    tol_cg::T
    tol_cg_up::T
    tol_cg_min::T
    eDIMACS::T
    preconditioner::Preconditioner
    erank::Int
    aamat::AMat
    fig_ev::Int
    verb::Verbosity
    initpoint::Initpoint
    maxit::Int

    const model::Model{T,I}

    tau::T
    expon::T

    const X::Vector{Matrix{T}} # symmetric
    const S::Vector{Matrix{T}} # symmetric
    const y::Vector{T}
    const delX::Vector{Matrix{T}} # symmetric
    const delS::Vector{Matrix{T}} # symmetric
    const dely::Vector{T}
    const Xn::Vector{Matrix{T}} # symmetric
    const Sn::Vector{Matrix{T}} # symmetric

    const X_lin::Vector{T}
    const S_lin::Vector{T}
    const Si_lin::Vector{T}
    const S_lin_inv::Vector{T}
    const delX_lin::Vector{T}
    const delS_lin::Vector{T}
    const Xn_lin::Vector{T}
    const Sn_lin::Vector{T}

    const D::Vector{Vector{T}}
    const G::Vector{Matrix{T}}
    const Gi::Vector{Matrix{T}}
    const W::Vector{Matrix{T}} # symmetric
    const Si::Vector{Matrix{T}} # symmetric
    const DDsi::Vector{Vector{T}}

    const Rp::Vector{T}
    const Rd::Vector{Matrix{T}} # symmetric
    const Rd_lin::Vector{T}
    const rhs::Vector{T}

    const RNT::Vector{Matrix{T}} # symmetric
    const RNT_lin::Vector{T}

    const alpha::Vector{T}
    const beta::Vector{T}
    alpha_lin::T
    beta_lin::T

    cg_iter_tot::Int
    cg_iter_pre::Int
    cg_iter_cor::Int

    sigma::T
    mu::T
    iter::Int
    DIMACS_error::T
    status::Status

    err1::T
    err2::T
    err3::T
    err4::T
    err5::T
    err6::T

    @doc """
        Solver(model; tol_cg=0.01, tol_cg_up=0.5, tol_cg_min=1e-7, eDIMACS=1e-7,
            preconditioner=PRECONDITIONER_HALPHA, erank=1, aamat=AMAT_DIAGAᵀA, fig_ev=0,
            verb=VERBOSITY_SHORT, initpoint=INITPOINT_LORAINE, maxit=100)

    Defines a solver for a previously defined model. Only the iterative conjugate gradient method is implemented.

    See also [`Model`](@ref), [`solve!`](@ref).
    """
    function Solver(model::Model{T,I}; tol_cg::T=T(1//100), tol_cg_up::T=T(1//2), tol_cg_min::T=T(1//10_000_000),
        eDIMACS::T=T(1//10_000_000), preconditioner::Preconditioner=PRECONDITIONER_HALPHA, erank::Int=1,
        aamat::AMat=AMAT_DIAGAᵀA, fig_ev::Int=0, verb::Verbosity=VERBOSITY_SHORT, initpoint::Initpoint=INITPOINT_LORAINE,
        maxit::Int=100) where {T<:Real,I<:Integer}
        iszero(model.nlmi) && error("No linear matrix inequalities")
        if tol_cg < tol_cg_min
            solver.tol_cg = tol_cg_min
            @printf(" ---Parameter tol_cg smaller than tol_cg_min, setting tol_cg = %7.1e\n", solver.tol_cg)
        end
        if tol_cg_min > eDIMACS
            solver.tol_cg_min = eDIMACS
            @printf(" ---Parameter tol_cg_min switched to eDIMACS = %7.1e\n", eDIMACS)
        end
        if erank < 0
            solver.erank = 1
            @printf(" ---Parameter erank negative, setting erank = %1d\n", solver.erank)
        end

        coneDims = model.coneDims
        solver = new{T,I}(
            tol_cg, tol_cg_up, tol_cg_min, eDIMACS, preconditioner, erank, aamat, fig_ev, verb, initpoint, maxit,
            model,

            T(95//100), # tau; lower value such as 0.9 leads to more iterations but more robust algo
            T(3),       # expon

            [Matrix{T}(undef, s, s) for s in coneDims], # X
            [Matrix{T}(undef, s, s) for s in coneDims], # S
            Vector{T}(undef, model.n),                  # y
            [Matrix{T}(undef, s, s) for s in coneDims], # delX
            [Matrix{T}(undef, s, s) for s in coneDims], # delS
            Vector{T}(undef, model.n),                  # dely
            [Matrix{T}(undef, s, s) for s in coneDims], # Xn
            [Matrix{T}(undef, s, s) for s in coneDims], # Sn

            Vector{T}(undef, model.nlin),               # X_lin
            Vector{T}(undef, model.nlin),               # S_lin
            Vector{T}(undef, model.nlin),               # Si_lin
            Vector{T}(undef, model.nlin),               # S_lin_inv
            Vector{T}(undef, model.nlin),               # delX_lin
            Vector{T}(undef, model.nlin),               # delS_lin
            Vector{T}(undef, model.nlin),               # Xn_lin
            Vector{T}(undef, model.nlin),               # Sn_lin

            [Vector{T}(undef, s)    for s in coneDims], # D
            [Matrix{T}(undef, s, s) for s in coneDims], # G
            [Matrix{T}(undef, s, s) for s in coneDims], # Gi
            [Matrix{T}(undef, s, s) for s in coneDims], # W
            [Matrix{T}(undef, s, s) for s in coneDims], # Si
            [Vector{T}(undef, s)    for s in coneDims], # DDsi

            similar(model.b),                           # Rp
            [Matrix{T}(undef, s, s) for s in coneDims], # Rd
            Vector{T}(undef, model.nlin),               # Rd_lin
            similar(model.b),                           # rhs

            [Matrix{T}(undef, s, s) for s in coneDims], # RNT
            Vector{T}(undef, model.nlin),               # RNT_lin

            Vector{T}(undef, model.nlmi),               # alpha
            Vector{T}(undef, model.nlmi),               # beta
        )

        if verb > VERBOSITY_NONE
            println("\n *** Loraine solver from PolynomialOptimization, based on Loraine.jl v0.2.5 ***\n *** Initialisation STARTS")
            @printf(" Number of variables: %5d\n", model.n)
            @printf(" LMI constraints    : %5d\n", model.nlmi)
            println(" Matrix size(s)     :")
            Printf.format.((stdout,), (Printf.Format("%6d"),), model.coneDims)
            @printf("\n Linear constraints : %5d\n", model.nlin)
            println(" Preconditioner     : ", preconditioner)
        end

        return solver
    end
end