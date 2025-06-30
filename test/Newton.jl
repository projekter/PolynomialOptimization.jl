include("./shared.jl")
using LinearAlgebra
using PolynomialOptimization.Newton: halfpolytope, newton_methods
using Combinatorics: combinations

function remove_control(s::AbstractString)
    while startswith(s, "\e[2K") # clear line control character
        pos = findfirst('\r', s)
        if isnothing(pos)
            s = s[5:end]
        else
            s = s[pos+1:end] # actually clear the line
        end
    end
    while endswith(s, "\r")
        s = s[1:end-1]
    end
    return s
end
readlines_removecontrol(io::IO) = remove_control.(readlines(io))

# Adapter from https://discourse.julialang.org/t/a-minimal-example-with-base-redirect-stdout/64245/8
function capture_stdout(f::Function, @nospecialize(readfn=readlines_removecontrol))
    pipe = Pipe()
    started = Base.Event()
    writer = @async redirect_stdout($pipe) do
        try
            notify($started)
            return $f()
        finally
            close(Base.pipe_writer($pipe))
        end
    end
    wait(started)
    result = readfn(pipe)
    return fetch(writer), result
end

macro testnewtons(fn)
    esc(quote
        for method in newton_methods
            @testset let method=method
                $fn
            end
        end
    end)
end

struct FixAll{F,T,K} <: Function
    f::F
    args::T
    kwargs::K

    FixAll(f::F, args...; kwargs...) where {F} = new{F,typeof(args),typeof(kwargs)}(f, args, kwargs)
end

(fa::FixAll)() = fa.f(fa.args...; fa.kwargs...)

# BPT: Blekhermann, Parrilo, Thomas - Semidefinite Optimization and Convex Algebraic Geometry
# Actually, there only the full Newton polytope is given, so for the examples, the half polytope was calculated manually.
@testset "Newton polytope (BPT Example 3.92)" begin
    DynamicPolynomials.@polyvar x y
    obj = 5 - x*y - x^2*y^2 + 3y^2 + x^4
    result = monomial_vector([1, y, x, x*y, x^2])
    @testnewtons begin
        # First check the quick algo
        data, output = capture_stdout(FixAll(halfpolytope, method, obj, verbose=true))
        @test data == result
        @test output[6] == "Removing redundancies from the convex hull - quick heuristic, 5 initial candidates"
        @test startswith(output[8], "Found 4 potential extremal points of the convex hull in")
        # Then the same for the fine algo
        data, output = capture_stdout(FixAll(
            halfpolytope, method, obj, preprocess_quick=false, preprocess_fine=true, verbose=true
        ))
        @test data == result
        @test output[6] == "Removing redundancies from the convex hull - fine, 5 initial candidates"
        @test startswith(output[7], "Found 4 extremal points of the convex hull in")
    end
end

@testset "Newton polytope (BPT Example 3.93)" begin
    DynamicPolynomials.@polyvar x y
    obj = 1 - x^2 + x*y + 4y^4
    result = monomial_vector([1, y, x, y^2])
    @testnewtons begin
        data, output = capture_stdout(FixAll(halfpolytope, method, obj, verbose=true))
        @test data == result
        @test output[6] == "Removing redundancies from the convex hull - quick heuristic, 4 initial candidates"
        @test startswith(output[8], "Found 3 potential extremal points of the convex hull")

        data, output = capture_stdout(FixAll(
            halfpolytope, method, obj, preprocess_quick=false, preprocess_fine=true, verbose=true
        ))
        @test data == result
        @test output[6] == "Removing redundancies from the convex hull - fine, 4 initial candidates"
        @test startswith(output[7], "Found 3 extremal points of the convex hull")
    end
end

@testset "Newton polytope (BPT, Example 3.95)" begin
    DynamicPolynomials.@polyvar w x y z
    obj = PolynomialOptimization.IntPolynomial((w^4 + 1) * (x^4 + 1) * (y^4 + 1) * (z^4 + 1) + 2w + 3x + 4y + 5z)
    result = monomials(4, 0, 0:8, minmultideg=fill(0, 4), maxmultideg=fill(2, 4))
    @testnewtons begin
        data, output = capture_stdout(FixAll(halfpolytope, method, obj, verbose=true))
        @test data == result
        @test output[6] == "Removing redundancies from the convex hull - quick heuristic, 20 initial candidates"
        @test startswith(output[8], "Found 16 potential extremal points of the convex hull in")

        data, output = capture_stdout(FixAll(
            halfpolytope, method, obj, preprocess_quick=false, preprocess_fine=true, verbose=true
        ))
        @test data == result
        @test output[6] == "Removing redundancies from the convex hull - fine, 20 initial candidates"
        @test startswith(output[7], "Found 16 extremal points of the convex hull in")
    end
end

# Motzkin (BPT Exercise 3.97) is already checked as part of the documentation

@testset "Newton polytope (BPT Exercise 4.5)" begin
    DynamicPolynomials.@polyvar x[1:3]
    @testnewtons @test halfpolytope(method, x[1]*x[2]^2 + x[2]^2 + prod(x)) == monomial_vector([x[2]])
end

@testset "Newton polytope (mutually unbiased bases)" begin
    function findmubs(d, n)
        DynamicPolynomials.@complex_polyvar x[1:d, 1:d, 1:n]
        obj = zero(polynomial_type(x[1,1,1], Float64))
        # make it an ONB
        @views for i in 1:n
            obj += sum(z -> real(z)^2 + imag(z)^2, x[:, :, i] * x[:, :, i]' - I)
        end
        fac = inv(Float64(d))
        for i in 1:n-1
            for j in 1:i-1
                obj += sum(z -> (real(z)^2 + imag(z)^2 - fac)^2, x[:, :, i] * x[:, :, j]')
            end
        end
        obj
    end

    obj = findmubs(2, 3)
    @testnewtons begin
        data, output = capture_stdout(FixAll(halfpolytope, method, obj, preprocess=true, verbose=true))
        @test length(data) == 1453
        @test output[6] == "Removing redundancies from the convex hull - quick heuristic, 689 initial candidates"
        @test startswith(output[8], "Found 517 potential extremal points of the convex hull in")
        @test output[9] == "Removing redundancies from the convex hull - randomized, 517 initial candidates"
        m = match(r"^Found (\d+) extremal points of the convex hull via randomization in", output[10])
        @test !isnothing(m)
        pts = parse(Int, m[1])
        @test 57 ≤ pts ≤ 517
        @test output[11] == "Removing redundancies from the convex hull - fine, $pts initial candidates"
        @test startswith(output[12], "Found 57 extremal points of the convex hull in")
    end
end

@testset "\"Newton\" polytope (complex-valued)" begin
    DynamicPolynomials.@complex_polyvar z[1:2]

    @testnewtons @test Newton.halfpolytope(:complex, z[1] + conj(z[1]) + conj(z[1])*z[2]^2 + z[1]*conj(z[2])^2) ==
        monomial_vector([1, z[1], z[2]^2])
end

# see https://github.com/qubit-ulm/ChannelLoss
@testset "Loss-protecting quantum error correction in defined subbasis" begin
    function symmetricPolyMatrix(name::AbstractString, size::Integer)
        vars = Vector{DynamicPolynomials.Variable{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder},
                                                  Graded{LexOrder}}}(
            undef, (size * (size + 1)) >> 1
        )
        mat = Matrix{(polynomial_type)(eltype(vars),Float64)}(undef, size, size)
        idx = 1
        for j in 1:size # we need to access both row-wise and col-wise, so order doesn't matter
            for i in 1:j
                var = DynamicPolynomials.Variable(
                    "$name[$i, $j]", DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}
                )
                @inbounds mat[i, j] = mat[j, i] = vars[idx] = var
                idx += 1
            end
        end
        return mat, vars
    end
    function optProgFull(d::Integer, r::Integer, basis::AbstractMatrix{Bool}, pdist::Float64)
        s = size(basis, 1) -1
        @assert d >= 2 && s >= r >= 1 && 0 <= pdist <= 1
        ρ, ρVars = symmetricPolyMatrix("ρ", size(basis, 2))
        PT = polynomial_type(eltype(ρ), Float64)
        total_trace = zero(PT)
        total_overlap = zero(PT)
        conditions = PT[]
        psd_conditions = Matrix{PT}[]
        push!(psd_conditions, ρ)
        nslots = binomial(s, r)
        @inbounds for (combi, arrival) in enumerate(combinations(2:s+1, r))
            # ρRedBasis is constructed from the first basis element (Alice) + all the ones indexed in arrival.
            # The basis elements are identified by integers with r +1 bits. First (highest) is Alice, rest are received slots.
            ρRedBasisElems = Set{UInt}()
            for basisᵢ in eachcol(basis)
                i = UInt(basisᵢ[1])
                for arrivalᵢ in arrival
                    i <<= 1
                    basisᵢ[arrivalᵢ] && (i |= true)
                end
                push!(ρRedBasisElems, i)
            end
            ρRedBasis = sort!(collect(ρRedBasisElems))
            # ρRed is constructed by tracing over all lost ones
            lost = setdiff(2:s+1, arrival)
            ρRed = zeros(PT, length(ρRedBasis), length(ρRedBasis))
            for (jj, basisⱼ) in enumerate(eachcol(basis))
                j = UInt(basisⱼ[1])
                for arrivalⱼ in arrival
                    j <<= 1
                    basisⱼ[arrivalⱼ] && (j |= true)
                end
                jpos = findfirst(isequal(j), ρRedBasis)
                for (ii, basisᵢ) in enumerate(eachcol(basis))
                    @views basisⱼ[lost] == basisᵢ[lost] || continue
                    i = UInt(basisᵢ[1])
                    for arrivalᵢ in arrival
                        i <<= 1
                        basisᵢ[arrivalᵢ] && (i |= true)
                    end
                    ipos = findfirst(isequal(i), ρRedBasis)
                    ρRed[ipos, jpos] += ρ[ii, jj]
                end
            end
            aliceBit = one(UInt) << r
            receiveMask = aliceBit - one(UInt)
            # choiBasis needs all reduced basis elements (received part only) plus 0 and 1 for output as the highest bit (which
            # is the same as where the Alice bit was coded).
            choiBasis = sort!(collect(Set((x & receiveMask) | receiveOut
                                          for x in ρRedBasisElems for receiveOut in (zero(UInt), aliceBit))))
            choi, choiVars = symmetricPolyMatrix("choi" * string(combi), length(choiBasis))
            push!(psd_conditions, choi)

            # we must impose the tr_out choi ⪯ 𝟙 condition.
            choiPTBasis = sort!(collect(Set(x & receiveMask for x in choiBasis)))
            choiPT = Matrix{PT}(undef, length(choiPTBasis), length(choiPTBasis))
            choiPT .= LinearAlgebra.I(length(choiPTBasis))
            for pt in (zero(UInt), aliceBit)
                for (jj, basisⱼ) in enumerate(choiPTBasis)
                    j = findfirst(isequal(basisⱼ | pt), choiBasis)
                    isnothing(j) && continue
                    for (ii, basisᵢ) in enumerate(choiPTBasis)
                        i = findfirst(isequal(basisᵢ | pt), choiBasis)
                        isnothing(i) && continue
                        choiPT[ii, jj] -= choi[i, j]
                    end
                end
            end
            push!(psd_conditions, choiPT)

            # we are only interested in the overlap with |Φ⁺⟩ in the end, so we don't need to calculate the whole ρFinal. Only
            # the |00⟩⟨00|, |00⟩⟨11|, hc, |11⟩⟨11| components are required.
            # However, we also need the trace, meaning |01⟩⟨01| and |10⟩⟨10|.
            for (cj, choiBasisⱼ) in enumerate(choiBasis)
                j = findfirst(isequal(choiBasisⱼ), ρRedBasis) # this checks both that Alice's bit matches the output bit (from
                                                              # the overlap and that the input matches the state)
                if !isnothing(j)
                    for (ci, choiBasisᵢ) in enumerate(choiBasis)
                        i = findfirst(isequal(choiBasisᵢ), ρRedBasis)
                        isnothing(i) && continue
                        total_overlap += choi[ci, cj] * ρRed[i, j]
                        if choiBasisᵢ & aliceBit == choiBasisⱼ & aliceBit
                            total_trace += choi[ci, cj] * ρRed[i, j]
                        end
                    end
                end
                # non-matching outputs
                choiBasisⱼ ⊻= aliceBit
                j = findfirst(isequal(choiBasisⱼ), ρRedBasis)
                isnothing(j) && continue
                for (ci, choiBasisᵢ) in enumerate(choiBasis)
                    choiBasisᵢ ⊻= aliceBit
                    choiBasisᵢ & aliceBit == choiBasisⱼ & aliceBit || continue
                    i = findfirst(isequal(choiBasisᵢ), ρRedBasis)
                    isnothing(i) && continue
                    total_trace += choi[ci, cj] * ρRed[i, j]
                end
            end
        end
        total_overlap /= 2nslots
        total_trace /= nslots
        fidelity = total_overlap / pdist
        return poly_problem(-fidelity, zero=push!(conditions, total_trace - pdist, 1 - tr(ρ)),
            nonneg=PT[fidelity], # we just add this here so that all the methods are hit in the test
            psd=psd_conditions)
    end

    prob = optProgFull(2, 2, Bool[0 0 1 1; 0 1 0 1; 0 1 0 1; 1 0 0 1], .8)
    @testnewtons begin
        method === :COPT && continue # no license for the large problem
        rel, output = capture_stdout(FixAll(Relaxation.Newton, prob, 2; verbose=true, method))
        @test strRep(rel) == "Relaxation.Newton of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20], x[21], x[22], x[23], x[24], x[25], x[26], x[27], x[28], x[29], x[30], x[31], x[32], x[33], x[34], x[35], x[36], x[37], x[38], x[39], x[40], x[41], x[42], x[43], x[44], x[45], x[46], x[47], x[48], x[49], x[50], x[51], x[52], x[53], x[54], x[55], x[56], x[57], x[58], x[59], x[60], x[61], x[62], x[63], x[64], x[65], x[66], x[67], x[68], x[69], x[70], x[71], x[72], x[73], x[74], x[75], x[76], x[77], x[78], x[79], x[80], x[81], x[82], x[83], x[84], x[85], x[86], x[87], x[88], x[89], x[90], x[91], x[92]
PSD block sizes:
  [744 => 2, 372 => 4, 186 => 1, 123 => 1, 93 => 1]
Free block sizes:
  [93 => 2]"
        @test output[8] == "└ Total number of coefficients: 260432"
        @test startswith(output[12], "Found 122110 potential extremal points of the convex hull in")
        @test output[13] == "Starting point selection among 4279 possible monomials"
        @test startswith(output[end-2], "Found 123 elements in the Newton halfpolytope in")
    end

end