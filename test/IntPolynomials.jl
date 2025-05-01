# This testsuite uses some of MP's tests (and some more), but does so explicitly, since IntPolynomials on purpose does not
# support a lot of the things that you would want to do with polynomials.
using Test
using LinearAlgebra, SparseArrays

using MultivariatePolynomials
const MP = MultivariatePolynomials
import DynamicPolynomials, Combinatorics

using PolynomialOptimization.IntPolynomials
using PolynomialOptimization.IntPolynomials.MultivariateExponents

testdir = dirname(pathof(MP)) * "/../test"
include("$testdir/utils.jl")

@testset "Exponents" begin
    @testset failfast=true "ExponentsAll" begin
        ea = ExponentsAll{5,Int32}()
        @test stack(Iterators.take(Iterators.drop(ea, 5), 7)) == Int[
            1  0  0  0  0  0  0
            0  0  0  0  0  0  0
            0  0  0  0  1  1  2
            0  0  1  2  0  1  0
            0  2  1  0  1  0  0
        ]
        @test ea[Int32(8)] == [0, 0, 0, 1, 1]
        @test_throws BoundsError ea[Int32(0)]
        counts, success = index_counts(ea, 8)
        @test success
        # note: in a fresh Julia session, this should be the whole counts. However, if something has already run (or the test
        # suite is executed more than once), the previous runs will have populated the cache some more, so counts might have
        # more rows.
        @test counts[1:9, :] == Int32[   1   1   1  1 1 1
                                         6   5   4  3 2 1
                                        21  15  10  6 3 1
                                        56  35  20 10 4 1
                                       126  70  35 15 5 1
                                       252 126  56 21 6 1
                                       462 210  84 28 7 1
                                       792 330 120 36 8 1
                                      1287 495 165 45 9 1]
        alloc_test(let ea=ea; () -> index_counts(ea, 8) end, 0)

        ei = Vector{Int}(undef, 5)
        eiter = similar(ei)
        @test_throws DimensionMismatch iterate!(@view(eiter[1:3]), ea)
        copyto!(eiter, first(ea))
        for (i, eitervec, eitervec2) in zip(Int32(1):Int32(500), veciter(ea), veciter(ea, similar(ei)))
            alloc_test(let ea=ea, ei=ei, i=i; () -> copyto!(ei, exponents_from_index(ea, i)) end, 0)
            alloc_test(let ea=ea, ei=ei; () -> exponents_to_index(ea, ei) end, 0)
            @test exponents_to_index(ea, ei) === i
            @test eiter == ei
            @test iterate!(eiter, ea)
            @test eitervec == ei
            @test eitervec2 == ei
        end
        @test degree_from_index(ea, Int32(2000)) == 9
    end
    @testset failfast=true "ExponentsDegree" begin
        ed = ExponentsDegree{5,Int32}(2:7)
        @test length(ed) == 786
        @test length(unsafe, ed) == 786
        @test stack(Iterators.take(Iterators.drop(ed, 5), 7)) == Int[
            0  0  0  0  0  1  1
            0  1  1  1  2  0  0
            2  0  0  1  0  0  0
            0  0  1  0  0  0  1
            0  1  0  0  0  1  0
        ]
        @test ed[Int32(8)] == [0, 1, 0, 1, 0]
        @test_throws BoundsError ed[Int32(0)]
        @test_throws BoundsError ed[Int32(787)]
        # Checking index_counts does not make much sense, it is shared with the generic verison. Just ensure that no higher
        # degrees are created than what we had before (which was one higher than for 2000 -> 10)
        counts, success = index_counts(ed, 11)
        @test !success

        ei = Vector{Int}(undef, 5)
        copyto!(ei, exponents_from_index(ed, one(Int32)))
        @test ei == [0, 0, 0, 0, 2]
        eiter = similar(ei)
        @test_throws DimensionMismatch iterate!(@view(eiter[1:3]), ed)
        copyto!(eiter, first(ed))
        for (i, eitervec, eitervec2) in zip(Int32(1):Int32(786), veciter(ed), veciter(ed, similar(ei)))
            alloc_test(let ed=ed, ei=ei, i=i; () -> copyto!(ei, exponents_from_index(ed, i)) end, 0)
            alloc_test(let ed=ed, ei=ei; () -> exponents_to_index(ed, ei) end, 0)
            @test exponents_to_index(ed, ei) === i
            @test eiter == ei
            @test iterate!(eiter, ed) == (i != 786)
            @test eitervec == ei
            @test eitervec2 == ei
        end
        @test degree_from_index(ed, Int32(786)) == 7
        @test degree_from_index(ed, Int32(787)) == 8
    end
    @testset failfast=true "ExponentsMultideg" begin
        emd = ExponentsMultideg{5,Int32}(3:7, [1, 0, 0, 1, 0], [5, 3, 6, 7, 2])
        @test length(emd) == 223
        @test length(unsafe, emd) == 223
        @test stack(Iterators.take(Iterators.drop(emd, 5), 7)) == Int[
            1  1  1  1  1  1  1
            0  0  0  0  0  0  1
            0  0  0  1  1  2  0
            1  2  3  1  2  1  1
            2  1  0  1  0  0  1
        ]
        @test emd[Int32(8)] == [1, 0, 0, 3, 0]
        @test_throws BoundsError emd[Int32(0)]
        @test_throws BoundsError emd[Int32(224)]
        @test !index_counts(emd, 8)[2]
        counts, success = index_counts(emd, 7)
        @test success
        @test counts == Int32[  0    0   0   0  1  0
                                0    1   1   1  2  0
                                1    5   4   3  3  0
                                6   15  10   6  3  0
                               21   34  19   9  3  0
                               55   64  31  12  3  0
                              119  106  46  15  3  0
                              224  160  64  18  3  0]

        ei = Vector{Int}(undef, 5)
        copyto!(ei, exponents_from_index(emd, one(Int32)))
        @test ei == [1, 0, 0, 1, 1]
        eiter = similar(ei)
        @test_throws DimensionMismatch iterate!(@view(eiter[1:3]), emd)
        copyto!(eiter, first(emd))
        for (i, eitervec, eitervec2) in zip(Int32(1):Int32(223), veciter(emd), veciter(emd, similar(ei)))
            alloc_test(let emd=emd, ei=ei, i=i; () -> copyto!(ei, exponents_from_index(emd, i)) end, 0)
            alloc_test(let emd=emd, ei=ei; () -> exponents_to_index(emd, ei) end, 0)
            @test exponents_to_index(emd, ei) === i
            @test eiter == ei
            @test iterate!(eiter, emd) == (i != 223)
            @test eitervec == ei
            @test eitervec2 == ei
        end
        @test degree_from_index(emd, Int32(223)) == 7
        @test degree_from_index(emd, Int32(224)) == 8
    end
    @testset "Conversion" begin
        testexps(target, source, exps) =
            @test @inferred(convert_index(target, source, exponents_to_index(source, exps))) ===
                exponents_to_index(target, exps)

        target = ExponentsAll{4,Int32}()
        testexps(target, ExponentsAll{4,Int}(), [1, 2, 3, 4])
        testexps(target, ExponentsDegree{4,Int}(0, 7), [1, 0, 2, 3])
        testexps(target, ExponentsDegree{4,Int}(2, 7), [1, 0, 2, 3])
        testexps(target, ExponentsMultideg{4,Int}(0, 7, [0, 0, 0, 0], [3, 4, 5, 6]), [1, 0, 2, 0])
        testexps(target, ExponentsMultideg{4,Int}(0, 7, [0, 0, 0, 0], [3, 4, 5, 6]), [1, 0, 2, 3])
        testexps(target, ExponentsMultideg{4,Int}(0, 7, [1, 0, 1, 0], [3, 4, 5, 6]), [1, 0, 2, 0])
        testexps(target, ExponentsMultideg{4,Int}(2, 7, [0, 0, 0, 0], [3, 4, 5, 6]), [1, 0, 2, 0])

        for target_mindeg in (0, 2)
            target = ExponentsDegree{4,Int32}(target_mindeg, 6)
            testexps(target, ExponentsAll{4,Int}(), [1, 2, 1, 0])
            testexps(target, ExponentsAll{4,Int}(), [1, 2, 3, 4])
            testexps(target, ExponentsDegree{4,Int}(0, 7), [1, 0, 2, 3])
            testexps(target, ExponentsDegree{4,Int}(2, 7), [1, 0, 2, 3])
            testexps(target, ExponentsMultideg{4,Int}(0, 7, [0, 0, 0, 0], [3, 4, 5, 6]), [1, 0, 2, 0])
            testexps(target, ExponentsMultideg{4,Int}(0, 7, [0, 0, 0, 0], [3, 4, 5, 6]), [1, 0, 2, 4])
            testexps(target, ExponentsMultideg{4,Int}(0, 7, [1, 0, 1, 0], [3, 4, 5, 6]), [1, 0, 2, 0])
            testexps(target, ExponentsMultideg{4,Int}(2, 7, [0, 0, 0, 0], [3, 4, 5, 6]), [1, 0, 2, 0])
        end

        for target_mindeg in (0, 2)
            target = ExponentsMultideg{4,Int32}(target_mindeg, 6, [0, 0, 0, 0], [3, 4, 5, 6])
            testexps(target, ExponentsAll{4,Int}(), [1, 2, 1, 0])
            testexps(target, ExponentsAll{4,Int}(), [1, 2, 3, 4])
            testexps(target, ExponentsAll{4,Int}(), [1, 2, 3, 0])
            testexps(target, ExponentsDegree{4,Int}(0, 7), [1, 0, 2, 3])
            testexps(target, ExponentsDegree{4,Int}(2, 7), [1, 0, 2, 3])
            testexps(target, ExponentsMultideg{4,Int}(0, 7, [0, 0, 0, 0], [3, 4, 5, 6]), [1, 0, 2, 0])
            testexps(target, ExponentsMultideg{4,Int}(0, 7, [0, 0, 0, 0], [3, 4, 5, 6]), [1, 0, 2, 4])
            testexps(target, ExponentsMultideg{4,Int}(0, 7, [1, 0, 1, 0], [3, 4, 5, 6]), [1, 0, 2, 0])
            testexps(target, ExponentsMultideg{4,Int}(2, 7, [0, 0, 0, 0], [3, 4, 5, 6]), [1, 0, 2, 0])
        end

        for target_mindeg in (2, 4)
            target = ExponentsMultideg{4,Int32}(target_mindeg, 6, [1, 0, 1, 0], [3, 4, 5, 6])
            testexps(target, ExponentsAll{4,Int}(), [1, 2, 1, 0])
            testexps(target, ExponentsDegree{4,Int}(0, 7), [1, 0, 2, 3])
            for source_mindeg in (2, 3)
                testexps(target, ExponentsMultideg{4,Int}(source_mindeg, 7, [1, 0, 1, 0], [3, 4, 5, 6]), [1, 0, 2, 0])
                testexps(target, ExponentsMultideg{4,Int}(source_mindeg, 7, [1, 0, 1, 0], [3, 4, 5, 6]), [1, 0, 2, 4])
                testexps(target, ExponentsMultideg{4,Int}(source_mindeg, 7, [1, 0, 1, 0], [3, 5, 5, 6]), [1, 0, 2, 0])
                testexps(target, ExponentsMultideg{4,Int}(source_mindeg, 7, [1, 0, 1, 0], [3, 5, 5, 6]), [1, 0, 2, 1])
                testexps(target, ExponentsMultideg{4,Int}(source_mindeg, 7, [0, 0, 0, 0], [3, 5, 5, 6]), [1, 0, 2, 1])
            end
        end
    end
    @testset failfast=true "Comparison" begin
        ea = ExponentsAll{4,Int}()
        function testexps(e₁, exps₁, e₂, exps₂)
            ie₁ = exponents_to_index(e₁, exps₁)
            iszero(ie₁) && return
            ie₂ = exponents_to_index(e₂, exps₂)
            iszero(ie₂) && return
            ia₁ = exponents_to_index(ea, exps₁)
            ia₂ = exponents_to_index(ea, exps₂)
            for op in (==, <, ≤, >, ≤)
                @test @inferred(compare_indices(e₁, ie₁, op, e₂, ie₂)) === op(ia₁, ia₂)
            end
        end

        for (e₁, e₂) in Combinatorics.multiset_combinations(
            (ExponentsAll{4,Int}(), ExponentsDegree{4,Int}(0, 7), ExponentsDegree{4,Int}(2, 7),
            ExponentsMultideg{4,Int}(0, 7, [0, 0, 0, 0], [3, 4, 5, 6]),
            ExponentsMultideg{4,Int}(0, 7, [1, 0, 1, 0], [3, 4, 5, 6]),
            ExponentsMultideg{4,Int}(2, 7, [0, 0, 0, 0], [3, 4, 5, 6]),
            ExponentsMultideg{4,Int}(3, 7, [1, 0, 1, 0], [3, 5, 5, 6])), 2
        )
            for (i₁, i₂) in Iterators.product(Iterators.repeated(([1, 2, 3, 4], [1, 0, 2, 3], [1, 0, 2, 0],
                                                                  [1, 2, 1, 0], [1, 0, 2, 4],
                                                                  [1, 0, 2, 1]), 2)...)
                testexps(e₁, i₁, e₂, i₂)
            end
        end
    end
end

@testset "Variable" begin
    # No polyvar macro index set: we don't provide a macro for variable construction.
    # You shouldn't build up IntPolynomials from variables anyway!

    @test_throws DomainError IntRealVariable{1,0}(-1)
    @test_throws DomainError IntRealVariable{1,0}(0)
    @test_throws DomainError IntRealVariable{1,5}(2)
    @test_throws DomainError IntComplexVariable{0,1}(-1)
    @test_throws DomainError IntComplexVariable{0,1}(0)
    @test_throws DomainError IntComplexVariable{5,1}(2)
    for (t, i) in ((UInt8, Int8(1)), (UInt16, Int16(1)), (UInt32, Int32(1)), (UInt64, Int64(1)))
        x = IntRealVariable{typemax(t),0}(i)
        @test t(x) === t(i)
        @test convert(typeof(x), x) === x # these are a few of the only allowed converts
        alloc_test(() -> convert(typeof(x), x), 0)
        @test convert(variable_union_type(x), x) === x
        alloc_test(() -> convert(variable_union_type(x), x), 0)
        @test isreal(x) # test for short-circuit versions
        @test !isconj(x) # test for short-circuit versions
        @test ordinary_variable(x) === x
        z = IntComplexVariable{0,typemax(t)>>1}(i)
        @test t(z) === t(i)
        @test convert(typeof(z), z) === z
        alloc_test(() -> convert(typeof(z), z), 0)
        @test convert(variable_union_type(z), z) === z
        alloc_test(() -> convert(variable_union_type(z), z), 0)
    end

    x = IntRealVariable{1,1}(0x1)
    z = IntComplexVariable{1,1}(0x1)

    @test isreal(x)
    @test !isconj(x)
    @test ordinary_variable(x) === x

    @test !isreal(z)
    @test !isconj(z)
    @test ordinary_variable(z) === z

    zc = IntComplexVariable{1,1}(0x1, true)
    @test !isreal(zc)
    @test isconj(zc)
    @test ordinary_variable(zc) === z

    @test zc != z

    @test 1 != x
    @test x != 0
    @test nvariables(x) == 3 # different from MP - our ring is fixed
    @test !isapproxzero(x)
    @test !iszero(x)
    # no tests that involve zero(x): this would involve a conversion

    typetests(x)
    # no tests with polynomial(x): don't construct them in this way

    @test nterms(x) == 1
    # no tests for terms(x): we cannot convert, as we don't know whether to a dense or sparse monomial and we also don't know
    # the exponent type.

    @test degree(x, x) == 1
    @test degree(x, z) == 0
    @test length(exponents(x)) == 3
    @test first(exponents(x)) == 1
    @test isconstant(x) == false

    # no Issue #82, we cannot convert

    @testset "Effective variables" begin
        @test [x] == @inferred effective_variables(x)
        @test [z] == @inferred effective_variables(z)
    end

    # no creation of similar variables, this reeks of bottom-up construction
end

@testset "Monomial (real)" begin
    # lots of tests to skip, we don't construct monomials by multiplying variables

    @test_throws ArgumentError IntMonomial{2,0}([1])
    m = IntMonomial{7,0}([1, 0, 1, 0, 1, 0, 1])
    alloc_test(let v=[1, 0, 1, 0, 1, 0, 1]; () -> IntMonomial{7,0}(v) end, 0)
    @test m === IntMonomial{7,0}([1, 0, 1, 0, 1, 0, 1])
    alloc_test(() -> convert(typeof(m), m), 0)
    @test_throws ArgumentError IntMonomial{1,0}([1, 2])

    @test nvariables(IntMonomialVector{4,0}([2 0; 0 1; 0 1; 0 1])) == 4

    @test nterms(IntMonomial{1,0}([2])) == 1
    @test collect(m) == [(IntRealVariable{7,0}(1), 1), (IntRealVariable{7,0}(3), 1),
                         (IntRealVariable{7,0}(5), 1), (IntRealVariable{7,0}(7), 1)]
    @test collect(IntMonomial{7,0}([0, 1, 0, 2, 0, 3, 0])) ==
        [(IntRealVariable{7,0}(2), 1), (IntRealVariable{7,0}(4), 2), (IntRealVariable{7,0}(6), 3)]

    @test degree(IntMonomial{3,0}([1, 0, 2]), IntRealVariable{3,0}(1)) == UInt(1)
    @test degree(IntMonomial{3,0}([1, 0, 2]), IntRealVariable{3,0}(2)) == UInt(0)
    @test degree(IntMonomial{3,0}([1, 0, 2]), IntRealVariable{3,0}(3)) == UInt(2)

    @test_throws InexactError variable(IntMonomial{1,0}([2]))
    @test_throws InexactError variable(IntMonomial{2,0}([1, 1]))
    @test_throws InexactError variable(constant_monomial(typeof(IntRealVariable{1,0}(1))))

    x = IntRealVariable{2,0}(1)
    @test x != constant_monomial(typeof(x))
    @test constant_monomial(typeof(x)) != x

    # no tests for squaring monomials

    typetests(m)

    xmon = IntMonomial{2,0}([1, 0])
    @test variable(xmon) === x
    # no arithmetics
    @test variable(Term(1., xmon)) === x
    @test_throws InexactError variable(Term(3, xmon)) === x

    @test transpose(x) === x
    @test adjoint(x) === x
    @test transpose(m) === m
    @test adjoint(m) === m
    @test isreal(m)

    @testset "Effective variables" begin
        x = IntRealVariable{8,0}(1)
        y = ntuple(i -> IntRealVariable{8,0}(i+1), Val(7))
        T = variable_union_type(x)
        @test x isa T
        @test y[2] isa T
        @test T[x, y[2]] == @inferred effective_variables(IntMonomial{8,0}([1, 0, 1, 0, 0, 0, 0, 0]))
        @test T[x] == @inferred effective_variables(      IntMonomial{8,0}([1, 0, 0, 0, 0, 0, 0, 0]))
        @test T[y[2]] == @inferred effective_variables(   IntMonomial{8,0}([0, 0, 1, 0, 0, 0, 0, 0]))
    end

    # no mapexponents
end

@testset "Monomial (complex)" begin
    # lots of tests to skip, we don't construct monomials by multiplying variables

    @test_throws ArgumentError IntMonomial{0,1}([1, 0], [0])
    @test_throws ArgumentError IntMonomial{0,2}([1, 0], [0])
    m = IntMonomial{0,7}([1, 0, 0, 0, 1, 0, 1], [0, 0, 1, 0, 0, 0, 0])
    alloc_test(let v1=[1, 0, 0, 0, 1, 0, 1], v2=[0, 0, 1, 0, 0, 0, 0]; () -> IntMonomial{0,7}(v1, v2) end, 0)
    @test m == IntMonomial{0,7}([1, 0, 0, 0, 1, 0, 1], [0, 0, 1, 0, 0, 0, 0])
    alloc_test(() -> convert(typeof(m), m), 0)
    @test_throws ArgumentError IntMonomial{0,1}([1, 2], [0, 0])
    @test_throws ArgumentError IntMonomial{0,2}([1, 2], [0])
    @test_throws ArgumentError IntMonomial{0,2}([0], [1, 2])

    @test nvariables(IntMonomialVector{0,4}([2 0; 0 1; 0 0; 0 1], [0 0; 0 0; 0 1; 0 0])) == 8

    @test nterms(IntMonomial{0,1}([2], [0])) == 1
    @test collect(m) == [(IntComplexVariable{0,7}(1), 1), (IntComplexVariable{0,7}(3, true), 1),
                         (IntComplexVariable{0,7}(5), 1), (IntComplexVariable{0,7}(7), 1)]
    @test collect(IntConjMonomial(m)) == [(IntComplexVariable{0,7}(1, true), 1), (IntComplexVariable{0,7}(3), 1),
                         (IntComplexVariable{0,7}(5, true), 1), (IntComplexVariable{0,7}(7, true), 1)]
    @test collect(IntMonomial{0,7}([0, 1, 0, 0, 0, 3, 0], [0, 0, 0, 2, 0, 0, 0])) ==
                        [(IntComplexVariable{0,7}(2), 1), (IntComplexVariable{0,7}(4, true), 2),
                         (IntComplexVariable{0,7}(6), 3)]

    @test degree(IntMonomial{0,3}([1, 0, 2], [0, 1, 3]), IntComplexVariable{0,3}(1)) === 1
    @test degree(IntMonomial{0,3}([1, 0, 2], [0, 1, 3]), IntComplexVariable{0,3}(1, true)) === 0
    @test degree(IntMonomial{0,3}([1, 0, 2], [0, 1, 3]), IntComplexVariable{0,3}(2)) === 0
    @test degree(IntMonomial{0,3}([1, 0, 2], [0, 1, 3]), IntComplexVariable{0,3}(2, true)) === 1
    @test degree(IntMonomial{0,3}([1, 0, 2], [0, 1, 3]), IntComplexVariable{0,3}(3)) === 2
    @test degree(IntMonomial{0,3}([1, 0, 2], [0, 1, 3]), IntComplexVariable{0,3}(3, true)) === 3
    @test degree_complex(IntMonomial{0,3}([1, 0, 2], [0, 1, 3]), IntComplexVariable{0,3}(1)) === 1
    @test degree_complex(IntMonomial{0,3}([1, 0, 2], [0, 1, 3]), IntComplexVariable{0,3}(2)) === 1
    @test degree_complex(IntMonomial{0,3}([1, 0, 2], [0, 1, 3]), IntComplexVariable{0,3}(3)) === 3

    @test_throws InexactError variable(IntMonomial{0,1}([2], [0]))
    @test_throws InexactError variable(IntMonomial{0,2}([1, 0], [0, 1]))
    @test_throws InexactError variable(constant_monomial(typeof(IntComplexVariable{0,1}(1))))

    z = IntComplexVariable{0,2}(1)
    @test z != constant_monomial(typeof(z))
    @test constant_monomial(typeof(z)) != z

    # no tests for squaring monomials

    typetests(m)

    zmon = IntMonomial{0,2}([1, 0], [0, 0])
    @test variable(zmon) === z
    # no arithmetics
    @test variable(Term(1., zmon)) === z
    @test_throws InexactError variable(Term(3, zmon)) === z

    @test transpose(z) === z
    @test adjoint(z) === conj(z)
    @test transpose(m) === m
    @test adjoint(m) === conj(m) == IntMonomial{0,7}(UInt8[0, 0, 1, 0, 0, 0, 0], UInt8[1, 0, 0, 0, 1, 0, 1])
    @test !isreal(m)
    @test isreal(IntMonomial{0,7}(UInt8[0, 1, 2, 3, 2, 1, 0], UInt8[0, 1, 2, 3, 2, 1, 0]))

    @testset "Effective variables" begin
        y = ntuple(i -> IntComplexVariable{0,8}(i, true), Val(7))
        z = IntComplexVariable{0,8}(8)
        T = variable_union_type(z)
        @test z isa T
        @test y[2] isa T
        @test T[y[2], z] == @inferred effective_variables(IntMonomial{0,8}([0, 0, 0, 0, 0, 0, 0, 1],
                                                                              [0, 1, 0, 0, 0, 0, 0, 0]))
        @test T[z] == @inferred effective_variables(      IntMonomial{0,8}([0, 0, 0, 0, 0, 0, 0, 1],
                                                                              [0, 0, 0, 0, 0, 0, 0, 0]))
        @test T[y[2]] == @inferred effective_variables(   IntMonomial{0,8}([0, 0, 0, 0, 0, 0, 0, 0],
                                                                              [0, 1, 0, 0, 0, 0, 0, 0]))
    end

    # no mapexponents
end

@testset "Monomial (mixed)" begin
    # lots of tests to skip, we don't construct monomials by multiplying variables

    @test_throws ArgumentError IntMonomial{2,2}([1], [1, 2], [0, 1])
    @test_throws ArgumentError IntMonomial{2,2}([1, 0], [1], [0, 1])
    @test_throws ArgumentError IntMonomial{2,2}([1, 0], [1, 0], [0])
    @test_throws ArgumentError IntMonomial{2,2}([1, 0], [1], [0])
    m = IntMonomial{2,7}(UInt8[2, 3], UInt8[1, 0, 0, 0, 1, 0, 1], UInt8[0, 0, 1, 0, 0, 0, 0])
    alloc_test(let v1=UInt8[2, 3], v2=UInt8[1, 0, 0, 0, 1, 0, 1], v3=UInt8[0, 0, 1, 0, 0, 0, 0]
        () -> IntMonomial{2,7}(v1, v2, v3)
    end, 0)
    @test m == IntMonomial{2,7}(Int8[2, 3], Int8[1, 0, 0, 0, 1, 0, 1], Int8[0, 0, 1, 0, 0, 0, 0])
    alloc_test(() -> convert(typeof(m), m), 0)
    @test_throws ArgumentError IntMonomial{1,1}([1, 2], [1], [1])
    @test_throws ArgumentError IntMonomial{1,1}([1], [1, 2], [0, 0])
    @test_throws ArgumentError IntMonomial{1,2}([1], [1, 2], [0])
    @test_throws ArgumentError IntMonomial{1,2}([1], [0], [1, 2])

    @test nvariables(IntMonomialVector{3,4}([0 2; 0 0; 1 1], [0 2; 1 0; 0 0; 1 0], [0 0; 0 0; 1 0; 0 0])) == 11

    @test nterms(IntMonomial{1,1}([3], [2], [0])) == 1
    @test collect(m) == [(IntRealVariable{2,7}(1), 2), (IntRealVariable{2,7}(2), 3),
                         (IntComplexVariable{2,7}(1), 1), (IntComplexVariable{2,7}(3, true), 1),
                         (IntComplexVariable{2,7}(5), 1), (IntComplexVariable{2,7}(7), 1)]
    @test collect(IntConjMonomial(m)) == [(IntRealVariable{2,7}(1), 2), (IntRealVariable{2,7}(2), 3),
                         (IntComplexVariable{2,7}(1, true), 1), (IntComplexVariable{2,7}(3), 1),
                         (IntComplexVariable{2,7}(5, true), 1), (IntComplexVariable{2,7}(7, true), 1)]
    @test collect(IntMonomial{2,7}([0, 4], [0, 1, 0, 0, 0, 3, 0], [0, 0, 0, 2, 0, 0, 0])) ==
                         [(IntRealVariable{2,7}(2), 4),
                          (IntComplexVariable{2,7}(2), 1), (IntComplexVariable{2,7}(4, true), 2),
                          (IntComplexVariable{2,7}(6), 3)]

    @test degree(IntMonomial{2,3}([2, 1], [1, 0, 2], [0, 1, 3]), IntRealVariable{2,3}(1)) === 2
    @test degree(IntMonomial{2,3}([2, 1], [1, 0, 2], [0, 1, 3]), IntRealVariable{2,3}(2)) === 1
    @test degree(IntMonomial{2,3}([2, 1], [1, 0, 2], [0, 1, 3]), IntComplexVariable{2,3}(1)) === 1
    @test degree(IntMonomial{2,3}([2, 1], [1, 0, 2], [0, 1, 3]), IntComplexVariable{2,3}(1, true)) === 0
    @test degree(IntMonomial{2,3}([2, 1], [1, 0, 2], [0, 1, 3]), IntComplexVariable{2,3}(2)) === 0
    @test degree(IntMonomial{2,3}([2, 1], [1, 0, 2], [0, 1, 3]), IntComplexVariable{2,3}(2, true)) === 1
    @test degree(IntMonomial{2,3}([2, 1], [1, 0, 2], [0, 1, 3]), IntComplexVariable{2,3}(3)) === 2
    @test degree(IntMonomial{2,3}([2, 1], [1, 0, 2], [0, 1, 3]), IntComplexVariable{2,3}(3, true)) === 3
    @test degree_complex(IntMonomial{2,3}([2, 1], [1, 0, 2], [0, 1, 3]), IntRealVariable{2,3}(1)) === 2
    @test degree_complex(IntMonomial{2,3}([2, 1], [1, 0, 2], [0, 1, 3]), IntRealVariable{2,3}(2)) === 1
    @test degree_complex(IntMonomial{2,3}([2, 1], [1, 0, 2], [0, 1, 3]), IntComplexVariable{2,3}(1)) === 1
    @test degree_complex(IntMonomial{2,3}([2, 1], [1, 0, 2], [0, 1, 3]), IntComplexVariable{2,3}(2)) === 1
    @test degree_complex(IntMonomial{2,3}([2, 1], [1, 0, 2], [0, 1, 3]), IntComplexVariable{2,3}(3)) === 3

    @test_throws InexactError variable(IntMonomial{1,1}([3], [2], [0]))
    @test_throws InexactError variable(IntMonomial{1,2}([1], [1, 0], [0, 1]))
    @test_throws InexactError variable(constant_monomial(typeof(IntComplexVariable{1,1}(1))))

    # no constant monomial tests
    x = IntRealVariable{3,2}(1)
    @test x != constant_monomial(typeof(x))
    @test constant_monomial(typeof(x)) != x
    z = IntComplexVariable{3,2}(1)
    @test z != constant_monomial(typeof(z))
    @test constant_monomial(typeof(z)) != z
    @test x != z

    # no tests for squaring monomials

    typetests(m)

    xmon = IntMonomial{3,2}([1, 0, 0], [0, 0], [0, 0])
    zmon = IntMonomial{3,2}([0, 0, 0], [1, 0], [0, 0])
    @test variable(xmon) === x
    @test variable(zmon) === z
    # no arithmetics
    @test variable(Term(1., xmon)) === x
    @test variable(Term(1., zmon)) === z
    @test_throws InexactError variable(Term(3, xmon)) === z
    @test_throws InexactError variable(Term(3, zmon)) === z

    @test transpose(x) === x
    @test transpose(z) === z
    @test adjoint(x) === x
    @test adjoint(z) === conj(z)
    @test transpose(m) === m
    @test adjoint(m) === conj(m) == IntMonomial{2,7}(UInt8[2, 3], UInt8[0, 0, 1, 0, 0, 0, 0], UInt8[1, 0, 0, 0, 1, 0, 1])
    @test !isreal(m)
    @test isreal(IntMonomial{2,7}(UInt8[2, 3], UInt8[0, 1, 2, 3, 2, 1, 0], UInt8[0, 1, 2, 3, 2, 1, 0]))

    @testset "Effective variables" begin
        x = IntRealVariable{3,5}(2)
        y = IntComplexVariable{3,5}(1, true)
        z = IntComplexVariable{3,5}(3)
        T = IntVariable{3,5}
        @test z isa T
        @test y isa T
        @test T[x, y, z] == @inferred effective_variables(IntMonomial{3,5}([0, 2, 0], [0, 0, 1, 0, 0],
                                                                                         [1, 0, 0, 0, 0]))
        @test T[x] == @inferred effective_variables(      IntMonomial{3,5}([0, 2, 0], [0, 0, 0, 0, 0],
                                                                                         [0, 0, 0, 0, 0]))
        @test T[y] == @inferred effective_variables(      IntMonomial{3,5}([0, 0, 0], [0, 0, 0, 0, 0],
                                                                                         [1, 0, 0, 0, 0]))
        @test T[z] == @inferred effective_variables(      IntMonomial{3,5}([0, 0, 0], [0, 0, 1, 0, 0],
                                                                                         [0, 0, 0, 0, 0]))
    end

    # no mapexponents
end

@testset "Term" begin
    # as our implementation uses the default Term, there's not much to test
    t = Term(3, IntMonomial{2,0}([2, 4]))
    alloc_test(() -> convert(typeof(t), t), 0)
    typetests(t)
    typetests([t, Term(2, IntMonomial{2,0}([2, 0]))])
end

@testset "Monomial vector" begin
    # by default, our monomial vector does not impose any order
    @test_throws ArgumentError monomials(2, 0, 0x1:0x0)
    X = [IntMonomial{2,0}(UInt8[0, 2]),
         IntMonomial{2,0}(UInt8[1, 1]),
         IntMonomial{2,0}(UInt8[2, 0])]
    monos = monomials(2, 0, 0x2:0x2)
    @test length(monos) == length(X)
    @test monos ⊆ monos
    for (x, m) in zip(X, monos)
        @test m == x
        @test [x] ⊆ monos
    end
    m2 = monomials(2, 0, 0x1:0x2, filter_exps=e -> sum(e) == 2)
    @test m2 == monos
    @test m2 ⊆ monos
    @test monos ⊆ m2
    @test monomials(2, 0, 0x1:0x2, filter_mons=m -> degree(m) == 2) == monos
    @test monomials(2, 0, 0x1:0x3, filter_exps=e -> sum(e) > 1, filter_mons=m -> degree(m) < 3) == monos
    # we don't provide monomial_vector_type

    X = IntMonomialVector{2,0}(UInt8[1 0 1; 0 0 1])
    @test !(monos ⊆ X)
    @test !(X ⊆ monos)
    @test X == collect(X)
    @test nvariables(X) == 2
    @test variables(X)[1] == IntRealVariable{2,0}(1)
    @test variables(X)[2] == IntRealVariable{2,0}(2)
    @test X[2:3][1] == IntMonomial{2,0}([0x1, 0x0])
    @test X[2:3][2] == IntMonomial{2,0}([0x1, 0x1])
    @test X[3:3] ⊆ monos
    @test !(monos ⊆ X[3:3])

    _checkindex(::IntMonomialVector{<:Any,<:Any,<:Integer,<:Tuple}, indexed::Bool) = @test indexed
    _checkindex(::IntMonomialVector, indexed::Bool) = @test !indexed
    function req_same(out, ref, indexed)
        @test out == ref
        _checkindex(out, indexed)
    end

    @testset "merge_monomial_vectors" begin
        m1 = IntMonomialVector{2,0}([1; 1;; 1; 0])
        m2 = IntMonomialVector{2,0}([2; 1;; 1; 0])
        req_same(@inferred(merge_monomial_vectors([m1, m2])), IntMonomialVector{2,0}([1; 0;; 1; 1;; 2; 1]), true)
        req_same(merge_monomial_vectors(Any[m1, m2]), IntMonomialVector{2,0}([1; 0;; 1; 1;; 2; 1]), true)
        req_same(merge_monomial_vectors([m1, m2, monos]), IntMonomialVector{2,0}([1; 0;; 0; 2;; 1; 1;; 2; 0;; 2; 1]), true)
        req_same(merge_monomial_vectors([m1, monos, m2]), IntMonomialVector{2,0}([1; 0;; 0; 2;; 1; 1;; 2; 0;; 2; 1]), true)
        req_same(merge_monomial_vectors([monos, m1, m2]), IntMonomialVector{2,0}([1; 0;; 0; 2;; 1; 1;; 2; 0;; 2; 1]), true)
        req_same(merge_monomial_vectors([monomials(2, 0, 0x1:0x1), monomials(2, 0, 0x2:0x2)]), monomials(2, 0, 0x1:0x2), false)
        req_same(merge_monomial_vectors([monomials(2, 0, 0x1:0x1), monomials(2, 0, 0x1:0x2)]), monomials(2, 0, 0x1:0x2), false)
        req_same(merge_monomial_vectors([monomials(2, 0, 0x1:0x2), monomials(2, 0, 0x2:0x2)]), monomials(2, 0, 0x1:0x2), false)
        req_same(merge_monomial_vectors([monomials(2, 0, 0x1:0x1), monomials(2, 0, 0x2:0x2)]), monomials(2, 0, 0x1:0x2), false)
        req_same(merge_monomial_vectors([monomials(2, 0, 0x1:0x3), monomials(2, 0, 0x2:0x2)]), monomials(2, 0, 0x1:0x3), false)
        req_same(merge_monomial_vectors([monomials(2, 0, 0x1:0x3), monomials(2, 0, 0x2:0x2, minmultideg=[1, 0])]),
            monomials(2, 0, 0x1:0x3), false)
        req_same(merge_monomial_vectors([monomials(2, 0, 0x1:0x3, minmultideg=[2, 0], maxmultideg=[3, 1]),
            monomials(2, 0, 0x1:0x3, minmultideg=[2, 0], maxmultideg=[3, 1])]),
            monomials(2, 0, 0x1:0x3, minmultideg=[2, 0], maxmultideg=[3, 1]), false)
        req_same(merge_monomial_vectors([monomials(2, 0, 0x2:0x3, maxmultideg=[2, 2]),
            monomials(2, 0, 0x2:0x5, maxmultideg=[3, 4])]), monomials(2, 0, 0x2:0x5, maxmultideg=[3, 4]), false)
        req_same(merge_monomial_vectors([monomials(2, 0, 0x2:0x5, maxmultideg=[3, 4]),
            monomials(2, 0, 0x2:0x3, maxmultideg=[2, 2])]), monomials(2, 0, 0x2:0x5, maxmultideg=[3, 4]), false)
        req_same(merge_monomial_vectors([monomials(2, 0, 0x1:0x3), monomials(2, 0, 0x2:0x5, maxmultideg=[3, 4])]),
            monomials(2, 0, 0x1:0x5, maxmultideg=[3, 4]), false)
        req_same(merge_monomial_vectors([monomials(2, 0, 0x1:0x3, maxmultideg=[2, 2]),
            monomials(2, 0, 0x2:0x5, maxmultideg=[3, 4])]), monomials(2, 0, 0x1:0x5, maxmultideg=[3, 4]), false)
        req_same(merge_monomial_vectors([monomials(2, 0, 0x2:0x3, maxmultideg=[2, 2]),
            monomials(2, 0, 0x2:0x4, minmultideg=[2, 0], maxmultideg=[2, 2])]),
            IntMonomialVector{2,0}([0; 2;; 1; 1;; 2; 0;; 1; 2;; 2; 1;; 2; 2]), true)
    end

    @test monomials(1, 0, 0x1:0x3) == IntMonomialVector{1,0}(UInt8[1 2 3])

    @testset "monomials" begin
        @test monomials(3, 0, 0x0:0x3) == IntMonomialVector{3,0}(
            copy(transpose(UInt8[
                0 0 0
                0 0 1
                0 1 0
                1 1 1
                1 2 0
                2 0 1
                0 2 1
                0 3 0
                2 0 0
                0 2 0
                1 0 1
                0 0 3
                1 0 0
                0 0 2
                0 1 1
                1 1 0
                0 1 2
                1 0 2
                2 1 0
                3 0 0
            ]))
        )
    end

    @test_throws ArgumentError monomials(1, 0, -1:0)

    # new tests
    @testset "intersect" begin
        m1 = IntMonomialVector{2,0}([1; 1;; 1; 0])
        m2 = IntMonomialVector{2,0}([2; 1;; 1; 0])
        empty = IntMonomialVector{2,0}(zeros(Int, 2, 0))
        req_same(@inferred(intersect(m1, m2)), IntMonomialVector{2,0}([1; 0;;]), true)
        req_same(@inferred(intersect(m1, monos)), IntMonomialVector{2,0}([1; 1;;]), true)
        req_same(@inferred(intersect(monos, m1)), IntMonomialVector{2,0}([1; 1;;]), true)
        req_same(@inferred(intersect(m2, monos)), empty, true)
        req_same(@inferred(intersect(monos, m2)), empty, true)
        req_same(intersect(monos, monos), monos, false)
        req_same(intersect(monomials(2, 0, 0x1:0x1), monomials(2, 0, 0x2:0x2)), empty, true)
        req_same(intersect(monomials(2, 0, 0x1:0x1), monomials(2, 0, 0x1:0x2)), monomials(2, 0, 0x1:0x1), false)
        req_same(intersect(monomials(2, 0, 0x1:0x2), monomials(2, 0, 0x1:0x1)), monomials(2, 0, 0x1:0x1), false)
        req_same(intersect(monomials(2, 0, 0x1:0x2), monomials(2, 0, 0x2:0x2)), monomials(2, 0, 0x2:0x2), false)
        req_same(intersect(monomials(2, 0, 0x2:0x2), monomials(2, 0, 0x1:0x2)), monomials(2, 0, 0x2:0x2), false)
        req_same(intersect(monomials(2, 0, 0x1:0x3), monomials(2, 0, 0x2:0x2)), monomials(2, 0, 0x2:0x2), false)
        req_same(intersect(monomials(2, 0, 0x1:0x3), monomials(2, 0, 0x2:0x2, minmultideg=[1, 0])),
            monomials(2, 0, 0x2:0x2, minmultideg=[1, 0]), false)
        req_same(intersect(monomials(2, 0, 0x1:0x3, minmultideg=[2, 0], maxmultideg=[3, 1]),
            monomials(2, 0, 0x1:0x3, minmultideg=[2, 0], maxmultideg=[3, 1])),
            monomials(2, 0, 0x1:0x3, minmultideg=[2, 0], maxmultideg=[3, 1]), false)
        req_same(intersect(monomials(2, 0, 0x2:0x3, maxmultideg=ConstantVector(2, 2)),
            monomials(2, 0, 0x2:0x5, maxmultideg=[3, 4])), monomials(2, 0, 0x2:0x3, maxmultideg=[2, 2]), false)
        req_same(intersect(monomials(2, 0, 0x2:0x5, maxmultideg=[3, 4]),
            monomials(2, 0, 0x2:0x3, maxmultideg=ConstantVector(2, 2))), monomials(2, 0, 0x2:0x3, maxmultideg=[2, 2]), false)
        req_same(intersect(monomials(2, 0, 0x1:0x3), monomials(2, 0, 0x2:0x5, maxmultideg=[3, 4])),
            monomials(2, 0, 0x2:0x3), false)
        req_same(intersect(monomials(2, 0, 0x1:0x3, maxmultideg=[2, 2]),
            monomials(2, 0, 0x2:0x5, maxmultideg=[3, 4])), monomials(2, 0, 0x2:0x3, maxmultideg=[2, 2]), false)
        req_same(intersect(monomials(2, 0, 0x2:0x3, maxmultideg=[2, 2]),
            monomials(2, 0, 0x2:0x4, minmultideg=[2, 0], maxmultideg=[2, 2])),
            monomials(2, 0, 0x2:0x3, minmultideg=[2, 0], maxmultideg=[2, 2]), false)
        m3 = IntMonomialVector{2,0}(ExponentsMultideg{2,UInt}(0, 3, [0, 0], [1, 1]))
        req_same(intersect(m1, m3), m1, true)
        req_same(intersect(m2, m3), IntMonomialVector{2,0}([1; 0;;]), true)
    end

    @testset failfast=true "effective_variables" begin
        function evt(m::IntMonomialVector, v)
            @test effective_nvariables(m) == length(v)
            @test collect(effective_variables(m)) == v
        end
        evt(monomials(2, 3, 0x0:0x2), IntVariable{2,3}.(1:8))
        @test effective_nvariables(monomials(2, 3, 0x0:0x1), monomials(2, 3, 0:2)) == 8
        evt(monomials(2, 3, 0x0:0x0), IntVariable{2,3}[])
        evt(IntMonomialVector{2,0}([0; 0;; 1; 0]), IntVariable{2,0}.([1]))
        evt(IntMonomialVector{2,0}([1; 0;; 2; 0]), IntVariable{2,0}.([1]))
        evt(IntMonomialVector{2,0}([0; 1;; 1; 0]), IntVariable{2,0}.([1, 2]))
        @test effective_nvariables(IntMonomialVector{2,0}([0; 0;; 1; 0]), IntMonomialVector{2,0}([0; 0;; 1; 0])) == 1
        @test effective_nvariables(IntMonomialVector{2,0}([0; 0;; 1; 0]), IntMonomialVector{2,0}([0; 0;; 1; 1])) == 2

        evt(monomials(3, 0, 0x2:0x5, minmultideg=[1, 0, 2], maxmultideg=[7, 4, 3]), IntVariable{3,0}.(1:3))
        evt(monomials(3, 0, 0x2:0x5, minmultideg=[1, 0, 4], maxmultideg=[7, 4, 8]), IntVariable{3,0}.([1, 3]))
        mons = monomials(4, 0, 0:5, minmultideg=[2, 0, 0, 0], maxmultideg=[4, 5, 6, 7])
        for rstart in 1:length(mons), rend in rstart:length(mons)
            sub = @view(mons[rstart:rend])
            correct_efvars = Set{IntVariable{4,0}}()
            for mon in sub, (var, pow) in mon
                push!(correct_efvars, var)
            end
            evt(sub, sort!(collect(correct_efvars)))
        end
        for subs in Combinatorics.powerset(1:length(mons), 0, 4)
            sub = @view(mons[subs])
            correct_efvars = Set{IntVariable{4,0}}()
            for mon in sub, (var, pow) in mon
                push!(correct_efvars, var)
            end
            evt(sub, sort!(collect(correct_efvars)))
        end
    end

    @testset "iteration" begin
        mv = IntMonomialVector{2,0}([2, 3, 5, 8, 13])
        e = mv.e
        @testset "standard (indexed)" begin
            vi = veciter(mv)
            it, state = iterate(vi)
            @test it == exponents(IntMonomial{2,0}(unsafe, e, 2))
            it, state = iterate(vi, state)
            @test it == exponents(IntMonomial{2,0}(unsafe, e, 3))
            it, state = iterate(vi, state)
            @test it == exponents(IntMonomial{2,0}(unsafe, e, 5))
            it, state = iterate(vi, state)
            @test it == exponents(IntMonomial{2,0}(unsafe, e, 8))
            it, state = iterate(vi, state)
            @test it == exponents(IntMonomial{2,0}(unsafe, e, 13))
            @test isnothing(iterate(vi, state))
        end
        @testset "enumerated (indexed)" begin
            vi = veciter(mv, true)
            it, state = iterate(vi)
            @test it == (2, exponents(IntMonomial{2,0}(unsafe, e, 2)))
            it, state = iterate(vi, state)
            @test it == (3, exponents(IntMonomial{2,0}(unsafe, e, 3)))
            it, state = iterate(vi, state)
            @test it == (5, exponents(IntMonomial{2,0}(unsafe, e, 5)))
            it, state = iterate(vi, state)
            @test it == (8, exponents(IntMonomial{2,0}(unsafe, e, 8)))
            it, state = iterate(vi, state)
            @test it == (13, exponents(IntMonomial{2,0}(unsafe, e, 13)))
            @test isnothing(iterate(vi, state))
        end

        mv = IntMonomialVector{2,0}([2, 3, 5, 8, 9])
        @testset "standard (iterated)" begin
            vi = veciter(mv)
            it, state = iterate(vi)
            @test it == exponents(IntMonomial{2,0}(unsafe, e, 2))
            it, state = iterate(vi, state)
            @test it == exponents(IntMonomial{2,0}(unsafe, e, 3))
            it, state = iterate(vi, state)
            @test it == exponents(IntMonomial{2,0}(unsafe, e, 5))
            it, state = iterate(vi, state)
            @test it == exponents(IntMonomial{2,0}(unsafe, e, 8))
            it, state = iterate(vi, state)
            @test it == exponents(IntMonomial{2,0}(unsafe, e, 9))
            @test isnothing(iterate(vi, state))
        end
        @testset "enumerated (iterated)" begin
            vi = veciter(mv, true)
            it, state = iterate(vi)
            @test it == (2, exponents(IntMonomial{2,0}(unsafe, e, 2)))
            it, state = iterate(vi, state)
            @test it == (3, exponents(IntMonomial{2,0}(unsafe, e, 3)))
            it, state = iterate(vi, state)
            @test it == (5, exponents(IntMonomial{2,0}(unsafe, e, 5)))
            it, state = iterate(vi, state)
            @test it == (8, exponents(IntMonomial{2,0}(unsafe, e, 8)))
            it, state = iterate(vi, state)
            @test it == (9, exponents(IntMonomial{2,0}(unsafe, e, 9)))
            @test isnothing(iterate(vi, state))
        end

        mv = IntMonomialVector{2,0}(2:6)
        @testset "standard (unit range)" begin
            vi = veciter(mv)
            it, state = iterate(vi)
            @test it == exponents(IntMonomial{2,0}(unsafe, e, 2))
            it, state = iterate(vi, state)
            @test it == exponents(IntMonomial{2,0}(unsafe, e, 3))
            it, state = iterate(vi, state)
            @test it == exponents(IntMonomial{2,0}(unsafe, e, 4))
            it, state = iterate(vi, state)
            @test it == exponents(IntMonomial{2,0}(unsafe, e, 5))
            it, state = iterate(vi, state)
            @test it == exponents(IntMonomial{2,0}(unsafe, e, 6))
            @test isnothing(iterate(vi, state))
        end
        @testset "enumerated (unit range)" begin
            vi = veciter(mv, true)
            it, state = iterate(vi)
            @test it == (2, exponents(IntMonomial{2,0}(unsafe, e, 2)))
            it, state = iterate(vi, state)
            @test it == (3, exponents(IntMonomial{2,0}(unsafe, e, 3)))
            it, state = iterate(vi, state)
            @test it == (4, exponents(IntMonomial{2,0}(unsafe, e, 4)))
            it, state = iterate(vi, state)
            @test it == (5, exponents(IntMonomial{2,0}(unsafe, e, 5)))
            it, state = iterate(vi, state)
            @test it == (6, exponents(IntMonomial{2,0}(unsafe, e, 6)))
            @test isnothing(iterate(vi, state))
        end
    end

    @testset "monomial multiplication" begin
        m = IntMonomial{2,0}([1, 3])
        @testset "all" begin
            # with index conversion - mv has Int, m has UInt
            mv = IntMonomialVector{2,0}([2, 3, 5, 8, 13])
            mmv = @inferred m * mv
            @test mmv.e isa ExponentsAll{2,UInt}
            @test mmv.indices == UInt[0x11, 0x12, 0x18, 0x1f, 0x28]
        end
        @testset "degree, large enough" begin
            # with index conversion
            mv = IntMonomialVector{2,0}(ExponentsDegree{2,Int}(0, 8), [2, 3, 5, 8, 13])
            mmv = @inferred m * mv
            @test mmv.e isa ExponentsDegree{2,UInt} && mmv.e.mindeg == 0 && mmv.e.maxdeg == 8
            @test mmv.indices == UInt[0x11, 0x12, 0x18, 0x1f, 0x28]

            # no index conversion
            mv = IntMonomialVector{2,0}(ExponentsDegree{2,UInt}(0, 8), UInt[2, 3, 5, 8, 13])
            mmv = @inferred m * mv
            @test mmv.e === mv.e
            @test mmv.indices == UInt[0x11, 0x12, 0x18, 0x1f, 0x28]
        end
        @testset "degree, must enlarge" begin
            mv = IntMonomialVector{2,0}(ExponentsDegree{2,Int}(2, 3))
            mmv = @inferred mv * m
            @test mmv.e isa ExponentsDegree{2,UInt} && mmv.e.mindeg == 2 && mmv.e.maxdeg == 7
            @test mmv.indices == UInt[0x14, 0x15, 0x16, 0x1b, 0x1c, 0x1d, 0x1e]
        end
        @testset "multidegree" begin
            mv = IntMonomialVector{2,0}(ExponentsMultideg{2,Int}(1, 4, [0, 1], [2, 2]))
            mmv = @inferred m * mv
            @test mmv.e isa ExponentsMultideg{2,UInt} && mmv.e.mindeg == 5 && mmv.e.maxdeg == 8 &&
                mmv.e.minmultideg == [1, 4] && mmv.e.maxmultideg == [3, 5]
            @test mmv isa IntPolynomials.IntMonomialVectorComplete

            mv = IntMonomialVector{2,0}(ExponentsMultideg{2,Int}(1, 4, [0, 1], [2, 2]), 3:5)
            mmv = @inferred mv * m
            @test mmv.e isa ExponentsMultideg{2,UInt} && mmv.e.mindeg == 5 && mmv.e.maxdeg == 8 &&
                mmv.e.minmultideg == [1, 4] && mmv.e.maxmultideg == [3, 5]
            @test mmv.indices === UInt(3):UInt(5)

            mv = IntMonomialVector{2,0}(ExponentsMultideg{2,UInt}(1, 4, [0, 1], [2, 2]), UInt[2, 6])
            mmv = @inferred mv * m
            @test mmv.e isa ExponentsMultideg{2,UInt} && mmv.e.mindeg == 5 && mmv.e.maxdeg == 8 &&
                mmv.e.minmultideg == [1, 4] && mmv.e.maxmultideg == [3, 5]
            @test mmv.indices === mv.indices
        end
    end
end

@testset "Polynomial" begin
    # polynomial is just a very thin wrapper, and its main importance is in converting other MP polynomials to IntPolynomial
    # (which is also mainly a feature of IntMonomialVector, but we'll test it here)
    DynamicPolynomials.@polyvar x

    @test terms(IntPolynomial(polynomial([1, x^2, x, 2x^2]))) ==
        [Term(1, IntMonomial{1,0}([0])),
         Term(1, IntMonomial{1,0}([1])),
         Term(3, IntMonomial{1,0}([2]))]

    @test term(IntPolynomial(x + x^2 - x)) isa AbstractTerm

    DynamicPolynomials.@complex_polyvar y
    p = IntPolynomial(3x^2 * y^4 + 2x)
    @test terms(p)[1] == Term(2, IntMonomial{1,1}([1], [0], [0]))
    @test terms(p)[end] == Term(3, IntMonomial{1,1}([2], [4], [0]))
    typetests(p)
    typetests([p, IntPolynomial(x + y)])

    @test coefficient(IntPolynomial(2x + 4y^2 + 3), IntMonomial{1,1}([0], [2], [0])) == 4
    @test coefficient(IntPolynomial(2x + 4y^2 + 3), IntMonomial{1,1}([2], [0], [0])) == 0

    # no coefficient with variable selection, requires zero(::IntPolynomial)
    p = IntPolynomial(x^2 - x^2)
    @test maxdegree(p) == 0
    @test maxdegree(p, IntRealVariable{1,1}(1)) == 0
    @test maxdegree(p, IntComplexVariable{1,1}(1)) == 0
    @test mindegree(p) == 0
    @test mindegree(p, IntRealVariable{1,1}(1)) == 0
    @test mindegree(p, IntComplexVariable{1,1}(1)) == 0
    @test extdegree(p) == (0, 0)
    @test extdegree(p, IntRealVariable{1,1}(1)) == (0, 0)
    @test extdegree(p, IntComplexVariable{1,1}(1)) == (0, 0)
    q = IntPolynomial(x * y + 2 + x^2 * y + x + y)
    @test maxdegree(q) == 3
    @test maxdegree(q, IntRealVariable{1,1}(1)) == 2
    @test maxdegree(q, IntComplexVariable{1,1}(1)) == 1
    @test mindegree(q) == 0
    @test mindegree(q, IntRealVariable{1,1}(1)) == 0
    @test mindegree(q, IntComplexVariable{1,1}(1)) == 0
    @test extdegree(q) == (0, 3)
    @test extdegree(q, IntRealVariable{1,1}(1)) == (0, 2)
    @test extdegree(q, IntComplexVariable{1,1}(1)) == (0, 1)

    @test coefficients(IntPolynomial(x * y + 2 + 3x^2 * y + 4x + 6y)) == [2, 6, 4, 1, 3]
    @test coefficients(IntPolynomial(x * y + 2 + 3x^2 * y + 4x + 6y),
        [IntMonomial{1,1}([1], [0], [0]),
         IntMonomial{1,1}([1], [2], [0]),
         IntMonomial{1,1}([1], [1], [0]),
         IntMonomial{1,1}([2], [1], [0]),
         IntMonomial{1,1}([0], [1], [0]),
         IntMonomial{1,1}([3], [0], [0])
        ]) == [4, 0, 1, 3, 6, 0]
    @test monomials(IntPolynomial(4x^2 * y + x * y + 2x + 3))[1:1] == [constant_monomial(IntPolynomial(x * y))]

    for p in [IntPolynomial(polynomial([4, 9], [x, x * x])), IntPolynomial(polynomial([9, 4], [x * x, x]))]
        @test coefficients(p) == [4, 9]
        @test monomials(p)[1] == IntMonomial{1,0}([1])
        @test monomials(p)[2] == IntMonomial{1,0}([2])
        @test monomials(p)[1:2][1] == IntMonomial{1,0}([1])
        @test monomials(p)[1:2][2] == IntMonomial{1,0}([2])
    end

    @test IntPolynomial(x + 2y)' == IntPolynomial(x + 2conj(y))
    @test transpose(IntPolynomial(x + y)) == IntPolynomial(x + y)
    @test transpose(IntPolynomial(Term([1 2; 3 4], x^1))) == IntPolynomial(Term([1 3; 2 4], x^1))

    # none of the rest
end

# maybe promotion?
@testset "Monomial index and iterator" begin
    @test length(ExponentsDegree{5,UInt}(0, 4)) == 126
    mons = monomials(3, 2, 0:4)
    @test length(mons) == 330
    @test mons[35].index == 35
    @test mons[257].index == 257
    # mons[17] = z₁z̄₁, mons[25] = x₂z₂, mons[205] = x₂z₁z̄₂z₂
    @test (mons[17] * mons[25]).index == 205
    alloc_test(let m1=mons[17], m2=mons[25]; () -> (m1 * m2).index end, 0)

    @test monomial_index(IntRealVariable{3,5}(2)) == 13
    @test monomial_index(IntComplexVariable{3,5}(4)) == 5
    @test monomial_index(IntComplexVariable{3,5}(3, true)) == 6

    # mons[35] = x₁x₂, mons[113] = x₁x₂²
    @test monomial_index(mons[35], IntRealVariable{3,2}(2)) == 113
    # mons[109] = x₁x₂z₂
    @test monomial_index(mons[35], IntComplexVariable{3,2}(2)) == 109
    # mons[108] = x₁x₂z̄₂
    @test monomial_index(mons[35], IntComplexVariable{3,2}(2, true)) == 108

    @testset failfast=true "Index consistency with DynamicPolynomials" begin
        DynamicPolynomials.@polyvar x[1:3]
        multiit = Iterators.product(0:4, 0:4, 0:4)
        for mindeg in 0:4, maxdeg in 0:4, minmultideg in multiit, maxmultideg in multiit
            minm, maxm = collect(minmultideg), collect(maxmultideg)
            if max(mindeg, sum(minm)) > min(maxdeg, sum(maxm)) || any(minmultideg .> maxmultideg)
                @test_throws ArgumentError ExponentsMultideg{3,UInt}(mindeg, maxdeg, minm, maxm)
            else
                mi = ExponentsMultideg{3,UInt}(mindeg, maxdeg, minm, maxm)
                exp = MP.exponents.(monomials(x, Int(mindeg):Int(maxdeg), m -> all(minm .≤ MP.exponents(m) .≤ maxm)))
                @test collect(mi) == exp
                @test length(mi) == length(exp)
                exponents = Vector{Int}(undef, 3)
                for (i, mipow) in zip(Iterators.countfrom(one(UInt)), mi)
                    @test_throws(sum(mipow) == mi.maxdeg ? BoundsError : ArgumentError,
                        exponents_from_index(mi, i, sum(mipow) +1))
                    e = exponents_from_index(unsafe, mi, i)
                    copyto!(exponents, e)
                    mind = exponents_to_index(mi, exponents)
                    mind2 = exponents_to_index(mi, e)
                    if exponents != mipow || mind !== i || mind !== mind2
                        println(mindeg, " - ", maxdeg, " - ", minmultideg, " - ", maxmultideg, ": index ", i,
                            " gives exponents ", exponents, " which claim index ", mind, " via back-calculation and ",
                            mind2, " via short-path")
                    end
                    @test exponents == mipow
                    @test mind === i
                end
                @test_throws BoundsError exponents_from_index(mi, length(mi) +1)
                @test iszero(exponents_to_index(mi, maxm .+ one(eltype(maxm))))
            end
        end
    end
end