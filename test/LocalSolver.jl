using Test
using PolynomialOptimization
using MultivariatePolynomials
import DynamicPolynomials

const startval = [0.7210476969933126, 0.6380699455835934, 0.8532153291949746, 0.3952826714480566, 0.3064291214270869,
                  0.9778462338230722]

@testset "Some real-valued examples (scalar)" begin
    DynamicPolynomials.@polyvar x[1:2]
    prob = poly_problem(-(x[1] - 1)^2 - (x[1] - x[2])^2 - (x[2] - 3)^2, 1,
                        nonneg=[1 - (x[1] - 1)^2, 1 - (x[1] - x[2])^2, 1 - (x[2] - 3)^2])
    optval, optpt = poly_optimize(:LANCELOT, prob)(startval[1:2])
    @test optval ≈ -2 atol = 2e-5
    @test optpt ≈ [1., 2.] atol = 2e-5
end

@testset "Some real-valued examples (matrix)" begin
    DynamicPolynomials.@polyvar x[1:2]
    prob = poly_problem(-x[1]^2 - x[2]^2, 1, psd=[[1-4x[1]*x[2] x[1]; x[1] 4-x[1]^2-x[2]^2]])
    optval, optpt = poly_optimize(:LANCELOT, prob)(startval[1:2])
    @test optval ≈ -4 atol = 2e-5
    @test optpt ≈ [0., 2.] atol = 2e-5
end

@testset "Some real-valued examples (matrix and equality)" begin
    DynamicPolynomials.@polyvar x[1:2]
    prob = poly_problem(-x[1]^2 - x[2]^2, 1, zero=[x[1] + x[2] - 1], psd=[[1-4x[1]*x[2] x[1]; x[1] 4-x[1]^2-x[2]^2]])
    optval, optpt = poly_optimize(:LANCELOT, prob)(startval[1:2])
    # global would be -3.90489157833684, [-0.804778061268820., 1.804778061268820]
    @test optval ≈ -3.50442134489527 atol = 2e-5
    @test optpt ≈ [1.725646827649273, -0.7256475529616774] atol = 2e-5
end

@testset "Example 3.1 from correlative term sparsity paper" begin
    DynamicPolynomials.@polyvar x[1:3]
    prob = poly_problem(1 + x[1]^2 + x[2]^2 + x[3]^2 + x[1] * x[2] + x[2] * x[3] + x[3], 2)
    optval, optpt = poly_optimize(:LANCELOT, prob)(startval[1:3])
    @test optval ≈ 0.625 atol = 2e-5
    @test optpt ≈ [-0.25, 0.5, -0.75] atol = 2e-5
end

@testset "Example 3.4 from correlative term sparsity paper" begin
    DynamicPolynomials.@polyvar x[1:6]
    prob = poly_problem(1 + sum(x[i]^4 for i in 1:6) + x[1] * x[2] * x[3] + x[3] * x[4] * x[5] + x[3] * x[4] * x[6] +
                      x[3] * x[5] * x[6] + x[4] * x[5] * x[6], 2)
    optval, optpt = poly_optimize(:LANCELOT, prob)(startval[1:6])
    @test optval ≈ 0.9412824077069458 atol = 2e-5
    @test optpt ≈ [0.33346423231291855, -0.3334642984669165, 0.4447934153234278, 0.3948134188295012, -0.5023550515807035,
                   0.3948134185292575] atol = 2e-5
end

@testset "Example 4.1 from Zhen, Fantuzzi, Papachristodoulou review" begin
    DynamicPolynomials.@polyvar x[1:50]
    prob = poly_problem(sum((x[i-1] + x[i] + x[i+1])^4 for i in 2:49), 2)
    optval, optpt = poly_optimize(:LANCELOT, prob)(Float64.(1:50))
    @test optval ≈ 0 atol = 2e-5
    @test optpt ≈ [-8.770717648147754, 19.195148269050982, -10.4241040465211, -8.770587097595792, 19.19527829208145,
        -10.424013740859381, -8.770471446696972, 19.195447412580137, -10.423833297327013, -8.770341568400156,
        19.195532741833514, -10.423801795214992, -8.77028772970991, 19.19565649001814, -10.423598397978823,
        -8.77006973403981, 19.195852184485876, -10.423418747027155, -8.76990685774942, 19.19599374686716,
        -10.423299334909986, -8.769810139085106, 19.19606425096866, -10.423257370403109, -8.769794147119343,
        19.196047558946383, -10.423190385206096, -8.769633902470884, 19.19631134219725, -10.422925128314402,
        -8.769370293882075, 19.196578024647415, -10.422659090267413, -8.769104331027545, 19.196844022433854,
        -10.422393015240473, -8.768836359477682, 19.197109529367175, -10.422127765225504, -8.768563560489481,
        19.197373392714155, -10.421868309560763, -8.768275353928003, 19.197633449202378, -10.421636448256224,
        -8.76793545402412, 19.19789189224112, -10.421521467050434, -8.76749234945018, 19.198052783841895] atol = 2e-5
end

@testset "Example 4.4 from Zhen, Fantuzzi, Papachristodoulou review" begin
    DynamicPolynomials.@polyvar x[1:4]
    prob = poly_problem(2 + x[1]^2 * x[4]^2 * (x[1]^2 * x[4]^2 - 1) - x[1]^2 + x[1]^4 +
                        sum((x[i]^2 * x[i-1]^2 * (x[i]^2 * x[i-1]^2 - 1) - x[i]^2 + x[i]^4 for i in 2:4)), 4)
    optval, optpt = poly_optimize(:LANCELOT, prob)(startval[1:4])
    @test optval ≈ 0.110118425157690 atol = 2e-5
    @test optpt ≈ [0.793700525984100, 0.793700525984100, 0.793700525984100, 0.793700525984100] atol = 2e-5
end