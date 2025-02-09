using Test
using PolynomialOptimization.Solvers.SketchyCGAL
using LinearAlgebra
using SparseArrays
using MAT

optvals = Dict(
    # these values were calculated exactly using Mosek (with default precision)
    "G1.mat" => -1.2083197636e+04,
    "G2.mat" => -1.2089429689e+04,
    "G3.mat" => -1.2084333014e+04,
    "G4.mat" => -1.2111451013e+04,
    "G5.mat" => -1.2099886848e+04,
    "G6.mat" => -2.6561593626e+03,
    "G7.mat" => -2.4892623944e+03,
    "G8.mat" => -2.5069333027e+03,
    "G9.mat" => -2.5287311681e+03,
    "G10.mat" => -2.4850630812e+03,
    "G11.mat" => -6.2916464750e+02,
    "G12.mat" => -6.2387386994e+02,
    "G13.mat" => -6.4713639164e+02,
    "G14.mat" => -3.1915667947e+03,
    "G15.mat" => -3.1715580463e+03,
    "G16.mat" => -3.1750183547e+03,
    "G17.mat" => -3.1713268354e+03,
    "G18.mat" => -1.1660097961e+03,
    "G19.mat" => -1.0820103523e+03,
    "G20.mat" => -1.1113918128e+03,
    "G21.mat" => -1.1042837463e+03,
    "G22.mat" => -1.4135945710e+04,
    "G23.mat" => -1.4142108531e+04,
    "G24.mat" => -1.4140855773e+04,
    "G25.mat" => -1.4144245347e+04,
    "G26.mat" => -1.4132870713e+04,
    "G27.mat" => -4.1416580579e+03,
    "G28.mat" => -4.1007888828e+03,
    "G29.mat" => -4.2088860264e+03,
    "G30.mat" => -4.2153786416e+03,
    "G31.mat" => -4.1166803521e+03,
    "G32.mat" => -1.5676374418e+03,
    "G33.mat" => -1.5443115348e+03,
    "G34.mat" => -1.5466868817e+03,
    "G35.mat" => -8.0147396982e+03,
    "G36.mat" => -8.0059637581e+03,
    "G37.mat" => -8.0186232994e+03,
    "G38.mat" => -8.0149742412e+03,
    "G39.mat" => -2.8776418333e+03,
    "G40.mat" => -2.8647891658e+03,
    "G41.mat" => -2.8652147253e+03,
    "G42.mat" => -2.9462503935e+03,
    "G43.mat" => -7.0322218081e+03,
    "G44.mat" => -7.0278846753e+03,
    "G45.mat" => -7.0247813536e+03,
    "G46.mat" => -7.0299330812e+03,
    "G47.mat" => -7.0366604610e+03,
    "G48.mat" => -5.9999999500e+03,
    "G49.mat" => -5.9999999932e+03,
    "G50.mat" => -5.9881720520e+03,
    "G51.mat" => -4.0062554981e+03,
    "G52.mat" => -4.0096387634e+03,
    "G53.mat" => -4.0097184892e+03,
    "G54.mat" => -4.0061941087e+03,
    "G55.mat" => -1.1039460271e+04,
    "G56.mat" => -4.7600044021e+03,
    "G57.mat" => -3.8854865345e+03,
    "G58.mat" => -2.0136189748e+04,
    "G59.mat" => -7.3123104477e+03,
    "G60.mat" => -1.5222267971e+04,
    "G61.mat" => -6.8281043216e+03,
    "G62.mat" => -5.4309028843e+03,
    "G63.mat" => -2.8244417802e+04,
)
for (mat, optval) in optvals
    @testset "$mat" begin
        A = matread("./Gset/$mat")["Problem"]["A"]
        n = size(A, 1)
        C = spdiagm(A * ones(n)) - A
        C = .5*(C + C')
        C = -.25 .* C
        scale_x = 1/n
        scale_c = 1/norm(C)
        for method in (mat == "G01.mat" ? (:lanczos_space, :lanczos_time, :lobpcg_fast, :lobpcg_accurate) : (:auto,))
            opt = sketchy_cgal(
                (v, u, idx, α, β) -> mul!(v, C, u, α * scale_c, β),
                (v, u, z, idx, α, β) -> v .= α .* u .* z .+ β .* v,
                (v, u, i, α, β) -> v .= α .* u .^ 2 .+ β .* v,
                n,
                fill(scale_x, n);
                α=(1., 1.),
                rank=10,
                rescale_C=scale_c,
                rescale_X=scale_x,
                #
                primitive3_normsquare=1.,
                verbose=true,
                ϵ=0.01,
                method
            )[2]
            @test opt ≈ optval rtol=1e-2
        end
    end
end