using Documenter
using PolynomialOptimization
using MultivariatePolynomials
using DynamicPolynomials

makedocs(sitename="PolynomialOptimization.jl", modules=[PolynomialOptimization],
    format=Documenter.HTML(prettyurls=false),
    pages=[
        "index.md",
        "guide.md",
        "reference.md",
        "background.md",
    ])