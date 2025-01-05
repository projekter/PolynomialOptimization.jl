using Documenter
using PolynomialOptimization, GALAHAD
using MultivariatePolynomials, DynamicPolynomials

makedocs(sitename="PolynomialOptimization.jl",
    modules=filter(m -> startswith(string(m), "PolynomialOptimization"), Docs.modules),
    format=Documenter.HTML(prettyurls=false, size_threshold=nothing),
    pages=[
        "index.md",
        "guide.md",
        "reference.md",
        "includedsolvers.md",
        "backend.md",
        "auxreference.md",
        "simplepolynomials.md"
    ],
    warnonly=:missing_docs,
    doctest=false
)