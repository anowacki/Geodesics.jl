using Documenter, Geodesics

makedocs(
    sitename = "Geodesics.jl documentation",
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md",
        "Function index" => "function-index.md",
        ]
    )

deploydocs(
    repo = "github.com/anowacki/Geodesics.jl.git",
)
