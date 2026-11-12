using Documenter, Geodesics

makedocs(
    sitename = "Geodesics.jl documentation",
    pages = [
        "Home" => "index.md",
        "Quickstart" => "quickstart.md",
        "Ellipsoids" => "ellipsoids.md",
        "Algorithms" => "algorithms.md",
        "Manual" => "manual.md",
        "Unitful" => "unitful.md",
        "Type index" => "type-index.md",
        "Function index" => "function-index.md",
        "References" => "references.md",
        ]
    )

deploydocs(
    repo = "github.com/anowacki/Geodesics.jl.git",
)
