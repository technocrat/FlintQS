using Documenter, MaterialDocs, FlintQS

DocMeta.setdocmeta!(FlintQS, :DocTestSetup, :(using FlintQS); recursive = true)

makedocs(;
    modules = [FlintQS],
    authors = "technocrat",
    sitename = "FlintQS.jl",
    format = Material3(theme = :ocean_depth, dark_mode = :toggle, edit_link = "julia-port"),
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "How it works" => "algorithm.md",
        "API Reference" => "api.md",
        "Benchmarks" => "benchmarks.md",
    ],
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(; repo = "github.com/technocrat/FlintQS", devbranch = "julia-port")
end
