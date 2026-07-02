using Documenter
using MolecularGaussians
using MolecularGraph

DocMeta.setdocmeta!(
    MolecularGaussians,
    :DocTestSetup,
    :(using MolecularGaussians, MolecularGraph);
    recursive = true,
)

makedocs(;
    modules = [MolecularGaussians],
    authors = "Tom McGrath <tmcgrath325@gmail.com> and contributors",
    sitename = "MolecularGaussians.jl",
    format = Documenter.HTML(;
        canonical = "https://tmcgrath325.github.io/MolecularGaussians.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Concepts" => "concepts.md",
        "API reference" => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/tmcgrath325/MolecularGaussians.jl",
    devbranch = "main",
)
