using ALPINEExplorer
using Documenter

DocMeta.setdocmeta!(
    ALPINEExplorer, :DocTestSetup, :(using ALPINEExplorer); recursive=true
)

makedocs(;
    modules=[ALPINEExplorer],
    authors="Nicholas Minor <nrminor@wisc.edu>",
    repo="https://github.com/nrminor/ALPINEExplorer.jl/blob/{commit}{path}#{line}",
    sitename="ALPINEExplorer.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://nrminor.github.io/ALPINEExplorer.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/nrminor/ALPINEExplorer.jl", devbranch="main")
