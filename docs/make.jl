using ReactionSensitivity
using Documenter

DocMeta.setdocmeta!(ReactionSensitivity, :DocTestSetup, :(using ReactionSensitivity); recursive=true)

makedocs(;
    modules=[ReactionSensitivity],
    authors="Vinod Janardhanan",
    repo="https://github.com/vinodjanardhanan/ReactionSensitivity.jl/blob/{commit}{path}#{line}",
    sitename="ReactionSensitivity.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vinodjanardhanan.github.io/ReactionSensitivity.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/vinodjanardhanan/ReactionSensitivity.jl",
    devbranch="main",
)
