using ROMEO
using Documenter

DocMeta.setdocmeta!(ROMEO, :DocTestSetup, :(using ROMEO); recursive=true)

makedocs(;
    modules=[ROMEO],
    authors="Korbinian Eckstein korbinian90@gmail.com",
    repo="https://github.com/korbinian90/ROMEO.jl/blob/{commit}{path}#{line}",
    sitename="ROMEO.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://korbinian90.github.io/ROMEO.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/korbinian90/ROMEO.jl",
    devbranch="master",
)
