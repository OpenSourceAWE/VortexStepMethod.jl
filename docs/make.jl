using VortexStepMethod
using Pkg
if ("TestEnv" ∈ keys(Pkg.project().dependencies))
    if ! ("Documents" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
end
using Documenter

DocMeta.setdocmeta!(VortexStepMethod, :DocTestSetup, :(using VortexStepMethod); recursive=true)

makedocs(;
    modules=[VortexStepMethod],
    authors="Uwe Fechner <uwe.fechner.msc@gmail.com>, Bart van de Lint <bart@vandelint.net> and contributors",
    sitename="VortexStepMethod.jl",
    checkdocs=:none,
    format = Documenter.HTML(prettyurls = haskey(ENV, "CI")),
    pages=[
        "Home" => "index.md",
        "Exported Functions" => "functions.md",
    ],
)

deploydocs(;
    repo="github.com/Albatross-Kite-Transport/VortexStepMethod.jl",
    devbranch="main",
)
