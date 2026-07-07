using Makie
using VortexStepMethod
using AirfoilAero
using ObjAdapter
using Documenter

DocMeta.setdocmeta!(VortexStepMethod, :DocTestSetup, :(using VortexStepMethod); recursive=true)

makedocs(;
    modules=[VortexStepMethod,
             Base.get_extension(VortexStepMethod, :VortexStepMethodMakieExt),
             AirfoilAero,
             ObjAdapter,
             Base.get_extension(ObjAdapter, :ObjAdapterMakieExt)],
    authors="Uwe Fechner <uwe.fechner.msc@gmail.com>, Bart van de Lint <bart@vandelint.net> and contributors",
    sitename="VortexStepMethod.jl",
    checkdocs=:none,
    format = Documenter.HTML(prettyurls = haskey(ENV, "CI")),
    pages=[
        "Home" => "index.md",
        "How it works" => "explanation.md",
        "Examples" => "examples.md",
        "Exported Functions" => "functions.md",
        "Exported Types" => "types.md",
        "Private Functions" => "private_functions.md",
        "Private Types" => "private_types.md",
        "Reference Frames" => "reference_frames.md",
        "Tips and tricks" => "tips_and_tricks.md",
        "Glossary" => "glossary.md"
    ],
)

deploydocs(;
    repo="github.com/OpenSourceAWE/VortexStepMethod.jl",
    devbranch="main",
)
