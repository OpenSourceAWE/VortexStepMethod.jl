using Makie
using VortexStepMethod
using VortexStepMethod.AirfoilAero
using VortexStepMethod.ObjAdapter
using Documenter

DocMeta.setdocmeta!(VortexStepMethod, :DocTestSetup, :(using VortexStepMethod); recursive=true)

makedocs(;
    modules=[VortexStepMethod,
             Base.get_extension(VortexStepMethod, :VortexStepMethodMakieExt),
             VortexStepMethod.AirfoilAero,
             VortexStepMethod.ObjAdapter],
    authors="Uwe Fechner <uwe.fechner.msc@gmail.com>, Bart van de Lint <bart@vandelint.net> and contributors",
    sitename="VortexStepMethod.jl",
    warnonly=[:cross_references],
    format = Documenter.HTML(prettyurls = haskey(ENV, "CI")),
    pages=[
        "Home" => "index.md",
        "How it works" => "explanation.md",
        "CAD mesh to model" => "airfoil_pipeline.md",
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
