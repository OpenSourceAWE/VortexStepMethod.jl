using Pkg
if !("ControlPlots" ∈ keys(Pkg.project().dependencies))
    Pkg.activate(@__DIR__)
end

using ControlPlots
using VortexStepMethod
using REPL.TerminalMenus

url = "https://opensourceawe.github.io/VortexStepMethod.jl/dev"

examples_dir = joinpath(@__DIR__, "..", "examples")

example_files = [
    "V3_kite.jl",
    "billowing.jl",
    "pyramid_model.jl",
    "rectangular_wing.jl",
    "ram_air_kite.jl",
    "stall_model.jl",
    "bench.jl",
    "cleanup.jl",
]

function run_all()
    for f in example_files
        f == "cleanup.jl" && continue
        println("\n" * "="^60)
        println("Running: $f")
        println("="^60)
        try
            include(joinpath(examples_dir, f))
        catch e
            @error "Failed: $f" exception=(e, catch_backtrace())
        end
    end
    println("\nAll examples completed.")
end

function example_menu()
    options = [
        [("$(splitext(f)[1]) = include(\"../examples/$f\")")
         for f in example_files];
        "help_me = VortexStepMethod.help(\"$url\")";
        "quit"
    ]
    active = true
    while active
        menu = RadioMenu(options, pagesize=8)
        choice = request(
            "\nChoose function to execute or `q` to quit: ",
            menu)
        if choice != -1 && choice != length(options)
            eval(Meta.parse(options[choice]))
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end
end

if "--run-all" in ARGS
    run_all()
else
    example_menu()
end
