using Pkg
Pkg.activate(@__DIR__)

using GLMakie
using CairoMakie
using VortexStepMethod
using REPL.TerminalMenus

url = "https://opensourceawe.github.io/VortexStepMethod.jl/dev"

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
    failed_examples = String[]
    for f in example_files
        f == "cleanup.jl" && continue
        println("\n" * "="^60)
        println("Running: $f")
        println("="^60)
        try
            include(joinpath(@__DIR__, f))
        catch e
            @error "Failed: $f" exception=(e, catch_backtrace())
            push!(failed_examples, f)
        end
    end
    if isempty(failed_examples)
        println("\nAll examples completed.")
    else
        failed_list = join(failed_examples, ", ")
        throw(ErrorException("run_all failed for $(length(failed_examples)) example(s): $failed_list"))
    end
end

function example_menu()
    options = [
        [("$( splitext(f)[1]) = include(\"$f\")") for f in example_files];
        "GLMakie.activate!()";
        "CairoMakie.activate!()";
        "help_me = VortexStepMethod.help(\"$url\")";
        "quit"
    ]
    active = true
    while active
        menu = RadioMenu(options, pagesize=11)
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
