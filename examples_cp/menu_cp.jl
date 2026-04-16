using Pkg
if !("ControlPlots" ∈ keys(Pkg.project().dependencies))
    Pkg.activate(@__DIR__)
end

using ControlPlots
using VortexStepMethod
using REPL.TerminalMenus

url = "https://opensourceawe.github.io/VortexStepMethod.jl/dev"

options = [
           "V3_kite = include(\"../examples/V3_kite.jl\")",
           "billowing = include(\"../examples/billowing.jl\")",
           "pyramid_model = include(\"../examples/pyramid_model.jl\")",
           "rectangular_wing = include(\"../examples/rectangular_wing.jl\")",
           "ram_air_kite = include(\"../examples/ram_air_kite.jl\")",
           "stall_model = include(\"../examples/stall_model.jl\")",
           "bench = include(\"../examples/bench.jl\")",
           "cleanup = include(\"../examples/cleanup.jl\")",
           "help_me = VortexStepMethod.help(url)",
           "quit"]

function example_menu()
    active = true
    while active
        menu = RadioMenu(options, pagesize=8)
        choice = request("\nChoose function to execute or `q` to quit: ", menu)

        if choice != -1 && choice != length(options)
            eval(Meta.parse(options[choice]))
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end
end

example_menu()
