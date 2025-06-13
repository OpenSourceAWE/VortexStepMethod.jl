using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end

using ControlPlots
using VortexStepMethod
using REPL.TerminalMenus

url = "https://opensourceawe.github.io/VortexStepMethod.jl/dev"

options = ["rectangular_wing = include(\"rectangular_wing.jl\")",
           "ram_air_kite = include(\"ram_air_kite.jl\")",
           "stall_model = include(\"stall_model.jl\")",
           "bench = include(\"bench.jl\")",
           "cleanup = include(\"cleanup.jl\")",
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