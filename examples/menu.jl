using ControlPlots
using VortexStepMethod
using REPL.TerminalMenus

url = "https://albatross-kite-transport.github.io/VortexStepMethod.jl/dev"

function help() 
    if Sys.iswindows()
        run(`cmd /c start $url`)
    else
        io = IOBuffer()
        run(pipeline(`xdg-open $url`, stderr = io))
        out_data = String(take!(io)) 
    end
    nothing
end

options = ["help_me = help()",
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