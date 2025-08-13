# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# List all methods of VortexStepMethod and report if they are documented;
# By default, also exported constants are reported
using VortexStepMethod

exported_names = names(VortexStepMethod)
exported_functions = filter(n -> isdefined(VortexStepMethod, n) && isa(getfield(VortexStepMethod, n), Function), exported_names)

function find_exported_functions(exported_functions)
    total = 0
    for fun in exported_functions
        f = getfield(VortexStepMethod, fun)
        mes = methods(f)
        for me in mes
            println(me.sig)
            total += 1
        end
    end
    return total
end

println("Exported methods:\n")
total = find_exported_functions(exported_functions)
println("\nTotal: $total")

function check_exported_docs(mod::Module; only_functions=false)
    exported_symbols = names(mod, all=false)
    doc_status = Dict{Symbol,Bool}()
    for sym in exported_symbols
        if only_functions && !isa(getfield(mod, sym), Function)
            continue
        end
        doc_status[sym] = Base.Docs.hasdoc(mod, sym)
    end
    return doc_status
end

# Main execution block:
results = check_exported_docs(VortexStepMethod)
undocumented = filter(kv -> !kv[2], results)
println("\nUndocumented exported symbols: \n", keys(undocumented))
