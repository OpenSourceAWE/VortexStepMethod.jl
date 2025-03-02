# build and display the html documentation locally
# you must have installed the package LiveServer in your global environment

using Pkg
if !("Documenter" ∈ keys(Pkg.project().dependencies))
    using TestEnv
    TestEnv.activate()
end
using LiveServer; servedocs(launch_browser=true)
