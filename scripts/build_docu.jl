# build and display the html documentation locally
# run with: julia --project=docs scripts/build_docu.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "docs"))
using LiveServer; servedocs(launch_browser=true)
