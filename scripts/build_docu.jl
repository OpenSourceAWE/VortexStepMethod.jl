# build and display the html documentation locally
# run with: julia --project=docs scripts/build_docu.jl

using Pkg
Pkg.develop(path=dirname(@__DIR__))
Pkg.instantiate()
using LiveServer; servedocs(launch_browser=true)
