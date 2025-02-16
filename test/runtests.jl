using Test

cd("..")
println("Running tests...")
@testset verbose = true "Testing VortexStepMethod..." begin
    include("test_bound_filament.jl")
end