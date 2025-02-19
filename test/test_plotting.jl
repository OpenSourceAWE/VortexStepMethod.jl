using VortexStepMethod
using ControlPlots
using Test

plt.ioff()
@testset "Plotting" begin
    fig = plt.plot([1,2,3])
    @test fig isa Vector{plt.PyObject}
    show_plot(fig)
end
plt.ion()