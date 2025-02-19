using VortexStepMethod
using ControlPlots
using Test

plt.ioff()
@testset "Plotting" begin
    fig = plt.figure(figsize=(14, 14))
    res = plt.plot([1,2,3])
    @test fig isa plt.PyPlot.Figure
    @test res isa Vector{plt.PyObject}
    @test isfile("/tmp/plot.pdf")
    rm("/tmp/plot.pdf")
    show_plot(fig)
    save_plot(fig, "/tmp", "plot")
end
plt.ion()