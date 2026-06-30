"""
Benchmark comparing NeuralFoil neural-network evaluation against polar
interpolation evaluation.

The NeuralFoil timing mirrors the methodology of the upstream Python NeuralFoil
speed benchmark (https://github.com/peterdsharpe/NeuralFoil): a single-case run
and a vectorized batch of 100_000 cases, reported per model size. The
interpolation timing evaluates the linear `Interpolations.jl` lookups that the
VSM solver actually uses at runtime (built exactly as in `panel.jl`).
"""

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
using VortexStepMethod
using VortexStepMethod: load_neuralfoil_model, neuralfoil_aero, fit_kulfan_parameters,
                        prepare_inputs, nn_forward
using BenchmarkTools
using Interpolations
using Printf
using Statistics: mean

read_dat(path) = begin
    x = Float64[]
    y = Float64[]
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || !(isdigit(s[1]) || s[1] == '-' || s[1] == '.')) && continue
        parts = split(s)
        length(parts) >= 2 || continue
        push!(x, parse(Float64, parts[1]))
        push!(y, parse(Float64, parts[2]))
    end
    return x, y
end

n_batch = 100_000
reynolds = 5e5
model_sizes = ["xsmall", "small", "medium", "large", "xlarge"]

dat_path = joinpath(@__DIR__, "..", "data", "ram_air_kite", "ram_air_kite_foil.dat")
xc, yc = read_dat(dat_path)
params = fit_kulfan_parameters(xc, yc)

println("Airfoil: $(basename(dat_path)),  Re = $(reynolds),  batch = $(n_batch)\n")

alpha_single = 5.0
alpha_batch = collect(range(-10.0, 20.0, n_batch))

println("=== NeuralFoil (CL, CD, CM per case) ===")
@printf("%-9s %14s %16s %14s\n", "model", "single [ms]", "batch 100k [s]", "per-run [µs]")
for ms in model_sizes
    load_neuralfoil_model(ms)
    t_single = @belapsed neuralfoil_aero($params, $alpha_single, $reynolds; model_size=$ms)
    t_batch = @belapsed neuralfoil_aero($params, $alpha_batch, $reynolds; model_size=$ms)
    @printf("%-9s %14.3f %16.3f %14.3f\n",
        ms, t_single * 1e3, t_batch, t_batch / n_batch * 1e6)
end

println("\n=== Polar interpolation (Interpolations.jl linear lookup) ===")
nf = neuralfoil_aero(params, collect(-180.0:1.0:180.0), reynolds; model_size="xlarge")
alphas_rad = deg2rad.(nf.alpha)
extrap_flat = Interpolations.Flat()
extrap_line = Interpolations.Line()
cl_interp = linear_interpolation(alphas_rad, nf.CL; extrapolation_bc=extrap_flat)
cd_interp = linear_interpolation(alphas_rad, nf.CD; extrapolation_bc=extrap_line)
cm_interp = linear_interpolation(alphas_rad, nf.CM; extrapolation_bc=extrap_flat)

eval_interp(a) = (cl_interp(a), cd_interp(a), cm_interp(a))
function eval_interp_batch(as)
    s = 0.0
    for a in as
        cl, cd, cm = eval_interp(a)
        s += cl + cd + cm
    end
    return s
end

alpha_single_rad = deg2rad(alpha_single)
alpha_batch_rad = deg2rad.(alpha_batch)

t_single = @belapsed eval_interp($alpha_single_rad)
t_batch = @belapsed eval_interp_batch($alpha_batch_rad)
@printf("%-9s %14s %16s %14s\n", "", "single [ns]", "batch 100k [ms]", "per-run [ns]")
@printf("%-9s %14.1f %16.3f %14.1f\n",
    "linear", t_single * 1e9, t_batch * 1e3, t_batch / n_batch * 1e9)

println("\n=== Speed ratio (NeuralFoil xlarge / interpolation, single case) ===")
t_nf = @belapsed neuralfoil_aero($params, $alpha_single, $reynolds; model_size="xlarge")
@printf("NeuralFoil is ~%.0fx slower per single evaluation than linear interpolation\n",
    t_nf / t_single)
nothing
