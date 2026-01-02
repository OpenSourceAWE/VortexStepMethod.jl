using LinearAlgebra
using VortexStepMethod
using GLMakie
using DelimitedFiles

PLOT = true
USE_TEX = false
DEFORM = false

project_dir = dirname(dirname(pathof(VortexStepMethod)))  # Go up one level from src to project root#
literature_paths = [
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results","CFD_RANS_Rey_5e5_Poland2025_alpha_sweep_beta_0_NoStruts.csv"),
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results","CFD_RANS_Rey_10e5_Poland2025_alpha_sweep_beta_0.csv"),
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results","Python_VSM_Rey_5e5_Poland2025_alpha_sweep_beta_0.csv"),
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results","WindTunnel_Re_5e5_Poland2025_alpha_sweep_beta_0.csv"),
    ]
labels= [
    "Julia VSM 2D CFD PCHIP",
    "CFD RANS Re=5e5",
    "CFD RANS Re=10e5 (With Struts)",
    "Python VSM 2D CFD PCHIP Re=5e5",
    "Wind Tunnel Re=5e5 (With Struts)"
    ]

# Load VSM settings from YAML configuration file
settings = VSMSettings("TUDELFT_V3_KITE/vsm_settings.yaml")

# Create wing, body_aero, and solver objects using settings
wing = Wing(settings)
refine!(wing)
body_aero = BodyAerodynamics([wing])
VortexStepMethod.reinit!(body_aero)

if DEFORM
    VortexStepMethod.unrefined_deform!(
        wing,
        deg2rad.(range(-10, 10, length=wing.n_unrefined_sections)),
        deg2rad.(range(0, 0, length=wing.n_unrefined_sections));
        smooth=true
    )
    VortexStepMethod.reinit!(body_aero; init_aero=false)
end

solver = Solver(body_aero, settings)

# Set flight conditions from settings
set_va!(body_aero, settings)

# Extract values for plotting (optional - for reference)
wind_speed = settings.condition.wind_speed
angle_of_attack_deg = settings.condition.alpha
sideslip_deg = settings.condition.beta
yaw_rate = settings.condition.yaw_rate

# Using plotting modules, to create more comprehensive plots
PLOT = true
USE_TEX = false

# Solve and plot combined analysis
results = VortexStepMethod.solve(solver, body_aero; log=true)
fig1 = plot_combined_analysis(
    solver,
    body_aero,
    results;
    solver_label="VSM",
    literature_path_list=literature_paths,
    angle_range=range(-5, 25, length=31),
    angle_type="angle_of_attack",
    angle_of_attack=angle_of_attack_deg,
    side_slip=sideslip_deg,
    v_a=wind_speed,
    title="TU Delft V3 Kite",
    is_show=false,
    use_tex=USE_TEX,
    angle_of_attack_for_spanwise_distribution=10.0,
)
scr1 = display(fig1)
wait(scr1)


# Polar sweep including CMy, using solve! and calculate_results
function compute_polar_with_cm(
    solver,
    body_aero,
    angle_range;
    angle_type::String="angle_of_attack",
    angle_of_attack::Float64=0.0,
    side_slip::Float64=0.0,
    v_a::Float64=10.0
)
    n_angles = length(angle_range)
    cl = zeros(n_angles)
    cd = zeros(n_angles)
    cs = zeros(n_angles)
    cmy = fill(NaN, n_angles)
    reynolds_number = zeros(n_angles)

    gamma_prev = solver.sol.gamma_distribution
    for (i, angle_i) in enumerate(angle_range)
        if angle_type == "angle_of_attack"
            α = deg2rad(angle_i)
            β = deg2rad(side_slip)
        elseif angle_type == "side_slip"
            α = deg2rad(angle_of_attack)
            β = deg2rad(angle_i)
        else
            throw(ArgumentError("angle_type must be 'angle_of_attack' or 'side_slip'"))
        end

        set_va!(body_aero, [cos(α) * cos(β), sin(β), sin(α)] * v_a)
        solve!(solver, body_aero, gamma_prev; log=false)
        gamma_prev = solver.sol.gamma_distribution

        results_local = calculate_results(
            body_aero,
            solver.lr.gamma_new,
            zeros(MVec3),
            solver.density,
            solver.aerodynamic_model_type,
            solver.core_radius_fraction,
            solver.mu,
            solver.lr.alpha_dist,
            solver.lr.v_a_dist,
            solver.sol._chord_dist,
            solver.sol._x_airf_dist,
            solver.sol._y_airf_dist,
            solver.sol._z_airf_dist,
            solver.sol._va_dist,
            solver.br.va_norm_dist,
            solver.br.va_unit_dist,
            body_aero.panels,
            solver.is_only_f_and_gamma_output;
            correct_aoa=solver.correct_aoa
        )

        cl[i] = results_local["cl"]
        cd[i] = results_local["cd"]
        cs[i] = results_local["cs"]
        cmy[i] = get(results_local, "cmy", NaN)
        reynolds_number[i] = results_local["Rey"]
    end

    return (angle=angle_range, cl=cl, cd=cd, cs=cs, cmy=cmy, rey=reynolds_number)
end

function plot_polars_with_cmy(
    solver_list,
    body_aero_list,
    label_list;
    literature_path_list::Vector{String}=String[],
    angle_range=range(-10, 40, step=1),
    angle_type::String="angle_of_attack",
    angle_of_attack::Float64=0.0,
    side_slip::Float64=0.0,
    v_a::Float64=10.0,
    title::String="polar_with_cm",
    fig_size::Tuple{Int,Int}=(1200, 800),
    angle_xlim::Tuple{Real,Real}=(-10, 40)
)
    total_cases = length(body_aero_list) + length(literature_path_list)
    length(label_list) == total_cases || throw(ArgumentError("labels length ($(length(label_list))) must match number of cases ($total_cases)"))
    length(solver_list) == length(body_aero_list) || throw(ArgumentError("solver_list length must match body_aero_list length"))

    polar_data_list = Vector{Any}()
    labels_full = String[]

    # Computational cases
    for (solver_i, body, lbl) in zip(solver_list, body_aero_list, label_list[1:length(solver_list)])
        pd = compute_polar_with_cm(
            solver_i,
            body,
            angle_range;
            angle_type=angle_type,
            angle_of_attack=angle_of_attack,
            side_slip=side_slip,
            v_a=v_a
        )
        @info "polar sample (solver)" label=lbl first_cl=pd.cl[1] first_cd=pd.cd[1] first_cs=pd.cs[1] first_cmy=pd.cmy[1]
        push!(polar_data_list, pd)
        re_tag = round(Int, first(pd.rey) * 1e-5)
        push!(labels_full, "$(lbl) Re=$(re_tag)e5")
    end

    # Literature cases
    for (path, lbl) in zip(literature_path_list, label_list[length(solver_list)+1:end])
        data = readdlm(path, ',')
        header_raw = string.(data[1, :])
        header = lowercase.(strip.(header_raw))
        alpha_idx = findfirst(x -> occursin("alpha", x) || occursin("aoa", x), header)
        cl_idx    = findfirst(x -> occursin("cl", x), header)
        cd_idx    = findfirst(x -> occursin("cd", x), header)
        cs_idx    = findfirst(x -> occursin("cs", x), header)
        cmy_idx   = findfirst(x -> occursin("cmy", x), header)

        parse_col(col) = begin
            vals = Float64[]
            for v in col
                if v isa Real
                    push!(vals, Float64(v))
                else
                    s = strip(String(v))
                    y = tryparse(Float64, s)
                    push!(vals, isnothing(y) ? NaN : y)
                end
            end
            vals
        end

        alpha_col = parse_col(data[2:end, alpha_idx])
        cl_col    = parse_col(data[2:end, cl_idx])
        cd_col    = parse_col(data[2:end, cd_idx])
        cs_col    = cs_idx === nothing ? zeros(size(data, 1)-1) : parse_col(data[2:end, cs_idx])
        cmy_col   = cmy_idx === nothing ? nothing : parse_col(data[2:end, cmy_idx])

        push!(polar_data_list, (angle=alpha_col, cl=cl_col, cd=cd_col, cs=cs_col, cmy=cmy_col, rey=fill(NaN, length(alpha_col))))
        push!(labels_full, lbl)
    end

    fig = Figure(size=fig_size)
    Label(fig[0, :], title; fontsize=20, font=:bold)

    ax_cl    = Axis(fig[1, 1], title="CL vs $angle_type [deg]",  xlabel="$angle_type [deg]", ylabel="CL [-]")
    ax_cd    = Axis(fig[1, 2], title="CD vs $angle_type [deg]",  xlabel="$angle_type [deg]", ylabel="CD [-]")
    ax_cmy   = Axis(fig[1, 3], title="CMy vs $angle_type [deg]", xlabel="$angle_type [deg]", ylabel="CMy [-]")
    ax_cs    = Axis(fig[2, 1], title="CS vs $angle_type [deg]",  xlabel="$angle_type [deg]", ylabel="CS [-]")
    ax_polar = Axis(fig[2, 2], title="CL vs CD", xlabel="CD [-]", ylabel="CL [-]")

    xlims!(ax_cl, angle_xlim...)
    xlims!(ax_cd, angle_xlim...)
    xlims!(ax_cmy, angle_xlim...)
    xlims!(ax_cs, angle_xlim...)

    colors = Makie.wong_colors()

    for (idx, (pd, lbl)) in enumerate(zip(polar_data_list, labels_full))
        color = colors[mod1(idx, length(colors))]

        lines!(ax_cl, pd.angle, pd.cl; color, label=lbl)
        scatter!(ax_cl, pd.angle, pd.cl; color, markersize=6)

        lines!(ax_cd, pd.angle, pd.cd; color, label=lbl)
        scatter!(ax_cd, pd.angle, pd.cd; color, markersize=6)

        lines!(ax_cs, pd.angle, pd.cs; color, label=lbl)
        scatter!(ax_cs, pd.angle, pd.cs; color, markersize=6)

        lines!(ax_polar, pd.cd, pd.cl; color, label=lbl)
        scatter!(ax_polar, pd.cd, pd.cl; color, markersize=6)

        if pd.cmy !== nothing
            cmy_vals = Float64.(pd.cmy)
            if !all(isnan, cmy_vals)
                lines!(ax_cmy, pd.angle, cmy_vals; color, label=lbl)
                scatter!(ax_cmy, pd.angle, cmy_vals; color, markersize=6)
            end
        end
    end

    Legend(fig[2, 3], ax_cl)
    return fig
end


fig2 = plot_polars_with_cmy(
    [solver],
    [body_aero],
    labels;
    literature_path_list=literature_paths,
    angle_range=range(-5, 40, step=1),
    angle_type="angle_of_attack",
    angle_of_attack=angle_of_attack_deg,
    side_slip=sideslip_deg,
    v_a=wind_speed,
    title="Testing solve! on V3 kite",
    angle_xlim=(-5, 40)
)
scr2 = display(fig2)
wait(scr2)


nothing
