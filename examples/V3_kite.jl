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
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results","python_alpha_sweep.csv"),
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results","windtunnel_alpha_sweep_beta_00_0_Poland_2025_Rey_5e5.csv"),
    ]
labels= [
    "Julia VSM 2D CFD PCHIP NF",
    "CFD RANS Re=5e5",
    "CFD RANS Re=10e5 (With Struts)",
    "Python VSM 2D CFD PCHIP NF Re=5e5",
    "Wind Tunnel Re=5e5 (With Struts)"
    ]
beta_literature_paths = [
    joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results", "windtunnel_beta_sweep_alpha_07_4_Poland_2025_Rey_5e5.csv"),
    # joinpath(project_dir, "data", "TUDELFT_V3_KITE", "literature_results", "windtunnel_beta_sweep_alpha_12_5_Poland_2025_Rey_5e5.csv"),
]
beta_labels = [
    labels[1],
    "Wind Tunnel Re=5e5 beta sweep alpha=7.4",
    # "Wind Tunnel Re=5e5 beta sweep alpha=12.5",
]

# Load YAML settings directly (avoids VSMSettings(filename) conversion issues in live sessions)
settings_path = joinpath(project_dir, "data", "TUDELFT_V3_KITE", "vsm_settings.yaml")
settings_data = VortexStepMethod.YAML.load_file(settings_path)
condition_cfg = settings_data["condition"]
wing_cfg = settings_data["wings"][1]
solver_cfg = settings_data["solver_settings"]

# Create wing, body_aero, and solver objects using settings

wing = Wing(
    joinpath(project_dir, wing_cfg["geometry_file"]);
    n_panels=wing_cfg["n_panels"],
    spanwise_distribution=getproperty(VortexStepMethod, Symbol(wing_cfg["spanwise_panel_distribution"])),
    spanwise_direction=Float64.(wing_cfg["spanwise_direction"]),
    remove_nan=wing_cfg["remove_nan"],
)
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


# Construct Solver using keyword arguments from solver settings
solver = Solver(body_aero;
    solver_type=(solver_cfg["solver_type"] == "NONLIN" ? NONLIN : LOOP),
    aerodynamic_model_type=getproperty(VortexStepMethod, Symbol(solver_cfg["aerodynamic_model_type"])),
    density=solver_cfg["density"],
    max_iterations=solver_cfg["max_iterations"],
    rtol=solver_cfg["rtol"],
    tol_reference_error=solver_cfg["tol_reference_error"],
    relaxation_factor=solver_cfg["relaxation_factor"],
    is_with_artificial_damping=solver_cfg["artificial_damping"],
    artificial_damping=(k2=solver_cfg["k2"], k4=solver_cfg["k4"]),
    type_initial_gamma_distribution=getproperty(VortexStepMethod, Symbol(solver_cfg["type_initial_gamma_distribution"])),
    use_gamma_prev=get(solver_cfg, "use_gamma_prev", get(solver_cfg, "use_gamme_prev", true)),
    core_radius_fraction=solver_cfg["core_radius_fraction"],
    mu=solver_cfg["mu"],
    is_only_f_and_gamma_output=get(solver_cfg, "calc_only_f_and_gamma", false),
    correct_aoa=get(solver_cfg, "correct_aoa", false),
    reference_point=get(solver_cfg, "reference_point", [0.422646, 0.0, 9.3667]), #ref-point from WT exp.
)

# Extract values for plotting (optional - for reference)
wind_speed = condition_cfg["wind_speed"]
angle_of_attack_deg = condition_cfg["alpha"]
sideslip_deg = condition_cfg["beta"]
yaw_rate = condition_cfg["yaw_rate"]

# Set flight conditions from settings
α0 = deg2rad(angle_of_attack_deg)
β0 = deg2rad(sideslip_deg)
set_va!(body_aero, wind_speed .* [cos(α0) * cos(β0), sin(β0), sin(α0) * cos(β0)])

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
isinteractive() && wait(scr1)


# Polar sweep including force and moment coefficients, with selectable solve or solve! path
function compute_polar_input(
    solver,
    body_aero,
    angle_range;
    angle_type::String="angle_of_attack",
    angle_of_attack::Float64=0.0,
    side_slip::Float64=0.0,
    v_a::Float64=10.0,
    use_solve::Bool=false
)
    n_angles = length(angle_range)
    cl = zeros(n_angles)
    cd = zeros(n_angles)
    cs = zeros(n_angles)
    cmx = fill(NaN, n_angles)
    cmy = fill(NaN, n_angles)
    cmz = fill(NaN, n_angles)
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
        if use_solve
            results_local = solve(solver, body_aero, gamma_prev; log=false)
            gamma_prev = copy(solver.lr.gamma_new)
        else
            solve!(solver, body_aero, gamma_prev; log=false)
            gamma_prev = solver.sol.gamma_distribution

            results_local = calculate_results(
                body_aero,
                solver.lr.gamma_new,
                solver.reference_point,
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
        end

        cl[i] = results_local["cl"]
        cd[i] = results_local["cd"]
        cs[i] = results_local["cs"]
        cmx[i] = get(results_local, "cmx", NaN)
        cmy[i] = get(results_local, "cmy", NaN)
        cmz[i] = get(results_local, "cmz", NaN)
        reynolds_number[i] = results_local["Rey"]
    end

    return (angle=angle_range, cl=cl, cd=cd, cs=cs, cmx=cmx, cmy=cmy, cmz=cmz, rey=reynolds_number)
end

function plot_polars(
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
    angle_xlim::Tuple{Real,Real}=(-5, 20),
    use_solve::Bool=false
)
    total_cases = length(body_aero_list) + length(literature_path_list)
    length(label_list) == total_cases || throw(ArgumentError("labels length ($(length(label_list))) must match number of cases ($total_cases)"))
    length(solver_list) == length(body_aero_list) || throw(ArgumentError("solver_list length must match body_aero_list length"))

    polar_data_list = Vector{Any}()
    labels_full = String[]

    # Computational cases
    for (solver_i, body, lbl) in zip(solver_list, body_aero_list, label_list[1:length(solver_list)])
        pd = compute_polar_input(
            solver_i,
            body,
            angle_range;
            angle_type=angle_type,
            angle_of_attack=angle_of_attack,
            side_slip=side_slip,
            v_a=v_a,
            use_solve=use_solve
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
        angle_idx = if angle_type == "angle_of_attack"
            findfirst(x -> occursin("alpha", x) || occursin("aoa", x), header)
        elseif angle_type == "side_slip"
            findfirst(x -> occursin("beta", x) || occursin("side_slip", x) || occursin("sideslip", x), header)
        else
            throw(ArgumentError("angle_type must be 'angle_of_attack' or 'side_slip'"))
        end
        isnothing(angle_idx) && throw(ArgumentError(
            "Could not find angle column for angle_type='$angle_type' in literature file: $path"
        ))
        cl_idx    = findfirst(x -> occursin("cl", x), header)
        cd_idx    = findfirst(x -> occursin("cd", x), header)
        cs_idx    = findfirst(x -> occursin("cs", x), header)
        cmx_idx   = findfirst(x -> occursin("cmx", x), header)
        cmy_idx   = findfirst(x -> occursin("cmy", x), header)
        cmz_idx   = findfirst(x -> occursin("cmz", x), header)

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

        angle_col = parse_col(data[2:end, angle_idx])
        cl_col    = parse_col(data[2:end, cl_idx])
        cd_col    = parse_col(data[2:end, cd_idx])
        cs_col    = cs_idx === nothing ? zeros(size(data, 1)-1) : parse_col(data[2:end, cs_idx])
        cmx_col   = cmx_idx === nothing ? fill(NaN, size(data, 1)-1) : parse_col(data[2:end, cmx_idx])
        cmy_col   = cmy_idx === nothing ? fill(NaN, size(data, 1)-1) : parse_col(data[2:end, cmy_idx])
        cmz_col   = cmz_idx === nothing ? fill(NaN, size(data, 1)-1) : parse_col(data[2:end, cmz_idx])

        push!(polar_data_list, (angle=angle_col, cl=cl_col, cd=cd_col, cs=cs_col, cmx=cmx_col, cmy=cmy_col, cmz=cmz_col, rey=fill(NaN, length(angle_col))))
        push!(labels_full, lbl)
    end

    fig = Figure(size=fig_size)
    Label(fig[0, :], title; fontsize=20, font=:bold)

    ax_cl  = Axis(fig[1, 1], title="CL vs $angle_type [deg]",  xlabel="$angle_type [deg]", ylabel="CL [-]")
    ax_cd  = Axis(fig[1, 2], title="CD vs $angle_type [deg]",  xlabel="$angle_type [deg]", ylabel="CD [-]")
    ax_cs  = Axis(fig[1, 3], title="CS vs $angle_type [deg]",  xlabel="$angle_type [deg]", ylabel="CS [-]")
    ax_cmx = Axis(fig[2, 1], title="CMx vs $angle_type [deg]", xlabel="$angle_type [deg]", ylabel="CMx [-]")
    ax_cmy = Axis(fig[2, 2], title="CMy vs $angle_type [deg]", xlabel="$angle_type [deg]", ylabel="CMy [-]")
    ax_cmz = Axis(fig[2, 3], title="CMz vs $angle_type [deg]", xlabel="$angle_type [deg]", ylabel="CMz [-]")

    # Force truly equal columns independent of axis decoration sizes.
    colsize!(fig.layout, 1, Relative(1/3))
    colsize!(fig.layout, 2, Relative(1/3))
    colsize!(fig.layout, 3, Relative(1/3))

    # Keep the inner plotting rectangles aligned as well.
    for ax in (ax_cl, ax_cd, ax_cs, ax_cmx, ax_cmy, ax_cmz)
        ax.yticklabelspace = 36.0
        ax.xticklabelspace = 24.0
    end

    xlims!(ax_cl, angle_xlim...)
    xlims!(ax_cd, angle_xlim...)
    xlims!(ax_cs, angle_xlim...)
    xlims!(ax_cmx, angle_xlim...)
    xlims!(ax_cmy, angle_xlim...)
    xlims!(ax_cmz, angle_xlim...)

    colors = Makie.wong_colors()

    for (idx, (pd, lbl)) in enumerate(zip(polar_data_list, labels_full))
        color = colors[mod1(idx, length(colors))]

        lines!(ax_cl, pd.angle, pd.cl; color, label=lbl)
        scatter!(ax_cl, pd.angle, pd.cl; color, markersize=6)

        lines!(ax_cd, pd.angle, pd.cd; color, label=lbl)
        scatter!(ax_cd, pd.angle, pd.cd; color, markersize=6)

        lines!(ax_cs, pd.angle, pd.cs; color, label=lbl)
        scatter!(ax_cs, pd.angle, pd.cs; color, markersize=6)

        cmx_vals = Float64.(pd.cmx)
        cmy_vals = Float64.(pd.cmy)
        cmz_vals = Float64.(pd.cmz)
        if !all(isnan, cmx_vals)
            lines!(ax_cmx, pd.angle, cmx_vals; color, label=lbl)
            scatter!(ax_cmx, pd.angle, cmx_vals; color, markersize=6)
        end
        if !all(isnan, cmy_vals)
            lines!(ax_cmy, pd.angle, cmy_vals; color, label=lbl)
            scatter!(ax_cmy, pd.angle, cmy_vals; color, markersize=6)
        end
        if !all(isnan, cmz_vals)
            lines!(ax_cmz, pd.angle, cmz_vals; color, label=lbl)
            scatter!(ax_cmz, pd.angle, cmz_vals; color, markersize=6)
        end
    end

    Legend(fig[3, 1:3], ax_cl; orientation=:horizontal, tellwidth=false)
    return fig
end

fig2 = plot_polars(
    [solver],
    [body_aero],
    labels;
    literature_path_list=literature_paths,
    angle_range=range(-5, 20, step=1),
    angle_type="angle_of_attack",
    angle_of_attack=angle_of_attack_deg,
    side_slip=sideslip_deg,
    v_a=wind_speed,
    title="solve!",
    angle_xlim=(-5, 20)
)
scr2 = display(fig2)
isinteractive() && wait(scr2)

# Construct Solver using keyword arguments from solver settings
solver_cfg["solver_type"] = "LOOP"
solver = Solver(body_aero;
    solver_type=(solver_cfg["solver_type"] == "NONLIN" ? NONLIN : LOOP),
    aerodynamic_model_type=getproperty(VortexStepMethod, Symbol(solver_cfg["aerodynamic_model_type"])),
    density=solver_cfg["density"],
    max_iterations=solver_cfg["max_iterations"],
    rtol=solver_cfg["rtol"],
    tol_reference_error=solver_cfg["tol_reference_error"],
    relaxation_factor=solver_cfg["relaxation_factor"],
    is_with_artificial_damping=solver_cfg["artificial_damping"],
    artificial_damping=(k2=solver_cfg["k2"], k4=solver_cfg["k4"]),
    type_initial_gamma_distribution=getproperty(VortexStepMethod, Symbol(solver_cfg["type_initial_gamma_distribution"])),
    use_gamma_prev=get(solver_cfg, "use_gamma_prev", get(solver_cfg, "use_gamme_prev", true)),
    core_radius_fraction=solver_cfg["core_radius_fraction"],
    mu=solver_cfg["mu"],
    is_only_f_and_gamma_output=get(solver_cfg, "calc_only_f_and_gamma", false),
    correct_aoa=get(solver_cfg, "correct_aoa", false),
    reference_point=get(solver_cfg, "reference_point", [0.422646, 0.0, 9.3667])
)

fig3 = plot_polars(
    [solver],
    [body_aero],
    labels;
    literature_path_list=literature_paths,
    angle_range=range(-5, 20, step=1),
    angle_type="angle_of_attack",
    angle_of_attack=angle_of_attack_deg,
    side_slip=sideslip_deg,
    v_a=wind_speed,
    title="solve",
    angle_xlim=(-5, 20),
    use_solve=true
)
scr3 = display(fig3)
isinteractive() && wait(scr3)

fig4 = plot_polars(
    [solver],
    [body_aero],
    beta_labels;
    literature_path_list=beta_literature_paths,
    angle_range=range(0, 12, step=1),
    angle_type="side_slip",
    angle_of_attack=angle_of_attack_deg,
    side_slip=sideslip_deg,
    v_a=wind_speed,
    title="solve! beta sweep",
    angle_xlim=(0, 12),
)
scr4 = display(fig4)
isinteractive() && wait(scr4)

fig5 = plot_polars(
    [solver],
    [body_aero],
    beta_labels;
    literature_path_list=beta_literature_paths,
    angle_range=range(0, 12, step=1),
    angle_type="side_slip",
    angle_of_attack=angle_of_attack_deg,
    side_slip=sideslip_deg,
    v_a=wind_speed,
    title="solve beta sweep",
    angle_xlim=(0, 12),
    use_solve=true
)
scr5 = display(fig5)
isinteractive() && wait(scr5)


nothing
