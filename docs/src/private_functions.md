```@meta
CurrentModule = VortexStepMethod
```

## Private Functions

### Solver, forces and circulation
```@docs
calculate_AIC_matrices!
gamma_loop!
build_spanwise_laplacian!
local_lift_slope!
apply_artificial_viscosity!
frozen_wake!
calc_forces!
calculate_cl
calculate_cd
calculate_cm
calculate_cd_cm
set_pitch_rate_dist!
calculate_relative_alpha_and_velocity
calculate_relative_alpha_and_relative_velocity
update_effective_angle_of_attack!
calculate_stall_angle_list
wing_span_flip
calculate_circulation_distribution_elliptical_wing
_compute_reference_velocity_from_distribution
smooth_circulation!
smooth_distribution!
make_dual_shadow
```

### Panel aerodynamics
The per-panel aerodynamics, written once as pure, branch-free functions of the
section geometry and the flow. `SymbolicAWEModels` traces the same functions with
symbolic arguments to build its equations, so both packages evaluate one
definition of the physics.
```@docs
SMOOTH_FLOOR
smooth_norm
panel_span_vector
panel_chord_weight
panel_axes
effective_alpha
panel_inflow
dynamic_pressure
flow_curvature_cm
panel_force_directions
panel_moment
panel_couple_force
panel_loads
```

### Induced velocities
```@docs
velocity_3D_bound_vortex!
velocity_3D_trailing_vortex!
velocity_3D_trailing_vortex_semiinfinite!
calculate_velocity_induced_bound_2D!
calculate_velocity_induced_single_ring_semiinfinite!
cross3!
```

### Panels, filaments and geometry updates
```@docs
update_panel_properties!
build_interps
panel_interp_types
reinit!(wing::AbstractWing)
reinit!(panel::Panel, section_1::Section, section_2::Section, aero_center, control_point, bound_point_1, bound_point_2, x_airf, y_airf, z_airf, delta, vec)
rotated_te
calculate_filaments_for_plotting
```

### Mesh refinement and billowing
```@docs
unrefined_deform!
deform!
compute_refined_panel_mapping!
compute_refined_section_interpolation!
copy_sections
copy_sections_to_refined!
_apply_refined_section_thetas!
_panel_thetas_to_section_thetas!
_interpolate_unrefined_to_refined
span_order_key
normalize_span_order!
can_reuse_prior_refined_surface_tables
refine_mesh_for_linear_cosine_distribution!
refine_mesh_by_splitting_provided_sections!
refine_mesh_with_billowing!
apply_billowing_to_pair!
billowing_arc_length
flip_created_coord_in_pairs_if_needed!
update_non_deformed_sections!
```

### Aerodynamic data and Cp
```@docs
calculate_new_aero_data
set_polar!
rebuild_polar
window_alpha
assemble_polar_matrix
load_matrix_polar_data
read_aero_matrix
read_dat
read_node_table
write_node_rows
convert_node_table
delta_suffix
interpolate_matrix_nans!
remove_vector_nans
generate_polar_data
extract_literature_polar_data
parse_literature_column
interpolate_section_aero_to_refined!
validate_section_aero
```

### Examples
```@docs
copy_examples
```

## Airfoil aerodynamics (AirfoilAero)
```@meta
CurrentModule = VortexStepMethod.AirfoilAero
```

### Kulfan CST parametrization
```@docs
bernstein_basis
class_function
leading_edge_basis
compute_optimal_x_points
normalize_airfoil
get_lower_upper
turn_trailing_edge!
```

### Shrink-wrap distance field
```@docs
distance_parabolas!
squared_distance_transform!
grid_sampler
flood_outside
trace_level_set
largest_linking_gap
resample_arc
smooth_turning!
```

### NeuralFoil network
```@docs
load_neuralfoil_model
neuralfoil_section
neuralfoil_fused_output
fused_output
decode_surface_velocity
decode_coefficients
nn_forward
prepare_inputs
fill_case_input!
flip_inputs
flip_outputs
squared_mahalanobis_distance
swish
sigmoid
```

### Polars and airfoil IO
```@docs
create_2d_polars
lei_poly_coeffs
resolve_airfoil
read_dat_coordinates
write_dat
write_polar_csv
write_polar_matrix_csv
write_aero_matrix
write_node_table
flat_plate_cf
neuralfoil_contour_solution
contour_arc
arc_at_chord
trailing_edge_speed
velocity_knots
contour_pressure
fill_node_nans!
live_xfoil_solver
xfoil_ramp
```

## OBJ mesh conversion (ObjAdapter)
```@meta
CurrentModule = VortexStepMethod.ObjAdapter
```
```@docs
read_faces
slice_mesh_at_plane
order_segments_to_contour
station_indices
airfoil_frame
build_section
contour_to_airfoil
plane_contour_to_airfoil
reorder_airfoil_selig
densify_contour
create_interpolations
find_circle_center_and_radius
march_edges
calculate_inertia_tensor
center_to_com!
airfoils_from_yaml
write_geometry_yaml
resolve_aero_geometry
table_path_prefix
prefix_table_paths!
plot_airfoil_fit
migrate_node_tables
```

## Makie plotting internals
```@meta
CurrentModule = Base.get_extension(VortexStepMethod, :VortexStepMethodMakieExt)
```
```@docs
display_named
span_axis
create_geometry_plot_makie
plot_line_segment_makie!
set_axes_equal_makie!
map_airfoil_3d
fitted_airfoil_3d
generated_slices
airfoil_skin_geometry
panel_contour
panel_normal
plate_hinge_local
panel_plate_geometry
PLATE_FACES
Makie.plot!(ax, panel::VortexStepMethod.Panel)
Makie.plot!(ax, body::VortexStepMethod.BodyAerodynamics)
Makie.plot!(body::VortexStepMethod.BodyAerodynamics)
```
```@meta
CurrentModule = VortexStepMethod
```
