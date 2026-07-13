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
calculate_relative_alpha_and_velocity
calculate_relative_alpha_and_relative_velocity
update_effective_angle_of_attack!
calculate_stall_angle_list
calculate_circulation_distribution_elliptical_wing
_compute_reference_velocity_from_distribution
smooth_circulation!
smooth_distribution!
make_dual_shadow
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
rotated_te
calculate_filaments_for_plotting
```

### Mesh refinement and billowing
```@docs
unrefined_deform!
deform!
compute_refined_panel_mapping!
compute_refined_section_interpolation!
copy_sections_to_refined!
_apply_refined_section_thetas!
_panel_thetas_to_section_thetas!
_interpolate_unrefined_to_refined
_section_sort_key
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
assemble_polar_matrix
load_matrix_polar_data
read_aero_matrix
write_aero_matrix
interpolate_matrix_nans!
remove_vector_nans
generate_polar_data
extract_literature_polar_data
parse_literature_column
interpolate_cp_to_refined!
prepare_cp_output!
validate_cp_sections
```

### Examples and plotting backend
```@docs
copy_examples
set_plot_backend!
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
resample_x
smooth_turning!
split_surfaces
```

### NeuralFoil network
```@docs
load_neuralfoil_model
neuralfoil_section
neuralfoil_fused_output
nn_forward
prepare_inputs
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
fill_slice_nans!
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
plot_airfoil_fit
```

## Makie plotting internals
```@meta
CurrentModule = Base.get_extension(VortexStepMethod, :VortexStepMethodMakieExt)
```
```@docs
create_geometry_plot_makie
plot_line_segment_makie!
set_axes_equal_makie!
map_airfoil_3d
fitted_airfoil_3d
generated_slices
Makie.plot!(ax, panel::VortexStepMethod.Panel)
Makie.plot!(ax, body::VortexStepMethod.BodyAerodynamics)
Makie.plot!(body::VortexStepMethod.BodyAerodynamics)
```
```@meta
CurrentModule = VortexStepMethod
```
