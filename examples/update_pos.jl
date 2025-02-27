using VortexStepMethod
using VortexStepMethod: BoundFilament, MVec3, SemiInfiniteFilament, update_pos!, Panel, Section


x1 = zeros(MVec3)
direction = zeros(MVec3)
vel_mag = zero(Float64)
filament_direction = zero(Int64)
@time filament = SemiInfiniteFilament(x1, direction, vel_mag, filament_direction)
@time update_pos!(filament, x1, direction, vel_mag, filament_direction) 

bound_point_1 = zeros(MVec3)
bound_point_2 = zeros(MVec3)
@time filament = BoundFilament(bound_point_1, bound_point_2)
@time update_pos!(filament, bound_point_1, bound_point_2)

LE_point = zeros(MVec3)
TE_point = zeros(MVec3)
aero_input = :inviscid
@time section = Section(LE_point, TE_point, aero_input)
@time update_pos!(section, LE_point, TE_point)

aero_center = zeros(MVec3)
control_point = zeros(MVec3)
bound_point_1 = zeros(MVec3)
bound_point_2 = zeros(MVec3)
x_airf = zeros(MVec3)
y_airf = zeros(MVec3)
z_airf = zeros(MVec3)
@time panel = Panel(section, section, aero_center, control_point, bound_point_1, bound_point_2, x_airf, y_airf, z_airf)
@time update_pos!(panel, section, section, aero_center, control_point, bound_point_1, bound_point_2, x_airf, y_airf, z_airf)

