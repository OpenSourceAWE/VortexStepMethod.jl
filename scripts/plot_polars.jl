function plot_values(alphas, d_trailing_edge_angles, matrix, interp, name)
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    X_data = collect(d_trailing_edge_angles) .+ zeros(length(alphas))'
    Y_data = collect(alphas)' .+ zeros(length(d_trailing_edge_angles))

    matrix = Matrix{Float64}(matrix)
    interp_matrix = zeros(size(matrix)...)
    int_alphas, int_d_trailing_edge_angles = alphas .+ deg2rad(0.5), d_trailing_edge_angles .+ deg2rad(0.5)
    interp_matrix .= [interp(alpha, d_trailing_edge_angle) for alpha in int_alphas, d_trailing_edge_angle in int_d_trailing_edge_angles]
    X_int = collect(int_d_trailing_edge_angles) .+ zeros(length(int_alphas))'
    Y_int = collect(int_alphas)' .+ zeros(length(int_d_trailing_edge_angles))

    ax.plot_wireframe(X_data, Y_data, matrix, edgecolor="royalblue", lw=0.5, rstride=5, cstride=5, alpha=0.6)
    ax.plot_wireframe(X_int, Y_int, interp_matrix, edgecolor="orange", lw=0.5, rstride=5, cstride=5, alpha=0.6)
    plt.xlabel("Alpha")
    plt.ylabel("Flap angle")
    plt.zlabel("$name values")
    plt.title("$name for different d_flap and angle")
    plt.legend()
    plt.grid(true)
    plt.show()
end

cl_interp = extrapolate(scale(interpolate(cl_matrix, BSpline(Linear())), alphas, d_trailing_edge_angles), NaN)
cd_interp = extrapolate(scale(interpolate(cd_matrix, BSpline(Linear())), alphas, d_trailing_edge_angles), NaN)
cm_interp = extrapolate(scale(interpolate(cm_matrix, BSpline(Linear())), alphas, d_trailing_edge_angles), NaN)

plot_values(alphas, d_trailing_edge_angles, cl_matrix, cl_interp, "Cl")
plot_values(alphas, d_trailing_edge_angles, cd_matrix, cd_interp, "Cd")
plot_values(alphas, d_trailing_edge_angles, cm_matrix, cm_interp, "Cm")