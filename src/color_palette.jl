"""
Elegant color palette with main and darker versions
"""
const PALETTE = Dict{String, RGB{Float64}}(
    "Blue" => parse(RGB{Float64}, "#0076C2"),
    "Green" => parse(RGB{Float64}, "#4CAF50"),
    "Orange" => parse(RGB{Float64}, "#E0A458"),
    "Red" => parse(RGB{Float64}, "#D32F2F"),
    "Dark Blue" => parse(RGB{Float64}, "#005691"),
    "Dark Green" => parse(RGB{Float64}, "#388E3C"),
    "Dark Orange" => parse(RGB{Float64}, "#D17C1B"),
    "Dark Red" => parse(RGB{Float64}, "#B71C1C"),
    "White" => parse(RGB{Float64}, "#FFFFFF"),
    "Light Gray" => parse(RGB{Float64}, "#A9A9A9"),
    "Dark Gray" => parse(RGB{Float64}, "#696969")
)

"""
    get_color(color_name::String, alpha::Float64=1.0)

Return the RGBA color for the given color name with specified transparency.
"""
function get_color(color_name::String, alpha::Float64=1.0)
    color = get(PALETTE, color_name, RGB(0,0,0))  # Default to black if not found
    return RGBA(color, alpha)
end

"""
    get_color_list()

Return a vector of colors from the palette.
"""
function get_color_list()
    return collect(values(PALETTE))
end

"""
    visualize_palette()

Display the color palette in a plot.
"""
function visualize_palette()
    p = plot(
        legend=false,
        showaxis=false,
        grid=false,
        size=(800, 80)
    )
    
    n_colors = length(PALETTE)
    for (i, (name, color)) in enumerate(PALETTE)
        # Plot color rectangle
        x = [(i-1)/n_colors, i/n_colors]
        rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
        plot!(p, rectangle(1/n_colors, 1, (i-1)/n_colors, 0), 
              color=color, fillalpha=1)
        
        # Add color name
        annotate!(p, (i-0.5)/n_colors, 0.5, 
                 text(name, 8, name == "White" ? :black : :white))
    end
    
    plot!(p, ylims=(0,1), xlims=(0,1))
    display(p)
end

"""
    set_plot_style()

Set the default style for plots using LaTeX.
"""
function set_plot_style()
    # plt.style.use('seaborn-whitegrid')
    # plt.style.use("seaborn-v0_8-whitegrid")
    rcParams = plt.PyDict(plt.matplotlib."rcParams")
    rcParams["text.usetex"] = true
    rcParams["font.family"] = "serif"
    rcParams["font.serif"] = ["Computer Modern Roman"]
    rcParams["axes.titlesize"] = 28
    # rcParams["axes.ymargin"] = 0.1
    rcParams["axes.labelsize"] = 12
    rcParams["axes.linewidth"] = 1
    rcParams["lines.linewidth"] = 1
    rcParams["lines.markersize"] = 6
    rcParams["xtick.labelsize"] = 10
    rcParams["ytick.labelsize"] = 10
    rcParams["legend.fontsize"] = 10
    rcParams["figure.titlesize"] = 16
    rcParams["pgf.texsystem"] = "pdflatex"  # Use pdflatex
    rcParams["pgf.rcfonts"] = false
    rcParams["figure.figsize"] = (10, 6)  # Default figure size
end

# """
#     apply_palette!(p::Plots.Plot, colors::Vector{String})

# Apply the color palette to a plot.
# """
# function apply_palette!(p::Plots.Plot, colors::Vector{String})
#     for (i, series) in enumerate(p.series_list)
#         color_name = colors[mod1(i, length(colors))]
#         series.plotattributes[:linecolor] = get_color(color_name)
#     end
#     return p
# end