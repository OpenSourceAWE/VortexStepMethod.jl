"""
    set_plot_style()

Set the default style for plots using LaTeX.
"""
function set_plot_style(titel_size=16)
    rcParams = plt.PyDict(plt.matplotlib."rcParams")
    rcParams["text.usetex"] = true
    rcParams["font.family"] = "serif"
    rcParams["font.serif"] = ["Computer Modern Roman"]
    rcParams["axes.titlesize"] = titel_size
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
