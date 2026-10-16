### Fixed

- Under GLMakie, a plot function called with `is_save=true, is_show=true` no longer throws `GLMakie can not display a scene in multiple Screens`: `save_plot` closes the hidden screen the save opened.
