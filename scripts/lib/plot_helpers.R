save_gg_panel <- function(plot_object, stem, width, height) {
  png_path <- file.path(output_dir, paste0(stem, ".png"))
  tiff_path <- file.path(output_dir, paste0(stem, ".tiff"))
  ggsave(
    png_path, plot = plot_object, width = width, height = height,
    units = "in", dpi = 400, bg = "white", limitsize = FALSE
  )
  ggsave(
    tiff_path, plot = plot_object, width = width, height = height,
    units = "in", dpi = 600, device = "tiff", compression = "lzw",
    bg = "white", limitsize = FALSE
  )
  c(png = png_path, tiff = tiff_path)
}

open_png_device <- function(filename, width, height) {
  png_args <- list(
    filename = filename, width = width, height = height,
    units = "in", res = 400, bg = "white"
  )
  if (capabilities("cairo")) png_args$type <- "cairo-png"
  do.call(grDevices::png, png_args)
}

save_grid_panel <- function(draw_function, stem, width, height) {
  png_path <- file.path(output_dir, paste0(stem, ".png"))
  tiff_path <- file.path(output_dir, paste0(stem, ".tiff"))
  open_png_device(png_path, width, height)
  draw_function()
  invisible(dev.off())
  grDevices::tiff(
    tiff_path, width = width, height = height, units = "in",
    res = 600, compression = "lzw", bg = "white"
  )
  draw_function()
  invisible(dev.off())
  c(png = png_path, tiff = tiff_path)
}

significance_symbol <- function(q_value) {
  ifelse(
    is.na(q_value), "ns",
    ifelse(q_value < 0.001, "***",
      ifelse(q_value < 0.01, "**", ifelse(q_value < 0.05, "*", "ns"))
    )
  )
}


