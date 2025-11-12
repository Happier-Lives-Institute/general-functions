# Combine cowplot and a white background
theme_hli_wbg <- function() {
  theme_cowplot() +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

# Function to save a plot as both PNG and SVG
hli_double_save <- function(filename_no_end, plot, width, height, dpi, units = "in"){
  ggsave(
    filename = paste0(filename_no_end, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    units = units
  )
  ggsave(
    filename = paste0(filename_no_end, ".svg"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    units = units
  )
}