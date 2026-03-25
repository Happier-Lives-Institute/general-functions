# Combine cowplot and a white background
theme_hli_wbg <- function() {
  theme_cowplot() +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

# Function to save a plot as both PNG and SVG
hli_double_save <- function(filename_no_end, plot, width, height, dpi, units = "in") {

  # --- PNG ---
  ggsave(
    filename = paste0(filename_no_end, ".png"),
    plot     = plot,
    width    = width,
    height   = height,
    dpi      = dpi,
    units    = units
  )

  # --- SVG ---
  svg_path <- paste0(filename_no_end, ".svg")

  ggsave(
    filename = svg_path,
    plot     = plot,
  )

  svg_string <- readChar(svg_path, file.info(svg_path)$size)

  # Remove CDATA wrapper (breaks inline SVG in WordPress)
  svg_string <- gsub("<![CDATA[", "", svg_string, fixed = TRUE)
  svg_string <- gsub("]]>",       "", svg_string, fixed = TRUE)

  # Remove XML declaration (breaks inline SVG embedding)
  svg_string <- gsub("^<\\?xml[^\\?]*\\?>\\s*", "", svg_string, perl = TRUE)

  # Make responsive: width=100%, drop fixed height (viewBox preserves aspect ratio)
  svg_string <- gsub("(<svg[^>]*) width='[^']*'",  "\\1 width='100%'", svg_string, perl = TRUE)
  svg_string <- gsub("(<svg[^>]*) height='[^']*'", "\\1",              svg_string, perl = TRUE)

  writeChar(svg_string, svg_path, eos = NULL)
}