# Combine cowplot and a white background
# The sizes are optional, only applied when given.
theme_hli_wbg <- function(axis_title_size = NULL, axis_text_size = NULL,
                          legend_title_size = NULL, legend_text_size = NULL) {
  my_theme <- theme_cowplot() +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  if (!is.null(axis_title_size)) {
    my_theme <- my_theme + theme(axis.title = element_text(size = axis_title_size))
  }
  if (!is.null(axis_text_size)) {
    my_theme <- my_theme + theme(axis.text = element_text(size = axis_text_size))
  }
  if (!is.null(legend_title_size)) {
    my_theme <- my_theme + theme(legend.title = element_text(size = legend_title_size))
  }
  if (!is.null(legend_text_size)) {
    my_theme <- my_theme + theme(legend.text = element_text(size = legend_text_size))
  }

  my_theme
}

# Function to save a plot as both PNG and SVG
hli_double_save <- function(filename_no_end, plot, width, height, dpi, 
                            set_svg_same_ratio = FALSE, units = "in",
                            svg_title = NULL, svg_desc = NULL) {

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

  if(set_svg_same_ratio) {
    ggsave(
      filename = svg_path,
      plot     = plot,
      width    = width,
      height   = height,
      units     = units
    )
  } else {
    ggsave(
      filename = svg_path,
      plot     = plot
    )
  }

  svg_string <- readChar(svg_path, file.info(svg_path)$size)

  # Remove CDATA wrapper (breaks inline SVG in WordPress)
  svg_string <- gsub("<![CDATA[", "", svg_string, fixed = TRUE)
  svg_string <- gsub("]]>",       "", svg_string, fixed = TRUE)

  # Remove XML declaration (breaks inline SVG embedding)
  svg_string <- gsub("^<\\?xml[^\\?]*\\?>\\s*", "", svg_string, perl = TRUE)

  # Make responsive: width=100%, drop fixed height (viewBox preserves aspect ratio)
  svg_string <- gsub("(<svg[^>]*) width='[^']*'",  "\\1 width='100%'", svg_string, perl = TRUE)
  svg_string <- gsub("(<svg[^>]*) height='[^']*'", "\\1",              svg_string, perl = TRUE)

  # --- Accessibility ---
  if (!is.null(svg_title)) {
    title_id <- paste0(basename(filename_no_end), "-title")
    desc_id  <- paste0(basename(filename_no_end), "-desc")

    # Build aria-labelledby value
    labelledby <- title_id
    
    # Build nodes to inject
    a11y_nodes <- sprintf('<title id="%s">%s</title>', title_id, svg_title)

    if (!is.null(svg_desc)) {
      labelledby  <- paste(title_id, desc_id)
      a11y_nodes  <- paste0(a11y_nodes, sprintf('<desc id="%s">%s</desc>', desc_id, svg_desc))
    }

    # Add role and aria-labelledby to opening <svg> tag
    svg_string <- gsub("(<svg)([^>]*>)", 
                       sprintf('\\1 role="img" aria-labelledby="%s"\\2', labelledby), 
                       svg_string, perl = TRUE)

    # Inject title (and desc) immediately after the opening <svg ...> tag
    svg_string <- gsub("(<svg[^>]*>)", paste0("\\1", a11y_nodes), svg_string, perl = TRUE)
  }

  writeChar(svg_string, svg_path, eos = NULL)
}