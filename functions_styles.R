# Combine cowplot and a white background
theme_hli_wbg <- function() {
  theme_cowplot() +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}
