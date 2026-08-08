#~=======================================================~=
## Installing ----
#~=======================================================~=
# Making sure to install everything that is needed (but without loading it)
# because something we only use one function so we want to avoid overwriting
# and taking up too much space.

# Here's a function that will install any missing library.
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (length(find.package(pkg, quiet = TRUE)) == 0) {
      install.packages(pkg)
    }
  }
}

#~=======================================================~=
## Loading the packages ----
#~=======================================================~=

# Load a package without its attach banner or "masked from" report.
# NOTE: order matters, it decides which package masks which.
library_quietly <- function(...) {
  for (pkg in c(...)) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

# The attach messages are suppressed, so use this to check which package
# masks which function when something behaves unexpectedly.
show_conflicts <- function() {
  conflicts(detail = TRUE)
}

#~=======================================================~=
## Models ----
#~=======================================================~=

# Function to safely catch a list of models, without breaking if one of them fails
safe_model_list <- function(...) {
  args <- as.list(match.call(expand.dots = TRUE))[-1]
  caller <- parent.frame()
  result <- list()
  for (nm in names(args)) {
    val <- tryCatch(eval(args[[nm]], envir = caller), error = function(e) NULL)
    if (!is.null(val)) result[[nm]] <- val
  }
  result
}

#~=======================================================~=
## Tables ----
#~=======================================================~=

# Default theme for flextables.
# Set it with set_flextable_defaults(theme_fun = flextable_apa_theme)
flextable_apa_theme <- function(x) {
  apa.border <- fp_border(color = "black", width = 1.5, style = "solid")

  x %>%
    hline_top(part = "head", border = apa.border) %>%
    hline_bottom(part = "head", border = apa.border) %>%
    hline_top(part = "body", border = apa.border) %>%
    hline_bottom(part = "body", border = apa.border) %>%
    padding(padding.top = 2, padding.bottom = 2, padding.left = 0, padding.right = 0, part = "all") %>%
    line_spacing(space = 1, part = "all") %>%
    valign(valign = "center", part = "all")
}
