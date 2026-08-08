#~############################################################################~#
# Naming ----
#~############################################################################~#

SPLINE_PATTERN <- "rcs\\(|ns\\(|bs\\(|pspline\\("

# a function to extract variable names from moderator formulas
# it only works for models with up to 1 interaction for now.
extract_vars_flex <- function(formula) {
  terms <- all.vars(formula)
  formula_str <- as.character(formula)[2]
  if (grepl("\\*", formula_str)) {
    my_interaction <- unlist(strsplit(formula_str, " *\\* *"))
    my_interactions_var <- paste0(my_interaction, collapse = ":")
    terms <- c(terms, my_interactions_var)
  }
  return(unique(terms))
}

# Variable a model's spline is fitted on ("" if the model has no spline term).
# trimws because rcs(follow_up_years , 3) leaves a trailing space.
spline_var_of <- function(m) {
  spline_rows <- rownames(m$b)[grepl(SPLINE_PATTERN, rownames(m$b))]
  if (length(spline_rows) == 0) return("")
  trimws(sub("^[a-z]+\\(\\s*([^,)]+).*", "\\1", spline_rows[1]))
}

# Raw name for each coefficient of a model, before labelling.
coefficient_names <- function(m) {

  # spline_seq[ii] = position within the spline block (1 = linear term,
  # 2 = first non-linear term, ...), or NA for non-spline coefficients.
  spline_rows  <- grepl(SPLINE_PATTERN, rownames(m$b))
  spline_seq   <- match(seq_len(nrow(m$b)), which(spline_rows))
  n_parameters <- nrow(m$b)

  # A cell-means model (`~ 0 + moderator`) has no intercept, and its
  # coefficients name a *level* rather than a contrast against a base level.
  # Marking them with " = " keeps the labelling function from labelling them as
  # contrasts ("CBT vs Trauma" for what is simply "Trauma").
  separator <- if("intrcpt" %in% rownames(m$b)) " " else " = "

  vapply(seq_len(n_parameters), function(ii) {

    if(!is.na(spline_seq[ii])) {

      spline_var <- trimws(sub("^[a-z]+\\(\\s*([^,)]+).*", "\\1", rownames(m$b)[ii]))
      if(spline_seq[ii] == 1) {
        paste0(spline_var, "_spline_linear")
      } else {
        paste0(spline_var, "_spline_nl", spline_seq[ii] - 1)
      }

    } else if(rownames(m$b)[ii] != "intrcpt") {
      # Keyed on the coefficient name rather than the position, so that
      # cell-means models (`~ 0 + moderator`), which have no intercept at
      # all, name their first coefficient too.

      variable_name <- rownames(m$b)[ii]

      # If it's an interaction
      if(str_detect(variable_name, ":")) {
        str_replace(variable_name, ":", " x ")
      } else {

        # Match the variable name and the moderator name
        temp_mods <- extract_vars_flex(m$formula.mods)
        temp_mod  <- temp_mods[which(str_detect(variable_name, temp_mods))]

        # If it detects more than one
        if(length(temp_mod) > 1) {
          if(length(temp_mods[which(variable_name == temp_mods)]) > 0) {
            temp_mod <- temp_mods[which(variable_name == temp_mods)]
          } else {
            # Take the longest one
            temp_mod <- temp_mods[which(str_count(temp_mods) == max(str_count(temp_mods)))]
          }
        }

        # If the mod name isn't the whole thing, add spaces
        if(variable_name != temp_mod) {
          variable_name <- str_replace(variable_name, temp_mod,
                                       paste0(temp_mod, separator))
        }
        variable_name
      }

    } else {
      "Intercept"
    }
  }, character(1))
}
