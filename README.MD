These are the general functions for HLI projects. This allows us to keep the functions we use over and over again in the same place, so we don't have to copy them every time. Thereby, only one script needs to be edited to update these functions.

Some of these involve functions to simplify the making of graphs.

**Note that these will often require tidyverse**

# Functions
To load functions in your R script use:
```{r}
# For basic functions
devtools::source_url("https://raw.githubusercontent.com/Happier-Lives-Institute/general-functions/main/functions_basic.R")
# For reporting functions
devtools::source_url("https://raw.githubusercontent.com/Happier-Lives-Institute/general-functions/main/functions_reporting.R")
# For Cohen's d functions
devtools::source_url("https://raw.githubusercontent.com/Happier-Lives-Institute/general-functions/main/functions_cohen_d.R")
# For meta-analysis functions
devtools::source_url("https://raw.githubusercontent.com/Happier-Lives-Institute/general-functions/main/functions_meta_analysis.R")
```

To run it without the SHA message use:
```{r}
# For basic functions
suppressMessages(devtools::source_url("https://raw.githubusercontent.com/Happier-Lives-Institute/general-functions/main/functions_basic.R"))
# For reporting functions
suppressMessages(devtools::source_url("https://raw.githubusercontent.com/Happier-Lives-Institute/general-functions/main/functions_reporting.R"))
# For Cohen's d functions
suppressMessages(devtools::source_url("https://raw.githubusercontent.com/Happier-Lives-Institute/general-functions/main/functions_cohen_d.R"))
# For meta-analysis functions
suppressMessages(devtools::source_url("https://raw.githubusercontent.com/Happier-Lives-Institute/general-functions/main/functions_meta_analysis.R"))
```

# Setting up a project
For `install_if_missing`, `library_quietly`, `show_conflicts`, `safe_model_list`,
and `flextable_apa_theme`.
```{r}
suppressMessages(devtools::source_url("https://raw.githubusercontent.com/Happier-Lives-Institute/general-functions/main/functions_setup.R"))
```
Note that `source_url` needs devtools already installed, so `install_if_missing`
still has to be copied inline if you want it before anything else.

# Meta-analysis extras
These need `functions_meta_analysis.R`, so source that first.
```{r}
# Naming a model's coefficients (needed by functions_meta_flextable.R)
suppressMessages(devtools::source_url("https://raw.githubusercontent.com/Happier-Lives-Institute/general-functions/main/functions_coefficients.R"))
# Presenting several models side by side as a flextable
suppressMessages(devtools::source_url("https://raw.githubusercontent.com/Happier-Lives-Institute/general-functions/main/functions_meta_flextable.R"))
# Fitting and predicting from restricted cubic splines
suppressMessages(devtools::source_url("https://raw.githubusercontent.com/Happier-Lives-Institute/general-functions/main/functions_splines.R"))
# A limitmeta that does not crash
suppressMessages(devtools::source_url("https://raw.githubusercontent.com/Happier-Lives-Institute/general-functions/main/functions_publication_bias.R"))
# z-curve densities for ggplot, and the selection-bias tests
suppressMessages(devtools::source_url("https://raw.githubusercontent.com/Happier-Lives-Institute/general-functions/main/functions_zcurve.R"))
```

`meta_flextable` labels its rows with the function passed as `label_fun`. It shows
the raw coefficient names by default, so pass the project's own labelling function.

`get_spline_ma` fits with the function passed as `fit_rmv`, called as
`fit_rmv(dat = dat, mods = <formula>, method = <"ML"|"REML">)`, so the random
structure stays in the project.

# Graphs
```{r}
# For themes and saving
suppressMessages(devtools::source_url("https://raw.githubusercontent.com/Happier-Lives-Institute/general-functions/main/functions_styles.R"))
# For bar description
suppressMessages(devtools::source_url("https://raw.githubusercontent.com/Happier-Lives-Institute/general-functions/main/graph_bar_description.R"))
```

# Package requirements
Beyond tidyverse:

| Package | Needed by |
| --- | --- |
| `digest` | `functions_meta_analysis.R` (the null-model cache) and `save_docx_if_changed()` |
| `metafor` | `functions_meta_analysis.R`, `functions_splines.R` |
| `cowplot` | `functions_styles.R`, and the funnel plot in `functions_meta_analysis.R` |
| `rms` (and `Hmisc`, which comes with it) | `functions_splines.R`, and `make_mean_newmods()` when passed `knots` |
| `flextable`, `officer` | `functions_meta_flextable.R`, `flextable_apa_theme()`, `save_docx_if_changed()` |
| `metasens` | `functions_publication_bias.R` |
| `zcurve` | `functions_zcurve.R` |
| `msm`, `pracma` | `get_total_effect()` in `functions_reporting.R` |
| `rstan`, `HDInterval` | `stan_rows()` only — the rest of `functions_meta_flextable.R` works without them |

`functions_publication_bias.R` and `functions_zcurve.R` use `metasens:::` and
`zcurve:::` internals, so they may need updating when those packages change.
Written against metasens 1.5.3 and zcurve 2.4.6.
