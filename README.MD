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

# Graphs
```{r}
# For bar description
suppressMessages(devtools::source_url("https://raw.githubusercontent.com/Happier-Lives-Institute/general-functions/main/graph_bar_description.R"))
```