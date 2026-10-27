packages <- c("forecast", "ggplot2", "lmtest", "tseries", "dplyr", "tibble", "purrr", "stringr", "fpp2")
missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
