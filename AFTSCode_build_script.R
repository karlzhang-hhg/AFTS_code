# Rscript ../AFTSCode_build_script.R
library(devtools, roxygen2)

roxygen2::roxygenize()
devtools::build()
devtools::install()
