# 0) Make sure we're in the project and activate renv
renv::activate(".")

# 1) Decide if a restore is actually needed
needs_restore <- tryCatch({
  out <- capture.output(renv::status())
  !any(grepl("No issues found|The project is synchronized", out))
}, error = function(e) TRUE)

# 2) If needed, restore WITHOUT restarting the session
if (needs_restore) {
  renv::settings$bioconductor.version("3.20")
  Sys.setenv(RENV_CONFIG_RESTART = "FALSE")  # belt & suspenders
  renv::restore(prompt = FALSE, restart = FALSE)
  renv::snapshot(prompt = FALSE)
}

# 3) Knit in a **fresh child R process** so nothing can restart mid-run
if (!requireNamespace("callr", quietly = TRUE)) install.packages("callr")
if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown")

rmd <- normalizePath("LCiCAP_RFM4generation.Rmd", mustWork = TRUE)
out <- callr::r(
  function(file) {
    renv::activate(".")  # activate project in the child process
    wd <- dirname(file)
    knitr::opts_knit$set(root.dir = wd)
    rmarkdown::render(
      input         = file,
      output_format = "html_document",
      output_file   = "LCiCAP_RFM4generation.html",
      output_dir    = wd,
      clean         = TRUE,
      quiet         = FALSE
      # NOTE: don't set a super-isolated envir; let base::format() etc. be found
    )
  },
  args = list(rmd),
  env = c(Sys.getenv(), RENV_CONFIG_RESTART = "FALSE")  # also pass to child
)

cat("\nRendered HTML:\n", normalizePath(out, mustWork = TRUE), "\n")
