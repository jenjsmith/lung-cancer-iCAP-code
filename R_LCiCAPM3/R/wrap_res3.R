# Wrap a legacy res_3 (list of runs) into a structured object with meta

wrap_res3 <- function(fp_in = "results/res_3.rds",
                      fp_out = "results/res_3_wrapped.rds",
                      seed_list = NULL,
                      overwrite_out = FALSE) {
  stopifnot(file.exists(fp_in))
  if (!overwrite_out && file.exists(fp_out)) {
    stop("Output file already exists: ", fp_out,
         "\nSet overwrite_out=TRUE to overwrite.")
  }
  
  # 1) Load legacy object
  res_legacy <- readRDS(fp_in)
  
  # If it's already wrapped, just copy or return it
  if (is.list(res_legacy) && !is.null(res_legacy$results) && !is.null(res_legacy$meta)) {
    message("Input appears to already be wrapped; writing a copy to: ", fp_out)
    saveRDS(res_legacy, fp_out)
    return(invisible(res_legacy))
  }
  
  # Basic sanity: expect a list of runs, each with $auc and $final_model
  if (!is.list(res_legacy) || length(res_legacy) == 0)
    stop("Unexpected structure: expected a non-empty list of runs.")
  
  # 2) Extract AUCs safely
  auc_safe <- function(run) {
    a <- try(run$auc, silent = TRUE)
    if (inherits(a, "try-error") || is.null(a)) NA_real_ else as.numeric(a)
  }
  aucs <- vapply(res_legacy, auc_safe, numeric(1))
  if (all(is.na(aucs))) stop("All AUCs are NA; cannot select a model.")
  
  # 3) 0.75 quantile and chosen index (deterministic tie-break = lowest index)
  q75   <- as.numeric(stats::quantile(aucs, 0.75, type = 7, na.rm = TRUE))
  diffs <- abs(aucs - q75)
  idx   <- min(which(diffs == min(diffs, na.rm = TRUE)))
  
  # 4) Seed for the chosen run
  # prefer run$seed if present; else take from provided seed_list; else NA
  run_seed <- try(res_legacy[[idx]]$seed, silent = TRUE)
  if (inherits(run_seed, "try-error") || is.null(run_seed)) {
    if (!is.null(seed_list) && length(seed_list) >= idx) {
      chosen_seed <- seed_list[idx]
    } else {
      chosen_seed <- NA_integer_
    }
  } else {
    chosen_seed <- as.integer(run_seed)
  }
  
  # 5) Build meta & wrapped object
  meta <- list(
    chosen_index = idx,
    chosen_seed  = chosen_seed,
    chosen_auc   = aucs[idx],
    q75_auc      = q75,
    aucs         = aucs,
    method_note  = "Chosen run = closest to 0.75 quantile (type=7); ties -> lowest index.",
    timestamp    = Sys.time(),
    session_info = utils::sessionInfo()
  )
  
  # Special note for res_3 (legacy runs): seed was not stored separately,
  # but run number == seed number, so chosen_seed = chosen_index
  if (all(is.na(meta$chosen_seed))) {
    meta$note <- "Seeds not stored in this run; run number = seed number (so chosen_seed = chosen_index)."
  }
  
  res_wrapped <- list(results = res_legacy, meta = meta)
  
  # 6) Write out (make a backup next to the original for safety)
  dir.create(dirname(fp_out), showWarnings = FALSE, recursive = TRUE)
  saveRDS(res_wrapped, fp_out)
  message("Wrapped object written: ", normalizePath(fp_out, mustWork = TRUE))
  
  # small, human-readable meta file
  meta_txt <- file.path(dirname(fp_out), sub("[.]rds$", "_meta.txt", basename(fp_out)))
  writeLines(c(
    paste("chosen_index:", meta$chosen_index),
    paste("chosen_seed:",  meta$chosen_seed), #not saved for M3 but chosen seed = chosen index
    paste("chosen_auc:",   sprintf("%.6f", meta$chosen_auc)),
    paste("q75_auc:",      sprintf("%.6f", meta$q75_auc)),
    paste("timestamp:",    format(meta$timestamp, "%Y-%m-%d %H:%M:%S %Z"))
  ), meta_txt)
  message("Meta summary written:  ", normalizePath(meta_txt, mustWork = TRUE))
  
  invisible(res_wrapped)
}

# --- Usage ---
# Minimal (writes results/res_3_wrapped.rds):
#wrapped <- wrap_res3(fp_in = "results/res_3.rds")

# Load later, consistently:
# res_3 <- readRDS("results/res_3_wrapped.rds")
# model_3 <- res_3$results[[ res_3$meta$chosen_index ]]$final_model
