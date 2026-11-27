####################### INLA model results tidy #######################
prepare_coef_summary <- function(r_output_dir, phe_type, on_error = NULL){
  pattern = paste0(phe_type,"\\.rds$")
  files <- list.files(r_output_dir, pattern = pattern, full.names = TRUE)
  coef_summary <- data.frame()
  failed <- FALSE
  # for one genus
  for (file in files) {
    message(paste0("Processing file: ", file))
    result <- tryCatch({
      model_data <- readRDS(file)
      fix_effect <- model_data$fix_effect
      if (!is.data.frame(fix_effect) ||
          !all(c("Factor", "Sig", "Estimate") %in% names(fix_effect))) {
        stop("Missing or invalid fix_effect")
      }
      fix_effect
    }, error = function(e) {
      failed <<- TRUE
      if (is.null(on_error)) warning(file, ": ", conditionMessage(e)) else on_error(file, e)
      NULL
    })
    coef_summary <- dplyr::bind_rows(coef_summary, result)
  }
  if (!nrow(coef_summary)) return(coef_summary)
  coef_summary <- coef_summary %>%
    mutate(SigAlpha = ifelse(Sig == "*", 1, 0.2)) %>%
    mutate(Color = ifelse(Estimate >0, "#0D0887FF" , "#FCA636FF"))
  summary_path <- file.path(r_output_dir, paste0(phe_type, "_summary.rds"))
  if (!failed && !file.exists(summary_path)) write_rds(coef_summary, summary_path)
  return(coef_summary)
}



############################### MCMC model result #########################

prepare_coef_summary_mcmc <- function(r_output_dir, phe_type, on_error = NULL){
  files <- list.files(r_output_dir, pattern = paste0(phe_type, "\\.rds$"), full.names = TRUE)
  failed <- FALSE
  combined_param_data <- purrr::map(files, function(file) {
    tryCatch({
      model_data <- readRDS(file)
      genus <- sub(paste0("(.*?)_", phe_type, "\\.rds"), "\\1", basename(file))
      ggmcmc::ggs(model_data$samples) |>
        dplyr::filter(stringr::str_starts(Parameter, stringr::fixed("temp"))) |>
        dplyr::rename(ParameterOriginal = Parameter) |>
        dplyr::mutate(Parameter = genus, converge = model_data$converge)
    }, error = function(e) {
      failed <<- TRUE
      if (is.null(on_error)) warning(file, ": ", conditionMessage(e)) else on_error(file, e)
      NULL
    })
  }) |>
    bind_rows()
  if (!nrow(combined_param_data)) return(data.frame())
  
  # Summarize
  summary_data <- combined_param_data |>
    group_by(Parameter) |>
    summarise(
      Estimate = median(value, na.rm = TRUE),
      Lower = quantile(value, 0.025, na.rm = TRUE),
      Upper = quantile(value, 0.975, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      # waic = first(waic),
      converge = first(converge),
      # milk = first(milk),
      .groups = "drop"
    ) |>
    mutate(Factor= Parameter,
      Sig = ifelse(Upper * Lower > 0, "*", ""),
           SigAlpha = ifelse(Sig == "*", 1, 0.2),
           Color = ifelse(Estimate >0, "#0D0887FF" , "#FCA636FF")) 
  summary_path <- file.path(r_output_dir, paste0(phe_type, "_summary.rds"))
  if (!failed && !file.exists(summary_path)) write_rds(summary_data, summary_path)
  return(summary_data)
}
