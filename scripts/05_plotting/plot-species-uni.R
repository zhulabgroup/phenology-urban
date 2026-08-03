source(here("scripts/04_stats_model/fit_pglmm.R"))


plot_species_uni <- function(clim, phe_type, species_example, df_coef) {
  
  x_lab <- if (rlang::as_name(clim) %in% c("bio1")) {
    expression("Mean annual temperature ("*degree*C*")")
  } else if (rlang::as_name(clim) %in% c("bio12")) {
    "Mean annual precipitation (mm)"
  } else {
    rlang::as_name(clim)
  }
  
  pred_df <- lapply(
    species_example,
    function(sp) {
      fit_uni(
        phe_type = phe_type, input_data = df_coef, species_name = sp, clim = clim, bootstrap = T)$pred_table
    }
  ) |>
    bind_rows()
  
  fit_input_df <- df_coef |>
    filter(Factor %in% species_example, phe_type == !!phe_type) |>
    select(Estimate, Lower, Upper, !!clim, Factor)
  x_vals <- fit_input_df[[clim]]
  x_min  <- min(x_vals, na.rm = TRUE)
  x_max  <- max(x_vals, na.rm = TRUE)
  x_rng  <- x_max - x_min
  
  p_phe <- ggplot(fit_input_df, aes(x = .data[[clim]], y = Estimate, color = Factor)) +
    geom_pointrange(aes(ymin = Lower, ymax = Upper), alpha = 0.1, size = 0.2) +
    geom_ribbon(data = pred_df,
                aes(x = .data[[clim]], ymin = y_lwr, ymax = y_upr, fill = species),
                inherit.aes = FALSE,
                alpha = 0.2) +
    geom_line(data = pred_df,
              aes(x = .data[[clim]], y = y_pred, color = species),
              inherit.aes = FALSE,
              linewidth = 0.5) +
    guides(alpha = "none") +
    scale_x_continuous(
      limits = c(x_min - 0.15 * x_rng, x_max)
    ) +
    scale_y_continuous(
      limits = c(-30, 30)
    ) +
    theme_classic() +
    labs(
      x = x_lab,
      y = bquote(beta[Temp]~" estimate for "*.(phe_type)),   
      color = "Species",
      fill  = "Species"
    ) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(face = "italic")
      # axis.title = element_blank()
    )
  return(p_phe)
}
