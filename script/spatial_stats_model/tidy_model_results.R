####################### INLA model results tidy #######################

tree_dir <- "~/lab-data/datasets/vegetation/StreetTrees/cleaned_ps_input/"
rds_files <- list.files(tree_dir, pattern = "\\.rds$", full.names = TRUE)
cities <- gsub("\\.rds$", "", basename(rds_files))
phe_types <- c("SOS", "EOS")

prepare_coef_summary <- function(r_output_dir, phe_type){
  pattern = paste0(phe_type,"\\.rds$")
  files <- list.files(r_output_dir, pattern = pattern, full.names = TRUE)
  coef_summary <- data.frame()
  # for one genus
  for (file in files) {
    model_data <- readRDS(file)
    species <- model_data$species
    # beta_names <- c("beta[0]", "beta[AvgTemp_mean]", "beta[Sum_mm_sum]", "beta[srad_mean]")
    # df_list <- lapply(seq_along(model_data$res$marginals.fixed), function(i) {
    #   data.frame(x = model_data$res$marginals.fixed[[i]][,1], 
    #              density = model_data$res$marginals.fixed[[i]][,2], 
    #              beta = beta_names[i])
    # })
    # df <- bind_rows(df_list)
    beta_avg_temp <- model_data$res$marginals.fixed[[2]]
    
    coef_summary <- coef_summary %>%
      rbind(model_data$fix_effect)
  }
  
  coef_summary <- coef_summary %>%
    mutate(SigAlpha = ifelse(Sig == "*", 1, 0.2)) %>%
    mutate(Color = ifelse(Estimate >0, "#0D0887FF" , "#FCA636FF"))
  write_rds(coef_summary, paste0(r_output_dir, phe_type, "_summary.rds"))
}

for (city in cities[2:10]) {
  message(paste0("Processing city: ", city))
  r_output_dir <- paste0("~/phenology-urban/data/proc/model_output/final_models_Grid/inla/spatial_spde/", city, "/canonicalname/")
  for (phe_type in phe_types){
    prepare_coef_summary(r_output_dir, phe_type)
  }
}




plot_coef_sp_city <- function(r_output_dir, phe_type){
  coef_summary <- read_rds(paste0(r_output_dir, phe_type, "_summary.rds"))
  p_coef <-  ggplot(coef_summary, aes(x = factor(Factor, levels = rev(Factor)), 
                                      y = Estimate,
                                      colour = Color)) +
    geom_point(aes(alpha = SigAlpha), position = position_dodge(w = 0.5), size = 3) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper, alpha = SigAlpha), 
                  position = position_dodge(w = 0.5), 
                  width = 0.2) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(x = NULL, y = "Estimate") +
    coord_flip() +
    scale_alpha_identity() +
    scale_color_identity() +
    theme(legend.position = "right") +
    geom_text(aes(label = Sig, y = Upper + 0.05  ), hjust = -0.5, show.legend = FALSE) +
    geom_text(aes(label = Estimate, y = Estimate  ), vjust = -0.7, show.legend = FALSE, size = 3) +
    guides(alpha = "none") +
    # add title
    ggtitle(bquote(beta[Temp]~"for"~.(phe_type))) +
    theme_minimal()
  
  return(list(
    coef_summary_plot = p_coef
  ))
}

city = "Berkeley"
r_output_dir <- paste0("~/phenology-urban/data/proc/model_output/final_models_Grid/inla/spatial_spde/", city, "/canonicalname/")
plots_coef <- list()
assess_species <- data.frame()
for (phe_type in phe_types) {
  p1 <- plot_coef_sp_city(r_output_dir, phe_type)
  plots_coef[[phe_type]] <- p1$coef_summary_plot
  assess_species <- rbind(assess_species, p1$assess)
}
wrap_plots(plots_coef, ncol = 2)

