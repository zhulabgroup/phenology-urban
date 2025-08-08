run_nonsp_inla = function(model_input_all, r_output_dir, phe_type, city, 
                 taxon_col = "canonicalname", model_type = "no spatial") {
  # options(repos = c(
  #   INLA = "https://inla.r-inla-download.org/R/testing",
  #   getOption("repos")
  # ))
  # install.packages("inlabru")
  # install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
  
  taxa_list <- unique(model_input_all[[taxon_col]])
  boundary <- switch(
    city,
    "NY" = st_read("~/urban-cooling/data/raw/WU/NY/Untitled/nybb_dissolved.shp"),
    "DT" = st_read("~/urban-cooling/data/raw/WU/DT/City_of_Detroit_Boundary/City_of_Detroit_Boundary.shp"),
    "ST" = st_read("~/urban-cooling/data/raw/WU/ST/seattle_boundary/CityUGA.shp") %>%
      filter(CITY_NM == "Seattle"),
    "HT" = st_read("~/urban-cooling/data/raw/WU/HT/HGAC_City_Boundaries/HGAC_City_Boundaries.shp") %>%
      filter(NAME == "Houston"),  # change to some distance
    "TP" = st_read("~/urban-cooling/data/raw/WU/TP/Municipal_Boundary/Municipal_Boundary.shp") %>%
      filter(NAME == "Tampa"),
    "DV" = st_read("~/urban-cooling/data/raw/WU/DV/DRCOG_Municipalities/DRCOG_Municipalities.shp") %>%
      filter(city == "Denver"),
    "AT" = st_read("~/urban-cooling/data/raw/WU/AT/TxDOT_City_Boundaries/Cities.shp") %>%
      filter(CITY_NM == "Austin"),
    stop("Invalid city abbreviation.")
  )
  
  for (taxa in taxa_list) {
    if (model_type == "no spatial") {
      folder_path <- paste0(r_output_dir, "no_sp/",city,"/", taxon_col)
      }else if (model_type == "random spatial") {
      folder_path <- paste0(r_output_dir, "random_group_sp/",city,"/", taxon_col)
      }else if (model_type == "fix spatial") {
        folder_path <- paste0(r_output_dir, "fix_group_sp/",city,"/", taxon_col)
      }
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path, recursive = TRUE)
    }
    model_path<- paste0(folder_path, "/", gsub(" ", "_", taxa), "_", phe_type, ".rds")
    if (file.exists(model_path)) {
      cat(sprintf("\n File already exists for taxa: %s, phe_type: %s, city: %s. Skipping INLA run.\n", 
                  taxa, phe_type, city))
      next  # Skip the current iteration of the loop
    } # if the file exists, skip to the next iteration  
    
    cat(sprintf("\n Running INLA for taxa: %s, phe_type: %s, city: %s\n", 
                taxa, phe_type, city))
    if (phe_type %in% c("SOS","EOS")) {
      phe_metric <- "doy"
    } else {
      phe_metric <- "doy_diff"
    }
    
    model_input <- model_input_all %>%
      filter(!!sym(taxon_col) == !!taxa)
    data = data.frame(y = model_input[[phe_metric]],
                      Intercept = rep(1, nrow(model_input)),
                      temp = (model_input %>% select(starts_with("temp")) %>% pull()),
                      prcp = (model_input %>% select(starts_with("prcp")) %>% pull()),
                      srad = (model_input %>% select(starts_with("srad")) %>% pull()),
                      year = model_input$year,
                      hex_id = model_input$hex_id)
    
    ######################################## Run/Save the model ########################################
    fixed.priors <- list(
      mean = list(Intercept = 180, temp = 0, prcp = 0, srad = 0),
      prec = list(Intercept = 1 / (40^2), temp = 1 / (5^2), prcp = 1 / (1^2), srad = 1 / (1^2))
    )
    
    ### Run the model
    if (model_type == "no spatial") {
      formula <- y ~ 0 + Intercept + temp + prcp + srad + 
        f(year, model = "iid",
          hyper = list(prec = list(
            prior = "loggamma",
            param = c(0.01, 0.01),
            initial = 0, 
            fixed = FALSE
            )))
    }else if (model_type == "random spatial") {
      formula <- y ~ 0 + Intercept + temp + prcp + srad + 
        f(year, model = "iid",
          hyper = list(prec = list(
            prior = "pc.prec",
            param = c(0.5, 0.1)
          ))) +
        f(hex_id, model = "iid",
          hyper = list(prec = list(
            prior = "pc.prec",
            param = c(0.5, 0.05)
          )))
    }else if (model_type == "fix spatial") {
      formula <- y ~ 0 + Intercept + temp + prcp + srad + hex_id +
        f(year, model = "iid",
          hyper = list(prec = list(
            prior = "loggamma",
            param = c(0.01, 0.01),
            initial = 0, 
            fixed = FALSE
          ))) 
    }
    
    res <- inla(formula,
                family = "gaussian",
                data = data,
                control.inla = list(
                  strategy = "laplace",interpolator = "gaussian",fast = FALSE,restart = 5
                ),
                control.fixed = fixed.priors,
                control.compute = list(config = TRUE, cpo = TRUE, dic = TRUE, waic = TRUE)
                )
    # res$waic$waic
    # summary(res)
    
    fix_effect <- as.data.frame(summary(res)$fixed)%>%
      rename(Estimate = mean,Lower = `0.025quant`,Upper = `0.975quant`) %>%
      mutate(Sig = ifelse(Lower*Upper>0, "*", "")) 
    fix_effect <- data.frame(
      Factor = rownames(fix_effect),  # Add Factor column
      fix_effect, 
      row.names = NULL  # Remove row names to avoid duplication
    ) %>%
      filter(Factor == "temp") %>%
      mutate(Factor = ifelse(Factor == "temp", taxa, Factor))
    
    ##### get mlik? #####
    
    
    ######################################## Save the model results ########################################
    list(res = res,
         fix_effect = fix_effect,
         taxa = taxa) %>%
      write_rds(model_path, compress = "gz")
  }
}


plot_coef_taxa_city <- function(r_output_dir, phe_type){
  pattern = paste0(phe_type,"\\.rds$")
  files <- list.files(r_output_dir, pattern = pattern, full.names = TRUE)
  plots_tavg <- list()
  
  coef_summary <- data.frame()
  assess <- data.frame()
  random_sp_summary <- data.frame()
  random_year_summary <- data.frame()
  fixed_sp_summary <- data.frame()
  
  # for one taxa
  for (file in files) {
    model_data <- readRDS(file)
    taxa <- model_data$taxa
    # beta_names <- c("beta[0]", "beta[AvgTemp_mean]", "beta[Sum_mm_sum]", "beta[srad_mean]")
    # df_list <- lapply(seq_along(model_data$res$marginals.fixed), function(i) {
    #   data.frame(x = model_data$res$marginals.fixed[[i]][,1], 
    #              density = model_data$res$marginals.fixed[[i]][,2], 
    #              beta = beta_names[i])
    # })
    # df <- bind_rows(df_list)
    #### Temperature
    beta_avg_temp <- model_data$res$marginals.fixed$temp
    df <- data.frame(
      x = beta_avg_temp[, 1],
      density = beta_avg_temp[, 2],
      beta = "beta[temp]"
    )
    plot_tavg <- ggplot(df, aes(x = x, y = density)) +
      geom_line(linewidth = 1, color = "blue") +
      # facet_wrap(~beta, scales = "free", labeller = label_parsed) + 
      theme_minimal() +
      labs(x = "Coefficient Value", y = "Density") +
      ggtitle(paste(taxa)) +
      theme(legend.position = "none",
            plot.title = element_text(size = 8))
    plots_tavg[[taxa]] <- plot_tavg
    
    ### random hex_id
    random_sp <- model_data$res$marginals.random$hex_id
    if(!is.null(random_sp)) {
      random_sp_taxa <- do.call(rbind, lapply(names(random_sp), function(nm) {
        df <- as.data.frame(random_sp[[nm]])
        df$index <- nm
        return(df)
      }))
      random_sp_summary <- random_sp_summary %>%
        rbind(random_sp_taxa |> mutate(taxa = taxa))
    }
    ### fixed hex_id
    fixed_sp <- model_data$res$marginals.fixed[grep("^hex_id", names(model_data$res$marginals.fixed))]
    fixed_sp_taxa <- do.call(rbind, lapply(names(fixed_sp), function(nm) {
      df <- as.data.frame(fixed_sp[[nm]])
      df$index <- nm
      return(df)
    }))
    if(!is.null(fixed_sp_taxa)) {
      fixed_sp_summary <- fixed_sp_summary %>%
        rbind(fixed_sp_taxa |> mutate(taxa = taxa))
    }
    
    ### random year
    random_year <- model_data$res$marginals.random$year
    if(!is.null(random_year)) {
      random_year_taxa <- do.call(rbind, lapply(names(random_year), function(nm) {
        df <- as.data.frame(random_year[[nm]])
        df$index <- nm
        return(df)
      }))
      random_year_summary <- random_year_summary %>%
        rbind(random_year_taxa |> mutate(taxa = taxa))
    }
    
    coef_summary <- coef_summary %>%
      rbind(model_data$fix_effect)
    
    summary_res <- summary(model_data$res)

    assess <- assess |>
      rbind(data.frame(
        milk_type = rownames(summary_res$mlik),
        milk = as.numeric(summary_res$mlik[, 1]),
        Parameter = taxa,
        phe_type = phe_type,
        waic = summary_res$waic$waic,
        row.names = NULL
      ))
  }
  
  combined_plot_tavg <- wrap_plots(plots_tavg, ncol = 5)
  
  coef_summary <- coef_summary %>%
    mutate(SigAlpha = ifelse(Sig == "*", 1, 0.2))%>%
    mutate(Color = ifelse(Estimate >0, "#0D0887FF" , "#FCA636FF")) %>%
    filter(Factor != "Alnus")
  
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
    scale_colour_identity() +
    theme(legend.position = "right") +
    geom_text(aes(label = Sig, y = Upper + 0.05  ), hjust = -0.5, show.legend = FALSE) +
    geom_text(aes(label = Estimate, y = Estimate  ), vjust = -0.7, show.legend = FALSE, size = 3) +
    guides(alpha = "none") +
    ggtitle(bquote(beta[Temp]~"for"~.(phe_type))) +
    theme_minimal()
  
  if (nrow(random_sp_summary) == 0) {
    p_random_sp <- NULL
  } else {
    p_random_sp <-  ggplot(random_sp_summary, aes(x = x, y = y, color = index)) +
      geom_line(linewidth = 0.2)+
      facet_wrap(~taxa, ncol = 1, strip.position = "left") +
      guides(color = "none") +
      ggtitle(bquote("Spatial effect for"~.(phe_type))) +
      theme_classic() +
      scale_x_continuous(limits = c(-10, 10)) +
      theme(
        axis.line.y  = element_blank(),      
        axis.text.y  = element_blank(), 
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),  
        strip.placement = "outside",  
        strip.text.y.left = element_text(angle = 0) 
      )
  }
  
  if (nrow(fixed_sp_summary) == 0) {
    p_fixed_sp <- NULL
  } else {
    p_fixed_sp <-  ggplot(fixed_sp_summary, aes(x = x, y = y, color = index)) +
      geom_line(linewidth = 0.2)+
      facet_wrap(~taxa, ncol = 1, strip.position = "left") +
      guides(color = "none") +
      ggtitle(bquote("Spatial effect for"~.(phe_type))) +
      theme_classic() +
      theme(
        axis.line.y  = element_blank(),      
        axis.text.y  = element_blank(), 
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),  
        strip.placement = "outside",  
        strip.text.y.left = element_text(angle = 0) 
      )
  }
  
  p_random_year <-  ggplot(random_year_summary, aes(x = x, y = y, color = index)) +
    geom_line(linewidth = 0.2)+
    facet_wrap(~taxa, ncol = 1, strip.position = "left") +
    guides(color = "none") +
    ggtitle(bquote("Year effect for"~.(phe_type))) +
    theme_classic() +
    theme(
      axis.line.y  = element_blank(),      
      axis.text.y  = element_blank(), 
      axis.ticks.y = element_blank(),
      strip.background = element_blank(),  
      strip.placement = "outside",  
      strip.text.y.left = element_text(angle = 0) 
    )

  
  return(list(
    combined_plot_tavg = combined_plot_tavg,
    coef_summary_plot = p_coef,
    random_sp_plot = p_random_sp,
    fixed_sp_plot = p_fixed_sp,
    random_year_plot = p_random_year,
    assess = assess
  ))
}

