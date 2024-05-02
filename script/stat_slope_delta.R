stat_slope_delta <- function(dataset) {
  # Define the model
  code <- nimbleCode({
    # Prior distributions
    beta0_1 ~ dnorm(0, 0.001) # Prior for intercept of first dataset
    beta0_2 ~ dnorm(0, 0.001) # Prior for intercept of second dataset
    beta1 ~ dnorm(0, 0.001) # Prior for slope of first dataset
    delta ~ dnorm(0, 0.001) # Prior for difference between slopes
    sigma_1_sq ~ dunif(0, 100) # Prior for variance of first dataset
    sigma_2_sq ~ dunif(0, 100) # Prior for variance of second dataset
    
    # Likelihood for first dataset
    for(i in 1:n1) {
      y1[i] ~ dnorm(beta0_1 + beta1 * x1[i], sqrt(sigma_1_sq))
    }
    
    # Likelihood for second dataset
    for(i in 1:n2) {
      y2[i] ~ dnorm(beta0_2 + (beta1 + delta) * x2[i], sqrt(sigma_2_sq))
    }
  })

  # Create NIMBLE model
  # Define data
  rural <- subset(dataset, group == "No ALAN")
  urban <- subset(dataset, group == "ALAN")
  
  data <- list(x1 = rural$Tmax_Summer, y1 = rural$First_Yes_DOY, n1 = nrow(rural),
               x2 = urban$Tmax_Summer, y2 = urban$First_Yes_DOY, n2 = nrow(urban))
  
  model <- nimbleModel(code, data)
  
  # Configure and build the MCMC
  mcmc_config <- configureMCMC(model)
  mcmc <- buildMCMC(mcmc_config)
    # Compile the MCMC
  compiled_model <- compileNimble(model)
  compiled_mcmc <- compileNimble(mcmc, project = compiled_model, showCompilerOutput = TRUE) # project should be model/complied model?

  # Run the MCMC
  samples <- runMCMC(compiled_mcmc, niter = 30000, nburnin = 5000, nchains = 5)
  # mcmc_diagnose(samples)
  # nimble:::mcmc_diagnose(samples)
  
  # plot delta
  delta_samples <- data.frame(delta = numeric(0))
  for (i in 1:length(samples)) {
    chain_name <- paste("chain", i, sep = "")
    delta_samples <- rbind(delta_samples, data.frame(delta = samples[i][[chain_name]][,4]))
  }
  
  delta_plot <- ggplot(delta_samples, aes(x = delta)) +
    geom_density(fill = "lightblue", color = "black") +  # Use geom_density for smoothed curve
    labs(title = dataset$Genus[1], x = "Delta", y = "Density") +
    theme_classic()+
    geom_vline(xintercept = 0, linetype = "solid", color = "red") 
  delta_plot
  return(delta_plot)
  
}