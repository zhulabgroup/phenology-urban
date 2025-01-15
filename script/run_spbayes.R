library(spBayes)
library(sp)
library(sf)
library(dplyr)
library(stringr)

# Read file paths from environment variables
# r_output_file <- Sys.getenv("R_OUTPUT_FILE")
# genus <- Sys.getenv("GENUS")
# direction <- Sys.getenv("DIRECTION")
# thres <- as.numeric(Sys.getenv("THRES"))
# phe_metric <- Sys.getenv("PHE_METRIC")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("5 arguments needed: r_output_file, genus, direction, thres, phe_metric")
}

r_output_file <- args[1]
genus <- args[2]
direction <- args[3]
thres <- as.numeric(args[4])
phe_metric <- args[5]

# r_output_file <- "~/phenology-urban/data/proc/urban/sp_model/"
# genus <- "Carya"
# direction <- "up"
# thres <- 0.5
# phe_metric <- "doy"


tavg_allyear_var <- readRDS("~/phenology-urban/data/proc/urban/NY/tavg_allyear_var.rds")
tavg_phe <- tavg_allyear_var %>%
  filter(direction == !!direction & thres == !!thres & genus == !!genus) %>%
  na.omit()
coordinates(tavg_phe) <- ~ lon + lat
print(length(tavg_phe))

set.seed(2025)
coords <- coordinates(tavg_phe)

duplicated_flags <- duplicated(coords) | duplicated(coords, fromLast = TRUE)
jitter_matrix <- matrix(0, nrow = nrow(coords), ncol = 2)
jitter_matrix[duplicated_flags, ] <- matrix(rnorm(sum(duplicated_flags) * 2, mean = 0, sd = 1e-6),  ncol = 2)
coords_jittered <- coords + jitter_matrix

tavg_phe_df <- as.data.frame(tavg_phe)
tavg_phe_df$lon <- coords_jittered[, 1]
tavg_phe_df$lat <- coords_jittered[, 2]


formula <- as.formula(paste0(phe_metric, " ~ AvgTemp_mean + Sum_mm_sum + srad_mean"))

num_predictors <- 4
D <- as.matrix(dist(coords_jittered))
n.samples <- 100
starting <- list("phi" = 25, "sigma.sq" = 40, "tau.sq" = 100)
tuning <- list("phi" = 0.1, "sigma.sq" = 0.1, "tau.sq" = 0.1)
priors <- list(
  "beta.Norm" = list(mean = rep(0, num_predictors), var = diag(1000, num_predictors)),
  "phi.Unif" = c(5 / 1, 5 / 0.1),
  "sigma.sq.IG" = c(3, 80),
  "tau.sq.IG" = c(3, 200)
)

print(formula)
model <- spLM(formula, data=tavg_phe_df, coords=coords_jittered,
              starting=starting, tuning=tuning, priors=priors,knots = c(20, 20, 0),
              cov.model="exponential", verbose=FALSE, n.samples=n.samples)


burn.in <- 0.5* n.samples
#recover beta and spatial random effects
m_re <- spRecover(model, start=burn.in, verbose=FALSE)

if (phe_metric == "doy" & direction == "up" & thres == 0.5) {
  type <- "SOS"
} else if (phe_metric == "doy" & direction == "down" & thres == 0.5){
  type <- "EOS"
} else if (phe_metric == "doy_diff" & direction == "up" & thres == 0.2){
  type <- "GreenUp"
} else {
  type <- "GreenDown"
}

model_path <- paste0(r_output_file, "ny_", genus, "_tavg_", type,  ".RData")
model_re_path <- paste0(r_output_file, "ny_", genus, "_tavg_", type, "_re", ".RData")
save(model, file = model_path)
save(m_re, file = model_re_path)
