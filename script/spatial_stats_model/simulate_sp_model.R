library(MASS) 
library(ggplot2) 
library(dplyr)
library(sp)

set.seed(2025)
lon_range <- c(-74.25156, -73.70623)
lat_range <- c(40.50458, 40.91093)

N <- 1000       # Total points
N_unique <- 500 # Unique locations
N_duplicate <- N - N_unique


lon_unique <- runif(N_unique, min = lon_range[1], max = lon_range[2])
lat_unique <- runif(N_unique, min = lat_range[1], max = lat_range[2])
duplicate_indices <- sample(1:N_unique, N_duplicate, replace = TRUE)
lon <- c(lon_unique, lon_unique[duplicate_indices])
lat <- c(lat_unique, lat_unique[duplicate_indices])
years <- 2017:2023
nyears <- length(years)

# Assign years to unique locations first
year_unique <- sample(years, N_unique, replace = TRUE)
year_duplicate <- sapply(duplicate_indices, function(idx) {
  sample(setdiff(years, year_unique[idx]), 1)  # Ensure different year
})
year <- c(year_unique, year_duplicate)

# fix effects: intercept, tavg, precp, srad
intercept <- 120
tavg <- runif(N, 0, 25)
precp <- runif(N, 0, 500)
srad <- runif(N, 300, 450)

beta <- c(intercept, -2, 0.05, 0.05)
X <- cbind(1, tavg, precp, srad)

# spatial effect
sigma.sq <- 80
phi <- 25
spatial_effect <- numeric(N)

for (y in years) {
  year_idx <- which(year == y)
  
  if (length(year_idx) > 1) {
    lon_year <- lon[year_idx]
    lat_year <- lat[year_idx]
    dist_matrix_year <- as.matrix(dist(cbind(lon_year, lat_year)))
    cov_matrix_year <- sigma.sq * exp(-phi * dist_matrix_year)
    spatial_effect[year_idx] <- mvrnorm(1, mu = rep(0, length(year_idx)), Sigma = cov_matrix_year)
  } else {
    spatial_effect[year_idx] <- 0
  }
}

# year effect
year_effect <- rnorm(nyears, mean = 0, sd = 1)

mu <- numeric(N)
y <- numeric(N)
tau <- 1 
for (i in 1:N) {
  year_idx <- which(years == year[i])
  mu[i] <- X[i,] %*% beta + spatial_effect[i] + year_effect[year_idx]
  y[i] <- rnorm(1, mean = mu[i], sd = sqrt(1/tau))
}

sim_data <- data.frame(
  lon = lon,
  lat = lat,
  year = year,
  temp = tavg,
  prcp = precp,
  srad = srad,
  spatial_effect = spatial_effect,
  doy = y
)

# saveRDS(sim_data, "~/phenology-urban/data/proc/urban/SIM/tavg_allyear_var.rds")
