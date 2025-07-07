library(INLA)
library(Matrix)
library(ggplot2)
library(dplyr)

set.seed(123)

simulate_covariate_with_target_correlation <- function(
    spatial_component, x_coord, 
    cor_spatial = 0.5, cor_trend = 0.5, 
    scale_output = TRUE, range_min = -10, range_max = 10) {
  
  n <- length(spatial_component)
  
  # 1. Standardize spatial field and trend
  spatial_std <- scale(spatial_component)[, 1]
  trend_std <- scale(x_coord)[, 1]
  
  # 2. Create orthogonal noise (independent of spatial and trend)
  X <- cbind(spatial_std, trend_std)
  Q <- diag(n) - X %*% solve(t(X) %*% X) %*% t(X)  # residual-maker matrix
  noise <- Q %*% rnorm(n)
  noise <- scale(noise)[, 1]  # re-standardize
  
  # 3. Solve for coefficients that yield desired correlations
  # We want: var(result) = 1 and cor(result, spatial) = cor_spatial, etc.
  # Final variable: temp_raw = a * spatial_std + b * trend_std + c * noise
  a <- cor_spatial
  b <- cor_trend
  c <- sqrt(1 - a^2 - b^2)
  
  temp_raw <- a * spatial_std + b * trend_std + c * noise
  
  # 4. Scale to desired range
  if (scale_output) {
    temp_scaled <- (temp_raw - min(temp_raw)) / (max(temp_raw) - min(temp_raw))
    temp_scaled <- range_min + temp_scaled * (range_max - range_min)
    return(temp_scaled)
  } else {
    return(temp_raw)
  }
}


## 1. Simulate spatial coordinates and mesh
n <- 500
coords <- cbind(runif(n, 0, 10), runif(n, 0, 10))
years <- sample(2010:2015, n, replace = TRUE)

mesh <- inla.mesh.2d(loc = coords, cutoff = 0.1, max.edge = c(0.5, 1))
spde <- inla.spde2.matern(mesh, alpha = 2)
Q <- inla.spde2.precision(spde, theta = c(log(1), log(0.2)))
u <- as.vector(t(chol(solve(Q))) %*% rnorm(spde$n.spde))
A <- inla.spde.make.A(mesh, loc = coords)
spatial_field <- as.vector(A %*% u)
x_coord <- coords[, 1]
y_coord <- coords[, 2]

## 2. Simulate covariates
# Temperature: linear trend + spatial field + noise
temp <- simulate_covariate_with_target_correlation(
  spatial_component = spatial_field,
  x_coord = x_coord,
  cor_spatial = 0,
  cor_trend = 0,
  range_min = -10,
  range_max = 10
)

# Precipitation: spatial field (smooth) + noise
# simulate separate spatial field for prcp
prcp <- simulate_covariate_with_target_correlation(
  spatial_component = spatial_field,
  x_coord = x_coord,
  cor_spatial = 0,
  cor_trend = 0,
  range_min = 0,
  range_max = 500
)

# Solar radiation: spatial trend (e.g., higher in south) + smooth field
srad <- simulate_covariate_with_target_correlation(
  spatial_component = spatial_field,
  x_coord = x_coord,
  cor_spatial = 0,
  cor_trend = 0,
  range_min = 300,
  range_max = 450
)


## 3. Simulate year random effect
year_levels <- sort(unique(years))
year_effects <- rnorm(length(year_levels), sd = 1)
names(year_effects) <- year_levels
year_re <- year_effects[as.character(years)]

## 4. Simulate response y
X <- cbind(1, temp, precp, srad)
beta <- c(120, -2, 0.5, 0.5)
mu <- numeric(n)
y <- numeric(n)
tau <- 1 
for (i in 1:n) {
  mu[i] <- X[i,] %*% beta + spatial_field[i] + year_re[i]
  y[i] <- rnorm(1, mean = mu[i], sd = sqrt(1/tau))
}

## 5. Data setup
data_all <- data.frame(
  y = y, temp = temp, prcp = prcp, srad = srad, year = years
)

## 6. Model 1: No spatial random effect
model1 <- inla(
  y ~ temp + prcp + srad + f(year, model = "iid"),
  data = data_all,
  family = "gaussian",
  control.compute = list(dic = TRUE)
)

## 7. Model 2: With spatial random effect
s.index <- inla.spde.make.index("spatial", spde$n.spde)
stack <- inla.stack(
  tag = "est",
  data = list(y = y),
  A = list(1, A),
  effects = list(
    data.frame(temp = temp, prcp = prcp, srad = srad, year = years),
    s = s.index
  )
)

inla.setOption(num.threads = 1)
model2 <- inla(
  y ~ temp + prcp + srad + f(year, model = "iid") + f(spatial, model = spde),
  data = inla.stack.data(stack),
  family = "gaussian",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.compute = list(dic = TRUE),
  control.inla = list(strategy = "simplified.laplace", int.strategy = "grid"),
  verbose = TRUE  
)

## 8. Compare temp coefficients
beta1 <- model1$summary.fixed["temp", ]
beta2 <- model2$summary.fixed["temp", ]

print("Temp coefficient estimates:")
print(data.frame(
  Model = c("No Spatial", "With Spatial"),
  Mean = c(beta1$mean, beta2$mean),
  Lower = c(beta1$`0.025quant`, beta2$`0.025quant`),
  Upper = c(beta1$`0.975quant`, beta2$`0.975quant`)
))

## 9. Optional: check spatial confounding
spatial_field_mean <- model2$summary.random$spatial$mean
fitted_spatial <- as.vector(A %*% spatial_field_mean)
cor_temp_spatial <- cor(temp, fitted_spatial)
print(paste("Correlation between temp and spatial field:", round(cor_temp_spatial, 3)))

