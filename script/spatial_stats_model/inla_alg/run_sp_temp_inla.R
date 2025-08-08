library(sp)
library(sf)
library(dplyr)
library(stringr)
library(car)
library(gstat)
library(INLA)
library(ggregplot)
library(tidyverse)
library(viridis)
library(cowplot)
library(fmesher)
library(inlabru)
library(patchwork)

# options(repos = c(
#   INLA = "https://inla.r-inla-download.org/R/testing",
#   getOption("repos")
# ))
# install.packages("inlabru")
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

r_input_dir <- "~/phenology-urban/data/proc/urban/"
r_output_dir <- "~/phenology-urban/data/proc/urban/sp_model/"
city <- "NY"
phe_type <- "SOS"
# genus = "sp_random_yr"
genus <- "Acer"


######################################## Read fit data ######################################## 

ny_boundary <- st_read("~/phenology-urban/data/raw/NY/boundary/nybb_dissolved.shp")

if (phe_type == "SOS") {
  phe_metric <- "doy"
  direction <- "up"
  thres <- 0.5
} else if (phe_type == "EOS") {
  phe_metric <- "doy"
  direction <- "down"
  thres <- 0.5
} else if (phe_type == "GreenUp") {
  phe_metric <- "doy_diff"
  direction <- "up"
  thres <- 0.2
  tavg_allyear_var$doy_diff <- as.numeric(tavg_allyear_var$doy_diff)
} else {
  phe_metric <- "doy_diff"
  direction <- "down"
  thres <- 0.2
  tavg_allyear_var$doy_diff <- as.numeric(tavg_allyear_var$doy_diff)
}


tavg_allyear_var_sf <- readRDS(paste0(r_input_dir, city, "/tavg_allyear_var.rds")) %>%
  filter(direction == !!direction & thres == !!thres & genus == !!genus) %>%
  na.omit() %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = st_crs(ny_boundary)) %>%
  mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2], year = as.numeric(as.factor(year))) 

tavg_allyear_var <- tavg_allyear_var_sf%>% 
  st_drop_geometry()

print(paste("City:",city, "Genus:",genus, "; metric:", nrow(tavg_allyear_var)))
print(paste("Type:", phe_type, "; metric:", phe_metric, "; direction:", direction, "; thres:", thres))



######################################## Read prediction data ########################################

### Prediction data (Acer in 2023 for WU location)

mean_doy_acer <- readRDS(paste0("~/phenology-urban/data/proc/urban/",city,"/tavg_allyear_var.rds")) %>%
  filter(genus == !!genus & thres == !!thres & direction ==!!direction ) %>%
  summarise(doy = mean(doy)) %>%
  pull(doy)
preseason <- readRDS(paste0("~/phenology-urban/data/proc/urban/",city,"/tavg_allyear_pcorr.rds")) %>%
  filter(genus == !!genus & type == "SOS" ) %>%
  summarise(pre = mean(optimal_pre_length)) %>%
  pull(pre)
start_date = make_date(2023) + days(round(mean_doy_acer) - 1) - days(preseason)
end_date = make_date(2023) + days(round(mean_doy_acer) - 1)

NY_weather <- readRDS(paste0("~/urban-cooling/data/raw/WU/",city,"/",city,"_wu.rds")) %>%
  filter(Date >= start_date & Date <= end_date) %>%
  group_by(name) %>%
  summarise(AvgTemp_mean = mean(AvgTemp),
            Sum_mm_sum = sum(Sum_mm),
            year = 7)

NY_srad <- readRDS(paste0("~/phenology-urban/data/raw/",city,"/Daymet/daily2016-2024/srad_wu_2023.rds")) %>%
  filter(date >= start_date & date <= end_date) %>%
  group_by(Site) %>%
  summarise(srad_mean = mean(value)) %>%
  inner_join(NY_weather, by = c("Site" = "name"))

prediction <- read_csv(paste0("~/urban-cooling/data/raw/WU/",city,"/location.csv"), show_col_types = FALSE) %>%
  inner_join(NY_srad, by = "Site") %>%
  na.omit()

######################################## Prepare the mesh and spde ######################################## 

### Build the mesh 
coord <- cbind(tavg_allyear_var$lon, tavg_allyear_var$lat) # coordinates of observations
domain <- inla.nonconvex.hull(as.matrix(coord),
                              concave = -.07,convex = -0.05, resolution=c(100,100))

mesh <- inla.mesh.2d(boundary = domain, loc = coord,
                     max.edge = c(10000, 20000), 
                     cutoff = 1000 # prevent close points from being connected
)
# mesh$n
ggplot() +
  geom_sf(data=ny_boundary,color='turquoise',fill='transparent') +
  gg(mesh) +
  geom_sf(data=tavg_allyear_var_sf,col='purple',size=1.7,alpha=0.5)

### Build the spde 
Groups = "year"
NGroups <- nrow(unique(tavg_allyear_var[,Groups])) 

A.phe <- inla.spde.make.A(mesh = mesh, loc = coord, group = tavg_allyear_var$year, n.group = NGroups)
spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
s.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde, n.group = NGroups)
dim(s.index)

coorp <- cbind(prediction$Lon, prediction$Lat) 
A.phe.p <- inla.spde.make.A(mesh = mesh, loc = coorp,group = prediction$year, n.group = NGroups)
dim(A.phe)

######################################## Stack the data ########################################

stk.e <- inla.stack(
  tag = "est",
  data = list(y = tavg_allyear_var[[phe_metric]]),
  A = list(A.phe, 1),
  # effects = list(data.frame(Intercept = rep(1, nrow(coord)),
  #                           AvgTemp_mean = tavg_allyear_var$AvgTemp_mean,
  #                           Sum_mm_sum = tavg_allyear_var$Sum_mm_sum,
  #                           srad_mean = tavg_allyear_var$srad_mean),
  #                s = s.index)
  effects = list(c(s.index,
                   list(Intercept=1)),
                 list(tavg_allyear_var[,c(17,21,23)]))
)


stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(A.phe.p,1),
  # effects = list(data.frame(Intercept = rep(1, nrow(coorp)),
  #                           AvgTemp_mean = prediction$AvgTemp_mean,
  #                           Sum_mm_sum = prediction$Sum_mm_sum,
  #                           srad_mean = prediction$srad_mean,
  #                           year = prediction$year), 
  #                s = s.index)
  effects=list(c(s.index, 
                 list(Intercept=1)),
               list(prediction[,c(5,6,4)]))
)

# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.p)

######################################## Run/Save the model ########################################

### Run the model
formula <- y ~ 0 + Intercept + AvgTemp_mean + Sum_mm_sum + srad_mean + f(spatial.field, model = spde, 
                                                                         group = spatial.field.group,
                                                                         control.group = list(model = "rw1"))

res <- inla(formula,
            family = "gaussian",
            data = inla.stack.data(stk.full),
            control.compute = list(return.marginals = TRUE, waic = TRUE),
            control.predictor = list(
              compute = TRUE,
              A = inla.stack.A(stk.full)
            )
)
res$waic$waic
summary(res)

######################################## Inference from model ########################################

########################################
# Spatial field
########################################
# projection matrix for 300 x 300 locations on a grid that covers the mesh region.

proj <- inla.mesh.projector(mesh,dims = c(300, 300))
mean_s <- inla.mesh.project(proj, res$summary.random$spatial.field$mean)
sd_s <- inla.mesh.project(proj, res$summary.random$spatial.field$sd)
df <- expand.grid(x = proj$x, y = proj$y)
df$mean_s <- as.vector(mean_s)
df$sd_s <- as.vector(sd_s)

gmean <- ggplot(df, aes(x = x, y = y, fill = mean_s)) +
  geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

gsd <- ggplot(df, aes(x = x, y = y, fill = sd_s)) +
  geom_raster() +
  scale_fill_viridis(na.value = "transparent") +
  coord_fixed(ratio = 1) + theme_bw()

gmean + gsd

#Plot the decay of spatial autocorrelation across a user-defined range

spde.est <- inla.spde2.result(inla = res, name = "spatial.field",
                              spde = spde, do.transf = TRUE)
Kappa = inla.zmarginal(spde.est$marginals.kappa[[1]])[1]$mean
MaxRange <- 50000
Resolution <- 100
d.vec <- seq(0, MaxRange, length = Resolution)
Cor.M <- (Kappa * d.vec) * besselK(Kappa * d.vec, 1)
Cor.M[1] <- 1
CorData <- data.frame(Distance = d.vec, Correlation = Cor.M)
ggplot(CorData, aes(x = Distance, y = Correlation)) +
  geom_line(linewidth = 1, alpha = 0.8, color = "blue") +
  coord_fixed(ratio = MaxRange) +
  labs(x = "Distance (m)", y = "Correlation", title = "Matern Correlation Function") +
  theme_minimal()

########################################
# Random year effect
########################################

res$summary.random$year

########################################
# Fixed effect
########################################

# for one genus

beta_names <- c("beta[0]", "beta[AvgTemp_mean]", "beta[Sum_mm_sum]", "beta[srad_mean]")
df_list <- lapply(seq_along(res$marginals.fixed), function(i) {
  data.frame(x = res$marginals.fixed[[i]][,1], 
             density = res$marginals.fixed[[i]][,2], 
             beta = beta_names[i])
})
df <- bind_rows(df_list)
ggplot(df, aes(x = x, y = density)) +
  geom_line(linewidth = 1, color = "blue") +
  facet_wrap(~beta, scales = "free", labeller = label_parsed) + 
  theme_minimal() +
  labs(x = "Coefficient Value", y = "Density") +
  theme(legend.position = "none")  # Remove legend since we use facets

# for comparison among genus

fix_effect <- as.data.frame(summary(res)$fixed)%>%
  rename(Estimate = mean,Lower = `0.025quant`,Upper = `0.975quant`) %>%
  mutate(Sig = ifelse(Lower*Upper>0, "*", "")) 
fix_effect <- data.frame(
  Factor = rownames(fix_effect),  # Add Factor column
  fix_effect, 
  row.names = NULL  # Remove row names to avoid duplication
) %>%
  filter(Factor != "Intercept")

SigAlpha <- ifelse(fix_effect$Sig == "*", 1,0.5)

ggplot(fix_effect, aes(x = factor(Factor, levels = rev(Factor)), 
                       y = Estimate, 
                       colour = Factor)) +
  geom_point(alpha = SigAlpha, position = position_dodge(w = 0.5), size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, alpha = SigAlpha), 
                position = position_dodge(w = 0.5), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = NULL, y = "Estimate") +
  coord_flip() +
  scale_alpha_identity() +
  theme(legend.position = "right") +
  geom_text(aes(label = Sig, y = fix_effect$Upper + 0.05  ), hjust = -0.5, show.legend = FALSE) +
  guides(alpha = "none") +
  theme_minimal()

######################################## Prediction from model ########################################

### Prediction

index <- inla.stack.index(stack = stk.full, tag = "pred")$data

pred_mean <- res$summary.fitted.values[index, "mean"]
pred_ll <- res$summary.fitted.values[index, "0.025quant"]
pred_ul <- res$summary.fitted.values[index, "0.975quant"]

pred_ph_2023 <- rbind(
  data.frame(
    Lon = coorp[, 1], Lat = coorp[, 2],
    value = pred_mean, variable = "pred_mean"
  ),
  data.frame(
    Lon = coorp[, 1], Lat = coorp[, 2],
    value = pred_ll, variable = "pred_ll"
  ),
  data.frame(
    Lon = coorp[, 1], Lat = coorp[, 2],
    value = pred_ul, variable = "pred_ul"
  )
)
pred_ph_2023$variable <- as.factor(pred_ph_2023$variable)

ggplot(pred_ph_2023) +
  geom_point(aes(Lon, Lat, color  = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_color_gradient(
    name = "SOS in 2023",
    low = "blue", high = "orange"
  ) +
  theme_bw()
