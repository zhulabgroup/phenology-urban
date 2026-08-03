options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)
options(repos = c(CRAN = "https://cloud.r-project.org"))
pacman::p_load(sf, patchwork, tidyverse, leaflet)

# create_symlink_turbo("release/input/tree-list", "proj-phenology-urban/raw/CONUS/StreetTrees/cleaned_ps_input")
# create_symlink_turbo("release/input/tree-climate", "proj-phenology-urban/proc/tree_climate")
# create_symlink_turbo("release/input/phenology/doy", "datasets/vegetation/PS/CONUS_urban/doy")
# create_symlink_turbo("release/input/phenology/clean", "datasets/vegetation/PS/CONUS_urban/clean")
# create_symlink_turbo("release/input/phenology/no_outlier", "datasets/vegetation/PS/CONUS_urban/no_outlier")
# create_symlink_turbo("release/input/boundary", "datasets/boundary/USAdministrative")
# create_symlink_turbo("release/input/climate/UHAT", "datasets/climate/UHAT")
# create_symlink_turbo("release/input/climate/Daymet", "proj-phenology-urban/raw/CONUS/Daymet/all_city_trees")
# create_symlink_turbo("release/input/climate/CHELSA", "proj-phenology-urban/proc/tree_origin/chelsa_map")
# create_symlink_turbo("release/intermediate/variance_decompose", "proj-phenology-urban/proc/model_output/variance_decompose")
# create_symlink_turbo("release/intermediate/spatial/inla", "proj-phenology-urban/proc/model_output/final_models_clean_outlier/inla")
# create_symlink_turbo("release/intermediate/spatial/mcmc", "proj-phenology-urban/proc/model_output/final_models_clean_outlier/mcmc")
# create_symlink_turbo("release/intermediate/spatial/freq", "proj-phenology-urban/proc/model_output/final_models_clean_outlier/freq")
# create_symlink_turbo("release/intermediate/spatial/inla_gsl", "proj-phenology-urban/proc/model_output/final_models_clean_outlier/inla_gsl")
# create_symlink_turbo("release/output", "proj-phenology-urban/output")

city_label <- c(
  "DT" = "Detroit, MI",
  "TP" = "Tampa, FL",
  "Austin" = "Austin, TX",
  "Denver" = "Denver, CO",
  "Houston" = "Houston, TX",
  "NewYork" = "New York, NY",
  "Seattle" = "Seattle, WA",
  "Albuquerque" = "Albuquerque, NM",
  "Anaheim" = "Anaheim, CA",
  "Arlington" = "Arlington, TX",
  "Atlanta" = "Atlanta, GA",
  "Aurora" = "Aurora, CO",
  "Baltimore" = "Baltimore, MD",
  "Buffalo" = "Buffalo, NY",
  "CapeCoral" = "Cape Coral, FL",
  "ColoradoSprings" = "Colorado Springs, CO",
  "Columbus" = "Columbus, OH",
  "DesMoines" = "Des Moines, IA",
  "Durham" = "Durham, NC",
  "Greensboro" = "Greensboro, NC",
  "Knoxville" = "Knoxville, TN",
  "LasVegas" = "Las Vegas, NV",
  "LosAngeles" = "Los Angeles, CA",
  "Louisville" = "Louisville, KY",
  "Madison" = "Madison, WI",
  "Milwaukee" = "Milwaukee, WI",
  "Minneapolis" = "Minneapolis, MN",
  "Oakland" = "Oakland, CA",
  "OklahomaCity" = "Oklahoma City, OK",
  "Pittsburgh" = "Pittsburgh, PA",
  "Portland" = "Portland, OR",
  "Providence" = "Providence, RI",
  "Rochester" = "Rochester, NY",
  "Sacramento" = "Sacramento, CA",
  "SanDiego" = "San Diego, CA",
  "SanFrancisco" = "San Francisco, CA",
  "SanJose" = "San Jose, CA",
  "SantaRosa" = "Santa Rosa, CA",
  "SiouxFalls" = "Sioux Falls, SD",
  "StLouis" = "St. Louis, MO",
  "Stockton" = "Stockton, CA",
  "WashingtonDC" = "Washington, DC",
  "Allentown" = "Allentown, PA",
  "Amherst" = "Amherst, MA",
  "Boulder" = "Boulder, CO",
  "Bozeman" = "Bozeman, MT",
  "Cambridge" = "Cambridge, MA",
  "Champaign" = "Champaign, IL",
  "Charlottesville" = "Charlottesville, VA",
  "Escondido" = "Escondido, CA",
  "Philadelphia" = "Philadelphia, PA",
  "StAugustine" = "St. Augustine, FL",
  "Davenport" = "Davenport, IA",
  "Spokane" = "Spokane, WA",
  "Bloomington" = "Bloomington, IN",
  "Boston" = "Boston, MA",
  "Mesa" = "Mesa, AZ",
  "NewOrleans" = "New Orleans, LA",
  "Norfolk" = "Norfolk, VA",
  "Tempe" = "Tempe, AZ",
  "Livingston" = "Livingston, MT",
  "Riverside" = "Riverside, CA",
  "Boise" = "Boise, ID",
  "Moab" = "Moab, UT",
  "Murray" = "Murray, UT",
  "Nephi" = "Nephi, UT",
  "StGeorge" = "St. George, UT",
  "Berkeley" = "Berkeley, CA",
  "Ithaca" = "Ithaca, NY",
  "FortCollins" = "Fort Collins, CO",
  "Binghamton" = "Binghamton, NY",
  "Harrisburg" = "Harrisburg, PA",
  "Lancaster" = "Lancaster, PA",
  "Reading" = "Reading, PA",
  "Erie" = "Erie, PA",
  "Newburgh" = "Newburgh, NY"
)
