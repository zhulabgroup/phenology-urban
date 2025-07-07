data_folder <- "~/phenology-urban/data/raw/NY/HUMID"
file_list <- list.files(data_folder, pattern = "^NY_\\d{6}\\.txt$", full.names = TRUE)

# Read and combine all files
all_data <- map_dfr(file_list, function(file_path) {
  # Extract date string like "2016-01" from file name
  file_name <- basename(file_path)
  date_str <- str_extract(file_name, "\\d{6}") %>%
    str_replace("(\\d{4})(\\d{2})", "\\1-\\2")
  
  # Read file and add date
  read_tsv(file_path, show_col_types = FALSE) %>%
    mutate(date = date_str)
})

# Reorganize into separate long tables
temp_monthly <- all_data %>% select(id, value = temp, date)
prcp_monthly <- all_data %>% select(id, value = prcp, date)
srad_monthly <- all_data %>% select(id, value = srad, date)

# Save to .rds
saveRDS(temp_monthly, file = file.path(data_folder, "temp_monthly_all_tree.rds"))
saveRDS(prcp_monthly, file = file.path(data_folder, "prcp_monthly_all_tree.rds"))
saveRDS(srad_monthly, file = file.path(data_folder, "srad_monthly_all_tree.rds"))
