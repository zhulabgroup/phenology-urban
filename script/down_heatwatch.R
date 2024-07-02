library(httr)
library(jsonlite)
library(osfr)
library(dplyr)
library(stringr)

id_list <- c()

id_name <- c()

for (i in 1:6){
  api_url <- paste0("https://api.osf.io/v2/users/9neka/nodes/?page=", i)
  response <- GET(api_url)
  
  if (response$status_code == 200) {
    raw_content <- content(response, "text", encoding = "UTF-8")
    json_content <- fromJSON(raw_content, flatten = TRUE)
    id_page <- json_content$data$id
    id_name_page <- json_content$data$attributes.title
    print(id_name_page)
    id_list <- c(id_list, id_page)
  } else {
    stop("No response")
  }
}

for (node in 55:59){
  target_file <- osf_retrieve_node(id_list[node]) %>%
    osf_ls_files()
  if (any(str_detect(target_file$name, "All Data"))){
    target_file <- target_file %>%
      filter(str_detect(name, "All Data"))
    for (id in target_file$id){
      osf_retrieve_file(id) %>% osf_download(path="~/my_turbo_data/CONUS/Heat_Watch")
    }
    
  }
  else{
    combined_table <- data.frame()
    target_file <- target_file %>% filter(!str_detect(name, "\\.txt$"))
    for (id in target_file$id){
      table <- osf_retrieve_file(id) %>% osf_ls_files() %>%
        select(name, id)
      combined_table <- bind_rows(combined_table, table)
    }
    target_file <- combined_table %>%
      filter(str_detect(name, "All Data"))
    
    for (id in target_file$id){
      osf_retrieve_file(id) %>% osf_download(path="~/my_turbo_data/CONUS/Heat_Watch")
    }
  }
}
