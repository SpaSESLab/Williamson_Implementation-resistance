library(googledrive)
options(
  gargle_oauth_cache = ".secrets",
  gargle_oauth_email = TRUE
)


# Download PADUS ----------------------------------------------------------

folder_url <- "https://drive.google.com/drive/u/0/folders/15Z1A96UW43uvfi3O3DZDDerfhgr-lqES"
folder <- drive_get(as_id(folder_url))

gdrive_files <- drive_ls(folder)
#have to treat the gdb as a folder and download it into a gdb directory in order to deal with the fact that gdb is multiple, linked files
lapply(gdrive_files$id, function(x) drive_download(as_id(x), 
                                                   path = paste0(here::here("data/OriginalData/PAD_US2_1.gdb/"), gdrive_files[gdrive_files$id==x,]$name), overwrite = TRUE))


# Download Nolte et al data -----------------------------------------------
folder_url <- "https://drive.google.com/drive/u/0/folders/1oaEWHL3RWrMHZQyUEa4Iz4pSn6HsO4ri"
folder <- drive_get(as_id(folder_url))

gdrive_files <- drive_ls(folder)
lapply(gdrive_files$id, function(x) drive_download(as_id(x), 
                                                   path = paste0(here::here("data/OriginalData/LandValue/"), gdrive_files[gdrive_files$id==x,]$name), overwrite = TRUE))


# Download Theobald hmi ---------------------------------------------------
folder_url <- "https://drive.google.com/drive/u/0/folders/1uBKCAeZ4VStT7N0HYy-yd0IbEyrtEdtE"
folder <- drive_get(as_id(folder_url))

gdrive_files <- drive_ls(folder)
lapply(gdrive_files$id, function(x) drive_download(as_id(x), 
                                                   path = paste0(here::here("data/OriginalData/HMI/hm_fsum3_270/"), gdrive_files[gdrive_files$id==x,]$name), overwrite = TRUE))


folder_url <- "https://drive.google.com/drive/u/0/folders/1rugX-fF_qG68ip6gCWW6UeIPe1TizjHI"
folder <- drive_get(as_id(folder_url))

gdrive_files <- drive_ls(folder)
#have to tdownload the info folder so that R knows how to deal with the ArcGIS based file
lapply(gdrive_files$id, function(x) drive_download(as_id(x), 
                                                   path = paste0(here::here("data/OriginalData/HMI/info/"), gdrive_files[gdrive_files$id==x,]$name), overwrite = TRUE))

folder_url <- "https://drive.google.com/drive/u/0/folders/1YqE8yKcY8WTXkUYK3qnKY-hCGeKOwmEd"
folder <- drive_get(as_id(folder_url))

gdrive_files <- drive_ls(folder)
#have to treat the gdb as a folder and download it into a gdb directory in order to deal with the fact that gdb is multiple, linked files
lapply(gdrive_files$id[3:4], function(x) drive_download(as_id(x), 
                                                        path = paste0(here::here("data/OriginalData/HMI/"), gdrive_files[gdrive_files$id==x,]$name), overwrite = TRUE))


# download wolf preference data -------------------------------------------
folder_url <- "https://drive.google.com/drive/u/0/folders/1fcOrG_CSVf3wS4k4iVjYTUXrAZhJj-G_"
folder <- drive_get(as_id(folder_url))

gdrive_files <- drive_ls(folder)
lapply(gdrive_files$id, function(x) drive_download(as_id(x), 
                                                   path = paste0(here::here("data/OriginalData/wolves/"), gdrive_files[gdrive_files$id==x,]$name), overwrite = TRUE))


