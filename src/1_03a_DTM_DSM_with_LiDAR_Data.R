#script to derive DTM (digital terrain model) and DSM (digital surface model)
# using LiDAR data

#-------------------------------------------------------------------------------
#1 SET WORKING DIRECTORY AND LOAD PACKAGES
#-------------------------------------------------------------------------------
library(rLiDAR)
library(magrittr)
library(lidR)

#alternative using envimaR
# library(envimaR)
# root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
#                               alt_env_id = "COMPUTERNAME",
#                               alt_env_value = "PCRZP", 
#                               alt_env_root_folder = "F://BEN//edu")
# source(file.path(root_folder, "mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))

workingDir <- "D:/BEN/edu/data/"

#-------------------------------------------------------------------------------
#2 CREATING A NORMALIZED LAS FILE
#-------------------------------------------------------------------------------
# Construct a DTM with a resolution of 1m
lidar_file <- paste0(workingDir, "lidar/test/test_lasdata.las") %>% readLAS()

#file for quicker processing of .las files
paste0(workingDir, "lidar/test/test_lasdata.las") %>% writelax()

plot(lidar_file, bg ="white", color = "Z")

# Use interpolation to normalize
las_norm <- lidR::normalize_height(lidar_file, tin())

plot(las_norm)

#-------------------------------------------------------------------------------
#2 CREATING DSM
#-------------------------------------------------------------------------------
#to use it with a catalog
#lidR::opt_output_files(lcat_file)<-paste0(envrmt$path,"/{ID}filename"
#change lidar_file with lcat_file

dsm <-  grid_canopy(lidar_file, 0.5, p2r(subcircle = 0.2))
plot(dsm)

#-------------------------------------------------------------------------------
#2 CREATING DTM
#-------------------------------------------------------------------------------
dtm <-  grid_terrain(lidar_file,
                   res = 0.5,
                   algorithm = kriging(k = 10))

plot(dtm)
