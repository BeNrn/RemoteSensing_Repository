# Set working directory and load required packages
library(envimaR)
root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F:\\BEN\\edu")
source(path.expand("~/edu/mpg-envinsys-plygrnd/mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))


# START HERE IF YOU ARE HERE THE FIRST TIME
# AOI RGB
list.files(envrmt$path_aerial_processed)
img_org <- brick(file.path(envrmt$path_aerial_processed, "AOIimage.tif" ))
R <- img_org[[1]]
G <- img_org[[2]]
B <- img_org[[3]]

# Rename R G B
names(R) <- "AOIimage_R"
names(G) <- "AOIimage_G"
names(B) <- "AOIimage_B"

# Stack indices
files <- list.files(envrmt$path_aerial_processed, pattern = "AOI_index", full.names = T)
files
ind_rasters <- list()

for (i in seq(1:length(files))){
  if (i %% 2 != 0){
    ind_rasters[i] <- raster(files[i])
  }else {
    print("file does not need to be loaded")
  }
}
rm(files)
# Clean ind_rasters & stack indices + R + G + B
library(rlist)
ind_r <- list.clean(ind_rasters)
ind_rgb <- c(ind_r, R, G, B)
ind_stack <- stack(ind_rgb)
ind_stack
rm(ind_rasters)
rm(ind_r)
rm(ind_rgb)
rm(R)
rm(G)
rm(B)

# Normalize raster using RStoolbox
library(RStoolbox)
norm <- normImage(ind_stack, norm = T)
##norm

# Rename layers of norm
names(norm) <- names(ind_stack)
rm(ind_stack)

# Save norm
#writeRaster(norm, filename = file.path(envrmt$path_norm, "AOI_indices_rgb_norm.tif"))
filenames <-  paste0(names(norm),"_norm.tif")
writeRaster(norm, filename = file.path(envrmt$path_norm, filenames), bylayer = T)






# START HERE IF THE UPPER PART HAS BEEN EXECUTED BEFORE
# Read norm
norm <- raster::stack(file.path(envrmt$path_norm, "AOI_indices_rgb_norm.tif"))

# Change resolution of norm to xres=0.5, yres=0.5
res(norm) <- 0.5



# PCA
library(RStoolbox)
pca <- RStoolbox::rasterPCA(norm, nComp = 4, spca = T)
saveRSTBX(pca, filename = file.path(envrmt$path_aerial_processed, "AOI_PCA.tif"))

# Load PCA (if already created)
list.files(envrmt$path_aerial_processed)
pca <- readRSTBX(file.path(envrmt$path_aerial_processed, "AOI_PCA"))
pca
plot(pca$map)
