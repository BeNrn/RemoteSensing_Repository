#Set working directory and load required packages
library(envimaR)
library(beepr)
root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F:\\BEN\\edu")

source(file.path(root_folder, "/mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))

# Load SHP's
list.files(envrmt$path_run)
All <- readOGR(file.path(envrmt$path_run, "Trainingsgebiete_Einzelbaumbestimmung.shp"),
               layer = ogrListLayers(paste(envrmt$path_run, "Trainingsgebiete_Einzelbaumbestimmung.shp", sep = "/")))

Buche <- readOGR(file.path(envrmt$path_run, "Trainingsgebiete_Buche.shp"),
              layer =  ogrListLayers(paste(envrmt$path_run, "Trainingsgebiete_Buche.shp", sep = "/")))
Fichte <- readOGR(file.path(envrmt$path_run, "Trainingsgebiete_Fichte.shp"),
              layer =  ogrListLayers(paste(envrmt$path_run, "Trainingsgebiete_Fichte.shp", sep = "/")))
Eiche <- readOGR(file.path(envrmt$path_run, "Trainingsgebiete_Eiche.shp"),
              layer =  ogrListLayers(paste(envrmt$path_run, "Trainingsgebiete_Eiche.shp", sep = "/")))
Douglasie <- readOGR(file.path(envrmt$path_run, "Trainingsgebiete_Douglasie.shp"),
               layer =  ogrListLayers(paste(envrmt$path_run, "Trainingsgebiete_Douglasie.shp", sep = "/")))
Baumarten <- list(Buche, Fichte, Eiche, Douglasie)


# Load rsts
## Load RGB
list.files(envrmt$path_aerial_processed)
img_org <- brick(file.path(envrmt$path_aerial_processed, "AOIimage.tif"))

### Split into the 3 different visible colour-layers
lyr <- img_org[[1]]
lyg <- img_org[[2]]
lyb <- img_org[[3]]
rm(img_org)

## Load Indices
CI <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_CI.tif"))
GLI <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_GLI.tif"))
IO <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_IO.tif"))
NGRDI <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_NGRDI.tif"))
RGR <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_red_green_ratio.tif"))
SCI <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_SCI.tif"))
TGI <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_TGI.tif"))
VARI <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_VARI.tif"))
VVI <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_VVI.tif"))

PCA_gauss_3x3 <- raster(file.path(envrmt$path_aerial_processed,"AOI_PCA_filter _3x3_gauss.grd"))
PCA_gauss_5x5 <- raster(file.path(envrmt$path_aerial_processed, "AOI_PCA_filter _5x5_gauss.grd"))
PCA_gauss_7x7 <- raster(file.path(envrmt$path_aerial_processed, "AOI_PCA_filter _7x7_gauss.grd"))

PCA_glcm_3x3 <- raster(file.path(envrmt$path_aerial_processed, "AOI_PCA_filter _3x3_glcm.tif"))
PCA_glcm_5x5 <- raster(file.path(envrmt$path_aerial_processed, "AOI_PCA_filter _7x7_glcm.tif"))
PCA_glcm_15x15 <- raster(file.path(envrmt$path_aerial_processed, "AOI_PCA_filter _15x15_glcm.tif"))
PCA_glcm_31x31 <- raster(file.path(envrmt$path_aerial_processed, "AOI_PCA_filter _31x31_glcm.tif"))

PCA_mean_3x3 <- raster(file.path(envrmt$path_aerial_processed, "AOI_PCA_filter _3x3_mean.tif" ))
PCA_mean_5x5 <- raster(file.path(envrmt$path_aerial_processed, "AOI_PCA_filter _5x5_mean.tif"))
PCA_mean_7x7 <- raster(file.path(envrmt$path_aerial_processed, "AOI_PCA_filter _7x7_mean.tif"))
PCA_mean_15x15 <- raster(file.path(envrmt$path_aerial_processed, "AOI_PCA_filter _15x15_mean.tif"))
PCA_mean_31x31 <- raster(file.path(envrmt$path_aerial_processed, "AOI_PCA_filter _31x31_mean.tif"))

PCA_sd_3x3 <- raster(file.path(envrmt$path_aerial_processed,"AOI_PCA_filter _3x3_sd.tif"))
PCA_sd_5x5 <- raster(file.path(envrmt$path_aerial_processed, "AOI_PCA_filter _5x5_sd.tif" ))
PCA_sd_7x7 <- raster(file.path(envrmt$path_aerial_processed, "AOI_PCA_filter _7x7_sd.tif"))
PCA_sd_15x15 <- raster(file.path(envrmt$path_aerial_processed, "AOI_PCA_filter _15x15_sd.tif"))







rsts <- stack(lyr, lyg, lyb, CI, GLI, IO, NGRDI, RGR, SCI, TGI, VARI, VVI, 
              PCA_gauss_3x3, PCA_gauss_5x5, PCA_gauss_7x7, 
              PCA_mean_3x3, PCA_mean_5x5, PCA_mean_7x7, PCA_mean_15x15, 
              PCA_mean_31x31,
              PCA_glcm_3x3, PCA_glcm_5x5, PCA_glcm_15x15, PCA_glcm_31x31,
              PCA_sd_3x3, PCA_sd_5x5, PCA_sd_7x7, PCA_sd_15x15)

# Reproject Baumarten
projection <- VVI@crs
Buche@proj4string <- projection
Fichte@proj4string <- projection
Eiche@proj4string <- projection
Douglasie@proj4string <- projection



Bu1 <- raster::crop(raster::mask(VVI, Buche[1,]), Buche[1,])
writeRaster(Bu1, filename = file.path(envrmt$path_run, "test"))
plot(Bu1)
unique(Bu1)


for (i in seq(1: length(Buche@polygons))){
    Bu <- raster::crop(raster::mask(rsts, Buche[i,]), Buche[i,])
    beepr::beep(1)
    writeRaster(Bu, filename = file.path(envrmt$path_aerial_processed, paste0("Testgebiet_Buche_Baum_", i)))
    beepr::beep(2)
  }
beepr::beep(3)

for (i in seq(1: length(Fichte@polygons))){
  Bu <- raster::crop(raster::mask(rsts, Fichte[i,]), Fichte[i,])
  beepr::beep(1)
  writeRaster(Bu, filename = file.path(envrmt$path_aerial_processed, paste0("Testgebiet_Fichte_Baum_", i)))
  beepr::beep(2)
}
beepr::beep(3)

for (i in seq(1: length(Eiche@polygons))){
  Bu <- raster::crop(raster::mask(rsts, Eiche[i,]), Eiche[i,])
  beepr::beep(1)
  writeRaster(Bu, filename = file.path(envrmt$path_aerial_processed, paste0("Testgebiet_Eiche_Baum_", i)))
  beepr::beep(2)
}
beepr::beep(3)

for (i in seq(1: length(Douglasie@polygons))){
  Bu <- raster::crop(raster::mask(rsts, Douglase[i,]), Douglasie[i,])
  beepr::beep(1)
  writeRaster(Bu, filename = file.path(envrmt$path_aerial_processed, paste0("Testgebiet_Douglasie_Baum_", i)))
  beepr::beep(2)
}
beepr::beep(3)



