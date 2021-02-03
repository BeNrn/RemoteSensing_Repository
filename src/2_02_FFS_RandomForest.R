#Set working directory and load required packages
#.libPaths("F:/lib")

require(uavRst)
require(envimaR)
require(caret)
require(CAST)
library(doParallel)
require(stringr)
root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F:\\BEN\\edu")
source(file.path(root_folder, "/mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))
#-------------------------------------------------------------------------------
#read training shapes
dgl <- readOGR(file.path(envrmt$path_tmp, "Douglasie_Segmente.shp"))
bu <- readOGR(file.path(envrmt$path_tmp, "buche.shp"))
fi <- readOGR(file.path(envrmt$path_tmp, "Fichte.shp"))
ei <- readOGR(file.path(envrmt$path_tmp, "Trainingsgebiete_Eiche.shp"))
names(bu)[3] <- "id"
bu$id <- as.factor(bu$id)
all_spec <- rbind(dgl, bu, ei, fi)
rm(dgl, bu, ei, fi)

writeOGR(all_spec, file.path(envrmt$path_species_det, "species_poly.shp"),
         driver = "ESRI Shapefile",
         layer = "species_poly")

#for later call (FLC)
#all_spec <- readOGR(file.path(envrmt$path_tmp, "species_poly.shp"))

df_dgl <- data.frame(id = dgl@data$id, art = dgl@data$Baumart, abt = dgl@data$Abteilung)
df_bu <- data.frame(id = bu@data$id, art = bu@data$Baumart, abt = bu@data$Abteilung)
df_ei <- data.frame(id = ei@data$id, art = ei@data$Baumart, abt = ei@data$Abteilung)
df_fi <- data.frame(id = fi@data$id, art = fi@data$Baumart, abt = fi@data$Abteilung)
df_bu$id <- as.factor(df_bu$id)

df <- rbind(df_dgl, df_bu, df_ei, df_fi)

#-------------------------------------------------------------------------------
#read RGB
rgb <- raster::brick(file.path(envrmt$path_aerial_processed, "AOIimage.tif"))

#read shanno index
shn_list <- list.files(envrmt$path_lidar_processed, pattern = "shannon", full.names = TRUE)

shn_list <- shn_list[which(stringr::str_detect(shn_list, "ovr") == FALSE)]
shn_list <- shn_list[which(stringr::str_detect(shn_list, "aux") == FALSE)]
shn <- stack(shn_list)
rm(shn_list)

#read vegetation indices
CI <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_CI.tif"))
GLI <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_GLI.tif"))
IO <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_IO.tif"))
NGRDI <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_NGRDI.tif"))
RGR <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_red_green_ratio.tif"))
SCI <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_SCI.tif"))
TGI <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_TGI.tif"))
VARI <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_VARI.tif"))
VVI <- raster(file.path(envrmt$path_aerial_processed, "AOI_index_VVI.tif"))

veg_ind <- stack(CI, GLI, IO, NGRDI, RGR, SCI, TGI, VARI, VVI)
rm(CI, GLI, IO, NGRDI, RGR, SCI, TGI, VARI, VVI)

#read PCA
PCA <- raster::brick(file.path(envrmt$path_aerial_processed,"AOI_PCA.grd"))

PCA_mean_3x3 <- raster::brick(file.path(envrmt$path_aerial_processed, "AOI_PCA_filter _3x3_mean.tif"))
PCA_mean_5x5 <- raster::brick(file.path(envrmt$path_aerial_processed, "AOI_PCA_filter _5x5_mean.tif"))

PCA_sd_3x3 <- raster::brick(file.path(envrmt$path_aerial_processed,"AOI_PCA_filter _3x3_sd.tif"))
PCA_sd_5x5 <- raster::brick(file.path(envrmt$path_aerial_processed, "AOI_PCA_filter _5x5_sd.tif" ))


PCA_ind <- stack(PCA[[1:2]], PCA_mean_3x3[[1:2]], PCA_mean_5x5[[1:2]], PCA_sd_3x3[[1:2]], PCA_sd_5x5[[1:2]])
rm(PCA, PCA_mean_3x3, PCA_mean_5x5, PCA_sd_3x3, PCA_sd_5x5)

#stack indices
indices <- stack(rgb, veg_ind, PCA_ind)

#ind_resampl <- raster::resample(indices, shn)
#indices <- stack(ind_resampl, shn)


# [1] "AOIimage.1"                 "AOIimage.2"                 "AOIimage.3"                 "AOI_index_CI"              
#[5] "AOI_index_GLI"              "AOI_index_IO"               "AOI_index_NGRDI"            "AOI_index_red_green_ratio" 
#[9] "AOI_index_SCI"              "AOI_index_TGI"              "AOI_index_VARI"             "AOI_index_VVI"             
#[13] "PC1"                        "PC2"                        "AOI_PCA_filter__3x3_mean.1" "AOI_PCA_filter__3x3_mean.2"
#[17] "AOI_PCA_filter__5x5_mean.1" "AOI_PCA_filter__5x5_mean.2" "AOI_PCA_filter__3x3_sd.1"   "AOI_PCA_filter__3x3_sd.2"  
#[21] "AOI_PCA_filter__5x5_sd.1"   "AOI_PCA_filter__5x5_sd.2"   "shannon_level1"             "shannon_level2"            
#[25] "shannon_level3"             "shannon_level4"             "shannon_level5" 

#writeRaster(indices, file.path(envrmt$path_indices, "indices_all.tif"))
#writeRaster(indices, file.path(envrmt$path_indices, "indices_all.tif"), overwrite = T)
names(indices)[1:22] <- c("red", "green", "blue","CI","GLI", "IO", "NGRDI", "RGratio", "SCI", "TGI", "VARI", "VVI", "PC1", "PC2",
                          "PC1_mean3x3", "PC2_mean3x3", "PC1_mean5x5", "PC2_mean5x5", "PC1_std3x3", "PC2_std3x3", "PC1_std5x5", "PC2_std5x5")
#FLC
#indices <- raster::stack(file.path(envrmt$path_indices, "indices_all.tif"))

#extract pixel values of training samples
#IDs of polys as raster
indices <- indices[[1:22]]

ext = extent(indices)
all_spec_ras = gdalUtils::gdal_rasterize(src_datasource = file.path(envrmt$path_species_det, "species_poly.shp"), dst_filename = file.path(envrmt$path_species_det, "all_spec.tif"),
                                         l = "species_poly", tr = c(xres(indices),yres(indices)), te = c(ext[1],ext[3],ext[2],ext[4]),a = "id", output_Raster = T)

vals = values(all_spec_ras)
index0 = which(vals == 0)

idPoly = vals[-index0]
df = data.frame(ID = idPoly)



for(i in 1:length(indices@layers)){
  valsInd = values(indices[[i]])
  valsInd = valsInd[-index0]
  name <- paste0("valsInd_", i)
  assign(name, valsInd)
}
rm(valsInd)

for(j in 1:22){
  name <- get(paste0("valsInd_", j))
  df = cbind(df, name)
}

names(df)[2:23] <- c("red", "green", "blue","CI","GLI", "IO", "NGRDI", "RGratio", "SCI", "TGI", "VARI", "VVI", "PC1", "PC2",
                     "PC1_mean3x3", "PC2_mean3x3", "PC1_mean5x5", "PC2_mean5x5", "PC1_std3x3", "PC2_std3x3", "PC1_std5x5", "PC2_std5x5")

#assign tree species
species <- data.frame(species = NA)
df <- cbind(df[1], species, df[2:length(df)])

for(i in seq(nrow(df))){
  if(substr(df$ID[i], 1,1) == 1){
    df$species[i] <- "Douglasie"
  }else if(substr(df$ID[i], 1,1) == 2){
    df$species[i] <- "Buche"
  }else if(substr(df$ID[i], 1,1) == 3){
    df$species[i] <- "Eiche"
  }else if(substr(df$ID[i], 1,1) == 4){
    df$species[i] <- "Fichte"
  }
}

#assign forest section to data
section <- data.frame(section = NA)
df <- cbind(df[1:2], section, df[3:length(df)])

for(j in seq(length(all_spec))){
  df[df$ID == all_spec$id[j], 3] <- as.numeric(as.character(all_spec$Abteilung[j])) 
  #df[df$ID == df[j,1],2] <- as.numeric(as.character(df$section[j]))
}

df$species <- as.factor(df$species)

saveRDS(df, file.path(envrmt$path_indices, "df_400trees.rds"))
df <- readRDS("F:/BEN/edu/data/indices/df_400trees.rds")
#flc
#df<- read.csv(file.path(envrmt$path_indices, "train_area.csv"))
#-----------------------
#cross-validation
#ffs
#generate samples
smp <- lapply(seq(length(unique(df$ID))), function(i){
  set.seed(1589)
  smp_rows <- sample(nrow(df[df$ID == unique(df$ID)[i],]), 20)
  tmp <- df[df$ID == unique(df$ID)[i],][smp_rows,]
  #tmp1 <- train_area[train_area$ID == unique(train_area$ID)[i],][smp_rows[2],]
  #tmp <- rbind(tmp0, tmp1)
  return(tmp)
})

smp <- do.call(rbind, smp)

set.seed(1341324)
folds <- CAST::CreateSpacetimeFolds(smp, spacevar = "section", k = 5)

train_c <- caret::trainControl(method = "cv", number = 5, classProbs = TRUE, index = folds$index, indexOut = folds$indexOut)

cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

rfmod <- CAST::ffs(smp[4:25],smp$species, 
                   metric = "Kappa", method = "rf", withinSE = FALSE, importance = TRUE, trControl = train_c, verbose = TRUE)
stopCluster(cl)

saveRDS(rfmod, "J:/BEN/edu/data/indices/random_forest_mod3.rds")
#rfmod <- readRDS(file.path(envrmt$path_indices, "random_forest_mod3.rds"))

###############################################################
df_rownames <- c()
for(i in seq(1:nrow(df))){
  rnms<- row.names(df)[i]
  a <- row.names(smp)[which(row.names(smp) == rnms)]
  if(length(a > 0)){
    df_rownames <- c(df_rownames, row.names(df)[i])
  }
}

df_rownames <- as.numeric(df_rownames)
  
df[-df_rownames,]
###############################################################



#externe prediction, ergebnis in confusion matrix
test <- predict(rfmod, df[-df_rownames,])

confm <- caret::confusionMatrix(test, df$species[-df_rownames])
#saveRDS(confm, file.path(envrmt$path_indices, "cofusion_matrix.rds"))

pred_stack <- raster::predict(indices, rfmod)
writeRaster(pred_stack, file.path(envrmt$path_species_det , "speciesPrediction.tif"))

#variable importance
caret::varImp(rfmod)
#########################################
#ffs zur Vorhersage, welche Indizes am besten zur Unterscheidung der 
#verschiedener Baumarten geeignet sind
#
#random forest modell zur Vorhersage. welche Baumart an der Stelle der verschiedenen
#Baumkronensegmente vorliegt
#########################################
# foliage height density
# Originally MacArthur & MacArthur (1961)
# Implemented after:
# Hashimoto, H., Imanishi, J., Hagiwara, A., Morimoto, Y., & Kitada, K. (2004). 
#Estimating forest structure indices for evaluation of forest bird habitats by 
#an airborne laser scanner. In M. Thies, B. Koch, H. Spiecker, & H. Weinacker 
#(Eds.), Laser scanners for forest and landscape assessment: Proceedings of the 
#ISPRS Working Group VIII/2, Freiburg, 3-6 October 2004, 254-258.

# http://www.isprs.org/proceedings/XXXVI/8-W2/HASHIMOTO.pdf

fun_fhd <- function(a = XXX) {#input = rasterstack of lidar data
  l <- raster::nlayers(a)#get the number of layers of a raster file
  a[a<=0]=1
  p_i <- a/a[[l]]#pi = horizontal vegetation coverage in the ith layer
  r <- p_i * log(a / p_i)# = fhd for layer a
  sum(r[[1:(l-1)]])# = fhd for the tree/forest stand
}
