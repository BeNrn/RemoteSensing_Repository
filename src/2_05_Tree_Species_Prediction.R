library(envimaR)
root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F:\\BEN\\edu")
source(file.path(root_folder, "/mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))
#-------------------------------------------------------------------------------
spec_pred <- raster::raster(file.path(envrmt$path_species_det , "speciesPrediction.tif"))

#ids for the crown polygons
#crowns <- readOGR(file.path(envrmt$path_lidar_processed, "aoi_crowns.shp"))
#crowns$id <- 1:length(crowns)
#writeOGR(crowns, file.path(envrmt$path_species_det, "aoi_crowns_id.shp"), driver = "ESRI Shapefile", layer = "aoi_crowns_id.shp")


crowns <- readOGR(file.path(envrmt$path_species_det, "aoi_crowns_id.shp"))

#rasterite crown polys
ext = extent(spec_pred)
crown_rast = gdalUtils::gdal_rasterize(src_datasource = file.path(envrmt$path_species_det, "aoi_crowns_id.shp"), 
                                       dst_filename = file.path(envrmt$path_species_det, "crowns.tif"),
                                         l = "aoi_crowns_id", tr = c(xres(spec_pred),yres(spec_pred)),
                                       te = c(ext[1],ext[3],ext[2],ext[4]),a = "id", output_Raster = T)

saveRDS(crown_rast, file.path(envrmt$path_species_det, "crowns_rasterized.tif"))

#identify pixels with no tree polygon
vals = values(crown_rast)
index0 = which(vals == 0)

#assign id to polys
idPoly = vals[-index0]
df = data.frame(ID = idPoly)
saveRDS(df, file.path(envrmt$path_species_det, "basic_df_spec_det.rds"))

#assign pixel values to polgon ids
valsPred <-  values(spec_pred)
valsPred <- valsPred[-index0]

df <- cbind(df, valsPred)
saveRDS(df, file.path(envrmt$path_species_det, "basic_df_spec_det.rds"))


#detect probability for a tree species
df$valsPred <- as.factor(df$valsPred)

spec_prob <- table(df)

df_prob <- data.frame(id = rownames(spec_prob), bu = spec_prob[,1], dg = spec_prob[,2], ei = spec_prob[,3], fi = spec_prob[,4])
df_prob$sum <- raster::rowSums(df_prob[,2:5])

#saveRDS(df, file.path(envrmt$path_species_det, "final_df_spec_det.rds"))

problty <- unlist(lapply(seq(nrow(df_prob)), function(i){
  return(max(df_prob[i,2:5])/df_prob$sum[i])
}))

df_prob$prblty <- problty

species <- c("bu", "dg", "ei", "fi")
treeSpecies <- lapply(seq(nrow(df_prob)), function(j){
  return(species[which(df_prob[j,2:5] == max(df_prob[j,2:5]))])
})

l <- unlist(lapply(treeSpecies, length))
treeSpecies[l != 1] <- "undefined"

treeSpecies <- do.call(rbind, treeSpecies)
df_prob$species <- treeSpecies[,1]
saveRDS(df_prob, file.path(envrmt$path_species_det, "final_df_spec_det.rds"))

crowns$Species <- NA

#assign species to polygons
for(s in seq(crowns$id)){
  crowns$Species[s] <- df_prob[df_prob$id == crowns$id[s],8]
  print(s)
}

writeOGR(crowns, file.path(envrmt$path_species_det, "species_map.shp"), driver = "ESRI Shapefile", layer = "species_map.shp")
