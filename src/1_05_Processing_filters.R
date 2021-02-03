# Set working directory and load required packages
library(envimaR)
root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F:\\BEN\\edu")
source(path.expand("~/edu/mpg-envinsys-plygrnd/mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))



# Load in original image to bee processed
list.files(envrmt$path_aerial_processed)
img_org <- stack(file.path(envrmt$path_aerial_processed, "AOI_PCA.grd"))

# Saving name
## The name under which the new data is saved
## savename should be close to the original inputdata + "_filter"
## Exmpl.: envrmt$path_processed, savename, "_3x3_mean.tif"
## meaning with savename = "VVI_index_filter" -> saved under "VVI_index_filter_3x3_mean.tif"
savename <- "AOI_PCA_filter"

# saving file
savefolder <- envrmt$path_aerial_processed

# The Matrix --------------------------------------------------------------------------------
trinity <- matrix(1/9, ncol = 3, nrow = 3)
morpheus <- matrix(1/25 , ncol = 5, nrow = 5)
mrsmith <- matrix(1/49, ncol = 7, nrow = 7)
whiterabbit <- matrix(1/225 , ncol = 15, nrow = 15)
neo <- matrix(1/961, ncol = 31, nrow = 31)


#Mean filters ---------------------------------------------------------------------------------------
##for loop for mean filter 3x3
mean3x3 <- brick()
mean3x3 = img_org
for (i in seq(1:nlayers(img_org))){
  mean3x3[[i]] <- focal(x = img_org[[i]], fun = mean, w = trinity)
}
###plot the filter
if (class(mean3x3) == "RasterLayer"){
  plot(mean3x3, main = paste0(savename, "_3x3_mean"))
} else if (class(mean3x3) == "RasterBrick"){
  if (nbands(mean3x3) < 3){
    plot(mean3x3, main = paste0(savename, "_3x3_mean"))
  } else if (nbands(mean3x3) == 3){
    plotRGB(mean3x3, main = paste0(savename, "_3x3_mean"))
  } else if (nbands(mean3x3) > 3){
    plot(mean3x3, main = paste0(savename, "_3x3_mean"))
  }
} else if (class(mean3x3) == "RasterStack"){
  plot(mean3x3, main = paste0(savename, "_3x3_mean"))
}
###save the filter
writeRaster(mean3x3, filename = file.path(savefolder, paste(savename, "_3x3_mean.tif")))
###remove the filter
rm(mean3x3)

##for loop for mean filter 5x5
mean5x5 <- brick()
mean5x5 = img_org
for (i in seq(1:nlayers(img_org))){
  mean5x5[[i]] <- focal(x = img_org[[i]], fun = mean, w = morpheus)
}
###plot the filter
if (class(mean5x5) == "RasterLayer"){
  plot(mean5x5, main = paste0(savename, "_5x5_mean"))
} else if (class(mean5x5) == "RasterBrick"){
  if (nbands(mean5x5) < 3){
    plot(mean5x5, main = paste0(savename, "_5x5_mean"))
  } else if (nbands(mean5x5) == 3){
    plotRGB(mean5x5, main = paste0(savename, "_5x5_mean"))
  } else if (nbands(mean5x5) > 3){
    plot(mean5x5, main = paste0(savename, "_5x5_mean"))
  }
} else if (class(mean5x5) == "RasterStack"){
  plot(mean5x5, main = paste0(savename, "_5x5_mean"))
}
###save the filter
writeRaster(mean5x5, filename = file.path(savefolder, paste(savename, "_5x5_mean.tif")))
###remove the filter
rm(mean5x5)

###for loop for mean filter 7x7
mean7x7 <- brick()
mean7x7 = img_org
for (i in seq(1:nlayers(img_org))){
  mean7x7[[i]] <- focal(x = img_org[[i]], fun = mean, w = mrsmith)
}
###plot the filter
if (class(mean7x7) == "RasterLayer"){
  plot(mean7x7, main = paste0(savename, "_7x7_mean"))
} else if (class(mean7x7) == "RasterBrick"){
  if (nbands(mean7x7) < 3){
    plot(mean7x7, main = paste0(savename, "_7x7_mean"))
  } else if (nbands(mean7x7) == 3){
    plotRGB(mean7x7, main = paste0(savename, "_7x7_mean"))
  } else if (nbands(mean7x7) > 3){
    plot(mean7x7, main = paste0(savename, "_7x7_mean"))
  }
} else if (class(mean7x7) == "RasterStack"){
  plot(mean7x7, main = paste0(savename, "_7x7_mean"))
}
###save the filter
writeRaster(mean7x7, filename = file.path(savefolder, paste(savename, "_7x7_mean.tif")))
###remove the filter
rm(mean7x7)

##for loop for mean filter 15x15
mean15x15 <- brick()
mean15x15 = img_org
for (i in seq(1:nlayers(img_org))){
  mean15x15[[i]] <- focal(x = img_org[[i]], fun = mean, w = whiterabbit)
}
###plot the filter
if (class(mean15x15) == "RasterLayer"){
  plot(mean15x15, main = paste0(savename, "_15x15_mean"))
} else if (class(mean15x15) == "RasterBrick"){
  if (nbands(mean15x15) < 3){
    plot(mean15x15, main = paste0(savename, "_15x15_mean"))
  } else if (nbands(mean15x15) == 3){
    plotRGB(mean15x15, main = paste0(savename, "_15x15_mean"))
  } else if (nbands(mean15x15) > 3){
    plot(mean15x15, main = paste0(savename, "_15x15_mean"))
  }
} else if (class(mean15x15) == "RasterStack"){
  plot(mean15x15, main = paste0(savename, "_15x15_mean"))
}
###save the filter
writeRaster(mean15x15, filename = file.path(savefolder, paste(savename, "_15x15_mean.tif")))
###remove the filter
rm(mean15x15)

##for loop for mean filter 31x31
mean31x31 <- brick()
mean31x31 = img_org
for (i in seq(1:nlayers(img_org))){
  mean31x31[[i]] <- focal(x = img_org[[i]], fun = mean, w = neo)
}
###plot the filter
if (class(mean31x31) == "RasterLayer"){
  plot(mean31x31, main = paste0(savename, "_31x31_mean"))
} else if (class(mean31x31) == "RasterBrick"){
  if (nbands(mean31x31) < 3){
    plot(mean31x31, main = paste0(savename, "_31x31_mean"))
  } else if (nbands(mean31x31) == 3){
    plotRGB(mean31x31, main = paste0(savename, "_31x31_mean"))
  } else if (nbands(mean31x31) > 3){
    plot(mean31x31, main = paste0(savename, "_31x31_mean"))
  }
} else if (class(mean31x31) == "RasterStack"){
  plot(mean31x31, main = paste0(savename, "_31x31_mean"))
}
###save the filter
writeRaster(mean31x31, filename = file.path(savefolder, paste(savename, "_31x31_mean.tif")))
###remove the filter
rm(mean31x31)

# Grey level cooccurence matrices filters ---------------------------------------------------------------
library(glcm)
n_grey <- 2 # number of grey levels to use in texture calculation


##for loop for glcm filter 3x3
glcm3x3 <- brick()
glcm3x3 = img_org
for (i in seq(1:nlayers(img_org))){
  glcm3x3[[i]] <- glcm::glcm(x = img_org[[i]],n_grey = n_grey, window = c(3,3))
}
###plot the filter
if (class(glcm3x3) == "RasterLayer"){
  plot(glcm3x3, main = paste0(savename, "_3x3_glcm"))
} else if (class(glcm3x3) == "RasterBrick"){
  if (nbands(glcm3x3) < 3){
    plot(glcm3x3, main = paste0(savename, "_3x3_glcm"))
  } else if (nbands(glcm3x3) == 3){
    plotRGB(glcm3x3, main = paste0(savename, "_3x3_glcm"))
  } else if (nbands(glcm3x3) > 3){
    plot(glcm3x3, main = paste0(savename, "_3x3_glcm"))
  }
} else if (class(glcm3x3) == "RasterStack"){
  plot(glcm3x3, main = paste0(savename, "_3x3_glcm"))
}
###save the filter
writeRaster(glcm3x3, filename = file.path(savefolder, paste(savename, "_3x3_glcm.tif")))
###remove the filter
rm(glcm3x3)

##for loop for glcm filter 5x5
glcm5x5 <- brick()
glcm5x5 = img_org
for (i in seq(1:nlayers(img_org))){
  glcm5x5[[i]] <- glcm::glcm(x = img_org[[i]],n_grey = n_grey, window = c(5, 5))
}
###plot the filter
if (class(glcm5x5) == "RasterLayer"){
  plot(glcm5x5, main = paste0(savename, "_5x5_glcm"))
} else if (class(glcm5x5) == "RasterBrick"){
  if (nbands(glcm5x5) < 3){
    plot(glcm5x5, main = paste0(savename, "_5x5_glcm"))
  } else if (nbands(glcm5x5) == 3){
    plotRGB(glcm5x5, main = paste0(savename, "_5x5_glcm"))
  } else if (nbands(glcm5x5) > 3){
    plot(glcm5x5, main = paste0(savename, "_5x5_glcm"))
  }
} else if (class(glcm5x5) == "RasterStack"){
  plot(glcm5x5, main = paste0(savename, "_5x5_glcm"))
}
###save the filter
writeRaster(glcm5x5, filename = file.path(savefolder, paste(savename, "_5x5_glcm.tif")))
###remove the filter
rm(glcm5x5)

###for loop for glcm filter 7x7
glcm7x7 <- brick()
glcm7x7 = img_org
for (i in seq(1:nlayers(img_org))){
  glcm7x7[[i]] <- glcm::glcm(x = img_org[[i]],n_grey = n_grey, window = c(7, 7))
}
###plot the filter
if (class(glcm7x7) == "RasterLayer"){
  plot(glcm7x7, main = paste0(savename, "_7x7_glcm"))
} else if (class(glcm7x7) == "RasterBrick"){
  if (nbands(glcm7x7) < 3){
    plot(glcm7x7, main = paste0(savename, "_7x7_glcm"))
  } else if (nbands(glcm7x7) == 3){
    plotRGB(glcm7x7, main = paste0(savename, "_7x7_glcm"))
  } else if (nbands(glcm7x7) > 3){
    plot(glcm7x7, main = paste0(savename, "_7x7_glcm"))
  }
} else if (class(glcm7x7) == "RasterStack"){
  plot(glcm7x7, main = paste0(savename, "_7x7_glcm"))
}
###save the filter
writeRaster(glcm7x7, filename = file.path(savefolder, paste(savename, "_7x7_glcm.tif")))
###remove the filter
rm(glcm7x7)

##for loop for glcm filter 15x15
glcm15x15 <- brick()
glcm15x15 = img_org
for (i in seq(1:nlayers(img_org))){
  glcm15x15[[i]] <- glcm::glcm(x = img_org[[i]],n_grey = n_grey, window = c(15, 15))
}
###plot the filter
if (class(glcm15x15) == "RasterLayer"){
  plot(glcm15x15, main = paste0(savename, "_15x15_glcm"))
} else if (class(glcm15x15) == "RasterBrick"){
  if (nbands(glcm15x15) < 3){
    plot(glcm15x15, main = paste0(savename, "_15x15_glcm"))
  } else if (nbands(glcm15x15) == 3){
    plotRGB(glcm15x15, main = paste0(savename, "_15x15_glcm"))
  } else if (nbands(glcm15x15) > 3){
    plot(glcm15x15, main = paste0(savename, "_15x15_glcm"))
  }
} else if (class(glcm15x15) == "RasterStack"){
  plot(glcm15x15, main = paste0(savename, "_15x15_glcm"))
}
###save the filter
writeRaster(glcm15x15, filename = file.path(savefolder, paste(savename, "_15x15_glcm.tif")))
###remove the filter
rm(glcm15x15)

##for loop for glcm filter 31x31
glcm31x31 <- brick()
glcm31x31 = img_org
for (i in seq(1:nlayers(img_org))){
  glcm31x31[[i]] <- glcm::glcm(x = img_org[[i]],n_grey = n_grey, window = c(31, 31))
}
###plot the filter
if (class(glcm31x31) == "RasterLayer"){
  plot(glcm31x31, main = paste0(savename, "_31x31_glcm"))
} else if (class(glcm31x31) == "RasterBrick"){
  if (nbands(glcm31x31) < 3){
    plot(glcm31x31, main = paste0(savename, "_31x31_glcm"))
  } else if (nbands(glcm31x31) == 3){
    plotRGB(glcm31x31, main = paste0(savename, "_31x31_glcm"))
  } else if (nbands(glcm31x31) > 3){
    plot(glcm31x31, main = paste0(savename, "_31x31_glcm"))
  }
} else if (class(glcm31x31) == "RasterStack"){
  plot(glcm31x31, main = paste0(savename, "_31x31_glcm"))
}
###save the filter
writeRaster(glcm31x31, filename = file.path(savefolder, paste(savename, "_31x31_glcm.tif")))
###remove the filter
rm(glcm31x31)

##remove n_grey
rm(n_grey)

# Standard deviation filter ---------------------------------------------------------------------------------------
##for loop for sd filter 3x3
sd3x3 <- brick()
sd3x3 = img_org
for (i in seq(1:nlayers(img_org))){
  sd3x3[[i]] <- focal(x = img_org[[i]], fun = sd, w = trinity)
}
###plot the filter
if (class(sd3x3) == "RasterLayer"){
  plot(sd3x3, main = paste0(savename, "_3x3_sd"))
} else if (class(sd3x3) == "RasterBrick"){
  if (nbands(sd3x3) < 3){
    plot(sd3x3, main = paste0(savename, "_3x3_sd"))
  } else if (nbands(sd3x3) == 3){
    plotRGB(sd3x3, main = paste0(savename, "_3x3_sd"))
  } else if (nbands(sd3x3) > 3){
    plot(sd3x3, main = paste0(savename, "_3x3_sd"))
  }
} else if (class(sd3x3) == "RasterStack"){
  plot(sd3x3, main = paste0(savename, "_3x3_sd"))
}
###save the filter
writeRaster(sd3x3, filename = file.path(savefolder, paste(savename, "_3x3_sd.tif")))
###remove the filter
rm(sd3x3)

##for loop for sd filter 5x5
sd5x5 <- brick()
sd5x5 = img_org
for (i in seq(1:nlayers(img_org))){
  sd5x5[[i]] <- focal(x = img_org[[i]], fun = sd, w = morpheus)
}
###plot the filter
if (class(sd5x5) == "RasterLayer"){
  plot(sd5x5, main = paste0(savename, "_5x5_sd"))
} else if (class(sd5x5) == "RasterBrick"){
  if (nbands(sd5x5) < 3){
    plot(sd5x5, main = paste0(savename, "_5x5_sd"))
  } else if (nbands(sd5x5) == 3){
    plotRGB(sd5x5, main = paste0(savename, "_5x5_sd"))
  } else if (nbands(sd5x5) > 3){
    plot(sd5x5, main = paste0(savename, "_5x5_sd"))
  }
} else if (class(sd5x5) == "RasterStack"){
  plot(sd5x5, main = paste0(savename, "_5x5_sd"))
}
###save the filter
writeRaster(sd5x5, filename = file.path(savefolder, paste(savename, "_5x5_sd.tif")))
###remove the filter
rm(sd5x5)

##for loop for sd filter 7x7
sd7x7 <- brick()
sd7x7 = img_org
for (i in seq(1:nlayers(img_org))){
  sd7x7[[i]] <- focal(x = img_org[[i]], fun = sd, w = mrsmith)
}
###plot the filter
if (class(sd7x7) == "RasterLayer"){
  plot(sd7x7, main = paste0(savename, "_7x7_sd"))
} else if (class(sd7x7) == "RasterBrick"){
  if (nbands(sd7x7) < 3){
    plot(sd7x7, main = paste0(savename, "_7x7_sd"))
  } else if (nbands(sd7x7) == 3){
    plotRGB(sd7x7, main = paste0(savename, "_7x7_sd"))
  } else if (nbands(sd7x7) > 3){
    plot(sd7x7, main = paste0(savename, "_7x7_sd"))
  }
} else if (class(sd7x7) == "RasterStack"){
  plot(sd7x7, main = paste0(savename, "_7x7_sd"))
}
###save the filter
writeRaster(sd7x7, filename = file.path(savefolder, paste(savename, "_7x7_sd.tif")))
###remove the filter
rm(sd7x7)

##for loop for sd filter 15x15
sd15x15 <- brick()
sd15x15 = img_org
for (i in seq(1:nlayers(img_org))){
  sd15x15[[i]] <- focal(x = img_org[[i]], fun = sd, w = whiterabbit)
}
###plot the filter
if (class(sd15x15) == "RasterLayer"){
  plot(sd15x15, main = paste0(savename, "_15x15_sd"))
} else if (class(sd15x15) == "RasterBrick"){
  if (nbands(sd15x15) < 3){
    plot(sd15x15, main = paste0(savename, "_15x15_sd"))
  } else if (nbands(sd15x15) == 3){
    plotRGB(sd15x15, main = paste0(savename, "_15x15_sd"))
  } else if (nbands(sd15x15) > 3){
    plot(sd15x15, main = paste0(savename, "_15x15_sd"))
  }
} else if (class(sd15x15) == "RasterStack"){
  plot(sd15x15, main = paste0(savename, "_15x15_sd"))
}
###save the filter
writeRaster(sd15x15, filename = file.path(savefolder, paste(savename, "_15x15_sd.tif")))
###remove the filter
rm(sd15x15)

##for loop for sd filter 31x31
sd31x31 <- brick()
sd31x31 = img_org
for (i in seq(1:nlayers(img_org))){
  sd31x31[[i]] <- focal(x = img_org[[i]], fun = sd, w = neo)
}
###plot the filter
if (class(sd31x31) == "RasterLayer"){
  plot(sd31x31, main = paste0(savename, "_31x31_sd"))
} else if (class(sd31x31) == "RasterBrick"){
  if (nbands(sd31x31) < 3){
    plot(sd31x31, main = paste0(savename, "_31x31_sd"))
  } else if (nbands(sd31x31) == 3){
    plotRGB(sd31x31, main = paste0(savename, "_31x31_sd"))
  } else if (nbands(sd31x31) > 3){
    plot(sd31x31, main = paste0(savename, "_31x31_sd"))
  }
} else if (class(sd31x31) == "RasterStack"){
  plot(sd31x31, main = paste0(savename, "_31x31_sd"))
}
###save the filter
writeRaster(sd31x31, filename = file.path(savefolder, paste(savename, "_31x31_sd.tif")))
###remove the filter
rm(sd31x31)


# Gaussian filter ------------------------------------------------------------------------------------------------
## for loop for d = c(3,3)
gf3x3 <- brick()
gf3x3 = img_org
for (i in seq(1:nlayers(img_org))){
  gf3x3[[i]] <- focal(img_org[[i]], w = focalWeight(img_org[[i]], d = c(3,3), type = "Gauss"))
  
}
###plot the filter
if (class(gf3x3) == "RasterLayer"){
  plot(gf3x3, main = paste0(savename, "_3x3_gf"))
} else if (class(gf3x3) == "RasterBrick"){
  if (nbands(gf3x3) < 3){
    plot(gf3x3, main = paste0(savename, "_3x3_gf"))
  } else if (nbands(gf3x3) == 3){
    plotRGB(gf3x3, main = paste0(savename, "_3x3_gf"))
  } else if (nbands(gf3x3) > 3){
    plot(gf3x3, main = paste0(savename, "_3x3_gf"))
  }
} else if (class(gf3x3) == "RasterStack"){
  plot(gf3x3, main = paste0(savename, "_3x3_gf"))
}
###save the filter
writeRaster(gf3x3, filename = file.path(savefolder, paste(savename, "_3x3_gauss")))
### remove the filter
rm(gf3x3)


## for loop for d = c(5,5)
gf5x5 <- brick()
gf5x5 = img_org
for (i in seq(1:nlayers(img_org))){
  gf5x5[[i]] <- focal(img_org[[i]], w = focalWeight(img_org[[i]], d = c(5,5), type = "Gauss"))
  
}
###plot the filter
if (class(gf5x5) == "RasterLayer"){
  plot(gf5x5, main = paste0(savename, "_5x5_gf"))
} else if (class(gf5x5) == "RasterBrick"){
  if (nbands(gf5x5) < 3){
    plot(gf5x5, main = paste0(savename, "_5x5_gf"))
  } else if (nbands(gf5x5) == 3){
    plotRGB(gf5x5, main = paste0(savename, "_5x5_gf"))
  } else if (nbands(gf5x5) > 3){
    plot(gf5x5, main = paste0(savename, "_5x5_gf"))
  }
} else if (class(gf5x5) == "RasterStack"){
  plot(gf5x5, main = paste0(savename, "_5x5_gf"))
}
###save the filter
writeRaster(gf5x5, filename = file.path(savefolder, paste(savename, "_5x5_gauss")))
### remove the filter
rm(gf5x5)

## for loop for d = c(7,7)
gf7x7 <- brick()
gf7x7 = img_org
for (i in seq(1:nlayers(img_org))){
  gf7x7[[i]] <- focal(img_org[[i]], w = focalWeight(img_org[[i]], d = c(7,7), type = "Gauss"))
  
}
###plot the filter
if (class(gf7x7) == "RasterLayer"){
  plot(gf7x7, main = paste0(savename, "_7x7_gf"))
} else if (class(gf7x7) == "RasterBrick"){
  if (nbands(gf7x7) < 3){
    plot(gf7x7, main = paste0(savename, "_7x7_gf"))
  } else if (nbands(gf7x7) == 3){
    plotRGB(gf7x7, main = paste0(savename, "_7x7_gf"))
  } else if (nbands(gf7x7) > 3){
    plot(gf7x7, main = paste0(savename, "_7x7_gf"))
  }
} else if (class(gf7x7) == "RasterStack"){
  plot(gf7x7, main = paste0(savename, "_7x7_gf"))
}
###save the filter
writeRaster(gf7x7, filename = file.path(savefolder, paste(savename, "_7x7_gauss")))
### remove the filter
rm(gf7x7)
rm(trinity)
rm(morpheus)
rm(mrsmith)
rm(whiterabbit)
rm(neo)

# Edge detection ------------------------------------------------------------------------------
#library(imager)
#img_org_cimg <- as.imlist(img_org)


#ed_grad <- imlist()
#ed_x <- imlist()
#ed_y <- imlist()
#for (i in seq(1:nlayers(img_org))){
  #ed_x[[i]] <- imager::imgradient(img_org_cimg[[i]], "x")
  #ed_y[[i]] <- imager::imgradient(img_org_cimg[[i]], "y")
#}

#for (j in seq(1:length(img_org_cimg))){
  #ed_grad[[j]] <- sqrt(ed_x[[j]]^2+ed_y[[j]]^2)
#}
### plot the filter
#plot(ed_grad)


## Sobel filter ------------------------------------------------------------------------------------
### 3x3
kx = matrix(c(-1,-2,-1,0,0,0,1,2,1), ncol=3)
ky = matrix(c(1,0,-1,2,0,-2,1,0,-1), ncol=3)
k = (kx**2 + ky**2)**0.5

sf3x3 <- brick()
sf3x3 = img_org
for (i in seq(1:nlayers(img_org))){
  sf3x3[[i]] <- focal(img_org[[i]], w=k)
}

## plot the filter
if (class(sf3x3) == "RasterLayer"){
  plot(sf3x3, main = paste0(savename, "_3x3_sf"))
} else if (class(sf3x3) == "RasterBrick"){
  if (nbands(sf3x3) < 3){
    plot(sf3x3, main = paste0(savename, "_3x3_sf"))
  } else if (nbands(sf3x3) == 3){
    plotRGB(sf3x3, main = paste0(savename, "_3x3_sf"))
  } else if (nbands(sf3x3) > 3){
    plot(sf3x3, main = paste0(savename, "_3x3_sf"))
  }
} else if (class(sf3x3) == "RasterStack"){
  plot(sf3x3, main = paste0(savename, "_3x3_sf"))
}
## save the filter
writeRaster(sf3x3, filename = file.path(savefolder, paste(savename, "_3x3_sf.tif")))
## remove the filter
rm(kx)
rm(ky)
rm(k)
rm(sf3x3)

### 5x5
kx = matrix(c(-1,-1,-2,-1,-1, 0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  1,1,2,1,1), ncol=5)
ky = matrix(c(1,0,0,0,-1,  2,0,0,0,-2,  2,0,0,0,-2,  2,0,0,0,-2,  1,0,0,0,-1), ncol=5)
k = (kx**2 + ky**2)**0.5

sf5x5 <- brick()
sf5x5 = img_org
for (i in seq(1:nlayers(img_org))){
  sf5x5[[i]] <- focal(img_org[[i]], w=k)
}
## plot the filter
if (class(sf5x5) == "RasterLayer"){
  plot(sf5x5, main = paste0(savename, "_5x5_sf"))
} else if (class(sf5x5) == "RasterBrick"){
  if (nbands(sf5x5) < 3){
    plot(sf5x5, main = paste0(savename, "_5x5_sf"))
  } else if (nbands(sf5x5) == 3){
    plotRGB(sf5x5, main = paste0(savename, "_5x5_sf"))
  } else if (nbands(sf5x5) > 3){
    plot(sf5x5, main = paste0(savename, "_5x5_sf"))
  }
} else if (class(sf5x5) == "RasterStack"){
  plot(sf5x5, main = paste0(savename, "_5x5_sf"))
}
## save the filter
writeRaster(sf5x5, filename = file.path(savefolder, paste(savename, "_5x5_sf.tif")))
## remove the filter
rm(kx)
rm(ky)
rm(k)
rm(sf5x5)

### 7x7
kx = matrix(c(-1,-1,-2,-2,-2,-1,-1, -1,-1,-2,-2,-2,-1,-1,  0,0,0,0,0,0,0,  0,0,0,0,0,0,0,  
              0,0,0,0,0,0,0,        1,1,2,2,2,1,1,  1,1,2,2,2,1,1), ncol=7)

ky = matrix(c(1,1,0,0,0,-1,-1,  2,2,0,0,0,-2,-2,  2,2,0,0,0,-2,-2,  2,2,0,0,0,-2,-2,  
              2,2,0,0,0,-2,-2,  2,2,0,0,0,-2,-2,  1,1,0,0,0,-1,-1), ncol=7)
k = (kx**2 + ky**2)**0.5

sf7x7 <- brick()
sf7x7 = img_org
for (i in seq(1:nlayers(img_org))){
  sf7x7[[i]] <- focal(img_org[[i]], w = k)
}
## plot the filter
if (class(sf7x7) == "RasterLayer"){
  plot(sf7x7, main = paste0(savename, "_7x7_sf"))
} else if (class(sf7x7) == "RasterBrick"){
  if (nbands(sf7x7) < 3){
    plot(sf7x7, main = paste0(savename, "_7x7_sf"))
  } else if (nbands(sf7x7) == 3){
    plotRGB(sf7x7, main = paste0(savename, "_7x7_sf"))
  } else if (nbands(sf7x7) > 3){
    plot(sf7x7, main = paste0(savename, "_7x7_sf"))
  }
} else if (class(sf7x7) == "RasterStack"){
  plot(sf7x7, main = paste0(savename, "_7x7_sf"))
}
## save the filter
writeRaster(sf7x7, filename = file.path(savefolder, paste(savename, "_7x7_sf.tif")))
## remove the filter
rm(kx)
rm(ky)
rm(k)
rm(sf7x7)

### 15x15
kx = matrix(c(-1,-1,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,  
              -1,-1,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,
              -1,-1,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,  
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  
              1,1,2,2,2,2,2,2,2,2,2,2,2,1,1,  1,1,2,2,2,2,2,2,2,2,2,2,2,1,1,  
              1,1,2,2,2,2,2,2,2,2,2,2,2,1,1), ncol=15)

ky = matrix(c(1,1,1,0,0,0,0,0,0,0,0,0,-1,-1,-1,  2,2,2,0,0,0,0,0,0,0,0,0,-2,-2,-2,  
              2,2,2,0,0,0,0,0,0,0,0,0,-2,-2,-2,  2,2,2,0,0,0,0,0,0,0,0,0,-2,-2,-2,  
              2,2,2,0,0,0,0,0,0,0,0,0,-2,-2,-2,  2,2,2,0,0,0,0,0,0,0,0,0,-2,-2,-2,  
              2,2,2,0,0,0,0,0,0,0,0,0,-2,-2,-2,  2,2,2,0,0,0,0,0,0,0,0,0,-2,-2,-2,  
              2,2,2,0,0,0,0,0,0,0,0,0,-2,-2,-2,  2,2,2,0,0,0,0,0,0,0,0,0,-2,-2,-2,
              2,2,2,0,0,0,0,0,0,0,0,0,-2,-2,-2,  2,2,2,0,0,0,0,0,0,0,0,0,-2,-2,-2,  
              2,2,2,0,0,0,0,0,0,0,0,0,-2,-2,-2,  2,2,2,0,0,0,0,0,0,0,0,0,-2,-2,-2,  
              1,1,1,0,0,0,0,0,0,0,0,0,-1,-1,-1), 
            ncol=15)
k = (kx**2 + ky**2)**0.5

sf15x15 <- brick()
sf15x15 = img_org
for (i in seq(1:nlayers(img_org))){
  sf15x15[[i]] <- focal(img_org[[i]], w = k)
}
## plot the filter
if (class(sf15x15) == "RasterLayer"){
  plot(sf15x15, main = paste0(savename, "_15x15_sf"))
} else if (class(sf15x15) == "RasterBrick"){
  if (nbands(sf15x15) < 3){
    plot(sf15x15, main = paste0(savename, "_15x15_sf"))
  } else if (nbands(sf15x15) == 3){
    plotRGB(sf15x15, main = paste0(savename, "_15x15_sf"))
  } else if (nbands(sf15x15) > 3){
    plot(sf15x15, main = paste0(savename, "_15x15_sf"))
  }
} else if (class(sf15x15) == "RasterStack"){
  plot(sf15x15, main = paste0(savename, "_15x15_sf"))
}
## save the filter
writeRaster(sf15x15, filename = file.path(savefolder, paste(savename, "_15x15_sf.tif")))
## remove the filter
rm(kx)
rm(ky)
rm(k)
rm(sf15x15)

### 31x31
kx = matrix(c(-1,-1,-1,-1,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,-1,-1, 
              -1,-1,-1,-1,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,-1,-1,  
              -1,-1,-1,-1,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,-1,-1,  
              -1,-1,-1,-1,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,-1,-1,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,
              1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,
              1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,
              1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1),
            ncol=31)
ky = matrix(c(1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,
              1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1, 
              1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,
              1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,  
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,
              1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,
              1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1, 
              1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1, 
              1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,
              1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1),
            ncol=31)
k = (kx**2 + ky**2)**0.5

sf31x31 <- brick()
sf31x31 = img_org
for (i in seq(1:nlayers(img_org))){
  sf31x31[[i]] <- focal(img_org[[i]], w = k)
}
## plot the filter
if (class(sf31x31) == "RasterLayer"){
  plot(sf31x31, main = paste0(savename, "_31x31_sf"))
} else if (class(sf31x31) == "RasterBrick"){
  if (nbands(sf31x31) < 3){
    plot(sf31x31, main = paste0(savename, "_31x31_sf"))
  } else if (nbands(sf31x31) == 3){
    plotRGB(sf31x31, main = paste0(savename, "_31x31_sf"))
  } else if (nbands(sf31x31) > 3){
    plot(sf31x31, main = paste0(savename, "_31x31_sf"))
  }
} else if (class(sf31x31) == "RasterStack"){
  plot(sf31x31, main = paste0(savename, "_31x31_sf"))
}
## save the filter
writeRaster(sf31x31, filename = file.path(savefolder, paste(savename, "_31x31_sf.tif")))
## remove the filter
rm(kx)
rm(ky)
rm(k)
rm(sf31x31)

