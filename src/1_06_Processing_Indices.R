#Set working directory and load required packages
library(envimaR)
root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F:\\BEN\\edu")
source(path.expand("~/edu/mpg-envinsys-plygrnd/mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))

# Load in AOI
list.files(envrmt$path_aerial_processed)
img_org <- brick(file.path(envrmt$path_aerial_processed, "Testgebiet_Mischwald.tif"))
savename <- "/Testgebiet_Mischwald_index" #IMPORTANT: set "/" before the real savename
savefile <- envrmt$path_tmp

# Split into the 3 different visible colour-layers
lyr <- img_org[[1]]
lyg <- img_org[[2]]
lyb <- img_org[[3]]
## Plot them for visualisation
#plot(lyr)
#plot(lyg)
#plot(lyb)
## Remove img_org for RAM issues
rm(img_org)


# Calculate different indices

## 01 Green Leaf Index
### GLI = [(G-R)+(G-B)]/[(2*G)+R+B]
GLI <- ((lyg-lyr)+(lyg-lyb))/((2*lyg)+lyr+lyb)
#plot(GLI)
#### Save as raster
writeRaster(GLI, file = paste0(savefile, savename, "_GLI.tif"))
#### Remove GLI for RAM issues
rm(GLI)

## 02 Red - Green - Ratio
### RGR = R / G
RGR <- lyr/lyg
#plot(RGR)
### Save as raster
writeRaster(RGR, file = paste0(savefile, savename, "_red_green_ratio.tif"))
### Remove RGR_mat for RAM issues
rm(RGR)

## 03 Red - Blue - Ratio (for Iron Oxides = IO)
### IO = R/B
IO <- lyr/lyb
#plot(IO)
### Save as raster
writeRaster(IO, file = paste0(savefile, savename, "_IO.tif"))
### Remove IO for RAM issues
rm(IO)

## 04 Coloration Index
### CI = (R-B)/R
CI <- (lyr-lyb)/lyr
#plot(CI)
### Save as Raster
writeRaster(CI, file = paste0(savefile, savename, "_CI.tif"))
### Remove CI for RAM Issues
rm(CI)

## 05 Normalized Green Red Difference
### NGRDI = (G-R)/(G+R)
NGRDI <- (lyg-lyr)/(lyg+lyr)
#plot(NGRDI)
### Save as Raster
writeRaster(NGRDI, file = paste0(savefile, savename, "_NGRDI.tif"))
### Remove NGRDI for RAM issues
rm(NGRDI)

## 06 Visible Vegetation Index
### VVI = abs((1-abs((R-R0)/(R+R0)))*(1-abs((G-G0)/(G+G0)))*(1-abs((B-B0)/B+B0)))
###R0 = 30; G0 = 50; B0 = 0 <-- by default -> so using a bit more -> R0, G0, B0 + 0,01
VVI <- ( (1 - abs( (lyr - 30.01) / (lyr + 30.01) ) ) *
         (1 - abs( (lyg - 50.01) / (lyg + 50.01) ) ) *
         (1 - abs( (lyb -  0.01) / (lyb +  0.01)) ) )
#plot(VVI)
### Save as Raster
writeRaster(VVI, file = paste0(savefile, savename, "_VVI.tif"))
### Remove VVI for RAM issues
rm(VVI)

## 07 Soil Colour Index
### SCI = (R-G)/(R+G)
SCI <- (lyr-lyg)/(lyr+lyg)
#plot(SCI)
### Save as Raster
writeRaster(SCI, file = paste0(savefile, savename, "_SCI.tif"))
### Remove SCi for RAM issues
rm(SCI)

## 08 Triangular Greeness Index
### TGI = G-0.39*R-0.61*B
TGI <- lyg-0.39*lyr-0.61*lyb
#plot(TGI)
### Save as Raster
writeRaster(TGI, file = paste0(savefile, savename, "_TGI.tif"))
### Remove TGI for RAM issues
rm(TGI)

## 09 Visible Atmospherically Resistant Index
### VARI = (G-R)/(G+R-B)
VARI <- (lyg-lyr)/(lyg+lyr-lyb)

#plot(VARI)
### Save as Raster
writeRaster(VARI, file = paste0(savefile, savename, "_VARI.tif"))
### Remove VARI for RAM issues
rm(VARI)
