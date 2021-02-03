# Set working directory and load required packages
library(envimaR)
root_folder = alternativeEnvi(root_folder = "~/edu/mpg-envinsys-plygrnd", 
                              alt_env_id = "COMPUTERNAME",
                              alt_env_value = "PCRZP", 
                              alt_env_root_folder = "F:\\BEN\\edu")
source(path.expand("~/edu/mpg-envinsys-plygrnd/mpg-envinfosys-teams-2018-bjcm_rs_18/src/000_setup.R"))

# Libs for processings
library(ggplot2)
library(factoextra)
library(RStoolbox)
# Load PCA
list.files(envrmt$path_aerial_processed)
pca <- readRSTBX(file.path(envrmt$path_aerial_processed, "AOI_PCA.grd"))
pca
plot(pca$map)

# Plots
eigenvalue_pca <- factoextra::get_eigenvalue(pca$model) # eigenvalue of PCA

## variances distributed over the components
barplot <- factoextra::fviz_eig(pca$model, choice = c("variance"), ylim = c(0,80), geom = c("bar", "line"), barfill = "steelblue",
                                linecolor = "#000000", hjust = 0, addlabels = T, main = "Variances", 
                                ggtheme = ggplot2::theme_classic()) +
          ggplot2::theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 15), 
                         axis.title.x = element_text(size = 20, margin(t = 0, r = 0, l = 0, b = 0)), 
                         axis.title.y = element_text(size = 20, margin(t = 0, r = 0, l = 0, b = 0)),
                         title = element_text(size = 35))
  

plot(barplot)


## variable contribution

pca_var <- factoextra::get_pca_var(pca$model)
plot(pca_var$contrib)

pca_res <- get_pca(pca$model)
