# RemoteSensing_Repository
 This repository contains code from the Remote Sensing Course I participated. The aim of the course was to develop a tree delineation algorithm using only cheap and easy to record  remote sensing data.

## 1 Data acquisition
 Building on RGB-recordings and LiDAR data only considerably high amount of information can be extracted from aerial forest scenes without using multispectral information. Derivable Digital Elevation Models (DEM) and Digital Terrain Models (DTM) allow calculating a Canopy Height Model (CHM). Additionally, using the RGB-channels a variety of vegetation indizes can be calculated like the *Visible Vegetation Index* or the *Green Leaf Index*. Furthermore, spatial filters and the primary component analysis provide deeper insights. 

## 2 Tree delineation algorithm
 The study covers three experimental tree delineation algorithms and present their results:
  - itcSegment
  - ForestTools
  - lidR
 After comparing the three delineation algorithms, the trees in the study area were segmented with itcSegment. 
## 3 Forest structure analysis
 The forest structure analysis incorporates not only the spatial information in direction of latitude and longitude but also in the direction of the vertical axis. This allows the determination of the tree height. Based on the three dimensional information, shannon entropy and descriptive statistics were calculated. The results allow the estimation of forest structure information. 

## 4 Tree species forcast
The previous analysis steps combined can be used to develop a tree species forcasting algorithm. Utilizing an unsupervised random forest classification based on ground truth data return tree species estimations. To providemore robust results, the classification process is connected with a cross-validation.

![Image tree species map](https://github.com/BeNrn/RemoteSensing_Repository/tree/main/images/Species_Map.jpg)

## Sources
 - Dalponte, M. (2018): Individual Tree Crowns Segmentation. https://cran.r-project.org/web/packages/itcSegment/index.html. [Access: 10.12.2018]
 - Meyer, H.; Reudenbach, C.; Ludwig, M.; Nauss, T.; Pebesma, E. (2021): CAST. caret Applications for Spatial-Temporal Models. https://cran.r-project.org/web/packages/CAST/index.html [Access: 19.03.2021]
 - Plowright, A. (2018): Analyzing Remotely Sensed Forest Data. https://cran.r-project.org/web/packages/ForestTools/index.html. [Access: 10.12.2018]
 - Roussel, J.R.; Auty , D.; De Boissieu , F. & A. S. Meador (2018): Airborne LiDAR Data Manipulation and Visualization for Forestry Applications . https://cran.r-project.org/web/packages/lidR/index.html. [Access: 10.12.2018]
