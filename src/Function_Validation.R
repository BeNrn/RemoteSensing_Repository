val_ft <- function(vp,x, corr) {
  
  stat <- ForestTools::sp_summarise(vp, x) # compute data points in polygons
  stat[is.na(stat$TreeCount)] <- 0 # na to 0
  
  #treecount: points within polygon
  polkp <- (sum(stat$TreeCount<"1") + corr) # amount polygon without tree (miss)
  polp <- (sum(stat$TreeCount=="1") - corr) # amount polygon with exact 1 tree (hit)
  polmp <-(sum(stat$TreeCount>"1")) # amount polygon with more than 1 tree (miss)
  pkpol <- length(vp) - (sum(stat$TreeCount))# amount validation points w/o polygons 
  
  
  hit.ratio = round(polp/length(stat$TreeCount), 4) # calc hit ratio in percent (amount of exact trees) ##von allen innerhalb eines 
  #Polygons liegenden Punkte Bäumen wurden x% richtig zugeordnet
  empty.ratio = round(polkp/length(x),4) #calc empty ration in percent (amount of polygon without trees) ##von allen 
  miss.ratio = round(polmp/length(x),4) #  miss rate in percent (amount of polygons with more than 1 Tree)
  
  #print ratio and missed points
  print(paste("hit ratio: ", hit.ratio))
  print(paste("empty ratio: ",empty.ratio))
  print(paste("miss ratio: ",miss.ratio))
  print(paste("Validation point without polygon: ",pkpol))
}
