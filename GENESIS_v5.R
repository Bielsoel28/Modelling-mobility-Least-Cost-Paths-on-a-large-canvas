###########################################
######### R Code: GENESIS_v5 ##############
###########################################
##### Designed by: Biel Soriano Elias #####
###########################################

# Function to compute ZPMs from a series of equally distributed limit points #
# Required packages: raster, sf, gdistance, doSNOW, parallel #
# OUTPUTS: Shapefile (.shp) with the vector of the ZPMs and a Raster (.tiff) with the ZPMs map
# INPUTS:
  # y = DEM raster with a projected CRS
  # q = Distance in kilometres between the edge points
  # z = Vector map with the rivers (Same CRS as y needed)
  # w = Cost function name, to choose: "TOBLER", "LLOBERA"

GENESIS_v5 <- function(y,q,z,w) { 

  #Loading of the input DEM or "y"
  elevacio <- y

  #Creation of the limit edge points
  
    #DEM's Limit values conversion (from package "raster")
    extent_raster<-  extent(elevacio) #reading of the values
    extent_matrix <- matrix(extent_raster) #values to matrix
    
    #Calculation of number edge points depending on "q"
      
      #Reading of the limit values of y
      dim_raster <- dim(elevacio) #reading of the limit values
      res_raster <- res(elevacio) #reading of the cell resolution 
    
        #Calculation of "y" dimensions in kilometers
        limit_sup <- dim_raster[1] * res_raster[1] / 1000 #limit Y calculation
        limit_inf <- dim_raster[2] * res_raster[1] / 1000 #limit X calculation
      
        #Division of the limits values of "y" in kilometers between "q"
        q1 <- limit_sup / q #limit Y calculation
        q2 <- limit_inf / q #limit X calculation 
      
        #Recovery of X and Y coordinates for the edge points depending on "q"
        limit_Y <- seq(from = extent_matrix[3, ], to = extent_matrix[4, ], length.out = q1) #limit Y calculation
        limit_X <- seq(from = extent_matrix[1, ], to = extent_matrix[2, ], length.out = q2) #limit X calculation
        
      #Recovery of edge points coordinates of each edge of the map to and SpatialPoint format (from package "sp")
      punts_limit_1 <- SpatialPoints(data.frame(X = limit_X[2:q2], Y = extent_raster@ymin)) #edge 1
      punts_limit_2 <- SpatialPoints(data.frame(X = limit_X[2:q2], Y = extent_raster@ymax)) #edge 2
      punts_limit_3 <- SpatialPoints(data.frame(X = extent_raster@xmin, Y = limit_Y[2:q1])) #edge 3
      punts_limit_4 <- SpatialPoints(data.frame(X = extent_raster@xmax, Y = limit_Y[2:q1])) #edge 4
      
      #Recovery of the objective points coordinates (LCPs) for every edge 
      punts_desti_1 <- rbind(punts_limit_2, punts_limit_3, punts_limit_4) #edge 1 objetive
      punts_desti_2 <- rbind(punts_limit_1, punts_limit_3, punts_limit_4) #edge 2 objetive
      punts_desti_3 <- rbind(punts_limit_1, punts_limit_2, punts_limit_4) #edge 3 objetive
      punts_desti_4 <- rbind(punts_limit_1, punts_limit_2, punts_limit_3) #edge 4 objetive
  
print("Step 1 completed")  
  
  
  #Application of the cost function (creation of the trasition layer from package "gdistance")
      
    #Slope calculation for "y"
    altDiff <- function(ph) {ph[2] - ph[1]} #Slope calculation function
    hd <- transition(elevacio, altDiff, 16, symm= FALSE) #Creation of transition layer for slope 
    slope <- geoCorrection(hd) #Geocorrection to apply to the values for slope's transition layer
    
    #Cost multipliers application
    
      #Setting up of "y"
      mde_base <- y #loading of "y"
      mde_base[mde_base] = 1 #Setting up all values of "y" to 1
    
      #Multipliers for rivers and their buffers
      
#Code to allow parallel computation (from packages "doSNOW" and "parallel")
ncores <- detectCores() - 1 #detection of number of cores on the pc
      
        #Creation of rivers' buffers based on "z"
        buff_rius <- buffer(z, dist = 1000) #buffer creation
        
cl <- makeCluster(ncores, type = "SOCK") #creation of the cluster for parallel computing
registerDoSNOW(cl) #loading of the cluster 
library("foreach") #loading of packages inside the cluster
      
        #Computing of the rivers' buffers cost multipliers  
        mde_base[buff_rius] <- foreach(mult = mde_base[buff_rius], .combine = "cbind", .packages = c("sf", "raster", "sp", "gdistance")) %dopar% {mult * 1.5}   
      
        #Computing of the rivers cost multipliers  
        mde_base[z] <- foreach(mult = mde_base[z], .combine = "cbind", .packages = c("sf", "raster", "sp", "gdistance")) %dopar% {mult *2}
      
print("Step 1,5 completed") 
      
      #Slope cost multiplier calculation
      pendent <- terrain(y, opt="slope", unit="degrees", neighbors= 8) #Slope's raster calculation from "y"
      pendent <- tan(pendent * pi / 180) * 100 #Conversion of slope's raster values to percentages
      
        #Computing of the >60% slopes  cost multipliers  
        mde_base[pendent > 60] <- foreach(mult = mde_base[pendent > 60] , .combine = "cbind", .packages = c("sf", "raster", "sp", "gdistance")) %dopar% {mult * 3} 
        
        #Computing of the >100% slopes  cost multipliers
        mde_base[pendent > 100] <- foreach(mult = mde_base[pendent > 100] , .combine = "cbind", .packages = c("sf", "raster", "sp", "gdistance")) %dopar% {mult * 10}
      
      #conversion of the cost multiplers' raster to factor raster
      mde_base_cat <- asFactor(mde_base)
      
stopCluster(cl) #finishing the cluster
      
      #Creation of the Transition Layer for cost multipliers
      multiplicador <- transition(elevacio_cat, transitionFunction = "barriers", directions = 16, symm = TRUE) #Creation of Transition Stack from categorical raster
      multiplicador_unif <- sum(multiplicador, na.rm = FALSE) #Creation of Transition Layer from Transition Stack
      
print("Step 2 completed") 


    #Cost function application

      #Setting up of Correction's process for values of adjacent cells
      adj <- adjacent(elevacio, cells = 1:ncell(elevacio), pairs = TRUE, directions = 16) 
      
      #Choosing of the cost function
      if (w == "TOBLER") {
      
        #Calculation of cost values for moving from one cell to another (Tolber, 1993 cost function)
        speed <- slope * multiplicador_unif #Speed calculation's setting up
        speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] * multiplicador_unif[adj] + 0.05)) #Equation for speed calculation
      
      } else {
      
        #Calculation of cost values for moving from one cell to another (Llober & Sluckin, 2007 cost function)
        
          #Calculation of values for directions
          aspect_mde <- terrain(elevacio, opt = "aspect") #Recovery of direction from slope
          direccio <- function(d) {cos((d[2]-d[1])-d[1])} #Setting up of direction's calculation function
          trans_direccio <- transition(aspect_mde, direccio, 16, symm= FALSE) #Creation of direction's transiotion layer
          direction <- geoCorrection(trans_direccio) #Geocorrection to apply to the values for directions transition layer
      
          #Calculation of cost values for moving from one cell to another
          speed <- slope * direction * multiplicador_unif #Speed calculation's setting up
          speed[adj] <- 2.635 + 17.37 * (slope[adj]*direction[adj]* multiplicador_unif[adj]) + 42.37 * ((slope[adj]*direction[adj])^2* multiplicador_unif[adj]) - 21.43 * ((slope[adj]*direction[adj])^3* multiplicador_unif[adj]) + 14.93 * ((slope[adj]*direction[adj])^4 * multiplicador_unif[adj])  ##Equation for speed calculation
          
          }
      
    Conductance <- geoCorrection(speed) #Geocorrection to apply to the values for speed's transition layer
    
print("Step 3 completed") 
    
  
  #Computation of LCPs for ZPMs

    #Creation of empty object to store LCPs
    xarxa_camins_total <- SpatialLinesDataFrame(sl<-SpatialLines(LinesList<-list(Lines(list(Line(coords<-matrix(c(0,1,0,1),2,2))),ID=1))),data<-data.frame(a=1,b=2,c=3))[-1,]
 
    #Loading of the number of each edge points for loop functionality
    num_elementos_1 <- length(punts_limit_1)
    num_elementos_2 <- length(punts_limit_2)
    num_elementos_3 <- length(punts_limit_3)
    num_elementos_4 <- length(punts_limit_4)

print("Step 3,1 completed") 
  
  #Creation of LCPs 
    
    #Loop for Edge 1
    for (i in 1:num_elementos_1) {
    
      #Recovery of coordinates of origin point
      punt_objectiu_1 <- punts_limit_1[i, ]
    
      #LCPs calculation to final points of edge 1
      xarxa_camins <- shortestPath(Conductance, punt_objectiu_1, punts_desti_1, output = "SpatialLines")
    
      #Sumation of all edge 1 LCPs
      xarxa_camins_total <- bind(xarxa_camins_total, xarxa_camins)
    
      }
  
print("Step 3,2 completed")
    
    #Loop for Edge 2
    for (i in 1:num_elementos_2) {
    
      punt_objectiu_2 <- punts_limit_2[i, ]
      xarxa_camins <- shortestPath(Conductance, punt_objectiu_2, punts_desti_2, output = "SpatialLines")
      xarxa_camins_total <- bind(xarxa_camins_total, xarxa_camins)
    
      }
  
print("Step 3,3 completed")
  
    #Loop for Edge 3
    for (i in 1:num_elementos_3) {
    
      punt_objectiu_3 <- punts_limit_3[i, ]
      xarxa_camins <- shortestPath(Conductance, punt_objectiu_3, punts_desti_3, output = "SpatialLines")
      xarxa_camins_total <- bind(xarxa_camins_total, xarxa_camins)
    
      }
  
print("Step 3,4 completed")
  
    #Loop for Edge 4
    for (i in 1:num_elementos_4) {
    
      punt_objectiu_4 <- punts_limit_4[i, ]
      xarxa_camins <- shortestPath(Conductance, punt_objectiu_4, punts_desti_4, output = "SpatialLines")
      xarxa_camins_total <- bind(xarxa_camins_total, xarxa_camins)
    
      }
  
print("Step 4 completed") 
    
    #Saving of the LCPs to Shapefile
    nom_arxiu <- paste(substitute(y),substitute(q),substitute(w),"xarxa_camins_total_v5.shp") #File name creation
    xarxa_camins_total_sf <- st_as_sf(xarxa_camins_total) #Conversion of LCPs object from sp to sf (from package "sf")
    st_write(xarxa_camins_total_sf, dsn = nom_arxiu) #Saving of LCPs object from sf object to Shapefile

print("Step 5 completed") 


  #Rasteritzation of the LCPs net to ZPMs
  
    #Loading of the Shapefile and filtration of the NULL values
    mapa_buit <- st_is_empty(xarxa_camins_total_sf) #Identification of NULL values
    camins <- xarxa_camins_total_sf[!mapa_buit, ] #Filtration of NULL values
  
      #Creation of a empty raster for rasteritzation
      raster_buit <- raster(ext = extent_raster, res = 200)
  
cl <- makeCluster(ncores, type = "SOCK")
registerDoSNOW(cl)
library("foreach")
  
    #Counting of the number of LCP in each cell
    raster_count <- foreach(i= camins,.packages = c("raster", "sf")) %dopar% {rasterize(i, raster_buit, fun = "count") }
  
    #Improving of the raster's resolution
    raster_count_2 <- foreach(i2= raster_count[1], .packages = c("raster","sf")) %dopar% {raster::resample(i2, y, method = "bilinear")}
  
stopCluster(cl)

    #Saving the raster in a readable format
    raster_count_3 <- raster(raster_count_2[[1]])
    raster_count_3[ ] <- getValues((raster_count_2[[1]]))
  
print("Step 6 completed")

  
    #Saving of the file to TIFF
    nom_arxiu2 <- paste(substitute(y),substitute(q),substitute(w),"rasterize.tif") #File name creation
    writeRaster(raster_count_3, filename = nom_arxiu2) #Saving ZPMs raster to TIFF

print("Process ended succesfully")

}

