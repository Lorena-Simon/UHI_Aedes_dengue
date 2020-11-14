setwd("/Users/lorena.mendes.simon/Documents/Doutorado/Projeto Sanduiche/Projeto:artigo")

load("Script_final_SEM_Kernel.Rdata")

#### Point density Analysis ####


### Data ###

## Shapefile of US Southeast regions ##

library(rgdal)
library(raster)
library(sf)


US_Shp <- readOGR("/Users/lorena.mendes.simon/Desktop/Projeto Sanduiche/Projeto:artigo/Dados/GIS_Curis/csa_tmp_dissolve2.shp") #Curtis data of EUA Southeast

rs_proj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # Reproject the shape 
US_Shp_proj <- spTransform(US_Shp, rs_proj)

South_US_Sf <- st_as_sf(US_Shp_proj, "sf") # Make the SpatialPolygon into "sf" compatible object


gc()
South_US_Union <- st_union(South_US_Sf)

South_US_Union <- as(South_US_Union, Class = "Spatial")


## Ocurrence data (South US) ##

Ae_Aegyp <- read.csv("/Users/lorena.mendes.simon/Desktop/Projeto Sanduiche/Projeto:artigo/Dados/Aedes/Aeg_Point_final.csv", header = TRUE, sep = ";")

Ae_albo <- read.csv("/Users/lorena.mendes.simon/Desktop/Projeto Sanduiche/Projeto:artigo/Dados/Aedes/Alb_Point_final.csv", header = TRUE, sep = ";")

library(data.table)
library(sf)

Aeg_as_Matrix <- as.data.table(Ae_Aegyp)
sf_Aeg <- st_as_sf(Aeg_as_Matrix, coords = c("X","Y"), crs = 4326, agr = "constant")

Alb_as_Matrix <- as.data.table(Ae_albo)
sf_Alb <- st_as_sf(Alb_as_Matrix, coords = c("X","Y"), crs = 4326, agr = "constant")

# Figure

South_fig <- tm_shape(US_Shp_proj) + tm_fill(alpha = 0.3) + tm_polygons(border.alpha = 0.5) + tm_layout(frame = FALSE)

Aeg_points <- South_fig + tm_shape(sf_Aeg) + tm_dots(col = "darkorange3", size = 0.3, alpha = 0.8) 

Alb_points <- South_fig + tm_shape(sf_Alb) + tm_dots(col = "darkorange3", size = 0.3, alpha = 0.8)


### Rarefaction ###

library(letsR)
library(iNEXT)

# Aegytpi

Sp_Aeg <- paste("Ae_pt", 1:217, sep = "") # As if each column referred to a different specie
Aeg_presabs <- lets.presab.points(Ae_Aegyp, species = Sp_Aeg, xmn = -93, xmx = -75, ymn = 24, ymx = 40, resol = 0.25) # Create this Presence/Absence object based on South US region
rich <- values(Aeg_presabs$Richness_Raster)

rich.all <- iNEXT(rich, q = 1, datatype = "abundance", nboot = 1000)

ggiNEXT(rich.all, type = 2, facet.var = "site")


# Albopictus

Sp_Alb <- paste("Al_pt", 1:1221, sep = "")
Alb_presabs <- lets.presab.points(Ae_albo, species = Sp_Alb, xmn = -93, xmx = -75, ymn = 24, ymx = 40, resol = 0.5)
rich_2 <- values(Alb_presabs$Richness_Raster)

rich.all_2 <- iNEXT(rich_2, q = 1, datatype = "abundance", nboot = 1000)

ggiNEXT(rich.all_2, type = 2, facet.var = "site")


### Hypothesis test ###

library(spatstat)


SouthOwin <- as.owin(South_US_Union) ## South Observation Window

PPP_Aeg <- ppp(Ae_Aegyp[,1], Ae_Aegyp[,2], window = SouthOwin) ## Aegipty point pattern object

PPP_Alb <- ppp(Ae_albo[,1], Ae_albo[,2], window = SouthOwin) ## Albopictus point pattern object


ann_p_Alb <- mean(nndist(PPP_Alb, k = 1))


## Null model ##

# Aegytpi

ann_p_Aeg <- mean(nndist(PPP_Aeg, k = 1))

n <- 599L # number of simulations

ann_aeg <- vector( length = n) # Create an empty object to store simulated ANN values
for (i in 1:n) {
  rand.aeg <- rpoint(n=PPP_Aeg$n) # Generate random point locations
  ann_aeg [i] <- mean(nndist(rand.aeg, k = 1)) # Tally the ANN values
}

plot(rand.aeg, pch = 16, main = "aegypti Null", cols = rgb(0,0,0,0.5) )

hist(ann_aeg, main= expression(bolditalic("Aedes aegypti")), las=1, breaks=50, col="bisque", xlim=range(ann_p_Aeg, ann_aeg), xlab="")
abline(v=ann_p_Aeg, col="blue")


# Albopictus

ann_alb <- vector( length = n) 
for (i in 1:n) {
  rand.alb <- rpoint(n=PPP_Alb$n) 
  ann_alb [i] <- mean(nndist(rand.alb, k = 1)) 
}

plot(rand.alb, pch = 16, main = "albopictus Null", cols = rgb(0,0,0,0.5) )

hist(ann_alb, main= expression(bolditalic("Aedes albopictus")), las=1, breaks=50, col="bisque", xlim=range(ann_p_Alb, ann_alb), xlab="", ylab="")
abline(v=ann_p_Alb, col="blue")


### Kernel Density ###

library("KernSmooth")

Kernel_est_alb <- bkde2D(Ae_albo, bandwidth = c(0.5,0.5), gridsize = c(32,36), range.x = list(c(-93,-75), c(24,40)))

kernel_est_aeg <- bkde2D(Ae_Aegyp, bandwidth = c(0.5,0.5), gridsize = c(32,36), range.x = list(c(-93,-75), c(24,40)))


## create raster ##

# Albopisctus 

Kernel_raster_alb <- raster(list(x = Kernel_est_alb$x1, y = Kernel_est_alb$x2, z = Kernel_est_alb$fhat))
projection(Kernel_raster_alb) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
names(Kernel_raster_alb) <- "Kernel Density"

xmin(Kernel_raster_alb) <- -93
xmax(Kernel_raster_alb) <- -75
ymin(Kernel_raster_alb) <- 24
ymax(Kernel_raster_alb) <- 40

plot(Kernel_raster_alb) # Simple plot
plot(US_Shp_proj, add = T, lwd = 0.2)


map_nz = tm_shape(Kernel_raster_alb) + # Fancy Plot
  tm_raster(style="cont", 
            palette = "Oranges", n = 5, contrast = c(0, 1))+
  tm_legend(outside=TRUE)

map_nz1 = map_nz + tm_shape(US_Shp_proj) + tm_fill(alpha = 0) + tm_polygons(border.alpha = 0.5) 


save_tmap(map_nz1, filename = "Kernel_Alb.pdf")


## Aegypti ##

Kernel_raster_aeg <- raster(list(x = kernel_est_aeg$x1, y = kernel_est_aeg$x2, z = kernel_est_aeg$fhat))
projection(Kernel_raster_aeg) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
names(Kernel_raster_aeg) <- "Kernel Density"

xmin(Kernel_raster_aeg) <- -93
xmax(Kernel_raster_aeg) <- -75
ymin(Kernel_raster_aeg) <- 24
ymax(Kernel_raster_aeg) <- 40

plot(Kernel_raster_aeg)
plot(US_Shp_proj, add = T, lwd = 1)



map_nx = tm_shape(Kernel_raster_aeg) + # Fancy Plot
  tm_raster(style="cont", 
            palette = "Oranges", n = 5, contrast = c(0, 1))+
  tm_legend(outside=TRUE)

map_nx1 = map_nx + tm_shape(US_Shp_proj) + tm_fill(alpha = 0) + tm_polygons(border.alpha = 0.5) 

save_tmap(map_nx1, filename = "Kernel_Aeg.pdf")


#### Structural Equation Model ####

### Data ###

## Heat Island ##

US_extent <- extent(South_US_Union)

Heat_island_US <- crop(Heat_island, US_extent)

Mat_final <- as.data.frame(cbind( Heat_island_US@data$LONGITUDE,  Heat_island_US@data$LATITUDE,
                                  Heat_island_US@data$SQKM_FINAL,  Heat_island_US@data$URB_D_MEAN, 
                                  Heat_island_US@data$BUF_D_MEAN,  Heat_island_US@data$URB_N_MEAN,
                                  Heat_island_US@data$BUF_N_MEAN, Heat_island_US@data$D_T_DIFF,
                                  Heat_island_US@data$N_T_DIFF),1623, 7)

colnames(Mat_final) <- c("longitude","latitude", "SQKM_FINAL", "URB_D_MEAN", "BUF_D_MEAN", "URB_N_MEAN", "BUF_N_MEAN", "D_T_DIFF", "N_T_DIFF")


## Land Cover ##

f <- file.path("/Users/lorena.mendes.simon/Desktop/Projeto Sanduiche/Projeto:artigo/Dados/Land_cover", c("Consensus_reduced_class_1_.tif", "Consensus_reduced_class_3.tif",
                                                                                                         "Consensus_reduced_class_4.tif", "Consensus_reduced_class_6.tif",
                                                                                                         "Consensus_reduced_class_7.tif",  "Consensus_reduced_class_8.tif",
                                                                                                         "Consensus_reduced_class_9.tif"))

Land_cover_data <- lapply(f, raster)

Land_extract <- function(x){
  
  Cut <- crop(x, South_US_Union)
  Land_values <- extract(Cut, Mat_final[,c(1,2)]) # Voltar ao original
  return(Land_values)
  
}

Land_values <- sapply(Land_cover_data, Land_extract)

Land_values <- as.data.frame(Land_values)
colnames(Land_values) <- c("class_1", "class_3", "class_4", "class_6", "class_7", "class_8", "class_9")

Mat_final_2 <- cbind(Mat_final,Land_values)


## Nasa Power ##

Wind_speed <- raster("/Users/lorena.mendes.simon/Desktop/Projeto Sanduiche/Projeto:artigo/Dados/Nasa_Power/POWER_Global_Climatology_wind_speed_2mt/POWER_Global_Climatology_52daacfc_ann_ws2m_range_lst.tif")

Precipitation <- raster("/Users/lorena.mendes.simon/Desktop/Projeto Sanduiche/Projeto:artigo/Dados/Nasa_Power/POWER_Global_Climatology_Precipitation/POWER_Global_Climatology_9924f633_ann_prectot_lst.tif")

Nasa_Power_data <- as.list(Wind_speed,Precipitation)

Nasa_values <- sapply(Nasa_Power_data,Land_extract)

Nasa_values <- as.data.frame(Nasa_values)

colnames(Nasa_values) <- c("Wind_speed", "Precipitation")

Mat_final_3 <- cbind(Mat_final_2, Nasa_values)


## Dengue Suitability ##

# Aegypti

Ae_persist_DEN_Ae<-raster("/Users/lorena.mendes.simon/Desktop/Projeto Sanduiche/Projeto:artigo/Dados/Dengue_Suit/Annual_AUC_Ae_INTRO_DEN.tif") #(Global raster)

Aeg_Dengue_US <- crop(Ae_persist_DEN_Ae, South_US_Union)

Values_Aeg_Deng_US <- extract(Aeg_Dengue_US, Mat_final[,c(1,2)])

Mat_final_Aeg <- cbind( Mat_final_3,Values_Aeg_Deng_US) # Add Aegypti values to the Data frame


# Albopictus

Al_persist_DEN_Al <- raster("/Users/lorena.mendes.simon/Desktop/Projeto Sanduiche/Projeto:artigo/Dados/Dengue_Suit/Annual_AUC_Al_INTRO_DEN.tif")

Alb_Dengue_US <- crop(Al_persist_DEN_Al, South_US_Union)

Values_Alb_Deng_US <- extract(Alb_Dengue_US, Mat_final[,c(1,2)])

Mat_final_Alb <- cbind(Mat_final_3, Values_Alb_Deng_US) # Add albopictus values to the Data frame


## Kernel Density ##

# Aegypti

Values_Kernel_densi_2 <- extract(Kernel_raster_aeg, Mat_final[,c(1,2)])

Mat_final_Aeg_2 <- cbind(Mat_final_Aeg, Values_Kernel_densi_2)

Mat_final_Aeg_2 <- na.omit(Mat_final_Aeg_2) # Remove missing values

write.csv2(Mat_final_Aeg_2, file = "Tab_Aegi_SEM.csv", sep = ",")

# Albopictus

Values_Kernel_densi <- extract(Kernel_raster_alb, Mat_final[,c(1,2)])

Mat_final_Alb_2 <- cbind(Mat_final_Alb, Values_Kernel_densi)

Mat_final_Alb_2 <- na.omit(Mat_final_Alb_2)

write.csv2(Mat_final_Alb_2, file = "Tab_Albo_SEM.csv", sep = ",")


### Spatial Filters ###

## Adjust an OLS model between the variables ##

# Aegypti

Deng_Aeg_lm <- lm(Values_Aeg_Deng_US ~ SQKM_FINAL + URB_D_MEAN  + BUF_D_MEAN + URB_N_MEAN  + BUF_N_MEAN + class_1 + class_3 + 
                    class_4 + class_6 + class_7 + class_8 + class_9 + Wind_speed + Values_Kernel_densi_2 + Precipitation , data = Mat_final_Aeg_2)

summary(Deng_Aeg_lm)

# Albopictus

Deng_Alb_lm <- lm(Values_Alb_Deng_US ~ SQKM_FINAL + URB_D_MEAN  + BUF_D_MEAN + URB_N_MEAN  + BUF_N_MEAN + class_1 + class_3 + 
                    class_4 + class_6 + class_7 + class_8 + class_9 + Wind_speed + Precipitation, data = Mat_final_Alb_2)

summary(Deng_Alb_lm)


## Check the residual autocorrelation ##

library(letsR)
library(ape)
library(ecodist)

# Aegypti

Deng_Aeg_coord_dist <- lets.distmat(Mat_final_Aeg_2[,c(1,2)], asdist = FALSE) # Distance matrix 

Deng_Aeg_lm_moran <- lets.correl(residuals(Deng_Aeg_lm), Deng_Aeg_coord_dist, 10, plot = T) # Calc the Moran's I for different distance classes and look for the correlogram

mat_dist_w <- ifelse(test = Deng_Aeg_coord_dist > 300, yes = 4*300, no = Deng_Aeg_coord_dist)# Create a truncate matrix based on the distance where there is no autocorelation (check in correlogram)

w <- 1/mat_dist_w # Invert the matrix in order to shorter distances have more weight
diag(w) <- 0

Moran.I(residuals(Deng_Aeg_lm), w) # Corroborate Moran's I using the new matrix w


# Albopictus

Deng_Alb_coord_dist <- lets.distmat(Mat_final_Alb_2[,c(1,2)], asdist = FALSE) 

Deng_Alb_lm_moran <- lets.correl(residuals(Deng_Alb_lm), Deng_Alb_coord_dist, 10, plot = T) 

mat_dist_w_Alb <- ifelse(test = Deng_Alb_coord_dist > 350, yes = 4*300, no = Deng_Alb_coord_dist)

w_2 <- 1/mat_dist_w_Alb 
diag(w_2) <- 0

Moran.I(residuals(Deng_Alb_lm), w_2) 


## Make an eigeanalysis of the distance matrix ##

# Aegypti

Aegip_pcord <- pco(as.dist(mat_dist_w))

Aegip_vecold <- Aegip_pcord$vectors

n <- ncol(Aegip_vecold)

# Albopictus

Albo_pcord <- pco(as.dist(mat_dist_w_Alb))

Albo_vecold <- Albo_pcord$vectors

n_2 <- ncol(Albo_vecold)


## Select the vectors based on the minimization of global Moran's I ##

# Aegypti

for (i in 1:200) {
  print(i)
  ml <- lm (Mat_final_Aeg_2$Values_Aeg_Deng_US ~ Aegip_vecold[,1:i] + Mat_final_Aeg_2$SQKM_FINAL + Mat_final_Aeg_2$URB_D_MEAN + 
              Mat_final_Aeg_2$BUF_D_MEAN + Mat_final_Aeg_2$URB_N_MEAN + Mat_final_Aeg_2$BUF_N_MEAN + Mat_final_Aeg_2$class_1 +
              Mat_final_Aeg_2$class_3 + Mat_final_Aeg_2$class_4 + Mat_final_Aeg_2$class_6 + Mat_final_Aeg_2$class_7 +
              Mat_final_Aeg_2$class_8 + Mat_final_Aeg_2$class_9 + Mat_final_Aeg_2$Wind_speed + Mat_final_Aeg_2$Values_Kernel_densi_2 +
              Mat_final_Aeg_2$Precipitation)
  res <- ml$residuals
  mor <- Moran.I(res,w)
  p.mor <- mor$p.value
  if (p.mor<0.05){vec.sel <-(i+2)}
  if (p.mor>0.05){break}
  
}

Aegyp_vecsel <- Aegip_vecold[,1:vec.sel] # Store the vectors


# Albopictus

for (i in 1:200) {
  print(i)
  ml <- lm (Mat_final_Alb_2$Values_Alb_Deng_US ~ Albo_vecold[,1:i] + Mat_final_Alb_2$SQKM_FINAL + Mat_final_Alb_2$URB_D_MEAN + 
              Mat_final_Alb_2$BUF_D_MEAN + Mat_final_Alb_2$URB_N_MEAN + Mat_final_Alb_2$BUF_N_MEAN + Mat_final_Alb_2$class_1 +
              Mat_final_Alb_2$class_3 + Mat_final_Alb_2$class_4 + Mat_final_Alb_2$class_6 + Mat_final_Alb_2$class_7 +
              Mat_final_Alb_2$class_8 + Mat_final_Alb_2$class_9 + Mat_final_Alb_2$Wind_speed + Mat_final_Alb_2$Values_Kernel_densi +
              Mat_final_Alb_2$Precipitation)
  res <- ml$residuals
  mor <- Moran.I(res,w_2)
  p.mor <- mor$p.value
  if (p.mor<0.05){vec.sel2 <-(i+2)}
  if (p.mor>0.05){break}
  
}

Albo_vecsel <- Albo_vecold[,1:vec.sel2] 


## Adjust another OLS now including the spacial vectors ##

# Aegypti

Deng_Aeg_lm_vecs <- lm(Mat_final_Aeg_2$Values_Aeg_Deng_US ~ Mat_final_Aeg_2$SQKM_FINAL + Mat_final_Aeg_2$URB_D_MEAN + 
                         Mat_final_Aeg_2$BUF_D_MEAN + Mat_final_Aeg_2$URB_N_MEAN + Mat_final_Aeg_2$BUF_N_MEAN + Mat_final_Aeg_2$class_1 +
                         Mat_final_Aeg_2$class_3 + Mat_final_Aeg_2$class_4 + Mat_final_Aeg_2$class_6 + Mat_final_Aeg_2$class_7 +
                         Mat_final_Aeg_2$class_8 + Mat_final_Aeg_2$class_9 + Mat_final_Aeg_2$Wind_speed + Mat_final_Aeg_2$Values_Kernel_densi_2 +
                         Mat_final_Aeg_2$Precipitation + Aegyp_vecsel [,1:20]) 

summary(Deng_Aeg_lm_vecs)

vec_mat <- as.data.frame(cbind(Aegyp_vecsel[,1:20]))


# Albopictus

Deng_Albo_lm_vecs <- lm(Mat_final_Alb_2$Values_Alb_Deng_US ~ Mat_final_Alb_2$SQKM_FINAL + Mat_final_Alb_2$URB_D_MEAN + 
                          Mat_final_Alb_2$BUF_D_MEAN + Mat_final_Alb_2$URB_N_MEAN + Mat_final_Alb_2$BUF_N_MEAN + Mat_final_Alb_2$class_1 +
                          Mat_final_Alb_2$class_3 + Mat_final_Alb_2$class_4 + Mat_final_Alb_2$class_6 + Mat_final_Alb_2$class_7 +
                          Mat_final_Alb_2$class_8 + Mat_final_Alb_2$class_9 + Mat_final_Alb_2$Wind_speed + Mat_final_Alb_2$Values_Kernel_densi +
                          Mat_final_Alb_2$Precipitation + Albo_vecsel [,1:30])

summary(Deng_Albo_lm_vecs)

vec_mat_2 <- as.data.frame(cbind(Albo_vecsel[,1:30]))


### SEM ###

library(lavaan)
library(semPlot)
library(DiagrammeR)


## Aegypti ##

# List structured equations for lavaan

colnames(vec_mat) <- paste("V", 1:20, sep = "") # Filter columns names

Mat_Aeg_SEM_Fil <- cbind(Mat_final_Aeg_2[,3:20], vec_mat)

Mat_Aeg_SEM_Fil <- scale(Mat_Aeg_SEM_Fil[,c(1:5, 8:38)]) # scale the data with mean 0 and standard deviation 1 minus UHI data

DIFF_abs <- (Mat_final_Aeg_2[,c("D_T_DIFF", "N_T_DIFF")])

Mat_final_Aeg_SEM_Fil <- cbind(Mat_Aeg_SEM_Fil, DIFF_abs)


Den_Aeg_Model <- '

D_T_DIFF ~ class_9 + class_8

N_T_DIFF ~ class_9 + class_8


Values_Kernel_densi_2 ~ D_T_DIFF + N_T_DIFF + Wind_speed + Precipitation + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + 

V10 +  V11 + V12  + V13 +  V15 + V16 + V17 + V18 + V19 + V20 


Values_Aeg_Deng_US ~ D_T_DIFF + N_T_DIFF +  Wind_speed + Precipitation + Values_Kernel_densi_2 + V1 + V2 + V3 + V4 + V5 + V6 + V7 + 

V8 + V9 + V10 +  V11 + V12  + V13 +  V15 + V16 + V17 + V18 + V19 + V20 

'

Den_Aeg_Model_SEM <- sem(Den_Aeg_Model, Mat_final_Aeg_SEM_Fil, estimator = "MLM") # Fit vcov SEM

summary(Den_Aeg_Model_SEM, standardize = TRUE, fit.measures=TRUE) 

inspect(Den_Aeg_Model_SEM, "rsquare") # Get R2 for models


# Plot results

AEG <- grViz(" 
      digraph SEM{
      
      graph [layout = neato,
      overlap = true,
      outputorder = edgesfirst]
      
      node [shape = rectangle]
      
      
      b [pos = '-3.2,0.1!', label = <UHI<BR/>Daytime<BR/>R=0.06>]
      c [pos = '-3.2,-1.3!', label = <UHI<BR/>Nighttime<BR/>R=0.003>]
      
      e [pos = '-1.2, -0.6!',label = <<I>Ae.&nbsp;aegypti</I> <BR/>Transmission<BR/>R=0.92>, style = 'bold']
      f [pos = '-4.6,-0.5!', label = <Land&nbsp;use>]
      g [pos = '0.9,0.2!', label = <Wind&nbsp;speed>]
      h [pos = '0.9,-1.6!', label = <Precipitation>]
      i [pos = '-1.2,0.8!', label = <Spatial<BR/>Filters>, shape ='oval']
      j [pos = '-1.2,-2!', label = <Vector<BR/>Density<BR/>R=0.49>, style = 'bold']
      
      
      b -> e [label = '-0.03']
      c -> e [style = 'dashed']
      
      
      g -> e [label = '0.07']
      h -> e [style = 'dashed']
      i -> e [label = '0.05']
      j -> e [style = 'dashed']
      
      f -> b [label = '0.07']
      f -> c [label = '0.05']
      
      
      
      b -> j [label = '0.05']
      c -> j [label = '-0.04']
      
      g -> j [label = '-0.19']
      h -> j [label = '0.31']

      }
      ")

## Albopictus ##

colnames(vec_mat_2) <- paste("V", 1:30, sep = "")

Mat_Albo_SEM_Fil <- cbind(Mat_final_Alb_2[,3:20], vec_mat_2)

Mat_Albo_SEM_Fil <- scale(Mat_Albo_SEM_Fil[c(1:5, 8:48)]) # scale the data with mean 0 and standard deviation 1

DIFF_abs_2 <- (Mat_final_Alb_2[,c("D_T_DIFF", "N_T_DIFF")])

Mat_final_Albo_SEM_Fil <- cbind(Mat_Albo_SEM_Fil, DIFF_abs_2)


Den_Alb_Model <- '

D_T_DIFF ~ class_9 + class_8

N_T_DIFF ~ class_9 + class_8


Values_Kernel_densi ~ D_T_DIFF + N_T_DIFF + Wind_speed + Precipitation  

Values_Alb_Deng_US ~ D_T_DIFF + N_T_DIFF + Wind_speed + Precipitation + Values_Kernel_densi 
'

Den_Alb_Model_SEM <- sem(Den_Alb_Model, Mat_final_Albo_SEM_Fil, estimator = "MLM")

summary(Den_Alb_Model_SEM, standardize = TRUE, fit.measures=TRUE) 

inspect(Den_Alb_Model_SEM, "rsquare")


ALB <- grViz(" 
      digraph SEM{
      
      graph [layout = neato,
      overlap = true,
      outputorder = edgesfirst]
      
      node [shape = rectangle]
      
      
      b [pos = '-3.2,0.1!', label = <UHI<BR/>Daytime<BR/>R=0.06>]
      c [pos = '-3.2,-1.3!', label = <UHI<BR/>Nighttime<BR/>R=0.003>]
      
      e [pos = '-1.2, -0.6!',label = <<I>Ae.&nbsp;albopictus</I> <BR/>Transmission<BR/>R=0.90>, style = 'bold']
      f [pos = '-4.6,-0.5!', label = <Land&nbsp;use>]
      g [pos = '0.9,0.2!', label = <Wind&nbsp;speed>]
      h [pos = '0.9,-1.6!', label = <Precipitation>]
      i [pos = '-1.2,0.8!', label = <Spatial<BR/>Filters>, shape ='oval']
      j [pos = '-1.2,-2!', label = <Vector<BR/>Density<BR/>R=0.54>, style = 'bold']
      
      
      b -> e [label = '-0.03']
      c -> e [style = 'dashed']
      
      
      g -> e [label = '0.13']
      h -> e [style = 'dashed']
      i -> e [label = '0.02']
      j -> e [style = 'dashed']
      
      f -> b [label = '0.07']
      f -> c [label = '0.05']
      
      b -> j [style = 'dashed']
      c -> j [label = '-0.04']
      
      g -> j [label = '-0.07']
      h -> j [label = '0.23']

      }
      
      ")


## UHI Histogram with ggplot2 ##

library("ggplot2")
library("gridExtra")

 Day <- ggplot(data = DIFF_abs_2, aes(DIFF_abs_2$D_T_DIFF)) +
  geom_histogram(col = "Black", 
                 alpha = 0.6 ,
                 aes(fill=..x..)) +
  labs(title = "UHI Daytime", x="Range", y="Count") +
  scale_fill_gradientn("Range", colours = rev(heat.colors(30))) +
  theme_minimal()
  

Night <-  ggplot(data = DIFF_abs_2, aes(DIFF_abs_2$N_T_DIFF)) +
  geom_histogram(col = "Black", 
                 alpha = 0.6 ,
                 aes(fill=..x..)) +
  labs(title = "UHI Nighttime", x="Range", y="Count") +
  scale_fill_gradientn("Range", colours = rev(heat.colors(50))) +
  theme_minimal()

grid.arrange(Day,Night, nrow = 1)


save.image("Script_final_SEM_Kernel.Rdata")









