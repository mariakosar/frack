## Load packages

library(data.table)
library(plyr)
library(dplyr)
library(bit64)
library(gdata)
library(usethis)
library(devtools)
library(evobiR)
library("dbscan")
library("fpc")
library(yaImpute)
library(intrinsicDimension)
library(plm)
library(lfe)
library(zoo)
library(httr)
library(jsonlite)
library("maptools")
library(ggplot2)
library(GISTools)
library(rgdal)
library(sf)
library(multiwayvcov)
library("stargazer")
library(geosphere)
library(tidyr)

options(scipen=999)
gc(1:1000)

#### add 6x6 miles identifier #### 
gc(1:1000)
shp <- readOGR('Countries_WGS84/Countries_WGS84.shp')
shp<-shp[shp@data$CNTRY_NAME=="United States",]
shp<-st_as_sf(shp)
shp<-shp %>% st_transform(32610)
# 
# # create 10km grid
grid_6 <- st_make_grid(shp, cellsize = c(10000, 10000)) %>%
   st_sf(grid_id = 1:length(.))

# # create labels for each grid_id
grid_lab <- st_centroid(grid_6) %>% cbind(st_coordinates(.))
save.image("WGS84.RData")

# coordinates oil
load("WGS84.RData")

colnames(ex_oil)[6]<-"Surface.Hole.Latitude"
colnames(ex_oil)[7]<-"Surface.Hole.Longitude"

ex_oil<-ex_oil[!is.na(ex_oil$Surface.Hole.Latitude),]
ex_oil<-ex_oil[!is.na(ex_oil$Surface.Hole.Longitude),]

pts_oil <- st_as_sf(ex_oil, coords = c("Surface.Hole.Longitude", "Surface.Hole.Latitude"), crs = 4326)
pts_oil<-st_transform(x = pts_oil, crs = 32610)

# which grid square is each point in?
pts_oil<-pts_oil %>% st_join(grid_6, join = st_intersects) %>% as.data.frame

pts_oil$tic_grid_id<-paste0(pts_oil$tic,"", pts_oil$grid_id)

pts_oil<-pts_oil[order(pts_oil$tic, pts_oil$month),]

# coordinates gas
colnames(ex_gas)[6]<-"Surface.Hole.Latitude"
colnames(ex_gas)[7]<-"Surface.Hole.Longitude"

ex_gas<-ex_gas[!is.na(ex_gas$Surface.Hole.Latitude),]
ex_gas<-ex_gas[!is.na(ex_gas$Surface.Hole.Longitude),]

pts_gas <- st_as_sf(ex_gas, coords = c("Surface.Hole.Longitude", "Surface.Hole.Latitude"), crs = 4326)
pts_gas<-st_transform(x = pts_gas, crs = 32610)

# which grid square is each point in?
pts_gas<-pts_gas %>% st_join(grid_6, join = st_intersects) %>% as.data.frame

pts_gas$tic_grid_id<-paste0(pts_gas$tic,"", pts_gas$grid_id)

pts_gas$tic_grid_id<-paste0(pts_gas$tic,"", pts_gas$grid_id)

pts_gas<-pts_gas[order(pts_gas$tic, pts_gas$month),]

rm(grid_6, grid_lab)
gc(1:1000)

