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
setwd("F:/MK/Prod_US/FracFocusMaria/OPTICS/PercentHFJob")
load("out_oil_full.RData")
load("out_gas_full.RData")
gc(1:1000)

rownames(out_gas)<-data_gas$APINumber[501:(nrow(data_gas)-1)]
rownames(out_oil)<-data_oil$APINumber[501:(nrow(data_oil)-1)]

setwd("F:/MK/Prod_US/FracFocusMaria/OPTICS/PercentHFJob/EXPERIMENTATION")

data_tics_1<-as.data.frame(unique(data_oil$`Operator Ticker`))
data_tics_2<-as.data.frame(unique(data_gas$`Operator Ticker`))
colnames(data_tics_1)<-"Operator Ticker"
colnames(data_tics_2)<-"Operator Ticker"

data_tics<-rbind(data_tics_1,data_tics_2)
data_tics<-as.data.frame(data_tics[!duplicated(data_tics$`Operator Ticker`),])
colnames(data_tics)<-"Operator Ticker"
data_tics$`Operator Ticker`<-as.character(data_tics$`Operator Ticker`)
rm(data_tics_1,data_tics_2)

#### Functions #### 

# experim
experim_f<-function(out, dat){
  dat<-dat[dat$APINumber %in% rownames(out),]
  dat<-dat[match(rownames(out), dat$APINumber),]
  experim<-cbind(out, dat[,c("API10", "JobStartDate","Operator (Reported)", "OperatorName","Operator Ticker",
                             "Surface Hole Latitude (WGS84)", "Surface Hole Longitude (WGS84)",
                             "Bottom Hole Latitude (WGS84)", "Bottom Hole Longitude (WGS84)","State/Province")])
  experim$month<-substr(experim$JobStartDate, 1, 7)
  colnames(experim)[1]<-"Outlier"
  experim<-experim[order(experim$'Operator (Reported)', experim$JobStartDate),]
  experim$Outlier<-ifelse(experim$Outlier==0,1,0)
  
  experim$delay<-ifelse(experim$`State/Province` %in% c("MT","TX","WY","MS"),30,
                        ifelse(experim$`State/Province` %in% c("MI","OK","PA","UT"),60,
                               ifelse(experim$`State/Province` %in% c("NM"),45,
                                      ifelse(experim$`State/Province` %in% c("LA"),20,
                                             120))))
  
  experim$DeclDate<-experim$JobStartDate + experim$delay
  experim<-experim[,c("API10", "month", "Operator Ticker", "Operator (Reported)","OperatorName","Surface Hole Latitude (WGS84)", "Surface Hole Longitude (WGS84)",
                      "Bottom Hole Latitude (WGS84)", "Bottom Hole Longitude (WGS84)","JobStartDate", "DeclDate",
                      "Outlier")]
  return(experim)
}

# debt data

#### Create experimentation variable ####
# Gas

ex_gas<-experim_f(out_gas, data_gas)
ex_gas<-ex_gas[order(ex_gas$'Operator Ticker', ex_gas$month),]


# OIL
ex_oil<-experim_f(out_oil, data_oil)
ex_oil<-ex_oil[order(ex_oil$'Operator Ticker', ex_oil$month),]

#### add child and experim well indicator ####

ex_gas<-cbind(ex_gas, midPoint(cbind(ex_gas$`Bottom Hole Longitude (WGS84)`,ex_gas$`Bottom Hole Latitude (WGS84)`),
                               cbind(ex_gas$`Surface Hole Longitude (WGS84)`, ex_gas$`Surface Hole Latitude (WGS84)`)))

ex_gas$interference<-NA
ex_gas$experim<-NA

for (i in 2:nrow(ex_gas)){
  coords = cbind(ex_gas$lon[i],ex_gas$lat[i])
  comparison_set = cbind(ex_gas$lon[1:(i-1)],ex_gas$lat[1:(i-1)])
  ex_gas$interference[i]<-ifelse(length(which(distHaversine(coords, comparison_set,r=6378137)<=321.869))==0,0,1)
  ex_gas$experim[i]<-ifelse(length(which(distHaversine(coords, comparison_set,r=6378137)>3218.69))==nrow(comparison_set),1,0)
  rm(coords, comparison_set)
}

ex_gas$interference<-ifelse(is.na(ex_gas$interference)==T,0,ex_gas$interference)
ex_gas$experim<-ifelse(is.na(ex_gas$experim)==T,0,ex_gas$experim)

#
ex_oil<-cbind(ex_oil, midPoint(cbind(ex_oil$`Bottom Hole Longitude (WGS84)`,ex_oil$`Bottom Hole Latitude (WGS84)`),
                               cbind(ex_oil$`Surface Hole Longitude (WGS84)`, ex_oil$`Surface Hole Latitude (WGS84)`)))

ex_oil$interference<-NA
ex_oil$experim<-NA

for (i in 2:nrow(ex_oil)){
  coords = cbind(ex_oil$lon[i],ex_oil$lat[i])
  comparison_set = cbind(ex_oil$lon[1:(i-1)],ex_oil$lat[1:(i-1)])
  ex_oil$interference[i]<-ifelse(length(which(distHaversine(coords, comparison_set,r=6378137)<=321.869))==0,0,1)
  ex_oil$experim[i]<-ifelse(length(which(distHaversine(coords, comparison_set,r=6378137)>3218.69))==nrow(comparison_set),1,0)
  rm(coords, comparison_set)
}

ex_oil$interference<-ifelse(is.na(ex_oil$interference)==T,0,ex_oil$interference)
ex_oil$experim<-ifelse(is.na(ex_oil$experim)==T,0,ex_oil$experim)

rm(i)
gc(1:1000)

