library(data.table)
library(dplyr)
library(bit64)
library(gdata)
library(dplyr)
library(usethis)
library(devtools)
library(evobiR)
library("dbscan")
library("fpc")
library(yaImpute)
library(intrinsicDimension)
library(plm)
library(lubridate)
library(foreach)
library(opticskxi)
library(parallel)
library(zoo)

########################################## OPTICS Algorithm ################################################################
#### OIL ####
options(scipen=999)
setwd("F:/MK/Prod_US/FracFocus/Jaccard/PercentHFJob")
load("jaccard_oil.RData")
rm(data_bin_oil, jaccard_f)
gc(1:1000)

jaccard_index_oil[is.na(jaccard_index_oil)]<-0
jaccard_index_oil[jaccard_index_oil>1]<-1
jaccard_index_oil<-1-jaccard_index_oil
diag(jaccard_index_oil)<-0
gc(1:1000)

#subset states
data_oil<-data_oil[data_oil$`State/Province` %in% c("OK","TX","UT","PA", "KS"),]

jaccard_index_oil<-jaccard_index_oil[rownames(jaccard_index_oil) %in% data_oil$APINumber,
                                     colnames(jaccard_index_oil) %in% data_oil$APINumber]

# optics algorithm

# SlidingWindow function
SlidingWindow <- function(FUN, data, window, step){
  total.x <- ncol(data)
  spots.x <- seq(from = 1, to = (total.x-window), by = step)
  total.y <- nrow(data)
  spots.y <- seq(from = 1, to = (total.y-window), by = step)
  result <- vector("list", length(spots.y))
  for(i in 1:length(spots.y)){
    result[[i]] <- match.fun(FUN)(as.dist(data[spots.y[i]:(spots.y[i] + window - 1),
                                                 spots.x[i]:(spots.x[i] + window - 1)]))
    
    gc(1:100)
  }
  return(result)
}

# optics function
optics_func<-function(x, mp, e) dbscan::optics(x, eps=e, minPts = mp, search = "dist")

# optics output

output_f<-function(jaccard, mp, e, window, step) {
  exper <- SlidingWindow(FUN=function(y) optics_func(y,mp,e), data=jaccard, window=window, step=step)
  gc(1:100)
  return(exper)
}


out_oil_res<-output_f(jaccard_index_oil,mp=6,e=0.3,window=501,step=1)
gc(1:1000)

out_oil_res_1<-output_f(jaccard_index_oil,mp=23,e=0.4,window=501,step=1)
gc(1:1000)

out_oil_res_2<-output_f(jaccard_index_oil,mp=60,e=0.5,window=501,step=1)
gc(1:1000)

setwd("F:/MK/Prod_US/FracFocus/OPTICS/PercentHFJob")
rm(jaccard_index_oil)
save.image("optics_oil_full.RData")
rm(list = ls())

#### GAS ####
options(scipen=999)
setwd("F:/MK/Prod_US/FracFocus/Jaccard/PercentHFJob")
load("jaccard_gas.RData")
rm(data_bin_gas, jaccard_f)
gc(1:1000)

jaccard_index_gas[is.na(jaccard_index_gas)]<-0
jaccard_index_gas[jaccard_index_gas>1]<-1
jaccard_index_gas<-1-jaccard_index_gas
diag(jaccard_index_gas)<-0

gc(1:1000)

#subset states
data_gas<-data_gas[data_gas$`State/Province` %in% c("OK","TX","UT","PA", "KS"),]

jaccard_index_gas<-jaccard_index_gas[rownames(jaccard_index_gas) %in% data_gas$APINumber,
                                     colnames(jaccard_index_gas) %in% data_gas$APINumber]

# optics algorithm
# SlidingWindow function
SlidingWindow <- function(FUN, data, window, step){
  total.x <- ncol(data)
  spots.x <- seq(from = 1, to = (total.x-window), by = step)
  total.y <- nrow(data)
  spots.y <- seq(from = 1, to = (total.y-window), by = step)
  result <- vector("list", length(spots.y))
  for(i in 1:length(spots.y)){
    result[[i]] <- match.fun(FUN)(as.dist(data[spots.y[i]:(spots.y[i] + window - 1),
                                               spots.x[i]:(spots.x[i] + window - 1)]))
    
    gc(1:100)
  }
  return(result)
}

# optics function
optics_func<-function(x, mp, e) dbscan::optics(x, eps=e, minPts = mp, search = "dist")

# optics output

output_f<-function(jaccard, mp, e, window, step) {
  exper <- SlidingWindow(FUN=function(y) optics_func(y,mp,e), data=jaccard, window=window, step=step)
  gc(1:100)
  return(exper)
}

out_gas_res<-output_f(jaccard_index_gas,mp=7,e=0.3,window=501,step=1)
gc(1:1000)

out_gas_res_1<-output_f(jaccard_index_gas,mp=28,e=0.4,window=501,step=1)
gc(1:1000)

out_gas_res_2<-output_f(jaccard_index_gas,mp=72,e=0.5,window=501,step=1)
gc(1:1000)

setwd("F:/MK/Prod_US/FracFocus/OPTICS/PercentHFJob")
rm(jaccard_index_gas)
save.image("optics_gas_full.RData")
rm(list = ls())
