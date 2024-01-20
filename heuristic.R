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
setwd("F:/MK/Prod_US/FracFocusMaria/Jaccard/PercentHFJob")
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
    result[[i]] <- match.fun(FUN)(data[spots.y[i]:(spots.y[i] + window - 1),
                                                 spots.x[i]:(spots.x[i] + window - 1)][2:window,1])
    
    gc(1:100)
  }
  return(result)
}

# optics function
cluster_func<-function(x, e) sum(x<=e)

# optics output

output_f<-function(jaccard, mp, e, window, step) {
  exper <- SlidingWindow(FUN=function(y) cluster_func(y,e), data=jaccard, window=window, step=step)
  gc(1:100)
  return(exper)
}


out_oil_res<-output_f(jaccard_index_oil,e=0.3,window=501,step=1)
gc(1:1000)

out_oil_res_1<-output_f(jaccard_index_oil,e=0.2,window=501,step=1)
gc(1:1000)

out_oil_res_2<-output_f(jaccard_index_oil,e=0.1,window=501,step=1)
gc(1:1000)

out_oil_res_3<-output_f(jaccard_index_oil,e=0.05,window=501,step=1)
gc(1:1000)

out_oil_res_4<-output_f(jaccard_index_oil,e=0.4,window=501,step=1)
gc(1:1000)

out_oil_res_5<-output_f(jaccard_index_oil,e=0.5,window=501,step=1)
gc(1:1000)

setwd("F:/MK/Prod_US/FracFocusMaria/SIMPLE_CLUSTERING")
rm(jaccard_index_oil)
save.image("cluster_oil_full.RData")
rm(list = ls())

#### GAS ####
options(scipen=999)
setwd("F:/MK/Prod_US/FracFocusMaria/Jaccard/PercentHFJob")
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
    result[[i]] <- match.fun(FUN)(data[spots.y[i]:(spots.y[i] + window - 1),
                                       spots.x[i]:(spots.x[i] + window - 1)][2:window,1])
    
    gc(1:100)
  }
  return(result)
}

# optics function
cluster_func<-function(x, e) sum(x<=e)

# optics output

output_f<-function(jaccard, mp, e, window, step) {
  exper <- SlidingWindow(FUN=function(y) cluster_func(y,e), data=jaccard, window=window, step=step)
  gc(1:100)
  return(exper)
}

out_gas_res<-output_f(jaccard_index_gas,e=0.3,window=501,step=1)
gc(1:1000)

out_gas_res_1<-output_f(jaccard_index_gas,e=0.2,window=501,step=1)
gc(1:1000)

out_gas_res_2<-output_f(jaccard_index_gas,e=0.1,window=501,step=1)
gc(1:1000)

out_gas_res_3<-output_f(jaccard_index_gas,e=0.05,window=501,step=1)
gc(1:1000)

out_gas_res_4<-output_f(jaccard_index_gas,e=0.4,window=501,step=1)
gc(1:1000)

out_gas_res_5<-output_f(jaccard_index_gas,e=0.5,window=501,step=1)
gc(1:1000)


setwd("F:/MK/Prod_US/FracFocusMaria/SIMPLE_CLUSTERING")
rm(jaccard_index_gas)
save.image("cluster_gas_full.RData")
rm(list = ls())

######################################### Extract Outliers of the best model ###################################################################
#### OIL ####
load("F:/MK/Prod_US/FracFocusMaria/SIMPLE_CLUSTERING/cluster_oil_full.RData")
options(scipen=999)
gc(1:1000)

out_oil_res<-do.call("rbind", out_oil_res)
rownames(out_oil_res)<-data_oil$APINumber[501:(nrow(data_oil)-1)]
out_oil<-as.data.frame(out_oil_res)
sum(out_oil==0)

out_oil_res_1<-do.call("rbind", out_oil_res_1)
rownames(out_oil_res_1)<-data_oil$APINumber[501:(nrow(data_oil)-1)]
out_oil_1<-as.data.frame(out_oil_res_1)
sum(out_oil_1==0)

out_oil_res_2<-do.call("rbind", out_oil_res_2)
rownames(out_oil_res_2)<-data_oil$APINumber[501:(nrow(data_oil)-1)]
out_oil_2<-as.data.frame(out_oil_res_2)
sum(out_oil_2==0)

rm(out_oil_res,out_oil_res_1, out_oil_res_2, info_oil, jaccard_index_oil,cluster_func, output_f, i, SlidingWindow)

save.image("F:/MK/Prod_US/FracFocusMaria/SIMPLE_CLUSTERING/out_oil_full.RData")
rm(list = ls())
gc(1:1000)

#### GAS ####
load("F:/MK/Prod_US/FracFocusMaria/SIMPLE_CLUSTERING/cluster_gas_full.RData")
options(scipen=999)
gc(1:1000)

out_gas_res<-do.call("rbind", out_gas_res)
rownames(out_gas_res)<-data_gas$APINumber[501:(nrow(data_gas)-1)]
out_gas<-as.data.frame(out_gas_res)
sum(out_gas==0)

out_gas_res_1<-do.call("rbind", out_gas_res_1)
rownames(out_gas_res_1)<-data_gas$APINumber[501:(nrow(data_gas)-1)]
out_gas_1<-as.data.frame(out_gas_res_1)
sum(out_gas_1==0)

out_gas_res_2<-do.call("rbind", out_gas_res_2)
rownames(out_gas_res_2)<-data_gas$APINumber[501:(nrow(data_gas)-1)]
out_gas_2<-as.data.frame(out_gas_res_2)
sum(out_gas_2==0)

rm(out_gas_res,out_gas_res_1, out_gas_res_2, info_gas, jaccard_index_gas, cluster_func,output_f, i, SlidingWindow)

save.image("F:/MK/Prod_US/FracFocusMaria/SIMPLE_CLUSTERING/out_gas_full.RData")
rm(list = ls())
gc(1:1000)
