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
library(plot.matrix)

########### Functions #########

## DBSCAN algorithm
# sliding window function
SlidingWindow <- function(FUN, data, window, step, last_value){
  total.x <- ncol(data)
  spots.x <- seq(from = 1, to = (total.x-window), by = step)
  total.y <- nrow(data)
  spots.y <- seq(from = 1, to = (total.y-window), by = step)
  result <- matrix(, nrow=length(spots.y), ncol=1)
  for(i in 1:length(spots.y)){
    result_prelim <- match.fun(FUN)(as.dist(data[spots.y[i]:(spots.y[i] + window - 1),
                                                 spots.x[i]:(spots.x[i] + window - 1)]))
    
    result[i,1] <-result_prelim[["cluster"]][last_value]
    gc(1:100)
  }
  return(result)
}

# dbscan function
dbscan_func<-function(x, mp, e) fpc::dbscan(x, eps = e, MinPts = mp, scale = FALSE, method = "dist")

# output

output<-function(jaccard, mp, e, window, last_value, step) {
  exper <- SlidingWindow(FUN=function(y) dbscan_func(y,mp,e), data=jaccard, window=window, step=step, last_value=last_value)
  gc(1:100)
  return(exper)
}


#### OIL ####

options(scipen=999)
gc(1:1000)

setwd("F:/MK/Prod_US/FracFocusMaria/Jaccard/PercentHFJob")
load("jaccard_oil.RData")
gc(1:1000)

### distance matrix 
jaccard_index_oil[is.na(jaccard_index_oil)]<-0
jaccard_index_oil[jaccard_index_oil>1]<-1
jaccard_index_oil<-1-jaccard_index_oil
diag(jaccard_index_oil) <- 0

rm(data_bin_oil)

### choose MinPts
median(colSums(ifelse(jaccard_index_oil<0.01,1,0)))
median(colSums(ifelse(jaccard_index_oil<0.05,1,0)))
median(colSums(ifelse(jaccard_index_oil<0.1,1,0)))

### choose epsilon
# all wells
dbscan::kNNdistplot(as.dist(jaccard_index_oil), k =  5)
abline(h = 0.3, lty = 2)

dbscan::kNNdistplot(as.dist(jaccard_index_oil), k =  20)
abline(h = 0.4, lty = 2)

dbscan::kNNdistplot(as.dist(jaccard_index_oil), k =  56)
abline(h = 0.5, lty = 2)

#excluding early wells
window<-1000
dbscan::kNNdistplot(as.dist(jaccard_index_oil[(window+1):nrow(jaccard_index_oil),
                                              (window+1):ngccol(jaccard_index_oil)]), k =  5)
abline(h = 0.3, lty = 2)

### dbscan
window<-1000
out_oil<-output(jaccard_index_oil[order(nrow(jaccard_index_oil):1),order(ncol(jaccard_index_oil):1)], mp=5, e=0.3, window=window, last_value=1, step=1)
gc(1:1000)

out_oil_1<-output(jaccard_index_oil[order(nrow(jaccard_index_oil):1),order(ncol(jaccard_index_oil):1)], mp=20, e=0.4, window=window, last_value=1, step=1)
gc(1:1000)

out_oil_2<-output(jaccard_index_oil[order(nrow(jaccard_index_oil):1),order(ncol(jaccard_index_oil):1)], mp=56, e=0.5, window=window, last_value=1, step=1)
gc(1:1000)

rm(jaccard_index_oil)
#setwd("//tsclient/D/Prod_US/FracFocusMaria/DBSCAN")
save.image("out_oil_full.Rdata")
rm(list = ls())
gc(1:1000)


#### GAS ####
options(scipen=999)
gc(1:1000)

setwd("F:/MK/Prod_US/FracFocusMaria/Jaccard/PercentHFJob")
load("jaccard_gas.RData")
jaccard_index_gas[is.na(jaccard_index_gas)]<-0
jaccard_index_gas[jaccard_index_gas>1]<-1
jaccard_index_gas<-1-jaccard_index_gas
diag(jaccard_index_gas) <- 0

rm(data_bin_gas)

### choose MinPts
median(colSums(ifelse(jaccard_index_gas<0.01,1,0)))
median(colSums(ifelse(jaccard_index_gas<0.05,1,0)))
median(colSums(ifelse(jaccard_index_gas<0.1,1,0)))

### choose epsilon
# all wells
dbscan::kNNdistplot(as.dist(jaccard_index_gas), k =  7)
abline(h = 0.3, lty = 2)

dbscan::kNNdistplot(as.dist(jaccard_index_gas), k =  30)
abline(h = 0.4, lty = 2)

dbscan::kNNdistplot(as.dist(jaccard_index_gas), k =  73)
abline(h = 0.5, lty = 2)

#excluding early wells
window<-1000
dbscan::kNNdistplot(as.dist(jaccard_index_gas[(window+1):nrow(jaccard_index_gas),
                                              (window+1):ncol(jaccard_index_gas)]), k =  30)
abline(h = 0.3, lty = 2)

### gas 
window<-1000
out_gas<-output(jaccard_index_gas[order(nrow(jaccard_index_gas):1),order(ncol(jaccard_index_gas):1)], mp=7, e=0.3, window=window, last_value=1, step=1)
gc(1:1000)

out_gas_1<-output(jaccard_index_gas[order(nrow(jaccard_index_gas):1),order(ncol(jaccard_index_gas):1)], mp=30, e=0.4, window=window, last_value=1, step=1)
gc(1:1000)

out_gas_2<-output(jaccard_index_gas[order(nrow(jaccard_index_gas):1),order(ncol(jaccard_index_gas):1)], mp=73, e=0.5, window=window, last_value=1, step=1)
gc(1:1000)


rm(jaccard_index_gas)
#setwd("//tsclient/D/Prod_US/FracFocusMaria/DBSCAN")
save.image("out_gas_full.RData")
rm(list = ls())



