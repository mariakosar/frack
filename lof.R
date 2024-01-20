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
library(DMwR)

options(scipen=999)


########################################## LOF Algorithm ################################################################
#### OIL ####
options(scipen=999)
setwd("F:/MK/Prod_US/FracFocusMaria/OPTICS/PercentHFJob/INFOSET")
load("oil_infoset.RData")
rm(info_set, jaccard_f,info_set_oil, jaccard_index_oil)
gc(1:1000)

### Functions

# optics function
lof_func<-function(x, mp, col){
  outlier.scores<-DMwR::lofactor(x[,col:ncol(x)], k = mp)
  outlier.scores<-ifelse(is.nan(outlier.scores)==T|is.infinite(outlier.scores)==T, NA, outlier.scores)
  outlier.scores<-ifelse(outlier.scores>quantile(outlier.scores, 0.9, na.rm = T),1,0)
  outlier.scores[is.na(outlier.scores)]<-0
  return(outlier.scores)
  }

# optics output
output_f<-function(x,z, mp,col){
  result<-apply(x,2,function(y)
    match.fun(FUN=function(y) lof_func(y,mp,col))(z[z$APINumber %in% y,])
  )
  return(result)
}

#
out_oil_res<-output_f(info_oil, data_oil, 6, 21)
gc(1:1000)

out_oil_res_1<-output_f(info_oil, data_oil, 23, 21)
gc(1:1000)

out_oil_res_2<-output_f(info_oil, data_oil, 60, 21)
gc(1:1000)

rm(info_oil, lof_func, output_f)

setwd("F:/MK/Prod_US/FracFocusMaria/LOF")
save.image("lof_oil_infoset.RData")
rm(list = ls())
gc(1:1000)

#
setwd("F:/MK/Prod_US/FracFocusMaria/LOF")
load("lof_oil_infoset.RData")

out_oil<-as.data.frame(out_oil_res[501,])
colnames(out_oil)<-"Outlier"

out_oil_1<-as.data.frame(out_oil_res_1[501,])
colnames(out_oil_1)<-"Outlier"

out_oil_2<-as.data.frame(out_oil_res_2[501,])
colnames(out_oil_2)<-"Outlier"

rm(out_oil_res, out_oil_res_1, out_oil_res_2)

save.image("out_oil_infoset.RData")


#### GAS ####
options(scipen=999)
setwd("F:/MK/Prod_US/FracFocusMaria/OPTICS/PercentHFJob/INFOSET")
load("gas_infoset.RData")
rm(info_set, jaccard_f,info_set_gas, jaccard_index_gas)
gc(1:1000)

### Functions

# optics function
lof_func<-function(x, mp, col){
  outlier.scores<-DMwR::lofactor(x[,col:ncol(x)], k = mp)
  outlier.scores<-ifelse(is.nan(outlier.scores)==T|is.infinite(outlier.scores)==T, NA, outlier.scores)
  outlier.scores<-ifelse(outlier.scores>quantile(outlier.scores, 0.9, na.rm = T),1,0)
  outlier.scores[is.na(outlier.scores)]<-0
  return(outlier.scores)
}

# optics output
output_f<-function(x,z, mp,col){
  result<-apply(x,2,function(y)
    match.fun(FUN=function(y) lof_func(y,mp,col))(z[z$APINumber %in% y,])
  )
  return(result)
}

#
out_gas_res<-output_f(info_gas, data_gas, 6, 21)
gc(1:1000)

out_gas_res_1<-output_f(info_gas, data_gas, 23, 21)
gc(1:1000)

out_gas_res_2<-output_f(info_gas, data_gas, 60, 21)
gc(1:1000)

rm(info_gas, lof_func, output_f)


setwd("F:/MK/Prod_US/FracFocusMaria/LOF")
save.image("lof_gas_infoset.RData")
rm(list = ls())
gc(1:1000)


#
setwd("F:/MK/Prod_US/FracFocusMaria/LOF")
load("lof_gas_infoset.RData")

out_gas<-as.data.frame(out_gas_res[501,])
colnames(out_gas)<-"Outlier"

out_gas_1<-as.data.frame(out_gas_res_1[501,])
colnames(out_gas_1)<-"Outlier"

out_gas_2<-as.data.frame(out_gas_res_2[501,])
colnames(out_gas_2)<-"Outlier"

rm(out_gas_res, out_gas_res_1, out_gas_res_2)

save.image("out_gas_infoset.RData")


