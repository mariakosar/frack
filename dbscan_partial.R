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

########### Functions ###########
# dbscan function
dbscan_func<-function(x, mp, e) fpc::dbscan(x, eps = e, MinPts = mp, scale = FALSE, method = "dist")

# count

count_f<-function (x){
  print(count(x[,"APINumber"]))
}

# dbscan output

output_f<-function(x,z, mp,e){
  result<-apply(x,2,function(y)
     match.fun(FUN=function(y) dbscan_func(y,mp,e))(as.dist(z[rownames(z) %in% y,
                                                                                     colnames(z) %in% y]))[["cluster"]][1001]
         )
  return(result)
}

# dimentions

dims<-function(x,h,n){
  print(pcaOtpmPointwiseDimEst(x[,h:ncol(x)], N=n))
}

# info set

info_set<-function(x){
  x<-x[order(x$DeclDate),]
  inf_set<-as.data.frame(outer(x$DeclDate, x$JobStartDate, `<`))
  colnames(inf_set)<-x$APINumber
  rownames(inf_set)<-x$APINumber
  inf_set[inf_set==FALSE]<-0

  gc(1:1000)
  inf_set_2<-as.data.frame(outer(x$JobStartDate, x$JobStartDate, `<`))
  colnames(inf_set_2)<-x$APINumber
  rownames(inf_set_2)<-x$APINumber
  inf_set_2[inf_set_2==FALSE]<-0
  gc(1:1000)

  inf_set_3<-as.data.frame(outer(x$`Operator Company Name`, x$`Operator Company Name`, `==`))
  colnames(inf_set_3)<-x$APINumber
  rownames(inf_set_3)<-x$APINumber
  inf_set_3[inf_set_3==FALSE]<-0
  gc(1:1000)

  inf_set_2<-inf_set_2+inf_set_3
  inf_set_2[inf_set_2==1]<-0
  inf_set_2[inf_set_2==2]<-1

  rm(inf_set_3)

  inf_set<-inf_set+inf_set_2
  inf_set[inf_set==2]<-1

  rm(inf_set_2)
  gc(1:1000)


  inf_set<-inf_set[,colSums(inf_set)>=1000,]
  gc(1:1000)

  return(inf_set)

}

#### OIL ####
options(scipen=999)
setwd("F:/MK/Prod_US/FracFocusMaria/Jaccard/PercentHFJob")
load("jaccard_oil.RData")
gc(1:1000)

jaccard_index_oil[is.na(jaccard_index_oil)]<-0
jaccard_index_oil[jaccard_index_oil>1]<-1
jaccard_index_oil<-1-jaccard_index_oil
diag(jaccard_index_oil)<-0

rm(data_bin_oil)

#subset states
data_oil<-data_oil[data_oil$`State/Province` %in% c("OK","TX","UT","PA"),]
jaccard_index_oil<-jaccard_index_oil[rownames(jaccard_index_oil) %in% data_oil$APINumber,
                                     colnames(jaccard_index_oil) %in% data_oil$APINumber]
# delay data
data_oil$delay<-ifelse(data_oil$`State/Province` %in% c("MT","TX","WY","MS"),30,
                        ifelse(data_oil$`State/Province` %in% c("MI","OK","PA","UT"),60,
                               ifelse(data_oil$`State/Province` %in% c("NM"),45,
                                      ifelse(data_oil$`State/Province` %in% c("LA"),20,
                                             120))))

data_oil$DeclDate<-data_oil$JobStartDate + data_oil$delay
data_oil<-cbind(data_oil[,1:18], data_oil[,839:840], data_oil[,19:838])
oil_d<-data_oil[,c("APINumber", "Operator Ticker", "Operator Company Name","JobStartDate", "DeclDate")]
gc(1:1000)

### find infoset of each well
info_set_oil<-info_set(oil_d)
gc(1:1000)

rm(oil_d)
gc(1:1000) 

info_oil<-apply(info_set_oil,2, function(y){
  g<-as.data.frame(cbind(y, rownames(info_set_oil)))
  g[,1]<-as.numeric(g[,1])
  g<-g[g[,1]==2,]
  g<-tail(g,1000) #keep most recent wells
  g<-as.character(g[,2])
})

gc(1:1000)
info_oil<-rbind(info_oil, colnames(info_oil))
gc(1:1000) 

info_oil<-as.data.frame(info_oil)
info_oil<-info_oil %>% mutate_all(as.character)
gc(1:1000) 

### choose MinPts
median(colSums(ifelse(jaccard_index_oil<0.01,1,0)))
median(colSums(ifelse(jaccard_index_oil<0.05,1,0)))
median(colSums(ifelse(jaccard_index_oil<0.1,1,0)))

### choose epsilon
# all wells
dbscan::kNNdistplot(as.dist(jaccard_index_oil), k =  5)
abline(h = 0.2, lty = 2)

dbscan::kNNdistplot(as.dist(jaccard_index_oil), k =  20)
abline(h = 0.3, lty = 2)

dbscan::kNNdistplot(as.dist(jaccard_index_oil), k =  50)
abline(h = 0.4, lty = 2)

#excluding early wells
window<-1000
dbscan::kNNdistplot(as.dist(jaccard_index_oil[(window+1):nrow(jaccard_index_oil),
                                              (window+1):ncol(jaccard_index_oil)]), k =  9)
abline(h = 0.3, lty = 2)

# dbscan algorithm
out_oil<-output_f(info_oil, jaccard_index_oil, 5, 0.3)
gc(1:1000)
#
out_oil_1<-output_f(info_oil, jaccard_index_oil, 20, 0.4)
gc(1:1000)

#
out_oil_2<-output_f(info_oil, jaccard_index_oil, 50, 0.6)
gc(1:1000)

rm(jaccard_index_oil, info_set_oil)
setwd("F:/MK/Prod_US/FracFocusMaria/DBSCAN")
save.image("out_oil_infoset.RData")
rm(list = ls())
gc(1:1000)

#### GAS ####
options(scipen=999)
setwd("F:/MK/Prod_US/FracFocusMaria/Jaccard/PercentHFJob")
load("jaccard_gas.RData")
gc(1:1000)

jaccard_index_gas[is.na(jaccard_index_gas)]<-0
jaccard_index_gas[jaccard_index_gas>1]<-1
jaccard_index_gas<-1-jaccard_index_gas
diag(jaccard_index_gas)<-0

rm(data_bin_gas)
gc(1:1000)

#subset states
data_gas<-data_gas[data_gas$`State/Province` %in% c("OK","TX","UT","PA","WV"),]
jaccard_index_gas<-jaccard_index_gas[rownames(jaccard_index_gas) %in% data_gas$APINumber,
                                     colnames(jaccard_index_gas) %in% data_gas$APINumber]
# delay data
data_gas$delay<-ifelse(data_gas$`State/Province` %in% c("MT","TX","WY","MS"),30,
                       ifelse(data_gas$`State/Province` %in% c("MI","OK","PA","UT"),60,
                              ifelse(data_gas$`State/Province` %in% c("NM"),45,
                                     ifelse(data_gas$`State/Province` %in% c("LA"),20,
                                            120))))

data_gas$DeclDate<-data_gas$JobStartDate + data_gas$delay
data_gas<-cbind(data_gas[,1:18], data_gas[,679:680], data_gas[,19:678])
gas_d<-data_gas[,c("APINumber", "Operator Ticker", "Operator Company Name","JobStartDate", "DeclDate")]
gc(1:1000)

### find infoset of each well

info_set_gas<-info_set(gas_d)
gc(1:1000)

rm(gas_d)
gc(1:1000)

info_gas<-apply(info_set_gas,2, function(y){
  g<-as.data.frame(cbind(y, rownames(info_set_gas)))
  g[,1]<-as.numeric(g[,1])
  g<-g[g[,1]==2,]
  g<-tail(g,1000) #keep most recent wells
  g<-as.character(g[,2])
})

gc(1:1000) 
info_gas<-rbind(info_gas, colnames(info_gas))
gc(1:1000) 

info_gas<-as.data.frame(info_gas)
info_gas<-info_gas %>% mutate_all(as.character)
gc(1:1000) 

### choose MinPts
median(colSums(ifelse(jaccard_index_gas<0.01,1,0)))
median(colSums(ifelse(jaccard_index_gas<0.05,1,0)))
median(colSums(ifelse(jaccard_index_gas<0.1,1,0)))

### choose epsilon
# all wells
dbscan::kNNdistplot(as.dist(jaccard_index_gas), k =  8)
abline(h = 0.2, lty = 2)

dbscan::kNNdistplot(as.dist(jaccard_index_gas), k =  30)
abline(h = 0.3, lty = 2)

dbscan::kNNdistplot(as.dist(jaccard_index_gas), k =  69)
abline(h = 0.4, lty = 2)

#excluding early wells
window<-1000
dbscan::kNNdistplot(as.dist(jaccard_index_gas[(window+1):nrow(jaccard_index_gas),
                                              (window+1):ncol(jaccard_index_gas)]), k =  15)
abline(h = 0.3, lty = 2)

# dbscan algorithm

out_gas<-output_f(info_gas, jaccard_index_gas, 8, 0.3)
gc(1:1000)
# 
out_gas_1<-output_f(info_gas, jaccard_index_gas, 30, 0.4)
gc(1:1000)
#
out_gas_2<-output_f(info_gas, jaccard_index_gas, 69, 0.6)
gc(1:1000)

rm(jaccard_index_gas, info_set_gas)
#save.image("out_gas_infoset.RData")
rm(list = ls())
gc(1:1000)
