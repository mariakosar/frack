library(data.table)
library(dplyr)
library(bit64)
library(gdata)
library(dplyr)
options(scipen=999)
gc(1:100) 

data<-data[order(data$JobStartDate),]

## 1. subset data to horizontal wells only

data_oil<-data[data$'Well Type'=="OIL" & data$'Drill Type'=="H",]
data_gas<-data[data$'Well Type'=="GAS" & data$'Drill Type'=="H",]
data_oil_gas<-data[data$'Well Type'=="OIL & GAS" & data$'Drill Type'=="H",]

# order data
data_oil<-data_oil[order(data_oil$JobStartDate),]
data_gas<-data_gas[order(data_gas$JobStartDate),]
data_oil_gas<-data_oil_gas[order(data_oil_gas$JobStartDate),]

# create numeric matrixes of percentages of chemical j (column) in well i (row)

d<-data_oil[,19:ncol(data_oil),]
d<-d[,colSums(d, na.rm = T) != 0]
data_oil<-cbind(data_oil[,1:18],d)

d<-data_gas[,19:ncol(data_gas),]
d<-d[,colSums(d, na.rm = T) != 0]
data_gas<-cbind(data_gas[,1:18],d)

d<-data_oil_gas[,19:ncol(data_oil_gas),]
d<-d[,colSums(d, na.rm = T) != 0]
data_oil_gas<-cbind(data_oil_gas[,1:18],d)

rm(d)

## 2. construct binary matrices with value of 1 if chemicals present in well i' are also present in well i

data_bin_oil<-cbind(data_oil[,1:18], as.data.frame(ifelse(data_oil[,19:ncol(data_oil)]>0,1,0)))
data_bin_gas<-cbind(data_gas[,1:18], as.data.frame(ifelse(data_gas[,19:ncol(data_gas)]>0,1,0)))
data_bin_oil_gas<-cbind(data_oil_gas[,1:18], as.data.frame(ifelse(data_oil_gas[,19:ncol(data_oil_gas)]>0,1,0)))

rm(data)

save.image("jaccard_segm.RData")

## 3. compute jaccard index
#gas
load("jaccard_segm.RData")
jaccard_f<-function(z,y,h){
  jaccard_matrix_A<-t(apply(z[,h:ncol(z)], 1, function (x) rowSums(sweep(y[,h:ncol(z)], MARGIN=2, x, `*`))))
  jaccard_matrix_B<-t(jaccard_matrix_A)
  
  jaccard_matrix_AB<-jaccard_matrix_A*jaccard_matrix_B
  jaccard_matrix_DEN<-jaccard_matrix_A+jaccard_matrix_B - jaccard_matrix_AB
  jaccard_index<-(jaccard_matrix_AB)/jaccard_matrix_DEN
  
  dimnames(jaccard_index)<-list(z$APINumber,z$APINumber)
  jaccard_index[is.na(jaccard_index)]<-0
  jaccard_index[jaccard_index>1]<-1
  jaccard_index[is.infinite(jaccard_index)]<-0
  jaccard_index[!is.numeric(jaccard_index)]<-0
  return(jaccard_index)
}

rm(data_bin_oil, data_bin_oil_gas, data_oil, data_oil_gas)
gc(1:1000)
jaccard_index_gas<-jaccard_f(data_gas, data_bin_gas, h=19)
save.image("jaccard_gas.RData")
rm(list = ls())
gc(1:1000)

# oil
load("jaccard_segm.RData")
jaccard_f<-function(z,y,h){
  jaccard_matrix_A<-t(apply(z[,h:ncol(z)], 1, function (x) rowSums(sweep(y[,h:ncol(z)], MARGIN=2, x, `*`))))
  jaccard_matrix_B<-t(jaccard_matrix_A)
  
  jaccard_matrix_AB<-jaccard_matrix_A*jaccard_matrix_B
  jaccard_matrix_DEN<-jaccard_matrix_A+jaccard_matrix_B - jaccard_matrix_AB
  jaccard_index<-(jaccard_matrix_AB)/jaccard_matrix_DEN
  
  dimnames(jaccard_index)<-list(z$APINumber,z$APINumber)
  jaccard_index[is.na(jaccard_index)]<-0
  jaccard_index[jaccard_index>1]<-1
  jaccard_index[is.infinite(jaccard_index)]<-0
  jaccard_index[!is.numeric(jaccard_index)]<-0
  return(jaccard_index)
}

rm(data_bin_gas, data_bin_oil_gas, data_gas, data_oil_gas)
gc(1:1000)
jaccard_index_oil<-jaccard_f(data_oil, data_bin_oil, h=19)
save.image("jaccard_oil.RData")
rm(list = ls())
gc(1:1000)

#ngls
load("jaccard_segm.RData")
jaccard_f<-function(z,y,h){
  jaccard_matrix_A<-t(apply(z[,h:ncol(z)], 1, function (x) rowSums(sweep(y[,h:ncol(z)], MARGIN=2, x, `*`))))
  jaccard_matrix_B<-t(jaccard_matrix_A)
  
  jaccard_matrix_AB<-jaccard_matrix_A*jaccard_matrix_B
  jaccard_matrix_DEN<-jaccard_matrix_A+jaccard_matrix_B - jaccard_matrix_AB
  jaccard_index<-(jaccard_matrix_AB)/jaccard_matrix_DEN
  
  dimnames(jaccard_index)<-list(z$APINumber,z$APINumber)
  jaccard_index[is.na(jaccard_index)]<-0
  jaccard_index[jaccard_index>1]<-1
  jaccard_index[is.infinite(jaccard_index)]<-0
  jaccard_index[!is.numeric(jaccard_index)]<-0
  return(jaccard_index)
}

rm(data_bin_oil, data_bin_gas, data_oil, data_gas)
gc(1:1000)
jaccard_index_oil_gas<-jaccard_f(data_oil_gas, data_bin_oil_gas, h=19)
save.image("jaccard_ngl.RData")
rm(list = ls())
gc(1:1000)
