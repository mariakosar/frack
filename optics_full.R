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

setwd("F:/MK/Prod_US/FracFocusMaria/OPTICS/PercentHFJob")
rm(jaccard_index_oil)
save.image("optics_oil_full.RData")
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

setwd("F:/MK/Prod_US/FracFocusMaria/OPTICS/PercentHFJob")
rm(jaccard_index_gas)
save.image("optics_gas_full.RData")
rm(list = ls())

########################################## OPTICS model selection ####################################################################
#### OIL ####

n_xi<-c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)

setwd("F:/MK/Prod_US/FracFocusMaria/OPTICS/PercentHFJob")
load("optics_oil_full.RData")
gc(1:1000)

names(out_oil_res)<-data_oil$APINumber[501:(nrow(data_oil)-1)]
names(out_oil_res_1)<-data_oil$APINumber[501:(nrow(data_oil)-1)]
names(out_oil_res_2)<-data_oil$APINumber[501:(nrow(data_oil)-1)]

# eps = 0.3
check_f<-function(y){
  # % of zeros
  a<-median(as.numeric(lapply(y, function (x) sum(as.numeric(x)==0)))/10)
  # number of clusters
  b<-median(as.numeric(lapply(y, function(x) max(as.numeric(x)))))
  # max %points in cluster
  c<-median(as.numeric(lapply(y, function (x) max(as.data.frame(as.numeric(x)) %>% count(as.data.frame(as.numeric(x))[,1]))))/10)
  outp<-as.data.frame(cbind(a,b,c))
  colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
  return(outp)
}

out_oil_check<-lapply(out_oil_res, function (x) extractDBSCAN(x, 0.3)[["cluster"]])
oil_check<-check_f(out_oil_check)
rownames(oil_check)<-"dbscan_03"
gc(1:1000)

rm(out_oil_check)

check_f<-function(y){
  tryCatch({
    y<-as.numeric(as.character(y))
    # % of zeros
    a<-as.numeric(sum(y==0, na.rm = T))/10
    # number of clusters
    b<-as.numeric(max(y, na.rm = T))
    # max %points in cluster
    suppl<-function (x){
      max(as.data.frame(x) %>% count(as.data.frame(as.numeric(x))[,1]))
    }
    c<-as.numeric(suppl(y))/10
    outp<-as.data.frame(cbind(a,b,c))
    colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
    return(outp)
  },
  error = function(i) {
    outp<-data.frame("% of zeros"=NA, "number of clusters"=NA, "max %points in cluster"=NA)
    colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
    return(outp)
  })}

copies_of_r <- 6

cl <- makeCluster(copies_of_r)
clusterExport(cl, list("n_xi", "out_oil_res", "check_f"))
clusterEvalQ(cl, library(parallel))
clusterEvalQ(cl, library(dbscan))
clusterEvalQ(cl, library(stats))
clusterEvalQ(cl, library(dplyr))
gc(1:1000)

oil_check_i<-apply(as.data.frame(n_xi),1, function(y) do.call("rbind",parLapply(cl,out_oil_res,
                                                                                function(x) check_f(as.numeric(extractXi(x, xi=y)[["cluster"]])))))


stopCluster(cl)
rm(cl)
gc(1:1000)

# eps = 0.4
check_f<-function(y){
  # % of zeros
  a<-median(as.numeric(lapply(y, function (x) sum(as.numeric(x)==0)))/10)
  # number of clusters
  b<-median(as.numeric(lapply(y, function(x) max(as.numeric(x)))))
  # max %points in cluster
  c<-median(as.numeric(lapply(y, function (x) max(as.data.frame(as.numeric(x)) %>% count(as.data.frame(as.numeric(x))[,1]))))/10)
  outp<-as.data.frame(cbind(a,b,c))
  colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
  return(outp)
}

out_oil_check_1<-lapply(out_oil_res_1, function (x) extractDBSCAN(x, 0.4)[["cluster"]])
oil_check_1<-check_f(out_oil_check_1)
rownames(oil_check_1)<-"dbscan_04"
gc(1:1000)

rm(out_oil_check_1)


check_f<-function(y){
  tryCatch({
    y<-as.numeric(as.character(y))
    # % of zeros
    a<-as.numeric(sum(y==0, na.rm = T))/10
    # number of clusters
    b<-as.numeric(max(y, na.rm = T))
    # max %points in cluster
    suppl<-function (x){
      max(as.data.frame(x) %>% count(as.data.frame(as.numeric(x))[,1]))
    }
    c<-as.numeric(suppl(y))/10
    outp<-as.data.frame(cbind(a,b,c))
    colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
    return(outp)
  },
  error = function(i) {
    outp<-data.frame("% of zeros"=NA, "number of clusters"=NA, "max %points in cluster"=NA)
    colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
    return(outp)
  })}

copies_of_r <- 6

cl <- makeCluster(copies_of_r)
clusterExport(cl, list("n_xi", "out_oil_res_1", "check_f"))
clusterEvalQ(cl, library(parallel))
clusterEvalQ(cl, library(dbscan))
clusterEvalQ(cl, library(stats))
clusterEvalQ(cl, library(dplyr))
gc(1:1000)

oil_check_1_i<-apply(as.data.frame(n_xi),1, function(y) do.call("rbind",parLapply(cl,out_oil_res_1,
                                                                                  function(x) check_f(as.numeric(extractXi(x, xi=y)[["cluster"]])))))

stopCluster(cl)
rm(cl)
gc(1:1000)

# eps = 0.5
check_f<-function(y){
  # % of zeros
  a<-median(as.numeric(lapply(y, function (x) sum(as.numeric(x)==0)))/10)
  # number of clusters
  b<-median(as.numeric(lapply(y, function(x) max(as.numeric(x)))))
  # max %points in cluster
  c<-median(as.numeric(lapply(y, function (x) max(as.data.frame(as.numeric(x)) %>% count(as.data.frame(as.numeric(x))[,1]))))/10)
  outp<-as.data.frame(cbind(a,b,c))
  colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
  return(outp)
}

out_oil_check_2<-lapply(out_oil_res_2, function (x) extractDBSCAN(x, 0.5)[["cluster"]])
oil_check_2<-check_f(out_oil_check_2)
rownames(oil_check_2)<-"dbscan_05"
gc(1:1000)

rm(out_oil_check_2)

check_f<-function(y){
  tryCatch({
    y<-as.numeric(as.character(y))
    # % of zeros
    a<-as.numeric(sum(y==0, na.rm = T))/10
    # number of clusters
    b<-as.numeric(max(y, na.rm = T))
    # max %points in cluster
    suppl<-function (x){
      max(as.data.frame(x) %>% count(as.data.frame(as.numeric(x))[,1]))
    }
    c<-as.numeric(suppl(y))/10
    outp<-as.data.frame(cbind(a,b,c))
    colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
    return(outp)
  },
  error = function(i) {
    outp<-data.frame("% of zeros"=NA, "number of clusters"=NA, "max %points in cluster"=NA)
    colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
    return(outp)
  })}

copies_of_r <- 6

cl <- makeCluster(copies_of_r)
clusterExport(cl, list("n_xi", "out_oil_res_2", "check_f"))
clusterEvalQ(cl, library(parallel))
clusterEvalQ(cl, library(dbscan))
clusterEvalQ(cl, library(stats))
clusterEvalQ(cl, library(dplyr))
gc(1:1000)

oil_check_2_i<-apply(as.data.frame(n_xi),1, function(y) do.call("rbind",parLapply(cl,out_oil_res_2,
                                                                                  function(x) check_f(as.numeric(extractXi(x, xi=y)[["cluster"]])))))

stopCluster(cl)
rm(cl)
gc(1:1000)

### Criteria of Schubert et al. (2017):

optics_03<-data.frame("1"=rep(NA,10), "2"=rep(NA,10), "3"=rep(NA,10), "4"=rep(NA,10), "5"=rep(NA,10), "6"=rep(NA,10))
colnames(optics_03)<-c("median % of zeros", "median number of clusters", "median max %points in cluster", "max % of zeros","min number of clusters","max %points in cluster")
rownames(optics_03)<-c("optics=0.3,xi=0.01", "optics=0.3,xi=0.02",
                       "optics=0.3,xi=0.03","optics=0.3,xi=0.04",
                       "optics=0.3,xi=0.05", "optics=0.3,xi=0.1",
                       "optics=0.3,xi=0.2", "optics=0.3,xi=0.3",
                       "optics=0.3,xi=0.4", "optics=0.3,xi=0.5")

optics_04<-data.frame("1"=rep(NA,10), "2"=rep(NA,10), "3"=rep(NA,10), "4"=rep(NA,10), "5"=rep(NA,10), "6"=rep(NA,10))
colnames(optics_04)<-c("median % of zeros", "median number of clusters", "median max %points in cluster", "max % of zeros","min number of clusters","max %points in cluster")
rownames(optics_04)<-c("optics=0.4,xi=0.01", "optics=0.4,xi=0.02",
                       "optics=0.4,xi=0.03","optics=0.4,xi=0.04",
                       "optics=0.4,xi=0.05", "optics=0.4,xi=0.1",
                       "optics=0.4,xi=0.2", "optics=0.4,xi=0.3",
                       "optics=0.4,xi=0.4", "optics=0.4,xi=0.5")

optics_05<-data.frame("1"=rep(NA,10), "2"=rep(NA,10), "3"=rep(NA,10), "4"=rep(NA,10), "5"=rep(NA,10), "6"=rep(NA,10))
colnames(optics_05)<-c("median % of zeros", "median number of clusters", "median max %points in cluster", "max % of zeros","min number of clusters","max %points in cluster")
rownames(optics_05)<-c("optics=0.5,xi=0.01", "optics=0.5,xi=0.02",
                       "optics=0.5,xi=0.03","optics=0.5,xi=0.04",
                       "optics=0.5,xi=0.05", "optics=0.5,xi=0.1",
                       "optics=0.5,xi=0.2", "optics=0.5,xi=0.3",
                       "optics=0.5,xi=0.4", "optics=0.5,xi=0.5")

# median % of outliers across all comparison sets, for each model

optics_03[,1]<-do.call("rbind",lapply(oil_check_i, function(x) median(x$`% of zeros`, na.rm = T)))
optics_04[,1]<-do.call("rbind",lapply(oil_check_1_i, function(x) median(x$`% of zeros`, na.rm = T)))
optics_05[,1]<-do.call("rbind",lapply(oil_check_2_i, function(x) median(x$`% of zeros`, na.rm = T)))

# median number of clusters across all comparison sets, for each model

optics_03[,2]<-do.call("rbind",lapply(oil_check_i, function(x) median(x$`number of clusters`, na.rm = T)))
optics_04[,2]<-do.call("rbind",lapply(oil_check_1_i, function(x) median(x$`number of clusters`, na.rm = T)))
optics_05[,2]<-do.call("rbind",lapply(oil_check_2_i, function(x) median(x$`number of clusters`, na.rm = T)))

# median % of points that fall in one cluster across all comparison sets, for each model

optics_03[,3]<-do.call("rbind",lapply(oil_check_i, function(x) median(x$`max %points in cluster`, na.rm = T)))
optics_04[,3]<-do.call("rbind",lapply(oil_check_1_i, function(x) median(x$`max %points in cluster`, na.rm = T)))
optics_05[,3]<-do.call("rbind",lapply(oil_check_2_i, function(x) median(x$`max %points in cluster`, na.rm = T)))

# max % of outliers across all comparison sets, for each model

optics_03[,4]<-do.call("rbind",lapply(oil_check_i, function(x) max(x$`% of zeros`, na.rm = T)))
optics_04[,4]<-do.call("rbind",lapply(oil_check_1_i, function(x) max(x$`% of zeros`, na.rm = T)))
optics_05[,4]<-do.call("rbind",lapply(oil_check_2_i, function(x) max(x$`% of zeros`, na.rm = T)))

# min number of clusters across all comparison sets, for each model

optics_03[,5]<-do.call("rbind",lapply(oil_check_i, function(x) min(x$`number of clusters`, na.rm = T)))
optics_04[,5]<-do.call("rbind",lapply(oil_check_1_i, function(x) min(x$`number of clusters`, na.rm = T)))
optics_05[,5]<-do.call("rbind",lapply(oil_check_2_i, function(x) min(x$`number of clusters`, na.rm = T)))

# max % of points that fall in one cluster across all comparison sets, for each model

optics_03[,6]<-do.call("rbind",lapply(oil_check_i, function(x) max(x$`max %points in cluster`, na.rm = T)))
optics_04[,6]<-do.call("rbind",lapply(oil_check_1_i, function(x) max(x$`max %points in cluster`, na.rm = T)))
optics_05[,6]<-do.call("rbind",lapply(oil_check_2_i, function(x) max(x$`max %points in cluster`, na.rm = T)))

### select clustering model by Schubert criteria

optics_03<-optics_03[optics_03$`median % of zeros`<20&optics_03$`median max %points in cluster`<50&optics_03$`median number of clusters`>1,]
optics_04<-optics_04[optics_04$`median % of zeros`<20&optics_04$`median max %points in cluster`<50&optics_04$`median number of clusters`>1,]
optics_05<-optics_05[optics_05$`median % of zeros`<20&optics_05$`median max %points in cluster`<50&optics_05$`median number of clusters`>1,]

optics_03<-optics_03[optics_03$`max % of zeros`<20&optics_03$`max %points in cluster`<50&optics_03$`min number of clusters`>1,]
optics_04<-optics_04[optics_04$`max % of zeros`<20&optics_04$`max %points in cluster`<50&optics_04$`min number of clusters`>1,]
optics_05<-optics_05[optics_05$`max % of zeros`<20&optics_05$`max %points in cluster`<50&optics_05$`min number of clusters`>1,]

options_oil<-rbind(optics_03, optics_04, optics_05)

rm(optics_03, optics_04, optics_05, oil_check, oil_check_1, oil_check_2, copies_of_r)


save.image("optics_oil_full.RData")
rm(list = ls())
gc(1:1000)

#### GAS ####

n_xi<-c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)

setwd("F:/MK/Prod_US/FracFocusMaria/OPTICS/PercentHFJob")
load("optics_gas_full.RData")
gc(1:1000)

names(out_gas_res)<-data_oil$APINumber[1001:(nrow(data_gas)-1)]
names(out_gas_res_1)<-data_oil$APINumber[1001:(nrow(data_gas)-1)]
names(out_gas_res_2)<-data_oil$APINumber[1001:(nrow(data_gas)-1)]

# eps = 0.3
check_f<-function(y){
  # % of zeros
  a<-median(as.numeric(lapply(y, function (x) sum(as.numeric(x)==0)))/10)
  # number of clusters
  b<-median(as.numeric(lapply(y, function(x) max(as.numeric(x)))))
  # max %points in cluster
  c<-median(as.numeric(lapply(y, function (x) max(as.data.frame(as.numeric(x)) %>% count(as.data.frame(as.numeric(x))[,1]))))/10)
  outp<-as.data.frame(cbind(a,b,c))
  colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
  return(outp)
}

out_gas_check<-lapply(out_gas_res, function (x) extractDBSCAN(x, 0.3)[["cluster"]])
gas_check<-check_f(out_gas_check)
rownames(gas_check)<-"dbscan_03"
gc(1:1000)

rm(out_gas_check)

check_f<-function(y){
  tryCatch({
    y<-as.numeric(as.character(y))
    # % of zeros
    a<-as.numeric(sum(y==0, na.rm = T))/10
    # number of clusters
    b<-as.numeric(max(y, na.rm = T))
    # max %points in cluster
    suppl<-function (x){
      max(as.data.frame(x) %>% count(as.data.frame(as.numeric(x))[,1]))
    }
    c<-as.numeric(suppl(y))/10
    outp<-as.data.frame(cbind(a,b,c))
    colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
    return(outp)
  },
  error = function(i) {
    outp<-data.frame("% of zeros"=NA, "number of clusters"=NA, "max %points in cluster"=NA)
    colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
    return(outp)
  })}

copies_of_r <- 4

cl <- makeCluster(copies_of_r)
clusterExport(cl, list("n_xi", "out_gas_res", "check_f"))
clusterEvalQ(cl, library(parallel))
clusterEvalQ(cl, library(dbscan))
clusterEvalQ(cl, library(stats))
clusterEvalQ(cl, library(dplyr))
gc(1:1000)

gas_check_i<-apply(as.data.frame(n_xi),1, function(y) do.call("rbind",parLapply(cl,out_gas_res,
                                                                                function(x) check_f(as.numeric(extractXi(x, xi=y)[["cluster"]])))))


stopCluster(cl)
rm(cl)
gc(1:1000)

# eps = 0.4
check_f<-function(y){
  # % of zeros
  a<-median(as.numeric(lapply(y, function (x) sum(as.numeric(x)==0)))/10)
  # number of clusters
  b<-median(as.numeric(lapply(y, function(x) max(as.numeric(x)))))
  # max %points in cluster
  c<-median(as.numeric(lapply(y, function (x) max(as.data.frame(as.numeric(x)) %>% count(as.data.frame(as.numeric(x))[,1]))))/10)
  outp<-as.data.frame(cbind(a,b,c))
  colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
  return(outp)
}

out_gas_check_1<-lapply(out_gas_res_1, function (x) extractDBSCAN(x, 0.4)[["cluster"]])
gas_check_1<-check_f(out_gas_check_1)
rownames(gas_check_1)<-"dbscan_04"
gc(1:1000)

rm(out_gas_check_1)


check_f<-function(y){
  tryCatch({
    y<-as.numeric(as.character(y))
    # % of zeros
    a<-as.numeric(sum(y==0, na.rm = T))/10
    # number of clusters
    b<-as.numeric(max(y, na.rm = T))
    # max %points in cluster
    suppl<-function (x){
      max(as.data.frame(x) %>% count(as.data.frame(as.numeric(x))[,1]))
    }
    c<-as.numeric(suppl(y))/10
    outp<-as.data.frame(cbind(a,b,c))
    colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
    return(outp)
  },
  error = function(i) {
    outp<-data.frame("% of zeros"=NA, "number of clusters"=NA, "max %points in cluster"=NA)
    colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
    return(outp)
  })}

copies_of_r <- 4

cl <- makeCluster(copies_of_r)
clusterExport(cl, list("n_xi", "out_gas_res_1", "check_f"))
clusterEvalQ(cl, library(parallel))
clusterEvalQ(cl, library(dbscan))
clusterEvalQ(cl, library(stats))
clusterEvalQ(cl, library(dplyr))
gc(1:1000)

gas_check_1_i<-apply(as.data.frame(n_xi),1, function(y) do.call("rbind",parLapply(cl,out_gas_res_1,
                                                                                  function(x) check_f(as.numeric(extractXi(x, xi=y)[["cluster"]])))))

stopCluster(cl)
rm(cl)
gc(1:1000)

# eps = 0.5
check_f<-function(y){
  # % of zeros
  a<-median(as.numeric(lapply(y, function (x) sum(as.numeric(x)==0)))/10)
  # number of clusters
  b<-median(as.numeric(lapply(y, function(x) max(as.numeric(x)))))
  # max %points in cluster
  c<-median(as.numeric(lapply(y, function (x) max(as.data.frame(as.numeric(x)) %>% count(as.data.frame(as.numeric(x))[,1]))))/10)
  outp<-as.data.frame(cbind(a,b,c))
  colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
  return(outp)
}

out_gas_check_2<-lapply(out_gas_res_2, function (x) extractDBSCAN(x, 0.5)[["cluster"]])
gas_check_2<-check_f(out_gas_check_2)
rownames(gas_check_2)<-"dbscan_05"
gc(1:1000)

rm(out_gas_check_2)

check_f<-function(y){
  tryCatch({
    y<-as.numeric(as.character(y))
    # % of zeros
    a<-as.numeric(sum(y==0, na.rm = T))/10
    # number of clusters
    b<-as.numeric(max(y, na.rm = T))
    # max %points in cluster
    suppl<-function (x){
      max(as.data.frame(x) %>% count(as.data.frame(as.numeric(x))[,1]))
    }
    c<-as.numeric(suppl(y))/10
    outp<-as.data.frame(cbind(a,b,c))
    colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
    return(outp)
  },
  error = function(i) {
    outp<-data.frame("% of zeros"=NA, "number of clusters"=NA, "max %points in cluster"=NA)
    colnames(outp)<-c("% of zeros", "number of clusters", "max %points in cluster")
    return(outp)
  })}

copies_of_r <- 4

cl <- makeCluster(copies_of_r)
clusterExport(cl, list("n_xi", "out_gas_res_2", "check_f"))
clusterEvalQ(cl, library(parallel))
clusterEvalQ(cl, library(dbscan))
clusterEvalQ(cl, library(stats))
clusterEvalQ(cl, library(dplyr))
gc(1:1000)

gas_check_2_i<-apply(as.data.frame(n_xi),1, function(y) do.call("rbind",parLapply(cl,out_gas_res_2,
                                                                                  function(x) check_f(as.numeric(extractXi(x, xi=y)[["cluster"]])))))

stopCluster(cl)
rm(cl)
gc(1:1000)

### Criteria of Schubert et al. (2017):

optics_03<-data.frame("1"=rep(NA,10), "2"=rep(NA,10), "3"=rep(NA,10), "4"=rep(NA,10), "5"=rep(NA,10), "6"=rep(NA,10))
colnames(optics_03)<-c("median % of zeros", "median number of clusters", "median max %points in cluster", "max % of zeros","min number of clusters","max %points in cluster")
rownames(optics_03)<-c("optics=0.3,xi=0.01", "optics=0.3,xi=0.02",
                       "optics=0.3,xi=0.03","optics=0.3,xi=0.04",
                       "optics=0.3,xi=0.05", "optics=0.3,xi=0.1",
                       "optics=0.3,xi=0.2", "optics=0.3,xi=0.3",
                       "optics=0.3,xi=0.4", "optics=0.3,xi=0.5")

optics_04<-data.frame("1"=rep(NA,10), "2"=rep(NA,10), "3"=rep(NA,10), "4"=rep(NA,10), "5"=rep(NA,10), "6"=rep(NA,10))
colnames(optics_04)<-c("median % of zeros", "median number of clusters", "median max %points in cluster", "max % of zeros","min number of clusters","max %points in cluster")
rownames(optics_04)<-c("optics=0.4,xi=0.01", "optics=0.4,xi=0.02",
                       "optics=0.4,xi=0.03","optics=0.4,xi=0.04",
                       "optics=0.4,xi=0.05", "optics=0.4,xi=0.1",
                       "optics=0.4,xi=0.2", "optics=0.4,xi=0.3",
                       "optics=0.4,xi=0.4", "optics=0.4,xi=0.5")

optics_05<-data.frame("1"=rep(NA,10), "2"=rep(NA,10), "3"=rep(NA,10), "4"=rep(NA,10), "5"=rep(NA,10), "6"=rep(NA,10))
colnames(optics_05)<-c("median % of zeros", "median number of clusters", "median max %points in cluster", "max % of zeros","min number of clusters","max %points in cluster")
rownames(optics_05)<-c("optics=0.5,xi=0.01", "optics=0.5,xi=0.02",
                       "optics=0.5,xi=0.03","optics=0.5,xi=0.04",
                       "optics=0.5,xi=0.05", "optics=0.5,xi=0.1",
                       "optics=0.5,xi=0.2", "optics=0.5,xi=0.3",
                       "optics=0.5,xi=0.4", "optics=0.5,xi=0.5")

# median % of outliers across all comparison sets, for each model

optics_03[,1]<-do.call("rbind",lapply(gas_check_i, function(x) median(x$`% of zeros`, na.rm = T)))
optics_04[,1]<-do.call("rbind",lapply(gas_check_1_i, function(x) median(x$`% of zeros`, na.rm = T)))
optics_05[,1]<-do.call("rbind",lapply(gas_check_2_i, function(x) median(x$`% of zeros`, na.rm = T)))

# median number of clusters across all comparison sets, for each model

optics_03[,2]<-do.call("rbind",lapply(gas_check_i, function(x) median(x$`number of clusters`, na.rm = T)))
optics_04[,2]<-do.call("rbind",lapply(gas_check_1_i, function(x) median(x$`number of clusters`, na.rm = T)))
optics_05[,2]<-do.call("rbind",lapply(gas_check_2_i, function(x) median(x$`number of clusters`, na.rm = T)))

# median % of points that fall in one cluster across all comparison sets, for each model

optics_03[,3]<-do.call("rbind",lapply(gas_check_i, function(x) median(x$`max %points in cluster`, na.rm = T)))
optics_04[,3]<-do.call("rbind",lapply(gas_check_1_i, function(x) median(x$`max %points in cluster`, na.rm = T)))
optics_05[,3]<-do.call("rbind",lapply(gas_check_2_i, function(x) median(x$`max %points in cluster`, na.rm = T)))

# max % of outliers across all comparison sets, for each model

optics_03[,4]<-do.call("rbind",lapply(gas_check_i, function(x) max(x$`% of zeros`, na.rm = T)))
optics_04[,4]<-do.call("rbind",lapply(gas_check_1_i, function(x) max(x$`% of zeros`, na.rm = T)))
optics_05[,4]<-do.call("rbind",lapply(gas_check_2_i, function(x) max(x$`% of zeros`, na.rm = T)))

# min number of clusters across all comparison sets, for each model

optics_03[,5]<-do.call("rbind",lapply(gas_check_i, function(x) min(x$`number of clusters`, na.rm = T)))
optics_04[,5]<-do.call("rbind",lapply(gas_check_1_i, function(x) min(x$`number of clusters`, na.rm = T)))
optics_05[,5]<-do.call("rbind",lapply(gas_check_2_i, function(x) min(x$`number of clusters`, na.rm = T)))

# max % of points that fall in one cluster across all comparison sets, for each model

optics_03[,6]<-do.call("rbind",lapply(gas_check_i, function(x) max(x$`max %points in cluster`, na.rm = T)))
optics_04[,6]<-do.call("rbind",lapply(gas_check_1_i, function(x) max(x$`max %points in cluster`, na.rm = T)))
optics_05[,6]<-do.call("rbind",lapply(gas_check_2_i, function(x) max(x$`max %points in cluster`, na.rm = T)))

### select clustering model by Schubert criteria

optics_03<-optics_03[optics_03$`median % of zeros`<20&optics_03$`median max %points in cluster`<50&optics_03$`median number of clusters`>1,]
optics_04<-optics_04[optics_04$`median % of zeros`<20&optics_04$`median max %points in cluster`<50&optics_04$`median number of clusters`>1,]
optics_05<-optics_05[optics_05$`median % of zeros`<20&optics_05$`median max %points in cluster`<50&optics_05$`median number of clusters`>1,]

optics_03<-optics_03[optics_03$`max % of zeros`<20&optics_03$`max %points in cluster`<50&optics_03$`min number of clusters`>1,]
optics_04<-optics_04[optics_04$`max % of zeros`<20&optics_04$`max %points in cluster`<50&optics_04$`min number of clusters`>1,]
optics_05<-optics_05[optics_05$`max % of zeros`<20&optics_05$`max %points in cluster`<50&optics_05$`min number of clusters`>1,]

options_gas<-rbind(optics_03, optics_04, optics_05)

rm(optics_03, optics_04, optics_05, gas_check, gas_check_1, gas_check_2, copies_of_r)


save.image("optics_gas_full.RData")
rm(list = ls())
gc(1:1000)


######################################### Extract Outliers of the best model ###################################################################
#### OIL ####
load("F:/MK/Prod_US/FracFocusMaria/OPTICS/PercentHFJob/optics_oil_full.RData")
options(scipen=999)
gc(1:1000)

out_oil<-lapply(out_oil_res, function (x) extractXi(x, xi=0.01)[["cluster"]][501])
out_oil<-do.call("rbind",out_oil)
rownames(out_oil)<-names(out_oil_res)
sum(out_oil==0)


rm(out_oil_res,out_oil_res_1, out_oil_res_2,
   oil_check_i, oil_check_1_i, oil_check_2_i, options_oil, i, n_xi, check_f, optics_func,
   output_f)

save.image("out_oil_full.RData")
rm(list = ls())
gc(1:1000)

#### GAS ####
load("F:/MK/Prod_US/FracFocusMaria/OPTICS/PercentHFJob/optics_gas_full.RData")
options(scipen=999)
gc(1:1000)

out_gas<-lapply(out_gas_res, function (x) extractXi(x, xi=0.01)[["cluster"]][501])
out_gas<-do.call("rbind",out_gas)
rownames(out_gas)<-data_gas$APINumber[502:nrow(data_gas)]
sum(out_gas==0)


rm(out_gas_res,out_gas_res_1, out_gas_res_2,
   gas_check_i, gas_check_1_i, gas_check_2_i, options_gas, i, n_xi, check_f, optics_func,
   output_f)

save.image("out_gas_full.RData")
rm(list = ls())
gc(1:1000)
