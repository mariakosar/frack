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

debt_f<-function(debt, data){
  
  tickers<-unique(data$'Operator Ticker')
  tickers<-tickers[tickers!=""]
  debt<-debt[debt$tic %in% tickers,c("tic","datafqtr","leverage", "tot_debt", "seqq")]
  debt<-as.data.frame(debt)
  
  debt$leverage<-ifelse(debt$leverage>1,1,debt$leverage)
  
  debt<-debt %>%
    group_by(datafqtr) %>% 
    mutate(highlev = ifelse(leverage>quantile(leverage,0.5,na.rm=T),1,0))   
  
  debt<-as.data.frame(debt)
  debt<-ddply(debt, "tic", na.locf)
  
  debt<-debt[order(debt$tic, debt$datafqtr),]
  debt<-debt[,c("tic", "datafqtr","leverage","highlev","tot_debt", "seqq")]
  
  ## add leverage quantiles as of 2013Q4
  lev_2013q4<-debt[debt$datafqtr==c("2013Q4"),]
  lev_2013q4<-as.data.frame(cbind(lev_2013q4[,1],lev_2013q4[,3],lev_2013q4[,4]))
  colnames(lev_2013q4)<-c("tic","leverage","lev_top")
  lev_2013q4$tic<-as.character(lev_2013q4$tic)
  lev_2013q4$lev_top<-as.numeric(as.character(lev_2013q4$lev_top))
  lev_2013q4$leverage<-as.numeric(as.character(lev_2013q4$leverage))
  
  return(lev_2013q4)
}

#book

book_f<-function(debt, data){
  
  tickers<-unique(data$'Operator Ticker')
  tickers<-tickers[tickers!=""]
  debt<-debt[debt$tic %in% tickers,c("tic","datafqtr","atq","ltq")]
  debt<-as.data.frame(debt)
  
  debt<-ddply(debt, "tic", na.locf)
  
  ## add leverage quantiles as of 2013Q4
  lev_2013q4<-debt[debt$datafqtr==c("2013Q4"),]
  lev_2013q4<-lev_2013q4[,c("tic","atq","ltq")]
  colnames(lev_2013q4)<-c("tic","assets_pre","liabilities_pre")
  lev_2013q4$tic<-as.character(lev_2013q4$tic)
  lev_2013q4$assets_pre<-as.numeric(as.character(lev_2013q4$assets_pre))
  lev_2013q4$liabilities_pre<-as.numeric(as.character(lev_2013q4$liabilities_pre))
  
  return(lev_2013q4)
}
# mcap 

mcap_f <-function(data, mcap){
  tickers<-unique(data$'Operator Ticker')
  tickers<-tickers[tickers!=""]
  
  mcap<-mcap[mcap$tic %in% tickers,c("tic","datafqtr","mkvaltq")]
  mcap<-as.data.frame(mcap)
  
  mcap$mkval<-ifelse(mcap$mkvaltq>quantile(mcap$mkvaltq, 0.5, na.rm = T),1,0)
  mcap<-mcap[mcap$datafqtr=="2013Q4",]
  
  mcap<-mcap[,c("tic", "mkval","mkvaltq")]
  colnames(mcap)<-c("tic","mkvaltq","mkval")
  
  return(mcap)
}
