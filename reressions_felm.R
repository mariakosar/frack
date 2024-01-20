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

#### 1. Regressions main ####

setwd("F:/MK/Prod_US/FracFocusMaria/OPTICS/PercentHFJob/EXPERIMENTATION")
load("pts_infoset.RData")
#load("pts_full.RData")

pts_gas$year<-ifelse(substr(pts_gas$month, 6,7) %in% c("01","02","03","04",
              "05", "06", "07", "08", "09")&pts_gas$year=="2014","2013",pts_gas$year)

pts_oil$year<-ifelse(substr(pts_oil$month, 6,7) %in% c("01","02","03","04",
              "05", "06", "07", "08", "09")&pts_oil$year=="2014","2013",pts_oil$year)

pts_gas<-pts_gas[pts_gas$year!="2020",]
pts_oil<-pts_oil[pts_oil$year!="2020",]

# sensitivity analysis
# pts_gas<-pts_gas[!(pts_gas$tic %in% c("EQT","RRC","AR","SWN","CHK")),]
# pts_gas<-pts_gas[!(pts_gas$tic %in% c("EQT")),]
# pts_gas<-pts_gas[!(pts_gas$tic %in% c("RRC")),]
# pts_gas<-pts_gas[!(pts_gas$tic %in% c("AR")),]
# pts_gas<-pts_gas[!(pts_gas$tic %in% c("SWN")),]
# pts_gas<-pts_gas[!(pts_gas$tic %in% c("CHK")),]
# 
# pts_oil<-pts_oil[!(pts_oil$tic %in% c("EOG","MRO","DVN","SD","XOM")),]
# pts_oil<-pts_oil[!(pts_oil$tic %in% c("EOG")),]
# pts_oil<-pts_oil[!(pts_oil$tic %in% c("MRO")),]
# pts_oil<-pts_oil[!(pts_oil$tic %in% c("DVN")),]
# pts_oil<-pts_oil[!(pts_oil$tic %in% c("SD")),]
# pts_oil<-pts_oil[!(pts_oil$tic %in% c("XOM")),]


#### gas ####

# median

summary(felm(Outlier ~ factor(treat_m)*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_gas[!(pts_gas$year %in% c("2014","2015")),]))

g_1<-felm(Outlier ~ factor(treat_m)*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_gas[!(pts_gas$year %in% c("2014","2015")),])

summary(felm(Outlier ~ factor(treat_m)*factor(year)
             + lev_top*factor(year) + mkvaltq*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_gas[!(pts_gas$year %in% c("2014","2015")),]))

g_2<-felm(Outlier ~ factor(treat_m)*factor(year)
          + lev_top*factor(year) + mkvaltq*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_gas[!(pts_gas$year %in% c("2014","2015")),])

summary(felm(Outlier ~ factor(treat_m)*factor(year)
             + lev_top*factor(year) + mkvaltq*factor(year) + rev_q*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_gas[!(pts_gas$year %in% c("2014","2015")),]))

g_3<-felm(Outlier ~ factor(treat_m)*factor(year)
          + lev_top*factor(year) + mkvaltq*factor(year) + rev_q*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_gas[!(pts_gas$year %in% c("2014","2015")),])

# top
summary(felm(Outlier ~ factor(treat_d)*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_gas[!(pts_gas$year %in% c("2014","2015")),]))

g_4<-felm(Outlier ~ factor(treat_d)*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_gas[!(pts_gas$year %in% c("2014","2015")),])

summary(felm(Outlier ~ factor(treat_d)*factor(year)
             + lev_top*factor(year) + mkvaltq*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_gas[!(pts_gas$year %in% c("2014","2015")),]))

g_5<-felm(Outlier ~ factor(treat_d)*factor(year)
          + lev_top*factor(year) + mkvaltq*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_gas[!(pts_gas$year %in% c("2014","2015")),])

summary(felm(Outlier ~ factor(treat_d)*factor(year)
             + lev_top*factor(year) + mkvaltq*factor(year) + rev_q*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_gas[!(pts_gas$year %in% c("2014","2015")),]))

g_6<-felm(Outlier ~ factor(treat_d)*factor(year)
          + lev_top*factor(year) + mkvaltq*factor(year) + rev_q*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_gas[!(pts_gas$year %in% c("2014","2015")),])

#### oil ####

# median

summary(felm(Outlier ~ factor(treat_m)*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_oil[!(pts_oil$year %in% c("2014","2015")),]))

o_1<-felm(Outlier ~ factor(treat_m)*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_oil[!(pts_oil$year %in% c("2014","2015")),])

summary(felm(Outlier ~ factor(treat_m)*factor(year)
             + lev_top*factor(year) + mkvaltq*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_oil[!(pts_oil$year %in% c("2014","2015")),]))

o_2<-felm(Outlier ~ factor(treat_m)*factor(year)
          + lev_top*factor(year) + mkvaltq*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_oil[!(pts_oil$year %in% c("2014","2015")),])

summary(felm(Outlier ~ factor(treat_m)*factor(year)
             + lev_top*factor(year) + mkvaltq*factor(year) + rev_q*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_oil[!(pts_oil$year %in% c("2014","2015")),]))

o_3<-felm(Outlier ~ factor(treat_m)*factor(year)
          + lev_top*factor(year) + mkvaltq*factor(year) + rev_q*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_oil[!(pts_oil$year %in% c("2014","2015")),])

# top
summary(felm(Outlier ~ factor(treat_d)*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_oil[!(pts_oil$year %in% c("2014","2015")),]))

o_4<-felm(Outlier ~ factor(treat_d)*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_oil[!(pts_oil$year %in% c("2014","2015")),])

summary(felm(Outlier ~ factor(treat_d)*factor(year)
             + lev_top*factor(year) + mkvaltq*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_oil[!(pts_oil$year %in% c("2014","2015")),]))

o_5<-felm(Outlier ~ factor(treat_d)*factor(year)
          + lev_top*factor(year) + mkvaltq*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_oil[!(pts_oil$year %in% c("2014","2015")),])

summary(felm(Outlier ~ factor(treat_d)*factor(year)
             + lev_top*factor(year) + mkvaltq*factor(year) + rev_q*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_oil[!(pts_oil$year %in% c("2014","2015")),]))

o_6<-felm(Outlier ~ factor(treat_d)*factor(year)
          + lev_top*factor(year) + mkvaltq*factor(year) + rev_q*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_oil[!(pts_oil$year %in% c("2014","2015")),])

#### output to Latex table ####

stargazer(o_1, o_2, o_3, g_1, g_2, g_3,summary=F, title = "Impact of revation of reserves on experimentation, OIL and GAS firms", align=TRUE)

#### 2. Regressions interference ####
#### gas ####

summary(felm(Outlier ~ interference*factor(treat_m)*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_gas[!(pts_gas$year %in% c("2014","2015")),]))

g_1<-felm(Outlier ~ interference*factor(treat_m)*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_gas[!(pts_gas$year %in% c("2014","2015")),])

summary(felm(Outlier ~ interference*factor(treat_m)*factor(year)
             + lev_top*factor(year) + mkvaltq*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_gas[!(pts_gas$year %in% c("2014","2015")),]))

g_2<-felm(Outlier ~ interference*factor(treat_m)*factor(year)
          + lev_top*factor(year) + mkvaltq*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_gas[!(pts_gas$year %in% c("2014","2015")),])

summary(felm(Outlier ~ interference*factor(treat_m)*factor(year)
             + lev_top*factor(year) + mkvaltq*factor(year) + rev_q*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_gas[!(pts_gas$year %in% c("2014","2015")),]))

g_3<-felm(Outlier ~ interference*factor(treat_m)*factor(year)
          + lev_top*factor(year) + mkvaltq*factor(year) + rev_q*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_gas[!(pts_gas$year %in% c("2014","2015")),])


#### oil ####

summary(felm(Outlier ~ interference*factor(treat_m)*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_oil[!(pts_oil$year %in% c("2014","2015")),]))

o_1<-felm(Outlier ~ interference*factor(treat_m)*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_oil[!(pts_oil$year %in% c("2014","2015")),])

summary(felm(Outlier ~ interference*factor(treat_m)*factor(year)
             + lev_top*factor(year) + mkvaltq*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_oil[!(pts_oil$year %in% c("2014","2015")),]))

o_2<-felm(Outlier ~ interference*factor(treat_m)*factor(year)
          + lev_top*factor(year) + mkvaltq*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_oil[!(pts_oil$year %in% c("2014","2015")),])

summary(felm(Outlier ~ interference*factor(treat_m)*factor(year)
             + lev_top*factor(year) + mkvaltq*factor(year) + rev_q*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_oil[!(pts_oil$year %in% c("2014","2015")),]))

o_3<-felm(Outlier ~ interference*factor(treat_m)*factor(year)
          + lev_top*factor(year) + mkvaltq*factor(year) + rev_q*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_oil[!(pts_oil$year %in% c("2014","2015")),])


#### output to Latex table ####

stargazer(o_1, o_2, o_3, g_1, g_2, g_3,summary=F, title = "Impact of revation of reserves on experimentation, OIL and GAS firms, interference wells", align=TRUE)


#### 3. Regressions experim ####
#### gas ####

summary(felm(Outlier ~ experim*factor(treat_m)*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_gas[!(pts_gas$year %in% c("2014","2015")),]))

g_1<-felm(Outlier ~ experim*factor(treat_m)*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_gas[!(pts_gas$year %in% c("2014","2015")),])

summary(felm(Outlier ~ experim*factor(treat_m)*factor(year)
             + lev_top*factor(year) + mkvaltq*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_gas[!(pts_gas$year %in% c("2014","2015")),]))

g_2<-felm(Outlier ~ experim*factor(treat_m)*factor(year)
          + lev_top*factor(year) + mkvaltq*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_gas[!(pts_gas$year %in% c("2014","2015")),])

summary(felm(Outlier ~ experim*factor(treat_m)*factor(year)
             + lev_top*factor(year) + mkvaltq*factor(year) + rev_q*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_gas[!(pts_gas$year %in% c("2014","2015")),]))

g_3<-felm(Outlier ~ experim*factor(treat_m)*factor(year)
          + lev_top*factor(year) + mkvaltq*factor(year) + rev_q*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_gas[!(pts_gas$year %in% c("2014","2015")),])


#### oil ####

summary(felm(Outlier ~ experim*factor(treat_m)*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_oil[!(pts_oil$year %in% c("2014","2015")),]))

o_1<-felm(Outlier ~ experim*factor(treat_m)*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_oil[!(pts_oil$year %in% c("2014","2015")),])

summary(felm(Outlier ~ experim*factor(treat_m)*factor(year)
             + lev_top*factor(year) + mkvaltq*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_oil[!(pts_oil$year %in% c("2014","2015")),]))

o_2<-felm(Outlier ~ experim*factor(treat_m)*factor(year)
          + lev_top*factor(year) + mkvaltq*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_oil[!(pts_oil$year %in% c("2014","2015")),])

summary(felm(Outlier ~ experim*factor(treat_m)*factor(year)
             + lev_top*factor(year) + mkvaltq*factor(year) + rev_q*factor(year)|tic_grid_id + datafqtr|0|tic,
             data = pts_oil[!(pts_oil$year %in% c("2014","2015")),]))

o_3<-felm(Outlier ~ experim*factor(treat_m)*factor(year)
          + lev_top*factor(year) + mkvaltq*factor(year) + rev_q*factor(year)|tic_grid_id + datafqtr|0|tic,
          data = pts_oil[!(pts_oil$year %in% c("2014","2015")),])


#### output to Latex table ####

stargazer(o_1, o_2, o_3, g_1, g_2, g_3,summary=F, title = "Impact of revation of reserves on experimentation, OIL and GAS firms, interference wells", align=TRUE)


