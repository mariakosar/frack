library(ggplot2)
library(dplyr)
library(data.table)

load("F:/MK/Prod_US/FracFocusMaria/OPTICS/PercentHFJob/EXPERIMENTATION/TOTAL WELLS/ex_infoset.RData")
states <- map_data("state")
counties <- map_data("county")

####  all wells by experimentation indicator #### 

labs_oil_0 <- ex_oil[ex_oil$Outlier==0,c("API10", "lon", "lat")]
labs_oil_1 <- ex_oil[ex_oil$Outlier==1,c("API10", "lon", "lat")]

labs_gas_0 <- ex_gas[ex_gas$Outlier==0,c("API10", "lon", "lat")]
labs_gas_1 <- ex_gas[ex_gas$Outlier==1,c("API10", "lon", "lat")]

ggplot() + 
  geom_polygon(data = counties, aes(x = long, y = lat, fill = region, group = group),fill="white",color = "grey") +
  geom_polygon(data = states, aes(x = long, y = lat, fill = region, group = group),fill=NA,color = "black") +
  coord_fixed(1.3) +
  guides(fill=FALSE) +  # do this to leave off the color legend
  geom_point(data = labs_oil_0, aes(x = lon, y = lat,color = "myline1"), size=1)+
  geom_point(data = labs_oil_1, aes(x = lon, y = lat,color = "myline2"), size=1)+
  geom_point(data = labs_gas_0, aes(x = lon, y = lat, color = "myline3"), size=1)+
  geom_point(data = labs_gas_1, aes(x = lon, y = lat,color = "myline4"), size=1)+
  scale_colour_manual(name="All wells:",
                      values=c(myline1="blue", myline2="red", myline3="green", myline4="orange"),
                      breaks=c("myline1", "myline2", "myline3","myline4"),
                      labels = c("Oil, non-exp", "Oil, exp",
                                 "Gas, non-exp", "Gas, exp"))

#### oil wells  #### 

# all

labs_oil_0 <- ex_oil[ex_oil$Outlier==0,c("API10", "lon", "lat")]
labs_oil_1 <- ex_oil[ex_oil$Outlier==1,c("API10", "lon", "lat")]

ggplot() + 
  geom_polygon(data = counties, aes(x = long, y = lat, fill = region, group = group),fill="white",color = "grey") +
  geom_polygon(data = states, aes(x = long, y = lat, fill = region, group = group),fill=NA,color = "black") +
  coord_fixed(1.3) +
  guides(fill=FALSE) +  # do this to leave off the color legend
  geom_point(data = labs_oil_0, aes(x = lon, y = lat,color = "myline1"), size=1)+
  geom_point(data = labs_oil_1, aes(x = lon, y = lat,color = "myline2"), size=1)+
  scale_colour_manual(name="Oil wells:",
                      values=c(myline1="blue", myline2="red"),
                      breaks=c("myline1", "myline2"),
                      labels = c("Oil, non-exp", "Oil, exp"))

# by year

plot_oil_by_year<-function(data, year){
  
  labs_oil_0 <- data[data$Outlier==0&data$year==year,c("API10", "lon", "lat")]
  labs_oil_1 <- data[data$Outlier==1&data$year==year,c("API10", "lon", "lat")]
  
  return(ggplot() + 
    geom_polygon(data = counties, aes(x = long, y = lat, fill = region, group = group),fill="white",color = "grey") +
    geom_polygon(data = states, aes(x = long, y = lat, fill = region, group = group),fill=NA,color = "black") +
    coord_fixed(1.3) +
    guides(fill=FALSE) +  # do this to leave off the color legend
    geom_point(data = labs_oil_0, aes(x = lon, y = lat,color = "myline1"), size=1)+
    geom_point(data = labs_oil_1, aes(x = lon, y = lat,color = "myline2"), size=1)+
    scale_colour_manual(name="Oil wells:",
                        values=c(myline1="blue", myline2="red"),
                        breaks=c("myline1", "myline2"),
                        labels = c("Oil, non-exp", "Oil, exp")))
  
}

plot_oil_by_year(ex_oil, "2013")
plot_oil_by_year(ex_oil, "2014")
plot_oil_by_year(ex_oil, "2015")
plot_oil_by_year(ex_oil, "2016")
plot_oil_by_year(ex_oil, "2017")
plot_oil_by_year(ex_oil, "2018")
plot_oil_by_year(ex_oil, "2019")
plot_oil_by_year(ex_oil, "2020")


####  gas wells  #### 

# all

labs_gas_0 <- ex_gas[ex_gas$Outlier==0,c("API10", "lon", "lat")]
labs_gas_1 <- ex_gas[ex_gas$Outlier==1,c("API10", "lon", "lat")]

ggplot() + 
  geom_polygon(data = counties, aes(x = long, y = lat, fill = region, group = group),fill="white",color = "grey") +
  geom_polygon(data = states, aes(x = long, y = lat, fill = region, group = group),fill=NA,color = "black") +
  coord_fixed(1.3) +
  guides(fill=FALSE) +  # do this to leave off the color legend
  geom_point(data = labs_gas_0, aes(x = lon, y = lat, color = "myline3"), size=1)+
  geom_point(data = labs_gas_1, aes(x = lon, y = lat,color = "myline4"), size=1)+
  scale_colour_manual(name="All wells:",
                      values=c(myline3="green", myline4="orange"),
                      breaks=c("myline3","myline4"),
                      labels = c("Gas, non-exp", "Gas, exp"))

# by year


plot_gas_by_year<-function(data, year){
  
  labs_gas_0 <- data[data$Outlier==0&data$year==year,c("API10", "lon", "lat")]
  labs_gas_1 <- data[data$Outlier==1&data$year==year,c("API10", "lon", "lat")]
  
  return(ggplot() + 
    geom_polygon(data = counties, aes(x = long, y = lat, fill = region, group = group),fill="white",color = "grey") +
    geom_polygon(data = states, aes(x = long, y = lat, fill = region, group = group),fill=NA,color = "black") +
    coord_fixed(1.3) +
    guides(fill=FALSE) +  # do this to leave off the color legend
    geom_point(data = labs_gas_0, aes(x = lon, y = lat, color = "myline3"), size=1)+
    geom_point(data = labs_gas_1, aes(x = lon, y = lat,color = "myline4"), size=1)+
    scale_colour_manual(name="All wells:",
                        values=c(myline3="green", myline4="orange"),
                        breaks=c("myline3","myline4"),
                        labels = c("Gas, non-exp", "Gas, exp")))
}

plot_gas_by_year(ex_gas, "2013")
plot_gas_by_year(ex_gas, "2014")
plot_gas_by_year(ex_gas, "2015")
plot_gas_by_year(ex_gas, "2016")
plot_gas_by_year(ex_gas, "2017")
plot_gas_by_year(ex_gas, "2018")
plot_gas_by_year(ex_gas, "2019")
plot_gas_by_year(ex_gas, "2020")
