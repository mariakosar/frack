# load("jaccard_ffr.RData") 
# permits
setwd("F:/MK/Prod_US/FracFocus/Jaccard")
data$OperatorName<-toupper(data$OperatorName)

perm<-fread("PermitsDI/PermitsTable.csv", header = T) #permits data from 2012 till end 2019, O,G,O&G,H,NewDrill

perm$`Approved Date`<-as.Date(perm$`Approved Date`, format="%Y-%m-%d")
perm$`Expired Date`<-as.Date(perm$`Expired Date`, format="%Y-%m-%d")

perm<-select(perm, -c("Operator Alias (Legacy)", "Amended Date", "Operator Company Name"))
colnames(perm)[1]<-c("API10")
perm$API10<-as.character(perm$API10)
perm$API10<-ifelse(nchar(perm$API10)>10, substr(perm$API10, 4, 13), perm$API10)

perm<-perm[order(perm$API10, perm$`Approved Date`),]
perm<-perm[!duplicated(perm$API10),]

## merge chemicals with permits
data<-merge(data, perm, by="API10")
rm(perm)

data<-cbind(data[,c("UploadKey","API10","JobStartDate","JobEndDate","APINumber","OperatorName", "State/Province","County/Parish",
                    "Operator (Reported)", "Operator Ticker","Permit Type","Well Type","Drill Type","Approved Date",
                    "Surface Hole Longitude (WGS84)","Bottom Hole Latitude (WGS84)","Surface Hole Latitude (WGS84)","Bottom Hole Longitude (WGS84)",
                    "Expired Date")], 
            data[,!(colnames(data) %in% 
                      c("UploadKey","API10","JobStartDate","JobEndDate","APINumber","OperatorName","State/Province","County/Parish",
                        "Operator (Reported)", "Operator Ticker","Permit Type","Well Type","Drill Type","Approved Date",
                        "Surface Hole Longitude (WGS84)","Bottom Hole Latitude (WGS84)","Surface Hole Latitude (WGS84)","Bottom Hole Longitude (WGS84)",
                        "Expired Date"))])

data<-select(data,-c("UploadKey"))
#drill type=NewDrill
data<-data[data$'Permit Type'=="NewDrill",]

# don't do it for now
# data_n$keep<-ifelse(data_n$JobStartDate<=data_n$`Expired Date`,1,0)
# data_n<-data_n[data_n$keep==1,]
# 
