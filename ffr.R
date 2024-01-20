##### Load packages #####
library(data.table)
library(plyr)
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
library(lfe)
library(zoo)
library(tidyr)
library(bit64)
library(webchem)
library(stringr)

options(scipen=999)

#read data from excel files downloaded from FracFocus
setwd("F:/MK/Prod_US/FracFocus/ffr/fracfocuscsv")
file_list<-list.files(getwd())
file_list<-c(file_list[1], file_list[12], file_list[15:21], file_list[2:11], file_list[13:14])

data<-vector("list",length(file_list))
names(data)<-file_list

for (i in file_list){
  data[[i]]<-fread(i,header = T)
}

names(data)<-1:21

# convert integer64 to character
for(i in 1:21){
  data[[i]]$APINumber<-as.character(data[[i]]$APINumber)
}

data<-do.call("rbind",data)
gc(1:1000)

###valid CAS number
#remove spaces from CAS number

data$CASLabel<-gsub(" ", "", data$CASNumber)

#replace invalid CAS with missing

missing_entries<-c("TradeSecret","TRADESECRET","Tradesecret","TradeSeceret","TradeSecret,disc.","tradesecret","TRADESECRETS",
                   "TradeName","TradeSecert","tracesecret","tradeSecret","TRADESECERET","TradeSecrte","TradeSecrer",
                   "tradeseccret","TradSecret","ts","TS","Proprietary","PROPRIETARY","3rdPartyProprietar",
                   "proprietary","Propietary","Prop","ProprietaryBlend","Propreitary","PROPRITARY",
                   "PRORIETARY","Proprietaryl","Proptietary","PROP","Proprietarty","Prop.","Proprietory",
                   "Propreitory","PROPRIETARY0.10","proprietry","propietary","Proprietatry","Propriatary",
                   "propriatary","Propritary","Proprietar","PROPRIERTARY","Priprietary","7732-18-5/propr",
                   "Proprietart","Proprieatary","Proprietery","PRIOPRIETARY","PROPRIERARY","Propriety","prop",
                   "proprietarty","7732-18-5proprietary","propriety","Properitary","Confidential",
                   "ConfidentialBusines","Confidnetial","Confidentail","confidential","CONFIDENTIAL",
                   "CONFIDENTIALBUSINES","BusinessConfidental","ConfidentialInfo","COnfidential",
                   "Confidental","Conf","ConfBusInfo","Confid.Bus.Info","Confidenial","Confidentia1",
                   "ConfidentialInf","Confinential","Confindential","CBI","N/A","NotListed","NONE","none",
                   "NotAvailable","N.A.","NOTPROVIDED","None","n/a","NOTASSIGNED","\"N/A\"","-","na",
                   "Notapplicable.", "Notavailable.","N/D","NotEstablished","Notavailable","NOTAVAILABLE",
                   "NULL","notlisted","NA","NDA","Notlisted","N\\A","non","NotApplicable","N/a","Undisclosed",
                   "UNK","unknown","Unknown","undisclosed","unk","Unavailable","","ListedBelow")


data$CASLabel<-ifelse(data$CASLabel %in% missing_entries,"missing",data$CASLabel)
gc(1:1000)

#mathemitically verify CAS number - check-digit verification
# is.cas_f: source - webcham package

is.cas_f <-  function(x, verbose = TRUE) {
   # x <- '64-17-5'
   if (length(x) > 1) {
      stop('Cannot handle multiple input strings.')
   }
   
   # cas must not have any alpha characters
   if(grepl(pattern = "[[:alpha:]]", x = x)){
      if(isTRUE(verbose)){message("String contains alpha characters")}
      return(FALSE)
   }
   
   # cas must have two hyphens
   nsep <- str_count(x, '-')
   if (nsep != 2) {
      if (isTRUE(verbose))
         message('Less than 2 hyphens in string.')
      return(FALSE)
   }
   
   # first part 2 to 7 digits
   fi <- gsub('^(.*)-(.*)-(.*)$', '\\1', x)
   if (nchar(fi) > 7 | nchar(fi) < 2) {
      if (isTRUE(verbose))
         message('First part with more than 7 digits!')
      return(FALSE)
   }
   
   # second part must be two digits
   se <- gsub('^(.*)-(.*)-(.*)$', '\\2', x)
   if (nchar(se) != 2) {
      if (isTRUE(verbose))
         message('Second part has not two digits!')
      return(FALSE)
   }
   
   # third part (checksum) must be 1 digit
   th <- gsub('^(.*)-(.*)-(.*)$', '\\3', x)
   if (nchar(th) != 1) {
      if (isTRUE(verbose))
         message('Third part has not 1 digit!')
      return(FALSE)
   }
   
   # check checksum
   di <-  as.numeric(strsplit(gsub('^(.*)-(.*)-(.*)$', '\\1\\2', x), split = '')[[1]])
   checksum <- sum(rev(seq_along(di)) * di)
   
   if ((is.na(checksum)==T)||(is.na(checksum)==F&checksum %% 10 != as.numeric(th))) {
      if (isTRUE(verbose))
         message('Checksum is not correct! ', checksum %% 10, ' vs. ', th)
      return(FALSE)
   }
   
   return(TRUE)
}

data$CASValid<-apply(as.data.frame(data$CASLabel),1,function (x) is.cas_f(as.character(x), verbose=F)) # only NA errors this function
data$CASLabelV<-ifelse(data$CASValid==TRUE,data$CASLabel,NA)

#remove leading zeros

data$CASNumberClean<- gsub("(^|[^0-9])0+", "\\1", data$CASLabelV, perl = TRUE)

#form datasets with only chemicals data in percentages

data_chem_percent<-data[,c("UploadKey","CASNumberClean","PercentHFJob")]
data_chem_percent<-data_chem_percent[order(data_chem_percent$UploadKey, data_chem_percent$CASNumberClean),]

### sum similar CAS by UploadKey, keep unique uploadkey and CAS

data_chem_percent_clean<-tapply(data_chem_percent$PercentHFJob, paste0(data_chem_percent$UploadKey,"_",data_chem_percent$CASNumberClean), function(x) sum(x, na.rm = T))
data_chem_percent_clean<-as.data.frame(data_chem_percent_clean)
colnames(data_chem_percent_clean)<-"PercentHFJob"
data_chem_percent_clean<-cbind(data_chem_percent_clean, unique(cbind(data_chem_percent$UploadKey, data_chem_percent$CASNumber)))
data_chem_percent_clean$`1`<-as.character(data_chem_percent_clean$`1`)
data_chem_percent_clean$`2`<-as.character(data_chem_percent_clean$`2`)
colnames(data_chem_percent_clean)<-c("PercentHFJob", "UploadKey", "CASNumberClean")
rownames(data_chem_percent_clean)<-1:nrow(data_chem_percent_clean)

## from long to wide format

data_wide_percent<-spread(data_chem_percent_clean, CASNumberClean, PercentHFJob)
data_wide_percent[is.na(data_wide_percent)]<-0

gc(1:1000)
                                
###### create registry table by keeping unique variables

registry<-data[!duplicated(data$UploadKey), c("UploadKey","JobStartDate","JobEndDate","APINumber","OperatorName")]
colnames(registry)<-c("UploadKey","JobStartDate","JobEndDate","APINumber","OperatorName")

####  format Date ####
                                
registry$JobStartDate<-substr(registry$JobStartDate,1,10)
registry$JobEndDate<-substr(registry$JobEndDate,1,10)

registry$JobStartDate<-ifelse(substr(registry$JobStartDate,9,9)==" ", paste0(0,substr(registry$JobStartDate,1,2),0,substr(registry$JobStartDate,3,8)),registry$JobStartDate)
registry$JobStartDate<-ifelse(substr(registry$JobStartDate,10,10)==" ", ifelse(substr(registry$JobStartDate,2,2)=="/",paste0(0,registry$JobStartDate),
                                                                      paste0(substr(registry$JobStartDate,1,3),0,substr(registry$JobStartDate,
                                                                     4,9))), registry$JobStartDate)

registry$JobEndDate<-ifelse(substr(registry$JobEndDate,9,9)==" ", paste0(0,substr(registry$JobEndDate,1,2),0,substr(registry$JobEndDate,3,8)),registry$JobEndDate)
registry$JobEndDate<-ifelse(substr(registry$JobEndDate,10,10)==" ", ifelse(substr(registry$JobEndDate,2,2)=="/",paste0(0,registry$JobEndDate),
                                                                      paste0(substr(registry$JobEndDate,1,3),0,substr(registry$JobEndDate,
                                                                                                                        4,9))), registry$JobEndDate)
registry<-registry[!is.na(registry$JobStartDate),]
registry<-registry[!is.na(registry$JobEndDate),]

registry$JobStartDate<-paste0(substr(registry$JobStartDate,7,10),"-", substr(registry$JobStartDate,1,2),"-",substr(registry$JobStartDate,4,5))
registry$JobEndDate<-paste0(substr(registry$JobEndDate,7,10),"-", substr(registry$JobEndDate,1,2),"-",substr(registry$JobEndDate,4,5))

registry$JobStartDate<-as.Date(registry$JobStartDate, format = "%Y-%m-%d")
registry$JobEndDate<-as.Date(registry$JobEndDate, format = "%Y-%m-%d")

registry<-registry[!is.na(registry$JobStartDate),]
registry<-registry[!is.na(registry$JobEndDate),]

#### remove duplicated APINumber and JobStartDate ####

registry$APINumber<-ifelse(nchar(registry$APINumber)==13, paste0(0,"",registry$APINumber), 
                       ifelse(nchar(registry$APINumber)==10, paste0(registry$APINumber,"","0000"), registry$APINumber))

registry<-registry[order(registry$APINumber, registry$JobStartDate),]
registry<-registry[!is.na(registry$JobStartDate),]
registry<-registry[!is.na(registry$JobEndDate),]

registry<-registry[!duplicated(registry$APINumber,registry$JobStartDate),]

#### remove all entries with JobStartDate<31 May 2013 & >31 Dec 2020

registry<-registry[registry$JobStartDate >= "2013-05-31"&registry$JobStartDate <= "2020-12-31", ]

#### merge registry with wide chemicals data ####

data_wide_percent<-merge(data_wide_percent, registry ,by= "UploadKey")

#save.image("ffr.RData")

rm(data,data_chem_percent,data_chem_percent_clean,registry,file_list,i,missing_entries,is.cas_f)

#keep unique wells
data_wide_percent$API10<-as.character(substr( data_wide_percent$APINumber,1,10))
data_wide_percent$API10<-as.character( data_wide_percent$API10)
data_wide_percent<- data_wide_percent[order( data_wide_percent$API10,  data_wide_percent$JobStartDate),]
data_wide_percent<- data_wide_percent[!duplicated( data_wide_percent$API10),]

data_wide_percent<-cbind( data_wide_percent[,c("UploadKey","API10","JobStartDate","JobEndDate","APINumber","OperatorName","<NA>","--")], 
                          data_wide_percent[,!(colnames( data_wide_percent) %in% 
                                                 c("UploadKey","API10","JobStartDate","JobEndDate","APINumber","OperatorName","<NA>","--"))])

# exclude water, sand, and unknowns

data_wide_percent<-select( data_wide_percent, -c("7732-18-5",
                                                 "7631-86-9","14808-60-7","7440-21-3","<NA>","--",
                                                 "14808-60-7","1309-37-1","1344-28-1","13463-67-7",
                                                 "9016-87-9","7631-86-9","1305-78-8","1309-48-4",
                                                 "57029-46-6","1318-16-7","12136-45-7","64742-30-9",
                                                 "1341-49-7", "1313-59-3", "7664-39-3", "1317-71-1",
                                                 "1343-88-", "27176-87-", "9025-56-3", "108-95-2",
                                                 "308075-7-2", "66402-68-4", "9003-35-4", "64742-94-5",
                                                 "68476-25-5", "34590-94-8", "1302-93-8", "14464-46-1",
                                                 "50-28-2","14608-60-7", "57018-52-7"))


data_wide_percent$JobStartDate<-as.character(data_wide_percent$JobStartDate)
data_wide_percent$JobEndDate<-as.character(data_wide_percent$JobEndDate)

# deal with missing values and negative values
                                
for (i in 1:6){
  data_wide_percent[,i]<-ifelse(is.na( data_wide_percent[,i])==T,"missing", data_wide_percent[,i])
}

data_wide_percent[is.na( data_wide_percent)]<-0
data_wide_percent[data_wide_percent<0]<-0
data_wide_percent[data_wide_percent=="missing"]<-NA

# clean data
data_wide_percent$Total<-rowSums( data_wide_percent[,!(colnames( data_wide_percent) %in% 
                                                         c("UploadKey","API10","JobStartDate","JobEndDate","APINumber","OperatorName"))])

data_wide_percent<- data_wide_percent[data_wide_percent$Total>0,]

data_wide_percent_ratio<-cbind( data_wide_percent[,c("UploadKey","API10","JobStartDate","JobEndDate","APINumber","OperatorName")],
                                data_wide_percent[,!(colnames( data_wide_percent) %in% 
                                                       c("UploadKey","API10","JobStartDate","JobEndDate","APINumber","OperatorName"))]/100)

data_wide_percent_ratio<- data_wide_percent_ratio[data_wide_percent_ratio$Total<=1,]

data_wide_percent_ratio<-cbind( data_wide_percent_ratio[,c("UploadKey","API10","JobStartDate","JobEndDate","APINumber","OperatorName")],
                                data_wide_percent_ratio[,!(colnames( data_wide_percent_ratio) %in% 
                                                         c("UploadKey","API10","JobStartDate","JobEndDate","APINumber","OperatorName"))]/ data_wide_percent_ratio$Total)

data_wide_percent_ratio<-select( data_wide_percent_ratio,-c("Total"))
data<- data_wide_percent_ratio

rm( data_wide_percent_ratio,  data_wide_percent)
gc(1:1000)

data<-data[,colSums(data[,!(colnames(data) %in% 
                              c("UploadKey","API10","JobStartDate","JobEndDate","APINumber","OperatorName"))], na.rm = T) != 0]

data$JobStartDate<-as.Date(data$JobStartDate, "%Y-%m-%d")
data$JobEndDate<-as.Date(data$JobEndDate, "%Y-%m-%d")

data<-data[order(data$APINumber, data$JobStartDate),]

#save.image("jaccard_ffr.RData")

