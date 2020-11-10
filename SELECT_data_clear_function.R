library(sas7bdat)
library(dplyr)

select_alldre<-read.sas7bdat('/Users/yiyaochen/Documents/Trails Data/SELECT Risk Calculator Data/select_alldre.sas7bdat')
select_allpsa<-read.sas7bdat('/Users/yiyaochen/Documents/Trails Data/SELECT Risk Calculator Data/select_allpsa.sas7bdat')
select_allbaseline<-read.sas7bdat('/Users/yiyaochen/Documents/Trails Data/SELECT Risk Calculator Data/select_baseline.sas7bdat')
select_allendpoint<-read.sas7bdat('/Users/yiyaochen/Documents/Trails Data/SELECT Risk Calculator Data/select_endpoint.sas7bdat')
select_allbiopsy <- read.sas7bdat('/Users/yiyaochen/Documents/Trails Data/SELECT negative biopsy data/select_negativebx.sas7bdat')

# exclude person has (prostate) cancer before baseline, baseline now has 34772 different ID. Man with ID 118992
# had cancer at enrolled
select_baseline <- subset(select_allbaseline, select_allbaseline$PROSTDX != 'Y')
select_baseline <- subset(select_allbaseline, select_allbaseline$CADX != 'Y') 
select_baseline <- select_baseline[which(select_baseline$ID != 118992),]
summary(levels(factor(select_baseline$ID)))
# 34772

# Pick out the PSA for pre-eligible participants
select_psa <- select_allpsa[which(select_allpsa$ID %in% select_baseline$ID),]
select_psa <- select_psa[complete.cases(select_psa$PSA),]
summary(levels(factor(select_allpsa$ID))) #35533
summary(levels(factor(select_psa$ID))) #34747

# Pick out the DRE for pre-eligible participants
select_dre <- select_alldre[which(select_alldre$ID %in% select_baseline$ID),]
select_dre <- select_dre[complete.cases(select_dre$dre),]
summary(levels(factor(select_alldre$ID))) #35506
summary(levels(factor(select_dre$ID))) # 34745

# Pick out the biopsy record for pre-eligible participants
select_biopsy <- select_allbiopsy[which(select_allbiopsy$ID %in% select_psa$ID),]
summary(levels(factor(select_allbiopsy$ID))) #3689
summary(levels(factor(select_biopsy$ID))) # 3610

# Pick out the cancer record for pre-eligible participants
select_endpoint <- select_allendpoint[which(select_allendpoint$ID %in% select_psa$ID), ]
summary(levels(factor(select_endpoint$ID))) #34747
table(select_endpoint$pca)
# 2489 with cancer

# ID merge cancer and biopsy since all cancer cases are biopsy-ascertained. 
# But for no PCa, only part of it are biopsy ascertained. 
ID_cancer = data.frame(ID = select_endpoint[select_endpoint$pca == 1, ]$ID) # 2489
ID_biopsy = data.frame(ID = unique(select_biopsy$ID)) # 3610
ID_biopsy = merge(ID_cancer, ID_biopsy, by = "ID", all = T) # 4975
ID_nocancer <- data.frame(ID = ID_biopsy[-which(ID_biopsy$ID %in% ID_cancer$ID),])
summary(levels(factor(ID_biopsy$ID))) # 4975 has biosy with 2489 with cancer


########################################################
# Functions use to filtering
########################################################
# dt1_pca: ID, pcadt
# dt2_psa: ID, PSA, psadt
# dt3_dre: ID, DRE, dredt
dt_cancer_latest_biopsy <- function(dt1_pca, dt2_psa, dt3_dre){
  tmp <- dt1_pca
  
  tmp$lastest_biopsy_date <- dt1_pca$pcadt
  
  # PSA: within 2 years before latest biopsy
  dt <- merge(dt2_psa, tmp[, c('ID','lastest_biopsy_date')], by = 'ID', all.x = TRUE)
  dt$date_gap = ifelse(dt$lastest_biopsy_date - dt$psadt < 0 | dt$lastest_biopsy_date - dt$psadt > 730,
                       NA, dt$lastest_biopsy_date - dt$psadt)
  dt = dt[complete.cases(dt$date_gap),]
  dt = dt %>% 
    group_by(ID) %>% 
    slice(which.min(date_gap))
  
  # DRE: within 2 years before latest biopsy
  dt2 <- merge(dt3_dre, tmp[, c('ID','lastest_biopsy_date')],  by = 'ID', all.x = TRUE)
  dt2$date_gap = ifelse(dt2$lastest_biopsy_date - dt2$dredt < 0 | dt2$lastest_biopsy_date - dt2$dredt > 730,
                       NA, dt2$lastest_biopsy_date - dt2$dredt)
  dt2 = dt2[complete.cases(dt2$date_gap),]
  dt2 = dt2 %>% 
    group_by(ID) %>% 
    slice(which.min(date_gap))
  
  
  # Merge dt, dt2, keep all ID in dt, remove the rest
  dt3 <- merge(dt[, c("ID","PSA", "psadt", "lastest_biopsy_date")],
               dt2[, c("ID" , "dre" ,"dredt")], by = 'ID', all.x = TRUE, all.y = FALSE)
  
  final <- merge(dt3, dt1_pca, by = 'ID', all.x=TRUE, all.y = FALSE)
  
  # Final contain: ID, PSA, psadt, dre, dredt, latest_biopsy_date, ......
  return(final)
}
dt1_pca <- select_endpoint[select_endpoint$ID %in% ID_cancer$ID, ]
dt2_psa <- select_psa[select_psa$ID %in% ID_cancer$ID, ]
dt3_dre <- select_alldre[select_alldre$ID %in% ID_cancer$ID, ]

cancer_select <- dt_cancer_latest_biopsy(dt1_pca, dt2_psa, dt3_dre)

dt2_psa <- select_psa[-which(select_psa$ID %in% ID_biopsy$ID), ]
dt3_dre <- select_alldre[-which(select_alldre$ID %in% ID_biopsy$ID), ]

dt_missing_latest_biopsy <- function(dt2_psa, dt3_dre){
  
  # PSA: latest PSA
  dt <- dt2_psa
  dt = dt %>% 
    group_by(ID) %>% 
    slice(which.max(psadt))
  
  # DRE: latest DRE
  dt2 <- dt3_dre
  dt2 = dt2 %>% 
    group_by(ID) %>% 
    slice(which.max(dredt))
  
  # Merge dt, dt2, keep all ID in dt, remove the rest
  dt3 <- merge(dt[, c("ID","PSA", "psadt")],
               dt2[, c("ID" , "dre" ,"dredt")], by = 'ID', all.x = TRUE, all.y = FALSE)
  
  final <- dt3
  
  # Final contain: ID, PSA, psadt, dre, dredt
  return(final)
}

dt2_psa <- select_psa[-which(select_psa$ID %in% ID_biopsy$ID), ]
dt3_dre <- select_alldre[-which(select_alldre$ID %in% ID_biopsy$ID), ]
nobiopp_select <- dt_missing_latest_biopsy(dt2_psa, dt3_dre)



dt_nocancer_latest_biopsy <- function(dt2_psa, dt3_dre, dt4_biopsy){
  
  # PSA
  eligible_psa <- data.frame(ID = rep(NA, nrow(ID_nocancer)),
                             lastest_biopsy_date =  rep(NA, nrow(ID_nocancer)),
                             PSA = rep(NA, nrow(ID_nocancer)),
                             psadt = rep(NA, nrow(ID_nocancer)))
  
  for (i in 1:nrow(ID_nocancer)) {
    tmp <- dt4_biopsy[which(dt4_biopsy$ID == ID_nocancer$ID[i]),]
    tmp_psa <- dt2_psa[which(dt2_psa$ID == ID_nocancer$ID[i]),]
    visitdt <- sort(tmp$visitdt, decreasing = TRUE)
    
    for(j in 1:length(visitdt)){
        lastest_biopsy_date <- visitdt[j]
        tmp_psa$gap <- ifelse(lastest_biopsy_date - tmp_psa$psadt >= 0 & lastest_biopsy_date - tmp_psa$psadt <= 730,
                              lastest_biopsy_date - tmp_psa$psadt, NA)
        na <- colSums(!is.na(tmp_psa))
        
        if(tail(na, n = 1) > 0){
          pick <-  tmp_psa[which(is.na(tmp_psa$gap) == FALSE & tmp_psa$gap == min(tmp_psa$gap, na.rm = TRUE)),c('ID','PSA','psadt')]
          eligible_psa$ID[i] <- ID_nocancer$ID[i]
          eligible_psa$lastest_biopsy_date[i] <- lastest_biopsy_date
          eligible_psa$PSA[i] <- pick$PSA
          eligible_psa$psadt[i] <- pick$psadt
          break
        } else if (tail(na, n = 1) == 0 & j == length(visitdt)) {
          pick <-  data.frame(ID = ID_nocancer$ID[i], PSA = NA, psadt = NA)
          eligible_psa$ID[i] <- ID_nocancer$ID[i]
          eligible_psa$lastest_biopsy_date[i] <- NA
          eligible_psa$PSA[i] <- NA
          eligible_psa$psadt[i] <- NA
          break
        }
        
    } # for loop for subject
  } # for loop for dataset
  
  eligible_psa <- eligible_psa[complete.cases(eligible_psa$lastest_biopsy_date),]
  
  # DRE: pick out the DRE based on the eligible_psa people, no more filtering
  dre<- merge(eligible_psa, dt3_dre, by = 'ID', all.x = TRUE, all.y = FALSE )
  dre$date_gap = ifelse(dre$lastest_biopsy_date - dre$dredt < 0 | 
                          dre$lastest_biopsy_date - dre$dredt > 730, NA, 
                        dre$lastest_biopsy_date- dre$dredt)
  dre = dre[complete.cases(dre$date_gap),] 
  dre = dre%>% 
    group_by(ID) %>% 
    slice(which.min(date_gap))
  
  
  final <- merge(eligible_psa, dre[, c('ID',"dre" ,"dredt")], by = 'ID', all.x = TRUE)
  
  
  return(final)
}


dt2_psa <- select_psa[select_psa$ID %in% ID_nocancer$ID, ]
dt3_dre <- select_alldre[select_alldre$ID %in% ID_nocancer$ID, ]
dt4_biopsy <- select_biopsy[which(select_biopsy$ID %in% ID_nocancer$ID),]

# dt <- as.data.frame(table(dt4_biopsy$ID))
# max(dt$Freq) = 5
# View(dt[which(dt$Freq == 5),])
nocancer_select <- dt_nocancer_latest_biopsy(dt2_psa,dt3_dre, dt4_biopsy)

###############################################
# Combined with baseline 
###############################################
nrow(cancer_select)
nrow(nocancer_select)
nocancer_select$pca <- rep(0, nrow(nocancer_select))
nrow(nobiopp_select)
nobiopp_select$pca <- rep(NA, nrow(nobiopp_select))
nobiopp_select$lastest_biopsy_date <- nobiopp_select$psadt
# 34541

psa_dre <- rbind(can<-cancer_select[,c("ID" ,"pca","PSA" , "psadt", "lastest_biopsy_date","dre","dredt")],
                 nocancer_select[, c("ID" ,"pca", "PSA" , "psadt", "lastest_biopsy_date","dre","dredt")],
                 nobiopp_select[, c("ID" ,"pca","PSA" , "psadt", "lastest_biopsy_date","dre","dredt")])
master_select<- merge(select_allbaseline, psa_dre, by = 'ID',all.x = FALSE, all.y = TRUE)
master_select$biopsyflag <- ifelse(is.na(master_select$pca) == FALSE, 1, 0)

master_select$fhflag <- ifelse(is.na(master_select$fampca) == T, NA, ifelse(master_select$fampca == 0, 0, 1))
master_select$PriorBioflag <- ifelse(master_select$PROSBX == 'Y',1, 0)
master_select$PSA4flag <- ifelse(master_select$PSA > 4, 1, 0)
master_select$raceflag <- ifelse(master_select$racesel == 'UNKNOWN', NA,ifelse(master_select$racesel == 'BLACK',1, 0))

# Calculate the age at chosen psa master_select
master_select$lamaster_selectAge <- round((master_select$lasmaster_select_biopsy_date - master_select$REGDT)/365, 2) + master_select$age
master_select$Age55flag <- ifelse(master_select$lamaster_selectAge < 55, 0, 1)
# Remove White with Age < 55
master_select <- master_select[- which(master_select$lamaster_selectAge <55 & master_select$raceflag == 0), ]
master_select$log2PSA <- ifelse(master_select$PSA == 0, log2(master_select$PSA + 0.01), log2(master_select$PSA))


dt_select <- subset(test, select = c("ID","latestAge", "PSA","log2PSA", "dre", "biopsyflag", "raceflag", "fhflag","PriorBioflag","PSA4flag",  'pca'))
colnames(dt_select) = c('ID','Age','PSA', "log2PSA", 'DRE_binary','biopsyflag','Race','FH','PB',"PSA4flag",'PCa')
SELECT <- dt_select[complete.cases(dt_select[, c(1:9)]),] # 33257
write.csv(SELECT, '/Users/yiyaochen/Documents/Causal Inference/code/master_SELECT_latestPSA_33257_Oct2020.csv')



