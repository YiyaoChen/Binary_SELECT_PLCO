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

#------ cancer------------------------------------------------
# For cancer case, pick the latest biopsy prior to cancer 
# Pick the latest PSA and/or DRE recorded 2 years prior than this biopsy 

cancer_psa = select_psa[select_psa$ID %in% ID_cancer$ID, ]
summary(levels(factor(cancer_psa$ID))) # 2489
cancer_psa_pca = merge(cancer_psa, select_endpoint[select_endpoint$pca == 1, ], by = "ID", all.x = T)
summary(levels(factor(cancer_psa_pca$ID))) #2489

# For cancer cases, the PSA records are within 2 years before the cancer diagnosis. 
# Criteria: remove psa test records with:
#             - pcadt (cancer diagnosis date) < psadt (PSA tests date): all those psa test records after prostate cancer diagnosis
#             - pcadt - psadt > 2*365 = 730: all those records taken more than 2 years ago.
cancer_psa_pca$date_gap = ifelse(cancer_psa_pca$pcadt - cancer_psa_pca$psadt < 0 |
                                  cancer_psa_pca$pcadt - cancer_psa_pca$psadt > 730, 
                                NA, cancer_psa_pca$pcadt - cancer_psa_pca$psadt)

# Remove those with cancer but without eligible PSA subject to above criteria records.
cancer_psa_pca = cancer_psa_pca[complete.cases(cancer_psa_pca$date_gap),]
summary(levels(factor(cancer_psa_pca$ID))) # 2364

# IMPORTANT Step: Pick out the records for those with cancer the PSA record, prior and cloest to the cancer dignosis
# i.e. the record with the smallest 'date_gap' among all his eligible psa test records.
cancer_psa2 = cancer_psa_pca %>% 
  group_by(ID) %>% 
  slice(which.min(date_gap))

# as.Date(15210, "1960-01-01") # 2001 Aug. 22
# Merge the eligible psa records with the baseline characteristic
cancer_psa = subset(cancer_psa2, select = c("ID","PSA","psadt", "pcadt","pca"))
cancer_f = merge(cancer_psa, select_baseline, by = "ID" )
# Calculated the age when the cancer is diagnosis. 
cancer_f$latestAge = round((cancer_f$pcadt - cancer_f$REGDT)/365, 2)+cancer_f$age
summary(levels(factor(cancer_f$ID))) # 2364

# Among the above eligible participants with cancer, pick their DRE records
cancer_dre = select_alldre[select_alldre$ID %in% cancer_f$ID, ]
summary(levels(factor(cancer_dre$ID))) #2364
# For cancer cases, the LAETEST! (no matter how many, pick the latest dre records within 2 years prior to their cancer dignosis)
# DRE records are within 2 years before the cancer diagnosis.
##    merge to have the "pcadt" for everyone
cancer_dre = merge(cancer_psa, cancer_dre, by = "ID", all.y = T)
##    remove the records that taken after diagnosis and those taken more than 2 years ago
cancer_dre$date_gap = ifelse(cancer_dre$pcadt - cancer_dre$dredt < 0 |
                              cancer_dre$pcadt - cancer_dre$dredt > 730,
                            NA, cancer_dre$pcadt - cancer_dre$dredt)
cancer_dre = cancer_dre[complete.cases(cancer_dre$date_gap),] 
summary(levels(factor(cancer_dre$ID))) # 2088
cancer_dre2 = cancer_dre%>% 
  group_by(ID) %>% 
  slice(which.min(date_gap))

cancer_dre2 = cancer_dre2[, c("ID","dre")] 
# As in SELECT, the dre = 0 stands for 'no/missing' so cannot distinguish those with missing and normal dre.
cancer_f = merge(cancer_f, cancer_dre2, by = "ID", all.x = T)
table(cancer_f$dre, exclude = NULL)
# 276 without eligible DRE records according to out criteria
# 566 with abnormal DRE
# 1522 with 'no/missing' DRE


# Pickout the latest biopsy date, i.e. the date when the PCa is diagnosis. 
cancer_f$latestBiopp = cancer_f$psadt
cancer_f$biopsyflag = rep(1, nrow(cancer_f))


################################################################################
# Below is a nesting steps to filter out the ID with non-missing PSA values taken with 2 years prior
# then the respective biopsy those without biopsy-ascertained cancer
# Step 1: Pick the latest biopsy date and filter for those with eligible PSA values
# Step 2: Remove the ID that already chosen in Step 1, filter for those ID with eligible PSA values w.r.t the 2nd latest biopsy
# Step 3: Remove the ID that already chosen in Step 1 & 2, filter for those ID with eligible PSA values w.r.t the 3rd latest biopsy
#############################
# Pick out those without biopsy-acsertained cancer (biopsy = yes, pca = 0)
nocancer_psa_all = select_psa[which(select_psa$ID %in% select_endpoint[select_endpoint$pca == 0, ]$ID & 
                                  select_psa$ID %in% ID_biopsy$ID), ]
summary(levels(factor(nocancer_psa_all$ID))) # 2486
ID_nocancer = data.frame(ID = unique(nocancer_psa_all$ID))
# Biopsy records for all IF without cancer
nocancer_biopsy = select_biopsy[which(select_biopsy$ID %in% ID_nocancer$ID),]
summary(levels(factor(nocancer_biopsy$ID))) # 2486



# Step 1:
   #    Find the latest biopsy date for no cancer cases since they do not have the cancer diagnosis date. 
   #     Latest biopsy date for each non-cancer case
nocancer_biopsy_1st= nocancer_biopsy %>% 
  group_by(ID) %>% 
  slice(which.max(visitdt))

   #    For the latest biopsy, check and pick the latest PSA that are taken within 2 years prior than the latest biopsy date.
   #    Merge the latest biopsy date for each ID, different PSA records (all.x = T). 
nocancer_psa_biopsy_1st = merge(nocancer_psa_all, nocancer_biopsy_1st, by = 'ID', all.x = T) 

nocancer_psa_biopsy_1st$date_gap = ifelse(nocancer_psa_biopsy_1st$visitdt - nocancer_psa_biopsy_1st$psadt < 0 |
                                              nocancer_psa_biopsy_1st$visitdt - nocancer_psa_biopsy_1st$psadt > 730, NA,
                                            nocancer_psa_biopsy_1st$visitdt - nocancer_psa_biopsy_1st$psadt)
nocancer_psa_biopsy_1st = nocancer_psa_biopsy_1st[complete.cases(nocancer_psa_biopsy_1st$date_gap),]
summary(levels(factor(nocancer_psa_biopsy_1st$ID))) # 2387
nocancer_psa_biopsy_1st <- nocancer_psa_biopsy_1st%>% 
  group_by(ID) %>% 
  slice(which.min(date_gap))

#----------------------------------------------
nocancer_latest <- nocancer_psa_biopsy_1st # 2387
#----------------------------------------------
# Step 2:
   # Pick the second latest biopsy date for ID:
   #  - With biopsy-ascertained no cancer
   #  - Without eligble PSA test (taken within 2 years prior than latest biopsy date) for latest biopsy date.
   #-----------------------------------
   # Second latest biopsy
   #   - remeber, we have 2486 with at least one biopsies during the study)
   #   - We have 2387 have the eligible psa test records with respect to their latest PSA test.

         # Remove the ID that already pick the latest PSA records 2 years prior than the latest biopsy
         # from the unfiltered nocancer_biopsy 
summary(levels(factor(nocancer_biopsy$ID)))
nocancer_biopsy_2nd_start = nocancer_biopsy[-which(nocancer_biopsy$ID %in% nocancer_latest$ID),]
summary(levels(factor(nocancer_biopsy_2nd_start$ID))) # 99 = 2486 - 2387

         # Remove rest ID with only 1 biopsy record as the biopsies for those ID does not have attached eligible PSA.
         # If they (the IDs) have eligible PSA attached to their only biopsy test, they should be picked out in the previous step already.
         # But they are not, so, the only biopsy record they have did not have eligible PSA (within 2 years prior than the biopsy) with respect to it.
n_occur <- data.frame(table(nocancer_biopsy_2nd_start$ID))
nocancer_biopsy_2nd_start = nocancer_biopsy_2nd_start[nocancer_biopsy_2nd_start$ID %in% n_occur$Var1[n_occur$Freq > 1],]
summary(levels(factor(nocancer_biopsy_2nd_start$ID))) #25: 74 remove without eligible PSA 

         # pick the second latest biopsy date from the 25
         #   - 1st remove the laste biopsy visidt
         #   - 2nd pick out the max biopsy visitdt after 1st step
nocancer_biopsy_2nd_date = nocancer_biopsy_2nd_start %>%
  group_by(ID) %>%
  filter(visitdt != max(visitdt))
nocancer_biopsy_2nd_date = nocancer_biopsy_2nd_date %>%
  group_by(ID) %>%
  filter(visitdt == max(visitdt))

         # merge the biopsy date with the psadt (for nocancer)
summary(levels(factor(nocancer_psa_all$ID))) # 2486
nocancer_psa_2nd = nocancer_psa_all[nocancer_psa_all$ID %in% nocancer_biopsy_2nd_date$ID, ]
nocancer_psa_biopsy_2nd = merge(nocancer_psa_2nd, nocancer_biopsy_2nd_date, by = 'ID', all.x = T)

nocancer_psa_biopsy_2nd$date_gap = ifelse(nocancer_psa_biopsy_2nd$visitdt - nocancer_psa_biopsy_2nd$psadt <0 |
                                        nocancer_psa_biopsy_2nd$visitdt - nocancer_psa_biopsy_2nd$psadt > 730, NA, 
                                      nocancer_psa_biopsy_2nd$visitdt - nocancer_psa_biopsy_2nd$psadt)
         # Pick out PSA records within 2 years prior than the 2nd latest biopsy
nocancer_psa_biopsy_2nd = nocancer_psa_biopsy_2nd[complete.cases(nocancer_psa_biopsy_2nd$date_gap), ]
         # Pick out the latest PSA records within 2 years prior than the 2nd latest biopsy
nocancer_psa_biopsy_2nd_latest <- nocancer_psa_biopsy_2nd %>% 
  group_by(ID) %>% 
  slice(which.min(date_gap)) # 18 eligible
#----------------------------------------------
nocancer_latest2 = nocancer_psa_biopsy_2nd_latest  # 15
#----------------------------------------------
summary(levels(factor(nocancer_biopsy_2nd_start$ID))) # 25 
# Step 3:
    # Go to third latest: in the started set of step 2, remove ID has been pick in step 2
length(unique(nocancer_biopsy_2nd_start$ID[!( nocancer_biopsy_2nd_start$ID %in% nocancer_latest2$ID)])) # 10 = 25 -15 
tmp_ls <- unique(nocancer_biopsy_2nd_start$ID[!( nocancer_biopsy_2nd_start$ID %in% nocancer_latest2$ID)])
    # The third latest biopsy for the rest 10 subjects
    # no_cancer_biopsy: unfiltered biopsy information for ID without cancer
nocancer_biopsy_3rd_start = nocancer_biopsy[which(nocancer_biopsy$ID %in% tmp_ls),]
summary(levels(factor(nocancer_biopsy_3rd_start$ID)))# 10 eligible 
    # Remove those with only 2 biopsies
n_occur <- data.frame(table(nocancer_biopsy_3rd_start$ID))
nocancer_biopsy_3rd_start = nocancer_biopsy_3rd_start[nocancer_biopsy_3rd_start$ID %in% n_occur$Var1[n_occur$Freq > 2],]
summary(levels(factor(nocancer_biopsy_3rd_start$ID))) # 3

    # pick the third latest biopsy date
nocancer_biopsy_3rd_date = nocancer_biopsy_3rd_start %>%
  group_by(ID) %>%
  filter(visitdt != max(visitdt))
nocancer_biopsy_3rd_date = nocancer_biopsy_3rd_date  %>%
  group_by(ID) %>%
  filter(visitdt != max(visitdt))
nocancer_biopsy_3rd_date = nocancer_biopsy_3rd_date %>%
  group_by(ID) %>%
  filter(visitdt == max(visitdt))

    # Merge with psadt
nocancer_psa3 = nocancer_psa_all[nocancer_psa_all$ID %in% nocancer_biopsy_3rd_date$ID, ]
nocancer_psa_biopsy_3rd = merge(nocancer_psa3, nocancer_biopsy_3rd_date, by = 'ID', all.x = T)
nocancer_psa_biopsy_3rd$date_gap = ifelse(nocancer_psa_biopsy_3rd$visitdt - nocancer_psa_biopsy_3rd$psadt <0 |
                                        nocancer_psa_biopsy_3rd$visitdt - nocancer_psa_biopsy_3rd$psadt > 730, NA, 
                                      nocancer_psa_biopsy_3rd$visitdt - nocancer_psa_biopsy_3rd$psadt)
summary(levels(factor(nocancer_psa_biopsy_3rd$ID))) # 3
nocancer_psa_biopsy_3rd = nocancer_psa_biopsy_3rd[complete.cases(nocancer_psa_biopsy_3rd$date_gap), ] # none

nocancer_psa_biopsy_3rd_latest <- nocancer_psa_biopsy_3rd %>% 
  group_by(ID) %>% 
  slice(which.min(date_gap)) # 18 eligible
#----------------------------------------------
nocancer_latest3 = nocancer_psa_biopsy_3rd_latest  # 3
#----------------------------------------------
#########################################################################################################
nocancer = rbind(nocancer_latest[, c('ID', "PSA", "psadt", "visitdt")], 
                 nocancer_latest2[, c('ID', "PSA", "psadt", "visitdt")],
                 nocancer_latest3[, c('ID', "PSA", "psadt", "visitdt")])

summary(levels(factor(nocancer$ID))) # 2405

nocancer_f = merge(nocancer, select_baseline, by = 'ID') # 2405

# Test with the 'one from function
# test2 <- merge(nocancer[, c('ID', 'psadt', 'visitdt')],
#                test[, c('ID', 'psadt','Biopsy_date')], by = 'ID')
# summary(test2$psadt.x - test2$psadt.y) # all 0
# summary(test2$visitdt - test2$Biopsy_date) # all 0

nocancer_f$latestAge = round((nocancer_f$visitdt - nocancer_f$REGDT)/365, 2) + nocancer_f$age


#----------------------
# Pick the DRE for no cancer  cases
nocancer_dre = select_alldre[select_alldre$ID %in% nocancer_f$ID, ]
summary(levels(factor(nocancer_dre$ID))) #2402
nocancer_dre = merge(nocancer_f[, c(1:4)], nocancer_dre, by = "ID", all.y = T)
nocancer_dre$date_gap = ifelse(nocancer_dre$visitdt - nocancer_dre$dredt < 0 | 
                                nocancer_dre$visitdt - nocancer_dre$dredt > 730, NA, 
                              nocancer_dre$visitdt - nocancer_dre$dredt)
summary(levels(factor(nocancer_dre$ID))) # 2402
nocancer_dre = nocancer_dre[complete.cases(nocancer_dre$date_gap),] 
summary(levels(factor(nocancer_dre$ID))) # 2215 with DRE, 187 missing
nocancer_dre2 = nocancer_dre%>% 
  group_by(ID) %>% 
  slice(which.min(date_gap))

##---------------------------------------
nocancer_dre2 = nocancer_dre2[, c("ID","dre")] 
nocancer_f = merge(nocancer_f, nocancer_dre2, by = "ID", all.x = T)
nocancer_f$latestBiopp = nocancer_f$visitdt
nocancer_f$pca = rep(0, nrow(nocancer_f))
nocancer_f$biopsyflag = rep(1, nrow(nocancer_f))

# Remove the redudant variable 
cancer_f <- cancer_f[,-which(colnames(cancer_f) == "pcadt")]
nocancer_f <- nocancer_f[,-which(colnames(nocancer_f) =="visitdt")]
##---------------------------------------
Biopp <- rbind(cancer_f, nocancer_f)
summary(levels(factor(Biopp$ID))) # 4766
##---------------------------------------

## No biopsy 
Nobiopp_psa = select_psa[-which(select_psa$ID %in% ID_biopsy$ID), ]
Nobiopp_psa_latest <- Nobiopp_psa%>% 
  group_by(ID) %>% 
  slice(which.max(psadt))

Nobiopp <- merge(Nobiopp_psa_latest[, c(1:3)],select_baseline, by = 'ID' )
Nobiopp$latestAge = round((Nobiopp$psadt- Nobiopp$REGDT)/365, 2)+ Nobiopp$age # Some have only the baseline PSA
summary(levels(factor(Nobiopp$ID))) # 29772

Nobiopp_dre = select_alldre[select_alldre$ID %in% Nobiopp$ID, ]
Nobiopp_dre = merge(Nobiopp[, c(1:3)], Nobiopp_dre, by = "ID", all.y = T)
summary(levels(factor(Nobiopp_dre$ID))) # 28770

Nobiopp_dre2 = Nobiopp_dre%>% 
  group_by(ID) %>% 
  slice(which.max(dredt))
summary(levels(factor(Nobiopp_dre2$ID))) # 29770
Nobiopp_dre2 = Nobiopp_dre2[, c("ID","dre")]
Nobiopp = merge(Nobiopp, Nobiopp_dre2, by = 'ID', all.x = T)

Nobiopp$latestBiopp <- rep(NA, nrow(Nobiopp))
Nobiopp$pca <- rep(NA, nrow(Nobiopp))
Nobiopp$biopsyflag = rep(0, nrow(Nobiopp))

select_f <- rbind(Biopp, Nobiopp) # 34538
master_select <- merge(select_f, select_endpoint[, c(1, 2)], by = 'ID')
write.csv(master_select, '/Users/yiyaochen/Documents/Causal Inference/code/master_SELECT_latestPSA_Oct2020.csv')

#---------------------------------------------------------------------
master_select$Age55flag <- ifelse(master_select$latestAge < 55, 0, 1)
table(master_select$Age55flag, master_select$racesel)
master_select$raceflag <- ifelse(master_select$racesel == 'UNKNOWN', NA,
                            ifelse(master_select$racesel == 'BLACK',1, 0))

master_select <- master_select[- which(master_select$latestAge <55 &
                                          master_select$raceflag == 0), ]

table(master_select$Age55flag, master_select$raceflag)

master_select$psaflag <- ifelse(master_select$PSA < 4, 1, 
                                ifelse(master_select$PSA >= 4 & master_select$PSA <= 10, 2, 3))

master_select$fhflag <- ifelse(is.na(master_select$fampca) == T, NA, ifelse(master_select$fampca == 0, 0, 1))
master_select$Gleasonflag <- ifelse(is.na(master_select$GLEASONS) == T, NA,
                                    ifelse(master_select$GLEASONS <= 6, 1, 
                                           ifelse(master_select$GLEASONS == 7, 2, 3)))
master_select$PriorBioflag <- ifelse(master_select$PROSBX == 'Y',1, 0)
master_select$PSA4flag <- ifelse(master_select$PSA > 4, 1, 0)


dt_select <- subset(master_select, select = c("ID","latestAge", "PSA", "dre", "biopsyflag", "raceflag", 
                                              "fhflag","PriorBioflag","psaflag","PSA4flag", "Gleasonflag", 'pca'))
colnames(dt_select) = c('ID','Age','PSA','DRE','biopsyflag','Race','FH','PB','psaflag',"PSA4flag",'Gleasonflag', 'PCa')

dt_select <- dt_select[complete.cases(dt_select[, -c(11, 12)]),-c(11)] # 33257
table(tmp$biopsyflag) 
tmp$PSA2 <- ifelse(tmp$PSA == 0, log2(tmp$PSA + 0.01), log2(tmp$PSA))
tmp$DRE2 <- tmp$DRE



rg <- glm(PCa ~  PSA2 + FH + PB + DRE2, 
                    family = binomial(link = "logit"), data = tmp[tmp$biopsyflag == 1, ])
round(exp(cbind(coef(rg), confint(rg))), 2)


table(tmp[tmp$biopsyflag == 1, ]$PCa, exclude = NULL)


