library(sas7bdat)

PersonalBase = read.sas7bdat("/Users/yiyaochen/Documents/Trails Data/PLCO Screening Study/package/Person based file/plco135_longitudinal_040315.sas7bdat")

subPersonalBase <- subset(PersonalBase, select = c('plco_id', 'age','arm','race7','pros_fh_cnt','pros_fh','trial_ph_any', 'trial_ph_pros','surg_biopsy','surg_any',
                                                   "psa_level0","psa_level1", "psa_level2","psa_level3", "psa_level4", "psa_level5", 
                                                   'psa_days0', 'psa_days1', 'psa_days2', 'psa_days3', 'psa_days4', 'psa_days5',
                                                   "dre_result0","dre_result1", "dre_result2", 'dre_result3',
                                                   'dre_days0', 'dre_days1', 'dre_days2', 'dre_days3',
                                                   "biopplink0", "biopplink1", "biopplink2", "biopplink3", "biopplink4", "biopplink5", 'numbiopp', 
                                                   'biopproc_daysp1', 'biopproc_daysp2','biopproc_daysp3','biopproc_daysp4','biopproc_daysp5','biopproc_daysp6',
                                                   'biopproc_daysp7', 'biopproc_daysp8', 'biopproc_daysp9',
                                                   'confirmed_pros','candxdaysp', 'dx_psa', 'dx_psa_gap', 'pros_gleason'))

subPersonalBase <- subPersonalBase[which(subPersonalBase$arm == 1), ]
subPersonalBase$biopsyflag <- ifelse(subPersonalBase$confirmed_pros == 1 | subPersonalBase$numbiopp > 0, 1, 0)


# No prior prostate(any) cancer: 37518
subPersonalBase <- subPersonalBase[which(subPersonalBase$trial_ph_any != 1),]
subPersonalBase <- subPersonalBase[which(subPersonalBase$trial_ph_pros != 1),]

# No PSA in the entire study: 35142
subPersonalBase <- subPersonalBase[which(is.na(subPersonalBase$psa_level0) == F |
                                           is.na(subPersonalBase$psa_level1) == F |
                                           is.na(subPersonalBase$psa_level2) == F |
                                           is.na(subPersonalBase$psa_level3) == F |
                                           is.na(subPersonalBase$psa_level4) == F |
                                           is.na(subPersonalBase$psa_level5) == F),]


subPersonalBase$dre_result0_m <- ifelse(subPersonalBase$dre_result0 == 1 , 'Normal', ifelse(subPersonalBase$dre_result0 == 2, 'Suspicious', ifelse(subPersonalBase$dre_result0 == 3, 'AbnormalNSP',NA)))
subPersonalBase$dre_result1_m <- ifelse(subPersonalBase$dre_result1 == 1 , 'Normal', ifelse(subPersonalBase$dre_result1 == 2, 'Suspicious', ifelse(subPersonalBase$dre_result1 == 3, 'AbnormalNSP',NA)))
subPersonalBase$dre_result2_m <- ifelse(subPersonalBase$dre_result2 == 1 , 'Normal', ifelse(subPersonalBase$dre_result2 == 2, 'Suspicious', ifelse(subPersonalBase$dre_result2 == 3, 'AbnormalNSP',NA)))
subPersonalBase$dre_result3_m <- ifelse(subPersonalBase$dre_result3 == 1 , 'Normal', ifelse(subPersonalBase$dre_result3 == 2, 'Suspicious', ifelse(subPersonalBase$dre_result3 == 3, 'AbnormalNSP',NA)))


# Without biopsy
NoBiopp <- subPersonalBase[which(subPersonalBase$numbiopp == 0 & subPersonalBase$confirmed_pros == 0),] #28515
NoBiopp$latestPSA <- rep(NA, nrow(NoBiopp))
NoBiopp$latestPSA_taken <- rep(NA, nrow(NoBiopp))
NoBiopp$latestAge <- rep(NA, nrow(NoBiopp))

NoBiopp$latestPSA[which(is.na(NoBiopp$psa_level5) == F)] <- NoBiopp[which(is.na(NoBiopp$psa_level5) == F), ]$psa_level5
NoBiopp$latestPSA_taken[which(is.na(NoBiopp$psa_level5) == F)] <- rep(5, nrow(NoBiopp[which(is.na(NoBiopp$psa_level5) == F), ]))
NoBiopp$latestAge[which(is.na(NoBiopp$psa_level5) == F)] <- round(NoBiopp[which(is.na(NoBiopp$psa_level5) == F), ]$psa_days5/365, 2)

NoBiopp$latestPSA[which(is.na(NoBiopp$psa_level4) == F & 
                          is.na(NoBiopp$psa_level5) == T)] <- NoBiopp[which(is.na(NoBiopp$psa_level4) == F & is.na(NoBiopp$psa_level5) == T), ]$psa_level4
NoBiopp$latestPSA_taken[which(is.na(NoBiopp$psa_level4) == F & 
                                is.na(NoBiopp$psa_level5) == T)] <- rep(4, nrow(NoBiopp[which(is.na(NoBiopp$psa_level4) == F & is.na(NoBiopp$psa_level5) == T), ])) 
NoBiopp$latestAge[which(is.na(NoBiopp$psa_level4) == F & is.na(NoBiopp$psa_level5) == T)] <- round(NoBiopp[which(is.na(NoBiopp$psa_level4) == F & is.na(NoBiopp$psa_level5) == T), ]$psa_days4/365, 2)

NoBiopp$latestPSA[which(is.na(NoBiopp$psa_level3) == F&
                          is.na(NoBiopp$psa_level4) == T &
                          is.na(NoBiopp$psa_level5) == T)] <- NoBiopp[which(is.na(NoBiopp$psa_level3) == F & 
                                                                              is.na(NoBiopp$psa_level4) == T & 
                                                                              is.na(NoBiopp$psa_level5) == T), ]$psa_level3
NoBiopp$latestPSA_taken[which(is.na(NoBiopp$psa_level3) == F&
                                is.na(NoBiopp$psa_level4) == T &
                                is.na(NoBiopp$psa_level5) == T)] <- rep(3, nrow(NoBiopp[which(is.na(NoBiopp$psa_level3) == F & 
                                                                                                is.na(NoBiopp$psa_level4) == T & 
                                                                                                is.na(NoBiopp$psa_level5) == T), ]))
NoBiopp$latestAge[which(is.na(NoBiopp$psa_level3) == F& 
                          is.na(NoBiopp$psa_level4) == T &
                          is.na(NoBiopp$psa_level5) == T)] <- round(NoBiopp[which(is.na(NoBiopp$psa_level3) == F& 
                                                                                    is.na(NoBiopp$psa_level4) == T &
                                                                                    is.na(NoBiopp$psa_level5) == T), ]$psa_days3/365, 2)

NoBiopp$latestPSA[which(is.na(NoBiopp$psa_level2) == F&
                          is.na(NoBiopp$psa_level3) == T &
                          is.na(NoBiopp$psa_level4) == T &
                          is.na(NoBiopp$psa_level5) == T)] <- NoBiopp[which(is.na(NoBiopp$psa_level2) == F & 
                                                                              is.na(NoBiopp$psa_level3) == T &
                                                                              is.na(NoBiopp$psa_level4) == T & 
                                                                              is.na(NoBiopp$psa_level5) == T), ]$psa_level2
NoBiopp$latestPSA_taken[which(is.na(NoBiopp$psa_level2) == F&
                                is.na(NoBiopp$psa_level3) == T &
                                is.na(NoBiopp$psa_level4) == T &
                                is.na(NoBiopp$psa_level5) == T)] <- rep(2, nrow(NoBiopp[which(is.na(NoBiopp$psa_level2) == F & 
                                                                                                is.na(NoBiopp$psa_level3) == T &
                                                                                                is.na(NoBiopp$psa_level4) == T & 
                                                                                                is.na(NoBiopp$psa_level5) == T), ]))
NoBiopp$latestAge[which(is.na(NoBiopp$psa_level2) == F& 
                          is.na(NoBiopp$psa_level3) == T &
                          is.na(NoBiopp$psa_level4) == T &
                          is.na(NoBiopp$psa_level5) == T)] <- round(NoBiopp[which(is.na(NoBiopp$psa_level2) == F& 
                                                                                    is.na(NoBiopp$psa_level3) == T &
                                                                                    is.na(NoBiopp$psa_level4) == T &
                                                                                    is.na(NoBiopp$psa_level5) == T), ]$psa_days2/365, 2)
NoBiopp$latestPSA[which(is.na(NoBiopp$psa_level1) == F&
                          is.na(NoBiopp$psa_level2) == T &
                          is.na(NoBiopp$psa_level3) == T &
                          is.na(NoBiopp$psa_level4) == T &
                          is.na(NoBiopp$psa_level5) == T)] <- NoBiopp[which(is.na(NoBiopp$psa_level1) == F & 
                                                                              is.na(NoBiopp$psa_level2) == T &
                                                                              is.na(NoBiopp$psa_level3) == T &
                                                                              is.na(NoBiopp$psa_level4) == T & 
                                                                              is.na(NoBiopp$psa_level5) == T), ]$psa_level1
NoBiopp$latestPSA_taken[which(is.na(NoBiopp$psa_level1) == F&
                                is.na(NoBiopp$psa_level2) == T &
                                is.na(NoBiopp$psa_level3) == T &
                                is.na(NoBiopp$psa_level4) == T &
                                is.na(NoBiopp$psa_level5) == T)] <- rep(1, nrow(NoBiopp[which(is.na(NoBiopp$psa_level1) == F & 
                                                                                                is.na(NoBiopp$psa_level2) == T &
                                                                                                is.na(NoBiopp$psa_level3) == T &
                                                                                                is.na(NoBiopp$psa_level4) == T & 
                                                                                                is.na(NoBiopp$psa_level5) == T), ]))
NoBiopp$latestAge[which(is.na(NoBiopp$psa_level1) == F&
                          is.na(NoBiopp$psa_level2) == T &
                          is.na(NoBiopp$psa_level3) == T &
                          is.na(NoBiopp$psa_level4) == T &
                          is.na(NoBiopp$psa_level5) == T)] <- round(NoBiopp[which(is.na(NoBiopp$psa_level1) == F & 
                                                                                    is.na(NoBiopp$psa_level2) == T &
                                                                                    is.na(NoBiopp$psa_level3) == T &
                                                                                    is.na(NoBiopp$psa_level4) == T & 
                                                                                    is.na(NoBiopp$psa_level5) == T), ]$psa_days1/365, 2)
NoBiopp$latestPSA[which(is.na(NoBiopp$psa_level0) == F&
                          is.na(NoBiopp$psa_level1) == T &
                          is.na(NoBiopp$psa_level2) == T &
                          is.na(NoBiopp$psa_level3) == T &
                          is.na(NoBiopp$psa_level4) == T &
                          is.na(NoBiopp$psa_level5) == T)] <- NoBiopp[which(is.na(NoBiopp$psa_level0) == F & 
                                                                              is.na(NoBiopp$psa_level1) == T &
                                                                              is.na(NoBiopp$psa_level2) == T &
                                                                              is.na(NoBiopp$psa_level3) == T &
                                                                              is.na(NoBiopp$psa_level4) == T & 
                                                                              is.na(NoBiopp$psa_level5) == T), ]$psa_level0
NoBiopp$latestPSA_taken[which(is.na(NoBiopp$psa_level0) == F&
                                is.na(NoBiopp$psa_level1) == T &
                                is.na(NoBiopp$psa_level2) == T &
                                is.na(NoBiopp$psa_level3) == T &
                                is.na(NoBiopp$psa_level4) == T &
                                is.na(NoBiopp$psa_level5) == T)] <- rep(0, nrow(NoBiopp[which(is.na(NoBiopp$psa_level0) == F & 
                                                                                                is.na(NoBiopp$psa_level1) == T &
                                                                                                is.na(NoBiopp$psa_level2) == T &
                                                                                                is.na(NoBiopp$psa_level3) == T &
                                                                                                is.na(NoBiopp$psa_level4) == T & 
                                                                                                is.na(NoBiopp$psa_level5) == T), ]))
NoBiopp$latestAge[which(is.na(NoBiopp$psa_level0) == F&
                          is.na(NoBiopp$psa_level1) == T &
                          is.na(NoBiopp$psa_level2) == T &
                          is.na(NoBiopp$psa_level3) == T &
                          is.na(NoBiopp$psa_level4) == T &
                          is.na(NoBiopp$psa_level5) == T)] <- round(NoBiopp[which(is.na(NoBiopp$psa_level0) == F & 
                                                                                    is.na(NoBiopp$psa_level1) == T &
                                                                                    is.na(NoBiopp$psa_level2) == T &
                                                                                    is.na(NoBiopp$psa_level3) == T &
                                                                                    is.na(NoBiopp$psa_level4) == T & 
                                                                                    is.na(NoBiopp$psa_level5) == T), ]$psa_days0/365, 2)

NoBiopp$latestAge <- NoBiopp$age+NoBiopp$latestAge
summary(NoBiopp)

NoBiopp$latestPSAday <- ifelse(NoBiopp$latestPSA_taken == 5, NoBiopp$psa_days5, 
                               ifelse(NoBiopp$latestPSA_taken == 4, NoBiopp$psa_days4, 
                                      ifelse(NoBiopp$latestPSA_taken== 3, NoBiopp$psa_days3, 
                                             ifelse(NoBiopp$latestPSA_taken == 2, NoBiopp$psa_days2, 
                                                    ifelse(NoBiopp$latestPSA_taken == 1, NoBiopp$psa_days1, NoBiopp$psa_days0)))))

NoBiopp$latestDRE <- ifelse(is.na(NoBiopp$dre_result3_m) == F, NoBiopp$dre_result3_m, 
                            ifelse(is.na(NoBiopp$dre_result2_m) == F,  NoBiopp$dre_result2_m, 
                                   ifelse(is.na(NoBiopp$dre_result1_m) == F,  NoBiopp$dre_result1_m,
                                          ifelse(is.na(NoBiopp$dre_result0_m) == F, NoBiopp$dre_result0_m, NA))))

NoBiopp$latestDRE_d <- ifelse(is.na(NoBiopp$dre_result3_m) == F, NoBiopp$dre_days3, 
                              ifelse(is.na(NoBiopp$dre_result2_m) == F,  NoBiopp$dre_days2, 
                                     ifelse(is.na(NoBiopp$dre_result1_m) == F,  NoBiopp$dre_days1,
                                            ifelse(is.na(NoBiopp$dre_result0_m) == F, NoBiopp$dre_days0, NA))))

#---------With Biopsy--------------------------------------------------------------------------------------
Biopp <- subPersonalBase[which(subPersonalBase$numbiopp > 0 | subPersonalBase$confirmed_pros == 1),] # 6627
Biopp$latestDRE <- ifelse(is.na(Biopp$dre_result3_m) == F, Biopp$dre_result3_m, 
                          ifelse(is.na(Biopp$dre_result2_m) == F,  Biopp$dre_result2_m, 
                                 ifelse(is.na(Biopp$dre_result1_m) == F,  Biopp$dre_result1_m,
                                        ifelse(is.na(Biopp$dre_result0_m) == F, Biopp$dre_result0_m, NA))))

Biopp$latestDRE_d <- ifelse(is.na(Biopp$dre_result3_m) == F, Biopp$dre_days3, 
                            ifelse(is.na(Biopp$dre_result2_m) == F,  Biopp$dre_days2, 
                                   ifelse(is.na(Biopp$dre_result1_m) == F,  Biopp$dre_days1,
                                          ifelse(is.na(Biopp$dre_result0_m) == F, Biopp$dre_days0, NA))))


cancer <- Biopp[which(Biopp$confirmed_pros == 1), ] #4144
nocancer <- Biopp[which(Biopp$confirmed_pros == 0), ] #2483

#---------Cancer-------------
# Remove those without PSA prior to the diagnosis biopsy within 2 years 3931
# dx_psa:PSA close to diagnosis
# dx_psa_gap: gap between diagnosis and closest PSA
cancer <- Biopp[which(Biopp$confirmed_pros == 1), ]
cancer$latestPSA <- cancer$dx_psa
cancer <- cancer[which(cancer$dx_psa_gap <= 730),]
# candxdaysp: Days from randomised until diagnosis of prostate cancer
cancer$latestPSA_d <- cancer$candxdaysp - cancer$dx_psa_gap
summary(cancer$latestPSA_d) # No negative values
# Take age at cancer diagnosis as new age
cancer$latestAge <- cancer$age + round(cancer$candxdaysp/365, 2)

# For DRE within 2 years prior to the cancer diagnosis (biopsy). 
# 3931
tmp <- cancer[, c('plco_id','dre_days0', 'dre_days1', 'dre_days2', 'dre_days3', 'candxdaysp')]
tmp$dre_days0 <- ifelse(tmp$candxdaysp - tmp$dre_days0 < 0 | tmp$candxdaysp - tmp$dre_days0> 731, NA, tmp$candxdaysp - tmp$dre_days0)
tmp$dre_days1 <- ifelse(tmp$candxdaysp - tmp$dre_days1< 0 |  tmp$candxdaysp - tmp$dre_days1> 731, NA, tmp$candxdaysp - tmp$dre_days1)
tmp$dre_days2 <- ifelse(tmp$candxdaysp - tmp$dre_days2 < 0 | tmp$candxdaysp - tmp$dre_days2> 731, NA, tmp$candxdaysp - tmp$dre_days2)
tmp$dre_days3 <- ifelse(tmp$candxdaysp - tmp$dre_days3 < 0 | tmp$candxdaysp - tmp$dre_days3> 731, NA, tmp$candxdaysp - tmp$dre_days3)
# At least one eligible DRE to the cancer diagnosis
# Remove those with all dre_day0-3 be NA, i.e.number of NA = number of column.
# 2022 has at least 1 eligible DRE
tmp <- tmp[which(rowSums(is.na(tmp[, c(2:5)])) != ncol(tmp[, c(2:5)])),]
tmp$minname_dre <- colnames(tmp[, c(2:5)])[apply(tmp[, c(2:5)],1,which.min)]

cancer <- merge(cancer, tmp[, c('plco_id','minname_dre')],by = 'plco_id', all.x = T) # 1909 missing DRE
# Latest DRE results 
# If there are several eligible, pick the latest no matter what DRE result it is
cancer$latestDRE <- ifelse(cancer$minname_dre == 'dre_days0', cancer$dre_result0_m, 
                           ifelse(cancer$minname_dre == 'dre_days1', cancer$dre_result1_m, 
                                  ifelse(cancer$minname_dre == 'dre_days2', cancer$dre_result2_m, 
                                         ifelse(cancer$minname_dre == 'dre_days3', cancer$dre_result3_m, NA))))
cancer$latestDRE_d <- ifelse(cancer$minname_dre == 'dre_days0', cancer$dre_days0, 
                             ifelse(cancer$minname_dre == 'dre_days1', cancer$dre_days1, 
                                    ifelse(cancer$minname_dre == 'dre_days2', cancer$dre_days2, 
                                           ifelse(cancer$minname_dre == 'dre_days3', cancer$dre_days3, NA))))



table(cancer$latestDRE, exclude = NULL)
# 1973 > 1909 --> there are 64 with latest DRE record be missing for corresponding DRE date. 


#---------No cancer-------------
nocancer <- Biopp[which(Biopp$confirmed_pros == 0), ]
table(nocancer$numbiopp) # 2483

# Take the latest biopproc_daysp. Renew the age
nocancer$latestBiopsyday <- ifelse(is.na(nocancer$biopproc_daysp7)== F, nocancer$biopproc_daysp7, 
                                   ifelse(is.na(nocancer$biopproc_daysp6) == F, nocancer$biopproc_daysp6,
                                          ifelse(is.na(nocancer$biopproc_daysp5) == F, nocancer$biopproc_daysp5,
                                                 ifelse(is.na(nocancer$biopproc_daysp4) == F, nocancer$biopproc_daysp4,
                                                        ifelse(is.na(nocancer$biopproc_daysp3) == F, nocancer$biopproc_daysp3,
                                                               ifelse(is.na(nocancer$biopproc_daysp2) == F, nocancer$biopproc_daysp2,
                                                                      ifelse(is.na(nocancer$biopproc_daysp1) == F, nocancer$biopproc_daysp1, NA)))))))
nocancer$latestAge <- nocancer$age + round(nocancer$latestBiopsyday/365, 2)

# Take the PSA prior within 2 years 
tmp <- subset(nocancer, select = c('plco_id','psa_days0', 'psa_days1', 'psa_days2', 'psa_days3', 'psa_days4', 'psa_days5','latestBiopsyday'))
tmp$psa_days0 <- ifelse(tmp$latestBiopsyday - tmp$psa_days0 < 0 | tmp$latestBiopsyday - tmp$psa_days0> 731, NA, tmp$latestBiopsyday - tmp$psa_days0)
tmp$psa_days1 <- ifelse(tmp$latestBiopsyday - tmp$psa_days1< 0 |  tmp$latestBiopsyday - tmp$psa_days1> 731, NA, tmp$latestBiopsyday - tmp$psa_days1)
tmp$psa_days2 <- ifelse(tmp$latestBiopsyday - tmp$psa_days2 < 0 | tmp$latestBiopsyday - tmp$psa_days2> 731, NA, tmp$latestBiopsyday - tmp$psa_days2)
tmp$psa_days3 <- ifelse(tmp$latestBiopsyday - tmp$psa_days3 < 0 | tmp$latestBiopsyday - tmp$psa_days3> 731, NA, tmp$latestBiopsyday - tmp$psa_days3)
tmp$psa_days4 <- ifelse(tmp$latestBiopsyday - tmp$psa_days4 < 0 | tmp$latestBiopsyday - tmp$psa_days4> 731, NA, tmp$latestBiopsyday - tmp$psa_days4)
tmp$psa_days5 <- ifelse(tmp$latestBiopsyday - tmp$psa_days5 < 0 | tmp$latestBiopsyday - tmp$psa_days5> 731, NA, tmp$latestBiopsyday - tmp$psa_days5)

##--  Find the most closest psa_day
tmp$min <- apply(tmp[, c(2:7)], 1, min, na.rm = T)
tmp$minname <- colnames(tmp[, c(2:7)])[apply(tmp[, c(2:7)],1,which.min)]
table(tmp$minname)

tmp2 <- merge(tmp[, c('plco_id','minname')],nocancer, by = 'plco_id')

##-- Take the eligible PSA value prior to the latest biopsy within 2 years. 
tmp2$latestPSA <- ifelse(tmp2$minname == 'psa_days0', tmp2$psa_level0, 
                         ifelse(tmp2$minname == 'psa_days1', tmp2$psa_level1,
                                ifelse(tmp2$minname == 'psa_days2', tmp2$psa_level2,
                                       ifelse(tmp2$minname == 'psa_days3', tmp2$psa_level3,
                                              ifelse(tmp2$minname == 'psa_days4', tmp2$psa_level4, 
                                                     ifelse(tmp2$minname == 'psa_days5', tmp2$psa_level5,NA))))))
tmp2$latestPSA_d <- ifelse(tmp2$minname == 'psa_days0', tmp2$psa_days0, 
                           ifelse(tmp2$minname == 'psa_days1', tmp2$psa_days1,
                                  ifelse(tmp2$minname == 'psa_days2', tmp2$psa_days2,
                                         ifelse(tmp2$minname == 'psa_days3', tmp2$psa_days3,
                                                ifelse(tmp2$minname == 'psa_days4', tmp2$psa_days4, 
                                                       ifelse(tmp2$minname == 'psa_days5', tmp2$psa_days5,NA))))))


##-- U-150491-0 second most cloest psa_day
## this subject has biopsied at 1st visit.
tmp2[which(tmp2$plco_id == 'U-150491-0'), 'latestPSA'] <- 0.39
tmp2[which(tmp2$plco_id == 'U-150491-0'), 'minname'] <- 'psa_days0'

#----------------Final------------
nocancer <- tmp2 #2483
#----------------Final------------
# Take the DRE prior within 2 years at the biopsy, otherwise missing.
tmp <- subset(nocancer, select = c('plco_id','dre_days0', 'dre_days1', 'dre_days2', 'dre_days3','latestBiopsyday')) #2483
tmp$dre_days0 <- ifelse(tmp$latestBiopsyday - tmp$dre_days0 < 0 | tmp$latestBiopsyday - tmp$dre_days0> 731, NA, tmp$latestBiopsyday - tmp$dre_days0)
tmp$dre_days1 <- ifelse(tmp$latestBiopsyday - tmp$dre_days1< 0 |  tmp$latestBiopsyday - tmp$dre_days1> 731, NA, tmp$latestBiopsyday - tmp$dre_days1)
tmp$dre_days2 <- ifelse(tmp$latestBiopsyday - tmp$dre_days2 < 0 | tmp$latestBiopsyday - tmp$dre_days2> 731, NA, tmp$latestBiopsyday - tmp$dre_days2)
tmp$dre_days3 <- ifelse(tmp$latestBiopsyday - tmp$dre_days3 < 0 | tmp$latestBiopsyday - tmp$dre_days3> 731, NA, tmp$latestBiopsyday - tmp$dre_days3)

tmp <- tmp[-which(is.na(tmp$dre_days0) == T &is.na(tmp$dre_days1) == T & is.na(tmp$dre_days2) == T & is.na(tmp$dre_days3) == T), ] #2158

##-- biopsy day fixed. Find the most closest dre_day within 2 years prior than biopsy. 
tmp$min <- apply(tmp[, c(2:5)], 1, min, na.rm = T)
ID2 <- tmp$plco_id # ID contain at less one eligible dre #2158
tmp$minname_dre <- colnames(tmp[, c(2:5)])[apply(tmp[, c(2:5)],1,which.min)]
table(tmp$minname_dre, exclude = NULL)

nocancer <- merge(nocancer, tmp[, c('plco_id','minname_dre')],by = 'plco_id', all.x = T)
nocancer$latestDRE <- ifelse(nocancer$minname_dre == 'dre_days0', nocancer$dre_result0_m, 
                           ifelse(nocancer$minname_dre == 'dre_days1', nocancer$dre_result1_m, 
                                  ifelse(nocancer$minname_dre == 'dre_days2', nocancer$dre_result2_m, 
                                         ifelse(nocancer$minname_dre == 'dre_days3', nocancer$dre_result3_m, NA))))
nocancer$latestDRE_d <- ifelse(nocancer$minname_dre == 'dre_days0', nocancer$dre_days0, 
                             ifelse(nocancer$minname_dre == 'dre_days1', nocancer$dre_days1, 
                                    ifelse(nocancer$minname_dre == 'dre_days2', nocancer$dre_days2, 
                                           ifelse(nocancer$minname_dre == 'dre_days3', nocancer$dre_days3, NA))))

#---------------------------------------------------------------
cancer <- cancer[, -which(colnames(cancer) %in% c("minname_dre"))] # 3931
nocancer <- nocancer[, -which(colnames(nocancer) %in% c("minname_dre","minname", "latestBiopsyday"))] # 2483

Biopp_f <- rbind(cancer, nocancer) # 6414
Biopp_f <- subset(Biopp_f, select = c("plco_id","age", "arm","race7", "pros_fh_cnt","pros_fh","trial_ph_any","trial_ph_pros","surg_biopsy","surg_any",
                                      "confirmed_pros","candxdaysp","dx_psa","dx_psa_gap",'pros_gleason',
                                      "latestPSA", "latestAge", "latestDRE","biopsyflag"))
rg <- glm(confirmed_pros~ latestDRE,family = binomial(link = "logit"), data = Biopp_f)
summary(rg)
table(Biopp_f$confirmed_pros, Biopp_f$latestDRE, exclude = NULL)

NoBiopp <- subset(NoBiopp, select = c("plco_id","age", "arm","race7", "pros_fh_cnt","pros_fh","trial_ph_any","trial_ph_pros","surg_biopsy","surg_any",
                                      "confirmed_pros","candxdaysp","dx_psa","dx_psa_gap",'pros_gleason',
                                      "latestPSA", "latestAge", "latestDRE", "biopsyflag"))

plco_f <- rbind(Biopp_f, NoBiopp)

# race7 = 7: missing; race7 = 2: black; 
plco_f$raceflag = ifelse(plco_f$race7 == 7, NA, ifelse(plco_f$race7 < 7 & plco_f$race7 == 2, 1, 0))
# pros_fh = 1: Yes,pros_fh = NA: missing
plco_f$fhflag = ifelse(plco_f$pros_fh == 1, 1, ifelse(is.na(plco_f$pros_fh) == T | plco_f$pros_fh == 9, NA, 0))

## 1: yes; 0: no, NA: missing
plco_f$PriorBioflag = ifelse(plco_f$surg_biopsy == 0, 0, ifelse(plco_f$surg_biopsy == 1 & is.na(plco_f$surg_biopsy) == F , 1, NA))

write.csv(plco_f, '/Users/yiyaochen/Documents/Causal Inference/code/master_PLCO_latestPSA_34929_Oct2020.csv')

# Test: 30246 have complete case
# dt_plco <- subset(plco_f, select = c("plco_id", "latestAge", "latestPSA", "latestDRE", "biopsyflag", "raceflag",
#                                           "fhflag","PriorBioflag", "confirmed_pros"))
# colnames(dt_plco) = c('ID','Age','PSA','DRE','biopsyflag','Race','FH','PB', 'PCa')
# test <- dt_plco[complete.cases(dt_plco),]
# 
# nrow(test[which(test$Age < 55),])




