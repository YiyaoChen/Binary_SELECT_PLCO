master_plco <- read.csv('/Users/yiyaochen/Documents/Causal Inference/code/Clear_data/master_PLCO_latestPSA_34929_Oct2020.csv')
master_select <- read.csv('/Users/yiyaochen/Documents/Causal Inference/code/Clear_data/master_SELECT_latestPSA_33257_Oct2020.csv')


master_plco$PSA4flag <- ifelse(master_plco$latestPSA <= 4, 0,1)
master_select$PSA4flag <- ifelse(master_select$PSA > 4, 1, 0)


dt_plco <- subset(master_plco, select = c("plco_id", "latestAge", "latestPSA", "latestDRE", "biopsyflag", "raceflag", 
                                          "fhflag","PriorBioflag", "PSA4flag", "confirmed_pros"))
colnames(dt_plco) = c('ID','Age','PSA','DRE','biopsyflag','Race','FH','PB','PSA4flag', 'PCa')

dt_select <- subset(master_select, select = c('ID','Age','PSA','DRE','biopsyflag','Race','FH','PB',"PSA4flag", 'PCa'))


rg_data_plco <- dt_plco
rg_data_select <- dt_select

rg_data_plco$DRE_binary <- ifelse(is.na(rg_data_plco$DRE) == T, 
                            NA, ifelse(rg_data_plco$DRE  == 'Suspicious',1, 0))
rg_data_select$DRE_binary <- rg_data_select$DRE

############################################################################################################################

# Remove NA from DRE, FH, PB, Race
rg_data_plco_c <- rg_data_plco[complete.cases(rg_data_plco[, c(1:11)]),] # 30246
rg_data_select_c <- rg_data_select[complete.cases(rg_data_select[, -c(10)]),] # 33257

# Create log2-scale PSA values
rg_data_plco_c$log2PSA <- ifelse(rg_data_plco_c$PSA == 0, log2(rg_data_plco_c$PSA + 0.01), log2(rg_data_plco_c$PSA))
rg_data_select_c$log2PSA <- ifelse(rg_data_select_c$PSA == 0, log2(rg_data_select_c$PSA + 0.01), log2(rg_data_select_c$PSA))

# Remove those with age less than 55 years old
rg_data_plco_c <- rg_data_plco_c[-which(rg_data_plco_c$Age < 55), ]
rg_data_select_c <- rg_data_select_c[-which(rg_data_select_c$Age < 55),]

# Output
PLCO <- rg_data_plco_c
SELECT <- rg_data_select_c
