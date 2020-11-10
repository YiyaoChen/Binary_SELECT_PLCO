# Table A1
## Probability of having cancer on PLCO, i.e. R
# rg <- glm(formula = PCa ~ (log2PSA + DRE_binary + FH + PB+ Age +Race)^2, family = binomial(link = "logit"), 
#                     data = PLCO[PLCO$biopsyflag == 1, ])
# rg <- step(rg, direction = 'both',k = log(3813))
# BIC(rg) # 4630.326
# rg_cancer_plco <- glm(formula = PCa ~ log2PSA +Race + FH+ PB+ log2PSA:DRE_binary + DRE_binary:Age,
#                       family = binomial(link = "logit"),
#                       data = PLCO[which(PLCO$biopsyflag == 1), ])
# BIC(rg_cancer_plco) # 4630.924
# round(exp(cbind(coef(rg_cancer_plco), confint(rg_cancer_plco))), 4)

rg_cancer_plco <- glm(formula = PCa ~ log2PSA +Race + FH+ PB+ log2PSA:DRE_binary,
                      family = binomial(link = "logit"), 
                      data = PLCO[which(PLCO$biopsyflag == 1), ])
round(exp(cbind(coef(rg_cancer_plco), confint(rg_cancer_plco))), 4)

# Table A2
## Probability of being in PLCO, i.e. p_s, on pooled data
## Only used the data for participants with age >= 55 in both.
tmp1 <- SELECT[which(SELECT$Age >= 55),c('Age','PSA','log2PSA', 'DRE_binary', 'biopsyflag','Race','FH','PB','PSA4flag','PCa')]
tmp2<- PLCO[which(PLCO$Age >= 55),c('Age','PSA','log2PSA', 'DRE_binary', 'biopsyflag','Race','FH','PB','PSA4flag','PCa')]

# Pooled data
tmp1$Group <- rep(0, nrow(tmp1))
tmp2$Group <- rep(1, nrow(tmp2))
pooled_data <- rbind(tmp1, tmp2) # 62874

rg_ps_pooled <- glm(formula = Group ~ PSA4flag + DRE_binary + Age + Race + FH + PB + 
                      PSA4flag:DRE_binary + PSA4flag:PB + DRE_binary:Race + DRE_binary:PB + Age:Race + 
                      Age:PB, family = binomial(link = "logit"), data = pooled_data)
round(exp(cbind(coef(rg_ps_pooled), confint(rg_ps_pooled))), 4)

# Table A3
## Probability of being biopsy in PLCO, i.e. p_v, on PLCO
rg_pv_plco <- glm(formula = biopsyflag ~ log2PSA + PSA4flag + DRE_binary + Age + log2PSA:DRE_binary + 
                    PSA4flag:DRE_binary, family = binomial(link = "logit"), data = PLCO)
round(exp(cbind(coef(rg_pv_plco), confint(rg_pv_plco))), 4)



# Table A4
## Probability of having cancer on SELECT
rg_cancer_select <- glm(formula = PCa ~ log2PSA + DRE_binary + FH + PB, family = binomial(link = "logit"), 
                        data = SELECT[which(SELECT$biopsyflag == 1), ])
round(exp(cbind(coef(rg_cancer_select), confint(rg_cancer_select))), 4)


