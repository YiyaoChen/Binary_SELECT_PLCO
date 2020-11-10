##########################################################################################
library(gridExtra)


## PSA
tb1 <- table(PLCO$biopsyflag, PLCO$PSA4flag)
tb2 <- table(SELECT$biopsyflag, SELECT$PSA4flag)
psa_OR_dt <- data.frame(OR = c(tb1[2,2]/tb1[2,1]/(tb1[1,2]/tb1[1,1]), tb2[2,2]/tb2[2,1]/(tb2[1,2]/tb2[1,1])),
                        Var = rep("PSA > 4 ng/ml", 2),
                        Group = c("PLCO", "SELECT"),
                        #CohortSize = c(sum(tb1), sum(tb2)),
                        Percent = c(round((tb1[2,2]+tb1[1,2])/sum(tb1) *100, 2), 
                                    round((tb2[2,2] + tb2[1, 2])/sum(tb2) *100, 2)))

xmin = 5
xmax = 15
jumpsize = (xmax-xmin)/10*1.5
psa_OR_dt$Position <- ifelse((xmax - psa_OR_dt$Percent) >= (psa_OR_dt$Percent - xmin), 
                             psa_OR_dt$Percent + jumpsize,psa_OR_dt$Percent - jumpsize )

p_PSA <- ggplot(psa_OR_dt, aes(x = Percent, y = OR))+
  geom_point(aes(color = Group, fill = Group, size = 5), shape = 21)+
  geom_text(aes(x = Position, label = Group),size= 5, check_overlap = F)+
  xlab('')+
  ylab('')+
  scale_y_log10(breaks = c(20, 25, 30,  50), limits = c(20, 50))+
  scale_x_continuous(breaks = c(5, 10, 15), limits = c(xmin, xmax),
                     labels = c('5%', '10%', '15%'))+
  theme_bw()+
  scale_fill_manual(values = c("lightgreen", "red"), guide = FALSE)+
  scale_color_manual(values = c("black", "black"), guide = FALSE)+
  theme(legend.position = 'none',
        legend.title = element_text(face="bold",size=12),
        legend.text=element_text(size=12),
        legend.key.size = unit(0.9, "cm"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12))

p_PSA <- p_PSA + facet_wrap(~ Var)+
  theme(strip.text.x = element_text(size = 14, face="bold"))

## DRE = abnormal
tb1 <- table(PLCO$biopsyflag, PLCO$DRE_binary, exclude = NULL)
tb2 <- table(SELECT$biopsyflag, SELECT$DRE_binary, exclude = NULL)
dre_OR_dt <- data.frame(OR = c(tb1[2,2]/tb1[2,1]/(tb1[1,2]/tb1[1,1]), tb2[2,2]/tb2[2,1]/(tb2[1,2]/tb2[1,1])),
                        Var = rep("Abnormal DRE", 2),
                        Group = c("PLCO", "SELECT"),
                        #CohortSize = c(sum(tb1), sum(tb2)),
                        Percent = c(round((tb1[2,2] + tb1[1, 2])/sum(tb1) *100, 2), round((tb2[2,2] + tb2[1, 2])/sum(tb2) *100, 2)))

xmin = 3
xmax = 10
jumpsize = (xmax-xmin)/10*1.7
dre_OR_dt$Position <- ifelse((xmax - dre_OR_dt$Percent) >= (dre_OR_dt$Percent - xmin), 
                             dre_OR_dt$Percent + jumpsize, dre_OR_dt$Percent - jumpsize )

p_DRE <- ggplot(dre_OR_dt, aes(x = Percent, y = OR))+
  geom_point(aes(color = Group, fill = Group, size = 5), shape = 21)+
  geom_text(aes(x = Position, label = Group),size= 5, check_overlap = F)+
  xlab('')+
  ylab('')+
  scale_y_log10(breaks = c(10, 13, 30, 33), limits = c(10, 33))+
  scale_x_continuous(breaks = c(3, 5, 10), limits = c(xmin, xmax+1), labels = c('3%','5%','10%'))+
  theme_bw()+
  scale_fill_manual(values = c("lightgreen", "red"), guide = FALSE)+
  scale_color_manual(values = c("black", "black"), guide = FALSE)+
  theme(legend.position = 'none',
        legend.title = element_text(face="bold",size=12),
        legend.text=element_text(size=12),
        legend.key.size = unit(0.9, "cm"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12))

p_DRE <- p_DRE + facet_wrap(~ Var)+
  theme(strip.text.x = element_text(size = 14, face="bold"))

## Race = AA
tb1 <- table(PLCO$biopsyflag, PLCO$Race, exclude = NULL)
tb2 <- table(SELECT$biopsyflag, SELECT$Race, exclude = NULL)
Race_OR_dt <- data.frame(OR = c(tb1[2,2]/tb1[2,1]/(tb1[1,2]/tb1[1,1]), tb2[2,2]/tb2[2,1]/(tb2[1,2]/tb2[1,1])),
                         Var = rep("African ancestry = yes", 2),
                         Group = c("PLCO", "SELECT"),
                         Percent = c(round((tb1[2,2] + tb1[1, 2])/sum(tb1) *100, 2), round((tb2[2,2] + tb2[1, 2])/sum(tb2) *100, 2)))

xmin = 1
xmax = 15
jumpsize = (xmax-xmin)/10*1.7
Race_OR_dt$Position <- ifelse((xmax - Race_OR_dt$Percent) >= (Race_OR_dt$Percent - xmin), 
                              Race_OR_dt$Percent + jumpsize, Race_OR_dt$Percent - jumpsize )

p_race <- ggplot(Race_OR_dt, aes(x = Percent, y = OR))+
  geom_point(aes(color = Group, fill = Group, size = 5), shape = 21)+
  geom_text(aes(x = Position, label = Group),size= 5, check_overlap = F)+
  xlab('')+
  ylab('')+
  scale_y_log10(breaks = c(0.85, 0.95, 1.3, 1.35), limits = c(0.85, 1.35))+
  scale_x_continuous(breaks = c(1, 5, 10, 15), limits = c(xmin, xmax+1), labels = c('1%','5%','10%','15%'))+
  theme_bw()+
  scale_fill_manual(values = c("lightgreen", "red"), guide = FALSE)+
  scale_color_manual(values = c("black", "black"), guide = FALSE)+
  theme(legend.position = 'none',
        legend.title = element_text(face="bold",size=12),
        legend.text=element_text(size=12),
        legend.key.size = unit(0.9, "cm"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12))

p_race <- p_race + facet_wrap(~ Var)+
  theme(strip.text.x = element_text(size = 14, face="bold"))

## FH = Yes
tb1 <- table(PLCO$biopsyflag, PLCO$FH, exclude = NULL)
tb2 <- table(SELECT$biopsyflag, SELECT$FH, exclude = NULL)
FH_OR_dt <- data.frame(OR = c(tb1[2,2]/tb1[2,1]/(tb1[1,2]/tb1[1,1]), tb2[2,2]/tb2[2,1]/(tb2[1,2]/tb2[1,1])),
                       Var = rep("Family history = yes", 2),
                       Group = c("PLCO", "SELECT"),
                       Percent = c(round((tb1[2,2] + tb1[1, 2])/sum(tb1) *100, 2), round((tb2[2,2] + tb2[1, 2])/sum(tb2) *100, 2)))

xmin = 5
xmax = 20
jumpsize = (xmax-xmin)/10*1.5
FH_OR_dt$Position <- ifelse((xmax - FH_OR_dt$Percent) >= (FH_OR_dt$Percent - xmin), 
                            FH_OR_dt$Percent + jumpsize, FH_OR_dt$Percent - jumpsize )

p_FH <- ggplot(FH_OR_dt, aes(x = Percent, y = OR))+
  geom_point(aes(color = Group, fill = Group, size = 5), shape = 21)+
  geom_text(aes(x = Position, label = Group),size= 5, check_overlap = F)+
  xlab('')+
  ylab('')+
  scale_y_log10(breaks = c(1.5, 1.7), limits = c(1.4, 1.72))+
  scale_x_continuous(breaks = c(5, 10, 15, 20), limits = c(5, 20), labels = c('5%','10%','15%', '20%'))+
  theme_bw()+
  scale_fill_manual(values = c("lightgreen", "red"), guide = FALSE)+
  scale_color_manual(values = c("black", "black"), guide = FALSE)+
  theme(legend.position = 'none',
        legend.title = element_text(face="bold",size=12),
        legend.text=element_text(size=12),
        legend.key.size = unit(0.9, "cm"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12))

p_FH <- p_FH + facet_wrap(~ Var)+
  theme(strip.text.x = element_text(size = 14, face="bold"))


## PB
tb1 <- table(PLCO$biopsyflag, PLCO$PB, exclude = NULL)
tb2 <- table(SELECT$biopsyflag, SELECT$PB, exclude = NULL)
PB_OR_dt <- data.frame(OR = c(tb1[2,2]/tb1[2,1]/(tb1[1,2]/tb1[1,1]), tb2[2,2]/tb2[2,1]/(tb2[1,2]/tb2[1,1])),
                       Var = rep("Prior negative biopsy = yes", 2),
                       Group = c("PLCO", "SELECT"),
                       Percent = c(round((tb1[2,2] + tb1[1, 2])/sum(tb1) *100, 2), round((tb2[2,2] + tb2[1, 2])/sum(tb2) *100, 2)))

xmin = 0
xmax = 10
jumpsize = (xmax-xmin)/10*1.5
PB_OR_dt$Position <- ifelse((xmax - PB_OR_dt$Percent) >= (PB_OR_dt$Percent - xmin), 
                            PB_OR_dt$Percent + jumpsize, PB_OR_dt$Percent - jumpsize )

p_PB <- ggplot(PB_OR_dt, aes(x = Percent, y = OR))+
  geom_point(aes(color = Group, fill = Group, size = 5), shape = 21)+
  geom_text(aes(x = Position, label = Group),size= 5, check_overlap = F)+
  xlab('')+
  ylab('')+
  scale_y_log10(breaks = c(1.5, 2,  3, 3.2), limits = c(1.5, 3.2))+
  scale_x_continuous(breaks = c(1, 5, 10), limits = c(xmin,xmax), labels = c('1%','5%','10%'))+
  theme_bw()+
  scale_fill_manual(values = c("lightgreen", "red"), guide = FALSE)+
  scale_color_manual(values = c("black", "black"), guide = FALSE)+
  theme(legend.position = 'none',
        legend.title = element_text(face="bold",size=12),
        legend.text=element_text(size=12),
        legend.key.size = unit(0.9, "cm"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12))

p_PB <- p_PB + facet_wrap(~ Var)+
  theme(strip.text.x = element_text(size = 14, face="bold"))



### Age > 75
PLCO$Age2flag <- as.factor(ifelse(PLCO$Age > 75, 1, 0))
SELECT$Age2flag <- as.factor(ifelse(SELECT$Age > 75, 1, 0))

tb1 <- table(PLCO$biopsyflag, PLCO$Age2flag, exclude = NULL)
tb2 <- table(SELECT$biopsyflag, SELECT$Age2flag, exclude = NULL)
Age_OR_dt <- data.frame(OR = c(tb1[2,2]/tb1[2,1]/(tb1[1,2]/tb1[1,1]), tb2[2,2]/tb2[2,1]/(tb2[1,2]/tb2[1,1])),
                        Var = rep("Age > 75", 2),
                        Group = c("PLCO", "SELECT"),
                        Percent = c(round((tb1[2,2] + tb1[1, 2])/sum(tb1) *100, 2), round((tb2[2,2] + tb2[1, 2])/sum(tb2) *100, 2)))

xmin = 5
xmax = 20
jumpsize = (xmax-xmin)/10*1.5
Age_OR_dt$Position <- ifelse((xmax - Age_OR_dt$Percent) >= (Age_OR_dt$Percent - xmin), 
                             Age_OR_dt$Percent + jumpsize, Age_OR_dt$Percent - jumpsize )

p_Age <- ggplot(Age_OR_dt, aes(x = Percent, y = OR))+
  geom_point(aes(color = Group, fill = Group, size = 5), shape = 21)+
  geom_text(aes(x = Position, label = Group),size= 5, check_overlap = F)+
  xlab('')+
  ylab('')+
  scale_y_log10(breaks = c(0.4, 0.45,0.5, 0.6), limits = c(0.4, 0.62))+
  scale_x_continuous(breaks = c(5, 10, 15, 20), limits = c(5, 20), labels = c('5%', '10%','15%', '20%'))+
  theme_bw()+
  scale_fill_manual(values = c("lightgreen", "red"), guide = FALSE)+
  scale_color_manual(values = c("black", "black"), guide = FALSE)+
  theme(legend.position = 'none',
        legend.title = element_text(face="bold",size=12),
        legend.text=element_text(size=12),
        legend.key.size = unit(0.9, "cm"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12))

p_Age <- p_Age + facet_wrap(~ Var)+
  theme(strip.text.x = element_text(size = 14, face="bold"))

p <- grid.arrange(p_PSA, p_DRE, p_PB, p_FH,p_Age, p_race, ncol=2)
ggsave(filename="/Users/yiyaochen/Binary_SELECT_PLCO/Biopsy_PBCG_OR_3.png", p, width = 8, height = 10.2, dpi = 500)

#-------------------------------------------------------------------------------------#
## For creating the x,y-axis. 
PBCG_OR_dt <- rbind(psa_OR_dt, dre_OR_dt, Race_OR_dt, FH_OR_dt, PB_OR_dt, Age_OR_dt)

p <- ggplot(PBCG_OR_dt, aes(x = Percent, y = OR, color = Group))+
  geom_point(aes(size = 5))+
  geom_text(aes(x = Position, label = Group),size= 5, check_overlap = F)+
  xlab('Prevalence of risk factor')+
  ylab('Odds ratio for having a biopsy')+
  # scale_y_log10(pretty_breaks(n = 2))+
  #scale_x_continuous(breaks = pretty_breaks(n = 3))+
  theme_bw()+
  scale_color_manual(values = c("#4393c3", "black"), guide = FALSE)+
  scale_size(limits = c(34500, 35500), guide = FALSE)+
  theme(legend.position = 'none',
        legend.title = element_text(face="bold",size=12),
        legend.text=element_text(size=12),
        legend.key.size = unit(0.9, "cm"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.title.x = element_text(face="bold",size=14),
        axis.title.y = element_text(face="bold",size=14))

p <- p + facet_wrap(~ Var, ncol = 2, scales = 'free')
ggsave(filename="/Users/yiyaochen/Binary_SELECT_PLCO/Biopsy_PBCG_OR_4.png", p, width = 8, height = 10.2, dpi = 500)

