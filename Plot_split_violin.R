# Functions from 
# https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

##########################################################
SELECT$pred_cancer_plco <- predict(rg_cancer_plco, newdata = SELECT, probability = TRUE,type = "response" )
SELECT$pred_pv <- predict(rg_pv_plco, newdata = SELECT, probability = TRUE,type = "response" )

SELECT_biopsy <- SELECT[which(SELECT$biopsyflag == 1), ]
SELECT_biopsy$w <- SELECT_biopsy$pred_pv/sum(SELECT_biopsy$pred_pv)
sum(SELECT_biopsy$w)

CIL_f <- function(dataset){
  dt <- dataset
  CIL <- sum(dt$pred_cancer_plco - dt$PCa)/nrow(dt)
  return(CIL)
}


plotdt_PSA4 <- data.frame(x = rep('PSA > 4 ng/ml', nrow(SELECT_biopsy)),
                          label = rep('PSA > 4 ng/ml', nrow(SELECT_biopsy)),
                          y = c(SELECT_biopsy[which(SELECT_biopsy$PSA4flag == 1), ]$w,
                                SELECT_biopsy[which(SELECT_biopsy$PSA4flag == 0), ]$w),
                          m = c(rep(1, nrow(SELECT_biopsy[which(SELECT_biopsy$PSA4flag == 1), ])),
                                rep(0, nrow(SELECT_biopsy[which(SELECT_biopsy$PSA4flag == 0), ]))),
                          l = c(rep('Yes', nrow(SELECT_biopsy[which(SELECT_biopsy$PSA4flag == 1), ])),
                                rep('No', nrow(SELECT_biopsy[which(SELECT_biopsy$PSA4flag == 0), ]))))
CIL_f(SELECT_biopsy[which(SELECT_biopsy$PSA4flag == 1), ])
CIL_f(SELECT_biopsy[which(SELECT_biopsy$PSA4flag == 0), ])


violin_p_PSA4flag <-   ggplot(plotdt_PSA4 , aes(x , y, group=m, fill = as.factor(m))) +
  geom_split_violin()+ 
  scale_y_continuous(limits = c(0, 0.0006), breaks = c(0, 0.0002,0.0004, 0.0006),
                     labels = c('0.0000', '0.0002','0.0004', '0.0006'))+
  theme_bw() +
  scale_x_discrete(labels=c("PSA > 4 ng/ml" = c('   No    Yes  ')))+
  xlab("") +
  ylab("") + 
  scale_fill_manual(name = "", values = c('lightgreen', 'red'))+
  theme(strip.text.x = element_text(size=12,face="bold"))+
  theme(legend.position = 'none',
        legend.title = element_text(face="bold",size=12),
        legend.text=element_text(size=12),
        legend.key.size = unit(0.9, "cm"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12)) 
violin_p_PSA4flag <- violin_p_PSA4flag + facet_wrap(~ label)+
  theme(strip.text.x = element_text(size = 14, face="bold"))+
 annotate(geom="text", x=1.4, y=0.0005,size = 5,fontface = "bold", label=paste0(round(CIL_f(SELECT_biopsy[which(SELECT_biopsy$PSA4flag == 1), ]),3)))+
  annotate(geom="text", x=0.65, y=0.0005,size = 5,fontface = "bold", label=paste0('CIL = ', round(CIL_f(SELECT_biopsy[which(SELECT_biopsy$PSA4flag == 0), ]), 3)))+
  annotate(geom="text", x=1.41, y=0.00045,size = 4,label=paste0(nrow(SELECT_biopsy[which(SELECT_biopsy$PSA4flag == 1), ])))+
  annotate(geom="text", x=0.655, y=0.00045, size = 4,label=paste0('n = ', nrow(SELECT_biopsy[which(SELECT_biopsy$PSA4flag == 0), ])))




plotdt_DRE <- data.frame(x = rep('Abnormal DRE', nrow(SELECT_biopsy)),
                         label = rep('Abnormal DRE', nrow(SELECT_biopsy)),
                          y = c(SELECT_biopsy[which(SELECT_biopsy$DRE_binary == 1), ]$w,
                                SELECT_biopsy[which(SELECT_biopsy$DRE_binary  == 0), ]$w),
                          m = c(rep(1, nrow(SELECT_biopsy[which(SELECT_biopsy$DRE_binary  == 1), ])),
                                rep(0, nrow(SELECT_biopsy[which(SELECT_biopsy$DRE_binary  == 0), ]))),
                          l = c(rep('Yes', nrow(SELECT_biopsy[which(SELECT_biopsy$DRE_binary  == 1), ])),
                                rep('No', nrow(SELECT_biopsy[which(SELECT_biopsy$DRE_binary  == 0), ]))))
CIL_f(SELECT_biopsy[which(SELECT_biopsy$DRE_binary  == 1), ])
CIL_f(SELECT_biopsy[which(SELECT_biopsy$DRE_binary  == 0), ])

violin_p_DRE <-   ggplot(plotdt_DRE  , aes(x , y, group=m, fill = as.factor(m))) +
  geom_split_violin()+ 
  scale_y_continuous(limits = c(0, 0.0006), breaks = c(0, 0.0002,0.0004, 0.0006),
                     labels = c('0.0000', '0.0002','0.0004', '0.0006'))+
  theme_bw() +
  scale_x_discrete(labels=c("Abnormal DRE" = c('   No    Yes  ')))+
  xlab("") +
  ylab("") + 
  scale_fill_manual(name = "", values = c('lightgreen', 'red'))+
  theme(strip.text.x = element_text(size=12,face="bold"))+
  theme(legend.position = 'none',
        legend.title = element_text(face="bold",size=12),
        legend.text=element_text(size=12),
        legend.key.size = unit(0.9, "cm"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12)) 
violin_p_DRE<- violin_p_DRE + facet_wrap(~ label)+
  theme(strip.text.x = element_text(size = 14, face="bold"))+
  annotate(geom="text", x=1.4, y=0.0005,size = 5,fontface = "bold", label=paste0( round(CIL_f(SELECT_biopsy[which(SELECT_biopsy$DRE_binary == 1), ]),3)))+
  annotate(geom="text", x=0.6, y=0.0005,size = 5,fontface = "bold", label=paste0( round(CIL_f(SELECT_biopsy[which(SELECT_biopsy$DRE_binary == 0), ]), 3)))+
  annotate(geom="text", x=1.41, y=0.00045,size = 4, label=paste0( nrow(SELECT_biopsy[which(SELECT_biopsy$DRE_binary == 1), ])))+
  annotate(geom="text", x=0.61, y=0.00045,size = 4,label=paste0( nrow(SELECT_biopsy[which(SELECT_biopsy$DRE_binary == 0), ])))


plotdt_Race <- data.frame(x = rep('African ancestry = yes', nrow(SELECT_biopsy)),
                          label = rep('African ancestry = yes', nrow(SELECT_biopsy)),
                          y = c(SELECT_biopsy[which(SELECT_biopsy$Race == 1), ]$w,
                                SELECT_biopsy[which(SELECT_biopsy$Race == 0), ]$w),
                          m = c(rep(1, nrow(SELECT_biopsy[which(SELECT_biopsy$Race == 1), ])),
                                rep(0, nrow(SELECT_biopsy[which(SELECT_biopsy$Race == 0), ]))),
                          l = c(rep('Yes', nrow(SELECT_biopsy[which(SELECT_biopsy$Race == 1), ])),
                                rep('No', nrow(SELECT_biopsy[which(SELECT_biopsy$Race == 0), ]))))
CIL_f(SELECT_biopsy[which(SELECT_biopsy$Race  == 1), ])
CIL_f(SELECT_biopsy[which(SELECT_biopsy$Race  == 0), ])


violin_p_Race <-   ggplot(plotdt_Race, aes(x , y, group=m, fill = as.factor(m))) +
  geom_split_violin()+ 
  scale_y_continuous(limits = c(0, 0.0006), breaks = c(0, 0.0002,0.0004, 0.0006),
                     labels = c('0.0000', '0.0002','0.0004', '0.0006'))+
  theme_bw() +
  scale_x_discrete(labels=c("African ancestry = yes" = c('   No    Yes  ')))+
  xlab("") +
  ylab("") + 
  scale_fill_manual(name = "", values = c('lightgreen', 'red'))+
  theme(strip.text.x = element_text(size=12,face="bold"))+
  theme(legend.position = 'none',
        legend.title = element_text(face="bold",size=12),
        legend.text=element_text(size=12),
        legend.key.size = unit(0.9, "cm"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12)) 
violin_p_Race<- violin_p_Race + facet_wrap(~ label)+
  theme(strip.text.x = element_text(size = 14, face="bold"))+
  annotate(geom="text", x=1.4, y=0.0005,size = 5,fontface = "bold", label=paste0( format(round(CIL_f(SELECT_biopsy[which(SELECT_biopsy$Race == 1), ]),3), nsmall = 3)))+
  annotate(geom="text", x=0.6, y=0.0005,size = 5,fontface = "bold", label=paste0( round(CIL_f(SELECT_biopsy[which(SELECT_biopsy$Race == 0), ]), 3)))+
  annotate(geom="text", x=1.41, y=0.00045,size = 4, label=paste0( nrow(SELECT_biopsy[which(SELECT_biopsy$Race == 1), ])))+
  annotate(geom="text", x=0.61, y=0.00045,size = 4, label=paste0( nrow(SELECT_biopsy[which(SELECT_biopsy$Race == 0), ])))




plotdt_FH <- data.frame(x = rep('Family history = yes', nrow(SELECT_biopsy)),
                        label = rep('Family history = yes', nrow(SELECT_biopsy)),
                        y = c(SELECT_biopsy[which(SELECT_biopsy$FH == 1), ]$w,
                                SELECT_biopsy[which(SELECT_biopsy$FH == 0), ]$w),
                        m = c(rep(1, nrow(SELECT_biopsy[which(SELECT_biopsy$FH== 1), ])),
                                rep(0, nrow(SELECT_biopsy[which(SELECT_biopsy$FH == 0), ]))),
                        l = c(rep('Yes', nrow(SELECT_biopsy[which(SELECT_biopsy$FH == 1), ])),
                              rep('No', nrow(SELECT_biopsy[which(SELECT_biopsy$FH == 0), ]))))
CIL_f(SELECT_biopsy[which(SELECT_biopsy$FH  == 1), ])
CIL_f(SELECT_biopsy[which(SELECT_biopsy$FH  == 0), ])

violin_p_FH <-   ggplot(plotdt_FH, aes(x , y, group=m, fill = as.factor(m))) +
  geom_split_violin()+ 
  scale_y_continuous(limits = c(0, 0.0006), breaks = c(0, 0.0002,0.0004, 0.0006),
                     labels = c('0.0000', '0.0002','0.0004', '0.0006'))+
  theme_bw() +
  scale_x_discrete(labels=c("Family history = yes" = c('   No    Yes  ')))+
  xlab("") +
  ylab("") + 
  scale_fill_manual(name = "", values = c('lightgreen', 'red'))+
  theme(strip.text.x = element_text(size=12,face="bold"))+
  theme(legend.position = 'none',
        legend.title = element_text(face="bold",size=12),
        legend.text=element_text(size=12),
        legend.key.size = unit(0.9, "cm"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12)) 
violin_p_FH<- violin_p_FH + facet_wrap(~ label)+
  theme(strip.text.x = element_text(size = 14, face="bold"))+
  annotate(geom="text", x=1.4, y=0.0005,size = 5,fontface = "bold", label=paste0(round(CIL_f(SELECT_biopsy[which(SELECT_biopsy$FH == 1), ]),3)))+
  annotate(geom="text", x=0.6, y=0.0005, size = 5,fontface = "bold",label=paste0( round(CIL_f(SELECT_biopsy[which(SELECT_biopsy$FH == 0), ]), 3)))+
  annotate(geom="text", x=1.41, y=0.00045,size = 4,label=paste0(nrow(SELECT_biopsy[which(SELECT_biopsy$FH == 1), ])))+
  annotate(geom="text", x=0.61, y=0.00045, size = 4,label=paste0(nrow(SELECT_biopsy[which(SELECT_biopsy$FH == 0), ])))




plotdt_PB <- data.frame(x = rep('Prior negative biopsy = yes', nrow(SELECT_biopsy)),
                        label = rep('Prior negative biopsy = yes', nrow(SELECT_biopsy)),
                        y = c(SELECT_biopsy[which(SELECT_biopsy$PB == 1), ]$w,
                                SELECT_biopsy[which(SELECT_biopsy$PB == 0), ]$w),
                        m = c(rep(1, nrow(SELECT_biopsy[which(SELECT_biopsy$PB == 1), ])),
                                rep(0, nrow(SELECT_biopsy[which(SELECT_biopsy$PB == 0), ]))),
                        l = c(rep('Yes', nrow(SELECT_biopsy[which(SELECT_biopsy$PB == 1), ])),
                              rep('No', nrow(SELECT_biopsy[which(SELECT_biopsy$PB == 0), ]))))
CIL_f(SELECT_biopsy[which(SELECT_biopsy$PB  == 1), ])
CIL_f(SELECT_biopsy[which(SELECT_biopsy$PB  == 0), ])


violin_p_PB <-   ggplot(plotdt_PB, aes(x , y, group=m, fill = as.factor(m))) +
  geom_split_violin()+ 
  scale_y_continuous(limits = c(0, 0.0006), breaks = c(0, 0.0002,0.0004, 0.0006),
                     labels = c('0.0000', '0.0002','0.0004', '0.0006'))+
  theme_bw() +
  scale_x_discrete(labels=c("Prior negative biopsy = yes" = c('   No    Yes  ')))+
  xlab("") +
  ylab("") + 
  scale_fill_manual(name = "", values = c('lightgreen', 'red'))+
  theme(strip.text.x = element_text(size=12,face="bold"))+
  theme(legend.position = 'none',
        legend.title = element_text(face="bold",size=12),
        legend.text=element_text(size=12),
        legend.key.size = unit(0.9, "cm"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12)) 
violin_p_PB <- violin_p_PB + facet_wrap(~ label)+
  theme(strip.text.x = element_text(size = 14, face="bold"))+
  annotate(geom="text", x=1.4, y=0.0005,size = 5,fontface = "bold", label=paste0( round(CIL_f(SELECT_biopsy[which(SELECT_biopsy$PB == 1), ]),3)))+
  annotate(geom="text", x=0.6, y=0.0005,size = 5,fontface = "bold", label=paste0(round(CIL_f(SELECT_biopsy[which(SELECT_biopsy$PB == 0), ]), 3)))+
  annotate(geom="text", x=1.41, y=0.00045, size = 4,label=paste0( nrow(SELECT_biopsy[which(SELECT_biopsy$PB == 1), ])))+
  annotate(geom="text", x=0.61, y=0.00045,size = 4,label=paste0(nrow(SELECT_biopsy[which(SELECT_biopsy$PB == 0), ])))



SELECT_biopsy$Age75 <- ifelse(SELECT_biopsy$Age <= 75, 0, 1)
plotdt_Age75 <- data.frame(x = rep('Age > 75', nrow(SELECT_biopsy)),
                           label = rep('Age > 75', nrow(SELECT_biopsy)),
                           y = c(SELECT_biopsy[which(SELECT_biopsy$Age75 == 1), ]$w,
                                 SELECT_biopsy[which(SELECT_biopsy$Age75 == 0), ]$w),
                           m = c(rep(1, nrow(SELECT_biopsy[which(SELECT_biopsy$Age75 == 1), ])),
                                 rep(0, nrow(SELECT_biopsy[which(SELECT_biopsy$Age75 == 0), ]))),
                           l = c(rep('Yes', nrow(SELECT_biopsy[which(SELECT_biopsy$Age75 == 1), ])),
                                 rep('No', nrow(SELECT_biopsy[which(SELECT_biopsy$Age75 == 0), ]))))
CIL_f(SELECT_biopsy[which(SELECT_biopsy$Age75  == 1), ])
CIL_f(SELECT_biopsy[which(SELECT_biopsy$Age75  == 0), ])

violin_p_Age75 <-   ggplot(plotdt_Age75, aes(x , y, group=m, fill = as.factor(m))) +
  geom_split_violin()+ 
  scale_y_continuous(limits = c(0, 0.0006), breaks = c(0, 0.0002,0.0004, 0.0006),
                     labels = c('0.0000', '0.0002','0.0004', '0.0006'))+
  theme_bw() +
  scale_x_discrete(labels=c("Age > 75" = c('   No    Yes  ')))+
  xlab("") +
  ylab("") + 
  scale_fill_manual(name = "", values = c('lightgreen', 'red'))+
  theme(strip.text.x = element_text(size=12,face="bold"))+
  theme(legend.position = 'none',
        legend.title = element_text(face="bold",size=12),
        legend.text=element_text(size=12),
        legend.key.size = unit(0.9, "cm"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12)) 
violin_p_Age75 <- violin_p_Age75 + facet_wrap(~ label)+
  theme(strip.text.x = element_text(size = 14, face="bold"))+
  annotate(geom="text", x=1.4, y=0.0005,size = 5,fontface = "bold",label=paste0( round(CIL_f(SELECT_biopsy[which(SELECT_biopsy$Age75 == 1), ]),3)))+
  annotate(geom="text", x=0.6, y=0.0005,size = 5,fontface = "bold",label=paste0(round(CIL_f(SELECT_biopsy[which(SELECT_biopsy$Age75 == 0), ]), 3)))+
  annotate(geom="text", x=1.41, y=0.00045,size = 4,label=paste0( nrow(SELECT_biopsy[which(SELECT_biopsy$Age75 == 1), ])))+
  annotate(geom="text", x=0.61, y=0.00045,size = 4,label=paste0(nrow(SELECT_biopsy[which(SELECT_biopsy$Age75 == 0), ])))



library(gridExtra)

p <- grid.arrange(violin_p_PSA4flag, violin_p_DRE,
                  violin_p_PB, violin_p_FH,
                  violin_p_Age75, violin_p_Race, ncol=2)
ggsave(filename="/Users/yiyaochen/Binary_SELECT_PLCO/Weights_CIL_by_factors_no_interact_1.png", p, width = 8, height = 10.2, dpi = 500)


plotdt <- rbind(plotdt_PSA4, plotdt_Race, plotdt_FH, plotdt_PB, plotdt_Age75, plotdt_DRE )


violin_p <-   ggplot(plotdt , aes(x , y, group=m, fill = as.factor(m))) +
    geom_split_violin()+ 
    scale_y_continuous(limits = c(0, 0.0006), breaks = c(0, 0.00025,0.0005, 0.0006))+
    theme_bw() +
    xlab("Levels of risk factor") +
    ylab("Density of weights") + 
    scale_fill_manual(name = "", values = c('lightgreen', 'red'))+
    theme(strip.text.x = element_text(size=12,face="bold"))+
    theme(legend.position = 'none',
        legend.title = element_text(face="bold",size=12),
        legend.text=element_text(size=12),
        legend.key.size = unit(0.9, "cm"),
        axis.text.x = element_text(face="bold",size=12),
        axis.text.y = element_text(face="bold",size=12),
        axis.title.x = element_text(face="bold",size=14),
        axis.title.y = element_text(face="bold",size=14))+
      facet_wrap(~ label, ncol = 2, scales = 'free')
ggsave(filename="/Users/yiyaochen/Binary_SELECT_PLCO/Weights_CIL_no_interact_2.png", violin_p, width = 8, height = 10.2, dpi = 500)

  
plotdt$l <-  factor(plotdt$l, levels = c( 'No', 'Yes'))
box_p<-   ggplot(plotdt , aes(l, y, group=m, fill = as.factor(m))) +
    geom_boxplot()+ 
    scale_y_continuous(limits = c(0, 0.0009), breaks = c(0, 0.00025,0.0005, 0.00075, 0.0009))+
    theme_bw() +
    facet_grid(. ~ x, scales="free_x", labeller= as_labeller(c("PSA4flag" = 'PSA > 4 ng/mL', 
                                                               "Race" = 'African ancestry',
                                                               "FH" = 'Family history', 
                                                               "PB" =  'Prior negative biopsy'))) +
    xlab("") +
    ylab("Weights") + 
    scale_fill_manual(name = "", values = c('lightgreen', 'red'))+
    theme(strip.text.x = element_text(size=12,face="bold"))+
    theme(legend.position = 'none',
          legend.title = element_text(face="bold",size=12),
          legend.text=element_text(size=12),
          legend.key.size = unit(0.9, "cm"),
          axis.text.x = element_text(face="bold",size=12),
          axis.text.y = element_text(face="bold",size=12),
          axis.title.x = element_text(face="bold",size=14),
          axis.title.y = element_text(face="bold",size=14))
  
ggsave(filename="/Users/yiyaochen/Binary_SELECT_PLCO/box_plot_weights_by_subgroup.png", box_p, width = 10, height = 6, dpi = 500)

                    