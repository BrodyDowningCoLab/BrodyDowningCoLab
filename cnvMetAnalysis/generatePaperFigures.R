library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggsci)

cnv_comb_auc1 <- readRDS("~/R/colab/BrodyDowningCoLab/cnvMetAnalysis/cnv_comb_auc1.rds")
cnv_comb_auc1 <- data.frame(cnv_comb_auc1)
cnv_comb_auc1 <- unique(cnv_comb_auc1)

temp <- cbind(cnv_comb_auc1[[1]], cnv_comb_auc1[[2]], round(cnv_comb_auc1[[3]], digits = 3))

result <- temp

temp <- cbind(cnv_comb_auc1[[1]], cnv_comb_auc1[[2]], round(cnv_comb_auc1[[4]], digits = 3))
result <- rbind(result, temp)
temp <- cbind(cnv_comb_auc1[[1]], cnv_comb_auc1[[2]], round(cnv_comb_auc1[[5]], digits = 3))
result <- rbind(result, temp)
temp <- cbind(cnv_comb_auc1[[1]], cnv_comb_auc1[[2]], round(cnv_comb_auc1[[6]], digits = 3))
result <- rbind(result, temp)
temp <- cbind(cnv_comb_auc1[[1]], cnv_comb_auc1[[2]], round(cnv_comb_auc1[[7]], digits = 3))
result <- rbind(result, temp)

result <- data.frame(result)

result$X3 <- as.numeric(result$X3)

cnvNames <- unique(result["X1"])
modelNames <- unique(result["X2"])

result <- result %>% 
  rename(
    CNV = X1,
    Model_Type = X2,
    AUC = X3
  )

####################
## Generate plots ##
####################
gp_cleanFull <- ggplot(subset(result, Model_Type %in% c("clean", "full")),
             aes(x = CNV, y = AUC,  colour = interaction(Model_Type, CNV), group = CNV)) + facet_wrap( ~ Model_Type) +
  geom_boxplot( fill="white", notch=FALSE) +
  geom_jitter( size=1, color="black", width=0.3) +
  #scale_y_discrete(breaks = c(0.50, 0.51,0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64)) +
  ggtitle("Comparison of Clean Versus Full CNV COAD Predictions") + 
  xlab("CNV Location") + ylab("AUC of Predictions") +
  theme(plot.title = element_text(size = 18, face = "bold"), 
        legend.position = "none",
        axis.text.x = element_text(angle = 90)) 

gp_cleanFull

gp_clean1up1down <- ggplot(subset(result, Model_Type %in% c("clean", "1 CNV upstream", "1 CNV downstream")),
                       aes(x = CNV, y = AUC,  colour = interaction(Model_Type, CNV), group = CNV)) + facet_wrap( ~ Model_Type) +
  geom_boxplot( fill="white", notch=FALSE) +
  geom_jitter( size=1, color="black", width=0.3) +
  #scale_y_discrete(breaks = c(0.50, 0.51,0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64)) +
  ggtitle("Comparison of Clean Versus Upstream and Downstream Methylation Beta Values COAD Predictions") + 
  xlab("CNV Location") + ylab("AUC of Predictions") +
  theme(plot.title = element_text(size = 18, face = "bold"), 
        legend.position = "none",
        axis.text.x = element_text(angle = 90))

gp_clean1up1down

gp_clean1up2up <- ggplot(subset(result, Model_Type %in% c("1 CNV upstream", "2 CNVs upstream")),
                           aes(x = CNV, y = AUC,  colour = interaction(Model_Type, CNV), group = CNV)) + facet_wrap( ~ Model_Type) +
  geom_boxplot( fill="white", notch=FALSE) +
  geom_jitter( size=1, color="black", width=0.3) +
  #scale_y_discrete(breaks = c(0.50, 0.51,0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64)) +
  ggtitle("Comparison of Clean Versus Upstream and Downstream Methylation Beta Values COAD Predictions") + 
  xlab("CNV Location") + ylab("AUC of Predictions") +
  theme(plot.title = element_text(size = 18, face = "bold"), 
        legend.position = "none",
        axis.text.x = element_text(angle = 90))

gp_clean1up2up

gp_clean1down2down <- ggplot(subset(result, Model_Type %in% c("1 CNV downstream", "2 CNVs downstream")),
                         aes(x = CNV, y = AUC,  colour = interaction(Model_Type, CNV), group = CNV)) + facet_wrap( ~ Model_Type) +
  geom_boxplot( fill="white", notch=FALSE) +
  geom_jitter( size=1, color="black", width=0.3) +
  #scale_y_discrete(breaks = c(0.50, 0.51,0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64)) +
  ggtitle("Comparison of Clean Versus Upstream and Downstream Methylation Beta Values COAD Predictions") + 
  xlab("CNV Location") + ylab("AUC of Predictions") +
  theme(plot.title = element_text(size = 18, face = "bold"), 
        legend.position = "none",
        axis.text.x = element_text(angle = 90))

gp_clean1down2down

gp_clean1up1down2up2down <- ggplot(subset(result, Model_Type %in% c("1 CNV downstream", "2 CNVs downstream", "1 CNV upstream", "2 CNVs upstream")),
                             aes(x = CNV, y = AUC,  colour = interaction(Model_Type, CNV), group = CNV)) + facet_wrap( ~ Model_Type) +
  geom_boxplot( fill="white", notch=FALSE) +
  geom_jitter( size=1, color="black", width=0.3) +
  #scale_y_discrete(breaks = c(0.50, 0.51,0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64)) +
  ggtitle("Comparison of Clean Versus Upstream and Downstream Methylation Beta Values COAD Predictions") + 
  xlab("CNV Location") + ylab("AUC of Predictions") +
  theme(plot.title = element_text(size = 18, face = "bold"), 
        legend.position = "none",
        axis.text.x = element_text(angle = 90))

gp_clean1up1down2up2down