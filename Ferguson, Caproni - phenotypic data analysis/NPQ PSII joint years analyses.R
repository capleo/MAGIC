
library(corrplot)
library(ggplot2)
library(forcats)
library(gridExtra)

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

#### Hertiability bar plot ####

CF_2022_H2 <- read.csv(
  "~/My Drive/Maize 2022 analyses/Chlorophyll Fluorescence 2022/chlFlu_2022_H2.csv",
  header = T, stringsAsFactors = F)
CF_2022_H2$Year <- as.factor(CF_2022_H2$Year)

CF_2021_H2 <- read.csv(
  "~/My Drive/Maize 2021 analyses/Chlorophyll Fluorescence 2021/chlFlu_2021_H2.csv", 
  header = T, stringsAsFactors = F
)
CF_2021_H2$Year <- as.factor(CF_2021_H2$Year)

CF_H2_bothYears <- rbind(CF_2021_H2, CF_2022_H2)

CF_joint_H2 <- read.csv(
  "~/My Drive/Maize joint years analyses/Chlorophyll fluorescence/chlFlu_joint_H2.csv",
  header = T, stringsAsFactors = F
)

CF_H2_allModel <- rbind(CF_H2_bothYears, CF_joint_H2)

ggplot(data=CF_H2_allModel, aes(x=fct_reorder(Trait,Hertiability), y=Hertiability, fill=Year)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_fill_manual(values=c('darkgoldenrod1','cornflowerblue','aquamarine')) + 
  scale_y_continuous(limits = c(0,0.75), breaks = seq(0,0.75,0.15)) + 
  theme_classic() 

rm(CF_2021_H2, CF_2022_H2, CF_H2_bothYears, CF_H2_allModel, CF_joint_H2)

mean(CF_2021_H2$Hertiability)
mean(CF_2022_H2$Hertiability)
mean(CF_joint_H2$Hertiability)


#### Corplots individual year basis ####

chlFlu_2022_BLUPs <- read.csv(
  "~/My Drive/Maize 2022 analyses/Chlorophyll Fluorescence 2022/chlFlu_2022_BLUPs.csv",
  header = T, stringsAsFactors = F)

cor.mat_2022 <- cor(chlFlu_2022_BLUPs[,2:ncol(chlFlu_2022_BLUPs)])
p.mat_2022 <- cor.mtest(chlFlu_2022_BLUPs[,2:ncol(chlFlu_2022_BLUPs)])

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(cor.mat_2022, method="color", col=col(200),  
         type="upper", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat_2022, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

chlFlu_2021_BLUPs <- read.csv(
  "~/My Drive/Maize 2021 analyses/Chlorophyll Fluorescence 2021/chlFlu_2021_BLUPs.csv",
  header = T, stringsAsFactors = F)

chlFlu_2021_BLUPs <- na.omit(chlFlu_2021_BLUPs) 

cor.mat_2021 <- cor(chlFlu_2021_BLUPs[,2:ncol(chlFlu_2021_BLUPs)])
p.mat_2021 <- cor.mtest(chlFlu_2021_BLUPs[,2:ncol(chlFlu_2021_BLUPs)])

corrplot(cor.mat_2021, method="color", col=col(200),  
         type="lower", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat_2021, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
) # export as 8x8 landscape

chlFlu_joint_BLUPs <- read.csv(
  "~/My Drive/Maize joint years analyses/Chlorophyll fluorescence/chlFlu_BLUPs_jointModel.csv",
  header = T, stringsAsFactors = F
)

cor.mat_joint <- cor(chlFlu_joint_BLUPs[,2:ncol(chlFlu_joint_BLUPs)])
p.mat_joint <- cor.mtest(chlFlu_joint_BLUPs[,2:ncol(chlFlu_joint_BLUPs)])

rm(cor.mat_2021, cor.mat_2022, p.mat_2021, p.mat_2022, cor.mat_joint, p.mat_joint)

#### Boxplots to show year effects on traits ####

chlFlu_2021_BLUPs$Year <- "2021"
chlFlu_2022_BLUPs$Year <- "2022"

chlFlu_2021_BLUPs$Line <- gsub("_", "_00", chlFlu_2021_BLUPs$Line)

identical(names(chlFlu_2021_BLUPs), names(chlFlu_2022_BLUPs))

chlFlu_bothYears_BLUPs <- rbind(chlFlu_2021_BLUPs, chlFlu_2022_BLUPs)

NPQ_ind_amp_boxplot <- ggplot(chlFlu_bothYears_BLUPs,
                         aes(x=Year, y=NPQ.induction.amplitude, fill=Year)) +
  geom_boxplot() + 
  ylab("NPQ induction (amplitude)") + 
  scale_fill_manual(values=c('#E69F00','#999999')) + 
  scale_y_continuous(limits = c(1,4), breaks = seq(1,4,0.5)) + 
  theme_classic() + 
  theme(legend.position = "none")

NPQ_ind_rate_boxplot <- ggplot(chlFlu_bothYears_BLUPs,
                              aes(x=Year, y=NPQ.induction.rate, fill=Year)) +
  geom_boxplot() + 
  ylab("NPQ induction (rate)") + 
  scale_fill_manual(values=c('#E69F00','#999999')) + 
  scale_y_continuous(limits = c(0,0.03), breaks = seq(0,0.03,0.01)) + 
  theme_classic() + 
  theme(legend.position = "none")

NPQ_rel_amp_boxplot <- ggplot(chlFlu_bothYears_BLUPs,
                               aes(x=Year, y=NPQ.relaxation.amplitude, fill=Year)) +
  geom_boxplot() + 
  ylab("NPQ relaxation (amplitude)") + 
  scale_fill_manual(values=c('#E69F00','#999999')) + 
  scale_y_continuous(limits = c(0,4), breaks = seq(0,4,1)) + 
  theme_classic() + 
  theme(legend.position = "none")

NPQ_rel_rate_boxplot <- ggplot(chlFlu_bothYears_BLUPs,
                              aes(x=Year, y=NPQ.relaxation.rate, fill=Year)) +
  geom_boxplot() + 
  ylab("NPQ relaxation (rate)") + 
  scale_fill_manual(values=c('#E69F00','#999999')) + 
  scale_y_continuous(limits = c(0,0.05), breaks = seq(0,0.05,0.01)) + 
  theme_classic() + 
  theme(legend.position = "none")

t.test(NPQ.relaxation.rate~Year, data = chlFlu_bothYears_BLUPs)


NPQ_rel_residual_boxplot <- ggplot(chlFlu_bothYears_BLUPs,
                               aes(x=Year, y=NPQ.relaxation.residual, fill=Year)) +
  geom_boxplot() + 
  ylab("NPQ relaxation (residual)") + 
  scale_fill_manual(values=c('#E69F00','#999999')) + 
  scale_y_continuous(limits = c(0,0.4), breaks = seq(0,0.4,0.1)) + 
  theme_classic() + 
  theme(legend.position = "none")

phiPSII_induction_amplitude_boxplot <- ggplot(chlFlu_bothYears_BLUPs,
                                   aes(x=Year, y=phiPSII.induction.amplitude, fill=Year)) +
  geom_boxplot() + 
  ylab("phiPSII induction (amplitude)") + 
  scale_fill_manual(values=c('#E69F00','#999999')) + 
  scale_y_continuous(limits = c(0,0.4), breaks = seq(0,0.4,0.1)) + 
  theme_classic() + 
  theme(legend.position = "none")

phiPSII_induction_rate_boxplot <- ggplot(chlFlu_bothYears_BLUPs,
                                              aes(x=Year, y=phiPSII.induction.rate, fill=Year)) +
  geom_boxplot() + 
  ylab("phiPSII induction (rate)") + 
  scale_fill_manual(values=c('#E69F00','#999999')) + 
  scale_y_continuous(limits = c(0,0.04), breaks = seq(0,0.04,0.01)) + 
  theme_classic() + 
  theme(legend.position = "none")

phiPSII_induction_residual_boxplot <- ggplot(chlFlu_bothYears_BLUPs,
                                         aes(x=Year, y=phiPSII.induction.residual, fill=Year)) +
  geom_boxplot() + 
  ylab("phiPSII induction (rate)") + 
  scale_fill_manual(values=c('#E69F00','#999999')) + 
  scale_y_continuous(limits = c(0.2,0.7), breaks = seq(0.2,0.7,0.1)) + 
  theme_classic() + 
  theme(legend.position = "none")

NPQ_induction_LM_boxplot <- ggplot(chlFlu_bothYears_BLUPs,
                                             aes(x=Year, y=NPQ.induction.linear, fill=Year)) +
  geom_boxplot() + 
  ylab("NPQ induction (LM slope)") + 
  scale_fill_manual(values=c('#E69F00','#999999')) + 
  scale_y_continuous(limits = c(0,0.05), breaks = seq(0,0.05,0.01)) + 
  theme_classic() + 
  theme(legend.position = "none")

maxNPQ_boxplot <- ggplot(chlFlu_bothYears_BLUPs,
                                   aes(x=Year, y=Maximum.NPQ, fill=Year)) +
  geom_boxplot() + 
  ylab("Maximum NPQ") + 
  scale_fill_manual(values=c('#E69F00','#999999')) + 
  scale_y_continuous(limits = c(0,4), breaks = seq(0,4,1)) + 
  theme_classic() + 
  theme(legend.position = "none")

finalNPQ_boxplot <- ggplot(chlFlu_bothYears_BLUPs,
                         aes(x=Year, y=Final.NPQ, fill=Year)) +
  geom_boxplot() + 
  ylab("Final NPQ") + 
  scale_fill_manual(values=c('#E69F00','#999999')) + 
  scale_y_continuous(limits = c(0,0.4), breaks = seq(0,0.4,0.1)) + 
  theme_classic() + 
  theme(legend.position = "none")

finalFvFm_boxplot <- ggplot(chlFlu_bothYears_BLUPs,
                           aes(x=Year, y=Final.FvFm, fill=Year)) +
  geom_boxplot() + 
  ylab("Final FvFm") + 
  scale_fill_manual(values=c('#E69F00','#999999')) + 
  scale_y_continuous(limits = c(0.6,0.85), breaks = seq(0.6,0.85,0.05)) + 
  theme_classic() + 
  theme(legend.position = "none")

initialFvFm_boxplot <- ggplot(chlFlu_bothYears_BLUPs,
                            aes(x=Year, y=Initial.FvFm, fill=Year)) +
  geom_boxplot() + 
  ylab("Initial FvFM") + 
  scale_fill_manual(values=c('#E69F00','#999999')) + 
  scale_y_continuous(limits = c(0.65,0.9), breaks = seq(0.65,0.9,0.05)) + 
  theme_classic() + 
  theme(legend.position = "none")

grid.arrange(initialFvFm_boxplot, finalFvFm_boxplot, 
             NPQ_induction_LM_boxplot,
             NPQ_ind_amp_boxplot, NPQ_ind_rate_boxplot, 
             maxNPQ_boxplot,
             NPQ_rel_amp_boxplot, NPQ_rel_rate_boxplot, NPQ_rel_residual_boxplot,
             finalNPQ_boxplot,
             phiPSII_induction_amplitude_boxplot, phiPSII_induction_rate_boxplot, phiPSII_induction_residual_boxplot,
             ncol = 5) 

#rm(initialFvFm_boxplot, finalFvFm_boxplot, NPQ_induction_LM_boxplot, NPQ_ind_amp_boxplot,
#   NPQ_ind_rate_boxplot, maxNPQ_boxplot, NPQ_rel_amp_boxplot, NPQ_rel_rate_boxplot, 
#   NPQ_rel_residual_boxplot, finalNPQ_boxplot, phiPSII_induction_amplitude_boxplot, 
#   phiPSII_induction_rate_boxplot, phiPSII_induction_residual_boxplot)


#### Year-on-year scatter plots ####

colnames(chlFlu_2022_BLUPs) <- paste(colnames(chlFlu_2022_BLUPs),"2022",sep="_")
names(chlFlu_2022_BLUPs)[names(chlFlu_2022_BLUPs) == 'Line_2022'] <- 'Line'
colnames(chlFlu_2021_BLUPs) <- paste(colnames(chlFlu_2021_BLUPs),"2021",sep="_")
names(chlFlu_2021_BLUPs)[names(chlFlu_2021_BLUPs) == 'Line_2021'] <- 'Line'

chlFlu_bothYears_BLUPs <- merge(chlFlu_2021_BLUPs, chlFlu_2022_BLUPs, by = "Line", all = T)

#write.csv(chlFlu_bothYears_BLUPs, 
#          "~/My Drive/Maize joint years analyses/Chlorophyll fluorescence/ChlFluor_BLUPs_bothYears.csv",
#          row.names = F)

jointModel_BLUPs <- read.csv("~/My Drive/Maize joint years analyses/Chlorophyll fluorescence/chlFlu_BLUPs_jointModel.csv", 
                             header = T, stringsAsFactors = F)
colnames(jointModel_BLUPs) <- paste(colnames(jointModel_BLUPs),"Joint",sep="_")
names(jointModel_BLUPs)[names(jointModel_BLUPs) == 'Line_Joint'] <- 'Line'

chlFlu_bothYears_JointModel_BLUPs <- merge(chlFlu_bothYears_BLUPs, jointModel_BLUPs, by = "Line", all = T)
chlFlu_bothYears_JointModel_BLUPs <- 
  chlFlu_bothYears_JointModel_BLUPs[,order(colnames(chlFlu_bothYears_JointModel_BLUPs))]

#write.csv(chlFlu_bothYears_JointModel_BLUPs, 
#          "~/My Drive/Maize joint years analyses/Chlorophyll fluorescence/CF_BLUPs_allModels.csv",
#          row.names = F)

initialFvFm_scatter <- ggplot(chlFlu_bothYears_BLUPs, 
                              aes(x=Initial.FvFm_2021, y=Initial.FvFm_2022)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0.7,0.85), breaks = seq(0.7,0.85,0.05)) + 
  scale_y_continuous(limits = c(0.7,0.85), breaks = seq(0.7,0.85,0.05)) + 
  geom_abline(slope=1, intercept=0) + 
  ylab("2022 Initial FvFm") + 
  xlab("2021 Initial FvFm") + 
  theme_replace() 

finalFvFm_scatter <- ggplot(chlFlu_bothYears_BLUPs, 
                              aes(x=Final.FvFm_2021, y=Final.FvFm_2022)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0.65,0.8), breaks = seq(0.65,0.8,0.05)) + 
  scale_y_continuous(limits = c(0.65,0.8), breaks = seq(0.65,0.8,0.05)) + 
  geom_abline(slope=1, intercept=0) + 
  ylab("2022 Final FvFm") + 
  xlab("2021 Final FvFm") +  
  theme_replace() 

NPQ_linearInduction_scatter <- ggplot(chlFlu_bothYears_BLUPs, 
                            aes(x=NPQ.induction.linear_2021, y=NPQ.induction.linear_2022)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0,0.05), breaks = seq(0,0.05,0.02)) + 
  scale_y_continuous(limits = c(0,0.05), breaks = seq(0,0.05,0.02)) + 
  geom_abline(slope=1, intercept=0) + 
  ylab("2022 NPQ induction (LM slope)") + 
  xlab("2021 NPQ induction (LM slope)") +  
  theme_replace() 

NPQ_ind_amp_scatter <- ggplot(chlFlu_bothYears_BLUPs, 
                                      aes(x=NPQ.induction.amplitude_2021, y=NPQ.induction.amplitude_2022)) + 
  geom_point() + 
  scale_x_continuous(limits = c(1.5,4), breaks = seq(1.5,4,0.8)) + 
  scale_y_continuous(limits = c(1.5,4), breaks = seq(1.5,4,0.8)) + 
  geom_abline(slope=1, intercept=0) + 
  ylab("2022 NPQ induction (amplitude)") + 
  xlab("2021 NPQ induction (amplitude)") +  
  theme_replace() 

NPQ_ind_rate_scatter <- ggplot(chlFlu_bothYears_BLUPs, 
                              aes(x=NPQ.induction.rate_2021, y=NPQ.induction.rate_2022)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0,0.025), breaks = seq(0,0.03,0.008)) + 
  scale_y_continuous(limits = c(0,0.025), breaks = seq(0,0.03,0.008)) + 
  geom_abline(slope=1, intercept=0) + 
  ylab("2022 NPQ induction (rate)") + 
  xlab("2021 NPQ induction (rate)") +  
  theme_replace() 

maxNPQ_scatter <- ggplot(chlFlu_bothYears_BLUPs, 
                         aes(x=Maximum.NPQ_2021, y=Maximum.NPQ_2022)) + 
  geom_point() + 
  scale_x_continuous(limits = c(1.5,4), breaks = seq(0,5,0.8)) + 
  scale_y_continuous(limits = c(1.5,4), breaks = seq(0,5,0.8)) + 
  geom_abline(slope=1, intercept=0) + 
  ylab("2022 Maximum NPQ") + 
  xlab("2021 Maximum NPQ") +   
  theme_replace() 

NPQ_rel_amp_scatter <- ggplot(chlFlu_bothYears_BLUPs, 
                         aes(x=NPQ.relaxation.amplitude_2021, y=NPQ.relaxation.amplitude_2022)) + 
  geom_point() + 
  scale_x_continuous(limits = c(1.4,3.6), breaks = seq(1,3.5,0.6)) + 
  scale_y_continuous(limits = c(1.4,3.6), breaks = seq(1,3.5,0.6)) + 
  geom_abline(slope=1, intercept=0) + 
  ylab("2022 NPQ relaxation (amplitude)") + 
  xlab("2021 NPQ relaxation (amplitude)") +    
  theme_replace() 

NPQ_rel_rate_scatter <- ggplot(chlFlu_bothYears_BLUPs, 
                              aes(x=NPQ.relaxation.rate_2021, y=NPQ.relaxation.rate_2022)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0.01,0.045), breaks = seq(0,1,0.009)) + 
  scale_y_continuous(limits = c(0.01,0.045), breaks = seq(0,1,0.009)) + 
  geom_abline(slope=1, intercept=0) + 
  ylab("2022 NPQ relaxation (rate)") + 
  xlab("2021 NPQ relaxation (reate)") +    
  theme_replace() 

NPQ_rel_residual_scatter <- ggplot(chlFlu_bothYears_BLUPs, 
                               aes(x=NPQ.relaxation.residual_2021, y=NPQ.relaxation.residual_2022)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0.15,0.36), breaks = seq(0,1,0.06)) + 
  scale_y_continuous(limits = c(0.15,0.36), breaks = seq(0,1,0.06)) + 
  geom_abline(slope=1, intercept=0) + 
  ylab("2022 NPQ relaxation (residual)") + 
  xlab("2021 NPQ relaxation (residual)") +    
  theme_replace() 

finalNPQ_scatter <- ggplot(chlFlu_bothYears_BLUPs, 
                                   aes(x=Final.NPQ_2021, y=Final.NPQ_2022)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0.1,0.35), breaks = seq(0.1,0.35,0.08)) + 
  scale_y_continuous(limits = c(0.1,0.35), breaks = seq(0.1,0.35,0.08)) + 
  geom_abline(slope=1, intercept=0) + 
  ylab("2022 Final NPQ") + 
  xlab("2021 Final NPQ") +    
  theme_replace() 

phiPSII_ind_amp_scatter <- ggplot(chlFlu_bothYears_BLUPs, 
                           aes(x=phiPSII.induction.amplitude_2021, y=phiPSII.induction.amplitude_2022)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0.15,0.35), breaks = seq(0.1,0.35,0.05)) + 
  scale_y_continuous(limits = c(0.15,0.35), breaks = seq(0.1,0.35,0.05)) + 
  geom_abline(slope=1, intercept=0) + 
  ylab("2022 phiPSII induction (amplitude)") + 
  xlab("2021 phiPSII induction (amplitude)") +    
  theme_replace() 

phiPSII_ind_rate_scatter <- ggplot(chlFlu_bothYears_BLUPs, 
                                  aes(x=phiPSII.induction.rate_2021, y=phiPSII.induction.rate_2022)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0.005,0.035), breaks = seq(0,0.35,0.008)) + 
  scale_y_continuous(limits = c(0.005,0.035), breaks = seq(0,0.35,0.008)) + 
  geom_abline(slope=1, intercept=0) + 
  ylab("2022 phiPSII rate (rate)") + 
  xlab("2021 phiPSII rate (rate)") +    
  theme_replace() 

phiPSII_ind_residual_scatter <- ggplot(chlFlu_bothYears_BLUPs, 
                                   aes(x=phiPSII.induction.residual_2021, y=phiPSII.induction.residual_2022)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0.35,0.6), breaks = seq(0,1,0.07)) + 
  scale_y_continuous(limits = c(0.35,0.6), breaks = seq(0,1,0.07)) + 
  geom_abline(slope=1, intercept=0) + 
  ylab("2022 phiPSII rate (residual)") + 
  xlab("2021 phiPSII rate (residual)") +    
  theme_replace() 

grid.arrange(initialFvFm_scatter, finalFvFm_scatter, 
             NPQ_linearInduction_scatter,
             NPQ_ind_amp_scatter, NPQ_ind_rate_scatter, 
             maxNPQ_scatter,
             NPQ_rel_amp_scatter, NPQ_rel_rate_scatter, NPQ_rel_residual_scatter,
             finalNPQ_scatter,
             phiPSII_ind_amp_scatter, phiPSII_ind_rate_scatter, phiPSII_ind_residual_scatter,
             ncol = 5) # export as landscape at 8 x 14

grid.arrange(initialFvFm_boxplot, NPQ_ind_rate_boxplot, maxNPQ_boxplot, NPQ_rel_rate_boxplot, phiPSII_induction_rate_boxplot,
             initialFvFm_scatter, NPQ_ind_rate_scatter, maxNPQ_scatter, NPQ_rel_rate_scatter, phiPSII_ind_rate_scatter,
             nrow = 2)

grid.arrange(finalFvFm_boxplot, NPQ_induction_LM_boxplot, NPQ_ind_amp_boxplot, NPQ_rel_amp_boxplot, NPQ_rel_residual_boxplot, finalNPQ_boxplot, phiPSII_induction_amplitude_boxplot, phiPSII_induction_residual_boxplot, 
             finalFvFm_scatter, NPQ_linearInduction_scatter, NPQ_ind_amp_scatter, NPQ_rel_amp_scatter, NPQ_rel_residual_scatter, finalNPQ_scatter, phiPSII_ind_amp_scatter, phiPSII_ind_residual_scatter,
             nrow = 2)

rm(initialFvFm_scatter, finalFvFm_scatter, NPQ_linearInduction_scatter, NPQ_ind_amp_scatter,
   NPQ_ind_rate_scatter, maxNPQ_scatter, NPQ_rel_amp_scatter, NPQ_rel_rate_scatter,
   NPQ_rel_residual_scatter, finalNPQ_scatter, phiPSII_ind_amp_scatter, phiPSII_ind_rate_scatter,
   phiPSII_ind_residual_scatter) 

cor.test(chlFlu_bothYears_BLUPs$Initial.FvFm_2021, 
         chlFlu_bothYears_BLUPs$Initial.FvFm_2022,
         method = "pearson")

