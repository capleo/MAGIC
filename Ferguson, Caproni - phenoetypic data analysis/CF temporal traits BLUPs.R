
library(lme4)
library(plyr)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(gridExtra)

temporal_data_21 <- read.csv(
  "~/My Drive/Maize 2021 analyses/Chlorophyll Fluorescence 2021/temporal_NPQ_phiPSII.csv", 
  header = T, stringsAsFactors = F)

temporal_data_21 <- subset(temporal_data_21, select = -c(PI))

temporal_data_21 <- aggregate(cbind(NPQ_L1, NPQ_L2, NPQ_L3, NPQ_L4, NPQ_L5, NPQ_L6, NPQ_L7, NPQ_L8, NPQ_L9, NPQ_L10,
                                    NPQ_L11, NPQ_Lss, NPQ_D1, NPQ_D2, NPQ_D3, NPQ_D4, NPQ_D5, NPQ_D6, NPQ_D7, NPQ_D8,
                                    Fv.Fm_L1, Fv.Fm_L2, Fv.Fm_L3, Fv.Fm_L4, Fv.Fm_L5, Fv.Fm_L6, Fv.Fm_L7, Fv.Fm_L8, Fv.Fm_L9,
                                    Fv.Fm_L10, Fv.Fm_L11, Fv.Fm_Lss, Fv.Fm_D1, Fv.Fm_D2, Fv.Fm_D3, Fv.Fm_D4, Fv.Fm_D5, 
                                    Fv.Fm_D6, Fv.Fm_D7, Fv.Fm_D8
) ~ PlotID, data = temporal_data_21, FUN = mean, na.rm = TRUE)

maize_ref_2021 <- read.csv(
  "~/My Drive/Maize 2021 analyses/Maize_reference.csv",
  header = T, stringsAsFactors = F)

temporal_data_21 <- merge(maize_ref_2021, temporal_data_21, by = "PlotID", all = T)
rm(maize_ref_2021)
temporal_data_21 <- na.omit(temporal_data_21)

temporal_data_21$Year <- 2021

temporal_data_21$PlotID <- as.factor(temporal_data_21$PlotID)
temporal_data_21$Block <- as.factor(temporal_data_21$Block)
temporal_data_21$Replicate <- as.factor(temporal_data_21$Replicate)
temporal_data_21$Accession <- as.factor(temporal_data_21$Accession)
temporal_data_21$Accession <-gsub("_", "_00", temporal_data_21$Accession)
temporal_data_21$Measurement_day <- as.factor(temporal_data_21$Measurement_day)
temporal_data_21$Measurement_silking_interval <- as.factor(temporal_data_21$Measurement_silking_interval)
temporal_data_21$Year <- as.factor(temporal_data_21$Year)
temporal_data_21 <- droplevels(temporal_data_21)

temporal_data_22 <- read.csv(
  "~/My Drive/Maize 2022 analyses/Chlorophyll Fluorescence 2022/temporal_NPQ_phiPSII.csv", 
  header = T, stringsAsFactors = F)

temporal_data_22 <- subset(temporal_data_22, select = -c(PI))

temporal_data_22 <- aggregate(cbind(NPQ_L1, NPQ_L2, NPQ_L3, NPQ_L4, NPQ_L5, NPQ_L6, NPQ_L7, NPQ_L8, NPQ_L9, NPQ_L10,
                                    NPQ_L11, NPQ_Lss, NPQ_D1, NPQ_D2, NPQ_D3, NPQ_D4, NPQ_D5, NPQ_D6, NPQ_D7, NPQ_D8,
                                    Fv.Fm_L1, Fv.Fm_L2, Fv.Fm_L3, Fv.Fm_L4, Fv.Fm_L5, Fv.Fm_L6, Fv.Fm_L7, Fv.Fm_L8, Fv.Fm_L9,
                                    Fv.Fm_L10, Fv.Fm_L11, Fv.Fm_Lss, Fv.Fm_D1, Fv.Fm_D2, Fv.Fm_D3, Fv.Fm_D4, Fv.Fm_D5, 
                                    Fv.Fm_D6, Fv.Fm_D7, Fv.Fm_D8
) ~ PlotID, data = temporal_data_22, FUN = mean, na.rm = TRUE)

maize_ref_2022 <- read.csv(
  "~/My Drive/Maize 2022 analyses/Maize_2022_ref.csv",
  header = T, stringsAsFactors = F)

temporal_data_22 <- merge(maize_ref_2022, temporal_data_22, by = "PlotID", all = T)
rm(maize_ref_2022)
temporal_data_22 <- na.omit(temporal_data_22)

temporal_data_22$Year <- 2022

temporal_data_22$PlotID <- as.factor(temporal_data_22$PlotID)
temporal_data_22$Block <- as.factor(temporal_data_22$Block)
temporal_data_22$Replicate <- as.factor(temporal_data_22$Replicate)
names(temporal_data_22)[names(temporal_data_22) == 'Line'] <- 'Accession'
temporal_data_22$Accession <- as.factor(temporal_data_22$Accession)
names(temporal_data_22)[names(temporal_data_22) == 'Measurement_date'] <- 'Measurement_day'
temporal_data_22$Measurement_day <- as.factor(temporal_data_22$Measurement_day)
names(temporal_data_22)[names(temporal_data_22) == 'Silking_measurement_interval'] <- 'Measurement_silking_interval'
temporal_data_22$Measurement_silking_interval <- as.factor(temporal_data_22$Measurement_silking_interval)
temporal_data_22$Year <- as.factor(temporal_data_22$Year) 
temporal_data_22 <- droplevels(temporal_data_22)

colnames(temporal_data_21)
colnames(temporal_data_22)

temporal_data_22 <- 
  temporal_data_22[, c("PlotID", "Accession", "Column","Block", "Replicate", "Measurement_day", 
                       "Measurement_silking_interval", "NPQ_L1", "NPQ_L2", "NPQ_L3", "NPQ_L4", "NPQ_L5", 
                       "NPQ_L6", "NPQ_L7", "NPQ_L8", "NPQ_L9", "NPQ_L10", "NPQ_L11", "NPQ_Lss", "NPQ_D1", 
                       "NPQ_D2", "NPQ_D3", "NPQ_D4", "NPQ_D5", "NPQ_D6", "NPQ_D7", "NPQ_D8", "Fv.Fm_L1", 
                       "Fv.Fm_L2", "Fv.Fm_L3",  "Fv.Fm_L4", "Fv.Fm_L5", "Fv.Fm_L6", "Fv.Fm_L7", "Fv.Fm_L8", 
                       "Fv.Fm_L9", "Fv.Fm_L10", "Fv.Fm_L11", "Fv.Fm_Lss", "Fv.Fm_D1", "Fv.Fm_D2", "Fv.Fm_D3", 
                       "Fv.Fm_D4", "Fv.Fm_D5", "Fv.Fm_D6", "Fv.Fm_D7", "Fv.Fm_D8", "Year")]

temporal_data_21_22 <- rbind(temporal_data_21, temporal_data_22)
temporal_data_21_22$Accession <- as.factor(temporal_data_21_22$Accession)

temporal_data_21_22 <- 
  temporal_data_21_22[, c("PlotID", "Accession", "Column","Block", "Replicate", "Measurement_day", 
                       "Measurement_silking_interval", "Year", "NPQ_L1", "NPQ_L2", "NPQ_L3", "NPQ_L4", "NPQ_L5", 
                       "NPQ_L6", "NPQ_L7", "NPQ_L8", "NPQ_L9", "NPQ_L10", "NPQ_L11", "NPQ_Lss", "NPQ_D1", 
                       "NPQ_D2", "NPQ_D3", "NPQ_D4", "NPQ_D5", "NPQ_D6", "NPQ_D7", "NPQ_D8", "Fv.Fm_L1", 
                       "Fv.Fm_L2", "Fv.Fm_L3",  "Fv.Fm_L4", "Fv.Fm_L5", "Fv.Fm_L6", "Fv.Fm_L7", "Fv.Fm_L8", 
                       "Fv.Fm_L9", "Fv.Fm_L10", "Fv.Fm_L11", "Fv.Fm_Lss", "Fv.Fm_D1", "Fv.Fm_D2", "Fv.Fm_D3", 
                       "Fv.Fm_D4", "Fv.Fm_D5", "Fv.Fm_D6", "Fv.Fm_D7", "Fv.Fm_D8")]


##### 2021 model ####

DataOutput_21 <- data.frame(matrix(vector(),316,1,dimnames = list(c(), c("Accession"))))

DataVarComp_21 <- data.frame()

drops <- c("var1","var2","sdcor")

#DataOutput$Accession <- unique(temporal_data_21[,2])

colnum <- c(8:ncol(temporal_data_21))

i = 1
for(i in  1:40){
  x=colnum[i]
  trait=colnames(temporal_data_21)[x]
  df1_21 <- temporal_data_21
  colnames(df1_21)[x]="y"
  
  model <- lmer(y ~ (1|Accession) + 
                  (1|Accession) + 
                  (1|Replicate) + (1|Block:Replicate) + 
                  (1|Measurement_day) + (1|Measurement_silking_interval), data = df1_21)
  
  summary(model)
  
  varComp_21 <- as.data.frame(VarCorr(model, comp="vcov"))
  
  blup = coef(model)$Accession
  
  colnames(blup) = trait
  
  DataOutput_21 <- cbind(DataOutput_21, blup) 
  
  varComp_21$Trait <- trait
  varComp_21 <- varComp_21[ , !(names(varComp_21) %in% drops)]
  
  DataVarComp_21 <- rbind(DataVarComp_21, varComp_21)
  
}

DataOutput_21$Accession <- rownames(DataOutput_21)
rownames(DataOutput_21) <- NULL

DataVarComp_21 <- dcast(DataVarComp_21, Trait ~ grp, value.var = "vcov")
DataVarComp_21$H2 <- DataVarComp_21$Accession/(DataVarComp_21$Accession+DataVarComp_21$Residual)

DataVarComp_21$TimePoint <- c(13,14,15,16,17,18,19,20,1,10,11,2,3,4,5,6,7,8,9,12,
                           13,14,15,16,17,18,19,20,1,10,11,2,3,4,5,6,7,8,9,12)

DataVarComp_21 <- DataVarComp_21 %>% 
  separate(Trait, into = c("Trait", "Point"), sep="_(?=[^_]+$)", remove = FALSE)

temporal_2021_H2 <- subset(DataVarComp_21, select = c(Trait, Point, TimePoint, H2))
write.csv(temporal_2021_H2, 
          "~/My Drive/Maize joint years analyses/Chlorophyll fluorescence/temporal_2021_H2.csv",
          row.names = F)

##### 2022 model ####

DataOutput_22 <- data.frame(matrix(vector(),312,1,dimnames = list(c(), c("Accession"))))

DataVarComp_22 <- data.frame()

drops <- c("var1","var2","sdcor")

#DataOutput$Accession <- unique(temporal_data_21[,2])

colnum <- c(8:ncol(temporal_data_22))

i = 1
for(i in  1:40){
  x=colnum[i]
  trait=colnames(temporal_data_22)[x]
  df1_22 <- temporal_data_22
  colnames(df1_22)[x]="y"
  
  model <- lmer(y ~ (1|Accession) + 
                  (1|Accession) + 
                  (1|Replicate) + (1|Block:Replicate) + 
                  (1|Measurement_day) + (1|Measurement_silking_interval), data = df1_22)
  
  summary(model)
  
  varComp_22 <- as.data.frame(VarCorr(model, comp="vcov"))
  
  blup = coef(model)$Accession
  
  colnames(blup) = trait
  
  DataOutput_22 <- cbind(DataOutput_22, blup) 
  
  varComp_22$Trait <- trait
  varComp_22 <- varComp_22[ , !(names(varComp_22) %in% drops)]
  
  DataVarComp_22 <- rbind(DataVarComp_22, varComp_22)
  
}

DataOutput_22$Accession <- rownames(DataOutput_22)
rownames(DataOutput_22) <- NULL

DataVarComp_22 <- dcast(DataVarComp_22, Trait ~ grp, value.var = "vcov")
DataVarComp_22$H2 <- DataVarComp_22$Accession/(DataVarComp_22$Accession+DataVarComp_22$Residual)

DataVarComp_22$TimePoint <- c(13,14,15,16,17,18,19,20,1,10,11,2,3,4,5,6,7,8,9,12,
                           13,14,15,16,17,18,19,20,1,10,11,2,3,4,5,6,7,8,9,12)

DataVarComp_22 <- DataVarComp_22 %>% 
  separate(Trait, into = c("Trait", "Point"), sep="_(?=[^_]+$)", remove = FALSE)

temporal_2022_H2 <- subset(DataVarComp_22, select = c(Trait, Point, TimePoint, H2))
write.csv(temporal_2022_H2, 
          "~/My Drive/Maize joint years analyses/Chlorophyll fluorescence/temporal_2022_H2.csv",
          row.names = F)


##### joint model #####

DataOutput_joint <- data.frame(matrix(vector(),326,1,dimnames = list(c(), c("Accession"))))

DataVarComp_joint <- data.frame()

drops <- c("var1","var2","sdcor")

#DataOutput$Accession <- unique(temporal_data_21[,2])

colnum_joint <- c(9:ncol(temporal_data_21_22))
i = 1
for(i in  1:40){
  x=colnum_joint[i]
  trait=colnames(temporal_data_21_22)[x]
  df1_joint <- temporal_data_21_22
  colnames(df1_joint)[x]="y"
  
  model <- lmer(y ~ (1|Accession) + (1|Year) + 
                  (1|Accession:Year) + 
                  (1|Replicate:Year) + (1|Block:Replicate:Year) + 
                  (1|Measurement_day:Year) + (1|Measurement_silking_interval:Year), data = df1_joint)
  
  summary(model)
  
  varComp_joint <- as.data.frame(VarCorr(model, comp="vcov"))
  
  blup = coef(model)$Accession
  
  colnames(blup) = trait
  
  DataOutput_joint <- cbind(DataOutput_joint, blup) 
  
  varComp_joint$Trait <- trait
  varComp_joint <- varComp_joint[ , !(names(varComp_joint) %in% drops)]
  
  DataVarComp_joint <- rbind(DataVarComp_joint, varComp_joint)
  
}

DataOutput_joint$Accession <- rownames(DataOutput_joint)
rownames(DataOutput_joint) <- NULL

DataVarComp_joint <- dcast(DataVarComp_joint, Trait ~ grp, value.var = "vcov")
DataVarComp_joint$H2 <- DataVarComp_joint$Accession/(DataVarComp_joint$Accession+DataVarComp_joint$Residual)
DataVarComp_joint$TimePoint <- c(13,14,15,16,17,18,19,20,1,10,11,2,3,4,5,6,7,8,9,12,
                           13,14,15,16,17,18,19,20,1,10,11,2,3,4,5,6,7,8,9,12)

DataVarComp_joint <- DataVarComp_joint %>% 
  separate(Trait, into = c("Trait", "Point"), sep="_(?=[^_]+$)", remove = FALSE)

phiPSII_H2 <- subset(DataVarComp_joint, Trait == "Fv.Fm")
phiPSII_H2$TimePoint <- as.factor(phiPSII_H2$TimePoint)
phiPSII_H2_plot <- ggplot(data=phiPSII_H2, aes(x=TimePoint, y=H2, fill = TimePoint)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c("yellow", "yellow", "yellow", "yellow", "yellow", "yellow", "yellow", "yellow", "yellow", "yellow", "yellow", "yellow", "black", "black", "black", "black", "black", "black", "black", "black")) + 
  scale_y_continuous(limits = c(0,0.8), breaks = seq(0,1,0.2)) + 
  ylab("phiPSII Heritability") + 
  theme_classic() + 
  theme(legend.position = "none")

NPQ_H2 <- subset(DataVarComp_joint, Trait == "NPQ")
NPQ_H2$TimePoint <- as.factor(NPQ_H2$TimePoint)
NPQ_H2_plot <- ggplot(data=NPQ_H2, aes(x=TimePoint, y=H2, fill = TimePoint)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c("yellow", "yellow", "yellow", "yellow", "yellow", "yellow", "yellow", "yellow", "yellow", "yellow", "yellow", "yellow", "black", "black", "black", "black", "black", "black", "black", "black")) + 
  scale_y_continuous(limits = c(0,0.8), breaks = seq(0,1,0.2)) + 
  ylab("NPQ Heritability") + 
  theme_classic() + 
  theme(legend.position = "none")

grid.arrange(phiPSII_H2_plot, NPQ_H2_plot, ncol = 1)

write.csv(DataOutput, "~/My Drive/Maize joint years analyses/Chlorophyll fluorescence/CF Joint temporal BLUPs.csv", 
          row.names = F)

temporal_BLUPs_21 <- read.csv(
  "~/My Drive/Maize 2021 analyses/Chlorophyll Fluorescence 2021/CF 2021 temporal BLUPs.csv", header = T
)
colnames(temporal_BLUPs_21) <- paste(colnames(temporal_BLUPs_21),"21",sep="_")
names(temporal_BLUPs_21)[names(temporal_BLUPs_21) == 'Accession_21'] <- 'Accession'
temporal_BLUPs_21$Accession <-gsub("_", "_00", temporal_BLUPs_21$Accession)

temporal_BLUPs_22 <- read.csv(
  "~/My Drive/Maize 2022 analyses/Chlorophyll Fluorescence 2022/CF 2022 temporal BLUPs.csv", header = T
)
colnames(temporal_BLUPs_22) <- paste(colnames(temporal_BLUPs_22),"22",sep="_")
names(temporal_BLUPs_22)[names(temporal_BLUPs_22) == 'Accession_22'] <- 'Accession'

colnames(DataOutput) <- paste(colnames(DataOutput), "Joint", sep = "_")
names(DataOutput)[names(DataOutput) == 'Accession_Joint'] <- 'Accession'

temporal_BLUPs_allModels <- merge(temporal_BLUPs_21, temporal_BLUPs_22, by = "Accession", all = TRUE)
temporal_BLUPs_allModels <- merge(temporal_BLUPs_allModels, DataOutput, by = "Accession", all = TRUE)

temporal_BLUPs_allModels <- temporal_BLUPs_allModels[,order(colnames(temporal_BLUPs_allModels))]

write.csv(temporal_BLUPs_allModels, 
          "~/My Drive/Maize joint years analyses/Chlorophyll fluorescence/CF allModels temporal BLUPs.csv", 
          row.names = F)

joint_temporal_H2 <- subset(DataVarComp_joint, select = c(Trait, Point, TimePoint, H2))
write.csv(joint_temporal_H2, 
          "~/My Drive/Maize joint years analyses/Chlorophyll fluorescence/temporal_joint_H2.csv",
          row.names = F)

