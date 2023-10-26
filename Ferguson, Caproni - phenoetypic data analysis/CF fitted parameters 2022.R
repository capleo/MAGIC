#set working directory 

setwd("~/My Drive/Maize 2022 analyses/Chlorophyll Fluorescence 2022//")

#Get packages ----

library(readxl)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(ggforce)
library(tidyverse)
library(dplyr)
library(writexl)
library(segmented)
library(minpack.lm)
library(data.table)
library(reshape2)

#get data ----

NPQdata <- read_excel("NPQdata.xlsx")
NPQdata$plot_id <- paste(NPQdata$Plot, NPQdata$Repeat)

#getting the phiPSII and NPQ values at each point

NPQ_L1_values <- subset(NPQdata, NPQ_names == "NPQ_L1")
NPQ_L1_values <- subset(NPQ_L1_values, select = c("plot_id", "NPQ_values"))
names(NPQ_L1_values)[names(NPQ_L1_values) == 'NPQ_values'] <- 'NPQ_L1'

NPQ_L2_values <- subset(NPQdata, NPQ_names == "NPQ_L2")
NPQ_L2_values <- subset(NPQ_L2_values, select = c("plot_id", "NPQ_values"))
names(NPQ_L2_values)[names(NPQ_L2_values) == 'NPQ_values'] <- 'NPQ_L2'

NPQ_L3_values <- subset(NPQdata, NPQ_names == "NPQ_L3")
NPQ_L3_values <- subset(NPQ_L3_values, select = c("plot_id", "NPQ_values"))
names(NPQ_L3_values)[names(NPQ_L3_values) == 'NPQ_values'] <- 'NPQ_L3'

NPQ_L4_values <- subset(NPQdata, NPQ_names == "NPQ_L4")
NPQ_L4_values <- subset(NPQ_L4_values, select = c("plot_id", "NPQ_values"))
names(NPQ_L4_values)[names(NPQ_L4_values) == 'NPQ_values'] <- 'NPQ_L4'

NPQ_L5_values <- subset(NPQdata, NPQ_names == "NPQ_L5")
NPQ_L5_values <- subset(NPQ_L5_values, select = c("plot_id", "NPQ_values"))
names(NPQ_L5_values)[names(NPQ_L5_values) == 'NPQ_values'] <- 'NPQ_L5'

NPQ_L6_values <- subset(NPQdata, NPQ_names == "NPQ_L6")
NPQ_L6_values <- subset(NPQ_L6_values, select = c("plot_id", "NPQ_values"))
names(NPQ_L6_values)[names(NPQ_L6_values) == 'NPQ_values'] <- 'NPQ_L6'

NPQ_L7_values <- subset(NPQdata, NPQ_names == "NPQ_L7")
NPQ_L7_values <- subset(NPQ_L7_values, select = c("plot_id", "NPQ_values"))
names(NPQ_L7_values)[names(NPQ_L7_values) == 'NPQ_values'] <- 'NPQ_L7'

NPQ_L8_values <- subset(NPQdata, NPQ_names == "NPQ_L8")
NPQ_L8_values <- subset(NPQ_L8_values, select = c("plot_id", "NPQ_values"))
names(NPQ_L8_values)[names(NPQ_L8_values) == 'NPQ_values'] <- 'NPQ_L8'

NPQ_L9_values <- subset(NPQdata, NPQ_names == "NPQ_L9")
NPQ_L9_values <- subset(NPQ_L9_values, select = c("plot_id", "NPQ_values"))
names(NPQ_L9_values)[names(NPQ_L9_values) == 'NPQ_values'] <- 'NPQ_L9'

NPQ_L10_values <- subset(NPQdata, NPQ_names == "NPQ_L10")
NPQ_L10_values <- subset(NPQ_L10_values, select = c("plot_id", "NPQ_values"))
names(NPQ_L10_values)[names(NPQ_L10_values) == 'NPQ_values'] <- 'NPQ_L10'

NPQ_L11_values <- subset(NPQdata, NPQ_names == "NPQ_L11")
NPQ_L11_values <- subset(NPQ_L11_values, select = c("plot_id", "NPQ_values"))
names(NPQ_L11_values)[names(NPQ_L11_values) == 'NPQ_values'] <- 'NPQ_L11'

NPQ_Lss_values <- subset(NPQdata, NPQ_names == "NPQ_Lss")
NPQ_Lss_values <- subset(NPQ_Lss_values, select = c("plot_id", "NPQ_values"))
names(NPQ_Lss_values)[names(NPQ_Lss_values) == 'NPQ_values'] <- 'NPQ_Lss'

NPQ_light_values <- cbind(NPQ_L1_values, NPQ_L2_values[,2],  NPQ_L3_values[,2],  NPQ_L4_values[,2],
                          NPQ_L5_values[,2],  NPQ_L6_values[,2],  NPQ_L7_values[,2],  NPQ_L8_values[,2],
                          NPQ_L9_values[,2],  NPQ_L10_values[,2],  NPQ_L11_values[,2],  NPQ_Lss_values[,2])

rm(NPQ_L1_values, NPQ_L2_values,  NPQ_L3_values,  NPQ_L4_values, NPQ_L5_values,  NPQ_L6_values,
   NPQ_L7_values,  NPQ_L8_values, NPQ_L9_values,  NPQ_L10_values,  NPQ_L11_values,  NPQ_Lss_values)

NPQ_light_values <- separate(data = NPQ_light_values, col = plot_id, into = c("PlotID", "Repeat"))

NPQ_D1_values <- subset(NPQdata, NPQ_names == "NPQ_D1")
NPQ_D1_values <- subset(NPQ_D1_values, select = c("plot_id", "NPQ_values"))
names(NPQ_D1_values)[names(NPQ_D1_values) == 'NPQ_values'] <- 'NPQ_D1'

NPQ_D2_values <- subset(NPQdata, NPQ_names == "NPQ_D2")
NPQ_D2_values <- subset(NPQ_D2_values, select = c("plot_id", "NPQ_values"))
names(NPQ_D2_values)[names(NPQ_D2_values) == 'NPQ_values'] <- 'NPQ_D2'

NPQ_D3_values <- subset(NPQdata, NPQ_names == "NPQ_D3")
NPQ_D3_values <- subset(NPQ_D3_values, select = c("plot_id", "NPQ_values"))
names(NPQ_D3_values)[names(NPQ_D3_values) == 'NPQ_values'] <- 'NPQ_D3'

NPQ_D4_values <- subset(NPQdata, NPQ_names == "NPQ_D4")
NPQ_D4_values <- subset(NPQ_D4_values, select = c("plot_id", "NPQ_values"))
names(NPQ_D4_values)[names(NPQ_D4_values) == 'NPQ_values'] <- 'NPQ_D4'

NPQ_D5_values <- subset(NPQdata, NPQ_names == "NPQ_D5")
NPQ_D5_values <- subset(NPQ_D5_values, select = c("plot_id", "NPQ_values"))
names(NPQ_D5_values)[names(NPQ_D5_values) == 'NPQ_values'] <- 'NPQ_D5'

NPQ_D6_values <- subset(NPQdata, NPQ_names == "NPQ_D6")
NPQ_D6_values <- subset(NPQ_D6_values, select = c("plot_id", "NPQ_values"))
names(NPQ_D6_values)[names(NPQ_D6_values) == 'NPQ_values'] <- 'NPQ_D6'

NPQ_D7_values <- subset(NPQdata, NPQ_names == "NPQ_D7")
NPQ_D7_values <- subset(NPQ_D7_values, select = c("plot_id", "NPQ_values"))
names(NPQ_D7_values)[names(NPQ_D7_values) == 'NPQ_values'] <- 'NPQ_D7'

NPQ_D8_values <- subset(NPQdata, NPQ_names == "NPQ_D8")
NPQ_D8_values <- subset(NPQ_D8_values, select = c("plot_id", "NPQ_values"))
names(NPQ_D8_values)[names(NPQ_D8_values) == 'NPQ_values'] <- 'NPQ_D8'

NPQ_dark_values <- cbind(NPQ_D1_values, NPQ_D2_values[,2],  NPQ_D3_values[,2],  NPQ_D4_values[,2],
                         NPQ_D5_values[,2],  NPQ_D6_values[,2],  NPQ_D7_values[,2],  NPQ_D8_values[,2])

rm(NPQ_D1_values, NPQ_D2_values,  NPQ_D3_values,  NPQ_D4_values, NPQ_D5_values,  NPQ_D6_values,
   NPQ_D7_values,  NPQ_D8_values)

NPQ_temporal_values <- cbind(NPQ_light_values, NPQ_dark_values[,-1])

phiPSII_L1_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_L1")
phiPSII_L1_values <- subset(phiPSII_L1_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_L1_values)[names(phiPSII_L1_values) == 'Fv/Fm_values'] <- 'Fv/Fm_L1'

phiPSII_L2_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_L2")
phiPSII_L2_values <- subset(phiPSII_L2_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_L2_values)[names(phiPSII_L2_values) == 'Fv/Fm_values'] <- 'Fv/Fm_L2'

phiPSII_L3_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_L3")
phiPSII_L3_values <- subset(phiPSII_L3_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_L3_values)[names(phiPSII_L3_values) == 'Fv/Fm_values'] <- 'Fv/Fm_L3'

phiPSII_L4_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_L4")
phiPSII_L4_values <- subset(phiPSII_L4_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_L4_values)[names(phiPSII_L4_values) == 'Fv/Fm_values'] <- 'Fv/Fm_L4'

phiPSII_L5_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_L5")
phiPSII_L5_values <- subset(phiPSII_L5_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_L5_values)[names(phiPSII_L5_values) == 'Fv/Fm_values'] <- 'Fv/Fm_L5'

phiPSII_L6_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_L6")
phiPSII_L6_values <- subset(phiPSII_L6_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_L6_values)[names(phiPSII_L6_values) == 'Fv/Fm_values'] <- 'Fv/Fm_L6'

phiPSII_L7_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_L7")
phiPSII_L7_values <- subset(phiPSII_L7_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_L7_values)[names(phiPSII_L7_values) == 'Fv/Fm_values'] <- 'Fv/Fm_L7'

phiPSII_L8_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_L8")
phiPSII_L8_values <- subset(phiPSII_L8_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_L8_values)[names(phiPSII_L8_values) == 'Fv/Fm_values'] <- 'Fv/Fm_L8'

phiPSII_L9_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_L9")
phiPSII_L9_values <- subset(phiPSII_L9_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_L9_values)[names(phiPSII_L9_values) == 'Fv/Fm_values'] <- 'Fv/Fm_L9'

phiPSII_L10_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_L10")
phiPSII_L10_values <- subset(phiPSII_L10_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_L10_values)[names(phiPSII_L10_values) == 'Fv/Fm_values'] <- 'Fv/Fm_L10'

phiPSII_L11_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_L11")
phiPSII_L11_values <- subset(phiPSII_L11_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_L11_values)[names(phiPSII_L11_values) == 'Fv/Fm_values'] <- 'Fv/Fm_L11'

phiPSII_Lss_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_Lss")
phiPSII_Lss_values <- subset(phiPSII_Lss_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_Lss_values)[names(phiPSII_Lss_values) == 'Fv/Fm_values'] <- 'Fv/Fm_Lss'

phiPSII_light_values <- cbind(phiPSII_L1_values, phiPSII_L2_values[,2],  phiPSII_L3_values[,2],
                              phiPSII_L4_values[,2], phiPSII_L5_values[,2],  phiPSII_L6_values[,2],
                              phiPSII_L7_values[,2],  phiPSII_L8_values[,2], phiPSII_L9_values[,2],
                              phiPSII_L10_values[,2],  phiPSII_L11_values[,2],  phiPSII_Lss_values[,2])

rm(phiPSII_L1_values, phiPSII_L2_values,  phiPSII_L3_values,  phiPSII_L4_values, phiPSII_L5_values,
   phiPSII_L6_values, phiPSII_L7_values,  phiPSII_L8_values, phiPSII_L9_values,  phiPSII_L10_values,
   phiPSII_L11_values,  phiPSII_Lss_values)

phiPSII_D1_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_D1")
phiPSII_D1_values <- subset(phiPSII_D1_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_D1_values)[names(phiPSII_D1_values) == 'Fv/Fm_values'] <- 'Fv/Fm_D1'

phiPSII_D2_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_D2")
phiPSII_D2_values <- subset(phiPSII_D2_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_D2_values)[names(phiPSII_D2_values) == 'Fv/Fm_values'] <- 'Fv/Fm_D2'

phiPSII_D3_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_D3")
phiPSII_D3_values <- subset(phiPSII_D3_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_D3_values)[names(phiPSII_D3_values) == 'Fv/Fm_values'] <- 'Fv/Fm_D3'

phiPSII_D4_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_D4")
phiPSII_D4_values <- subset(phiPSII_D4_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_D4_values)[names(phiPSII_D4_values) == 'Fv/Fm_values'] <- 'Fv/Fm_D4'

phiPSII_D5_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_D5")
phiPSII_D5_values <- subset(phiPSII_D5_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_D5_values)[names(phiPSII_D5_values) == 'Fv/Fm_values'] <- 'Fv/Fm_D5'

phiPSII_D6_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_D6")
phiPSII_D6_values <- subset(phiPSII_D6_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_D6_values)[names(phiPSII_D6_values) == 'Fv/Fm_values'] <- 'Fv/Fm_D6'

phiPSII_D7_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_D7")
phiPSII_D7_values <- subset(phiPSII_D7_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_D7_values)[names(phiPSII_D7_values) == 'Fv/Fm_values'] <- 'Fv/Fm_D7'

phiPSII_D8_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_D8")
phiPSII_D8_values <- subset(phiPSII_D8_values, select = c("plot_id", "Fv/Fm_values"))
names(phiPSII_D8_values)[names(phiPSII_D8_values) == 'Fv/Fm_values'] <- 'Fv/Fm_D8'

phiPSII_dark_values <- cbind(phiPSII_D1_values, phiPSII_D2_values[,2],  phiPSII_D3_values[,2],  
                             phiPSII_D4_values[,2], phiPSII_D5_values[,2],  phiPSII_D6_values[,2],
                             phiPSII_D7_values[,2],  phiPSII_D8_values[,2])

rm(phiPSII_D1_values, phiPSII_D2_values,  phiPSII_D3_values,  phiPSII_D4_values, phiPSII_D5_values,
   phiPSII_D6_values, phiPSII_D7_values,  phiPSII_D8_values)

phiPSII_temporal_values <- cbind(phiPSII_light_values, phiPSII_dark_values[,-1])

temporal_NPQ_phiPSII <- cbind(NPQ_temporal_values, phiPSII_temporal_values[,2:ncol(phiPSII_temporal_values)])

FvFm_values <- subset(NPQdata, `Fv/Fm` == "Fv/Fm_L0")
FvFm_values <- subset(FvFm_values, select = c("plot_id", "QY_max"))
FvFm_values <- separate(data = FvFm_values, col = plot_id, into = c("PlotID", "Repeat"))

temporal_NPQ_phiPSII$QY_max <- FvFm_values$QY_max

rm(FvFm_values)

temporal_NPQ_phiPSII$NPQ_ratio <- 1/(1+temporal_NPQ_phiPSII$NPQ_D8)

attach(temporal_NPQ_phiPSII)

temporal_NPQ_phiPSII$PI <- `Fv/Fm_D8` / 1 - ((1-QY_max)/(QY_max + ((1 - QY_max)/NPQ_ratio))/NPQ_ratio)

detach(temporal_NPQ_phiPSII)

rm(NPQ_dark_values, NPQ_light_values, phiPSII_dark_values, phiPSII_light_values,
   NPQ_temporal_values, phiPSII_temporal_values)

temporal_NPQ_phiPSII <- subset(temporal_NPQ_phiPSII, select = -c(QY_max, NPQ_ratio))

write.csv(temporal_NPQ_phiPSII, "temporal_NPQ_phiPSII.csv", row.names = F)

#getting initial Fv/Fm values ----

FvFm_values <- subset(NPQdata, select = c("plot_id", "QY_max"))

FvFm_values <- FvFm_values[!duplicated(FvFm_values), ]

#exponential fitting for induction test ----

induction_data <- filter(NPQdata, cumulative_time <= 600)
 
test <- induction_data[1:13,]
 
induction_model <- nlsLM(NPQ_values ~ a*(1-exp(-b*cumulative_time)),data=test, 
                          start = list(a=3.5, b=0.005))
 
a_fit <- induction_model$m$getAllPars()[1]
b_fit <- induction_model$m$getAllPars()[2]
 
plot(test$time_post_light_on, test$NPQ_values)
x <- seq(0,600, by=6)
y <- predict(induction_model, list(cumulative_time=x))
points(x,y, type="l")

#make function to use on all data for induction NPQ----

induction_data <- filter(NPQdata, cumulative_time <= 600)

NPQ_induction_fitting <- function(dataset,plot,a=3.5, b=0.005){
  
  subset <- subset(dataset, plot_id == plot) #subsets for each unique plot/repeat
  induction_model <- nlsLM(NPQ_values ~ a*(1-exp(-b*cumulative_time)),data=subset, 
                           start = list(a=a, b=b))
  a_fit <- induction_model$m$getAllPars()[1]
  a_confint <- confint(induction_model)[1,2] - confint(induction_model)[1,1]
  b_fit <- induction_model$m$getAllPars()[2]
  b_confint <- confint(induction_model)[2,2] - confint(induction_model)[2,1]
  induction_fitting <- data.frame(plot = plot, 
                                 a_fit = a_fit,
                                 a_confint = a_confint,
                                 b_fit = b_fit,
                                 b_confint=b_confint)
  
}

induction_plot_ids <-as.character(unique(induction_data$plot_id))
names(induction_plot_ids) <- induction_plot_ids

all_induction_fits <- lapply(induction_plot_ids, function(x) NPQ_induction_fitting(induction_data,x))
all_induction_fits <- rbindlist(all_induction_fits)

ggplot(all_induction_fits, aes(x=plot, y=a_confint)) + geom_point()
ggplot(all_induction_fits, aes(x=plot, y=b_confint)) + geom_point()

#plot outputs ----

# PDFpath <- "~/OneDrive - University of Cambridge/MPhil/NPQinduction.pdf"
# pdf(file=PDFpath)
# 
# for (value in unique(induction_data$plot_id)){
#   
#   subset <- subset(induction_data, plot_id == value)
#   plot(subset$cumulative_time, subset$NPQ_values, 
#        main=paste("Plot of", value))
#   x <- seq(0,600,by=6)
#   induction_model <- nlsLM(NPQ_values ~ a*(1-exp(-b*cumulative_time)),data=subset, 
#                            start = list(a=3, b=0.005))
#   y <- predict(induction_model, list(cumulative_time=x))
#   points(x,y, type="l")
# }
# 
# dev.off()
  
#exponential fitting for relaxation ----
# 
# relaxation_data <- filter(NPQdata, cumulative_time >= 600)
# 
# test2 <- relaxation_data[1:9,]
# 
# relaxation_model <- nlsLM(NPQ_values ~ c*(exp(-d*time_post_light_off)) + e, data=test2,
#                           start=list(c=3,d=0.005,e=0.5))
# 
# c_fit <- relaxation_model$m$getAllPars()[1]
# d_fit <- relaxation_model$m$getAllPars()[2]
# e_fit <- relaxation_model$m$getAllPars()[3]
# 
# plot(test2$time_post_light_off, test2$NPQ_values)
# x <- seq(0,720,by=5)
# y <- predict(relaxation_model, list(time_post_light_off=x))
# points(x,y, type="l")

#make function to use on all data for relaxation NPQ----

relaxation_data <- filter(NPQdata, cumulative_time >= 600)

NPQ_relaxation_fitting <- function(dataset,plot,c=3,d=0.005,e=0.5){
  
  subset <- subset(dataset, plot_id == plot) #subsets for each unique plot/repeat
  relaxation_model <- nlsLM(NPQ_values ~ c*(exp(-d*time_post_light_off)) + e, data=subset,
                            start=list(c=c,d=d,e=e))
  c_fit <- relaxation_model$m$getAllPars()[1]
  d_fit <- relaxation_model$m$getAllPars()[2]
  e_fit <- relaxation_model$m$getAllPars()[3]
  relaxation_fitting <- data.frame(plot = plot, 
                                 c_fit = c_fit,
                                 d_fit = d_fit,
                                 e_fit=e_fit)
  
}

relaxation_plot_ids <-as.character(unique(relaxation_data$plot_id))
names(relaxation_plot_ids) <- relaxation_plot_ids

all_relaxation_fits <- lapply(relaxation_plot_ids, function(x) NPQ_relaxation_fitting(relaxation_data,x))
all_relaxation_fits <- rbindlist(all_relaxation_fits)

#Plot outputs ----

# PDFpath2 <- "~/OneDrive - University of Cambridge/MPhil/NPQrelaxation.pdf"
# pdf(file=PDFpath2)
# 
# for (value in unique(relaxation_data$plot_id)){
#   
#   subset <- subset(relaxation_data, plot_id == value)
#   plot(subset$time_post_light_off, subset$NPQ_values, 
#        main=paste("Plot of", value))
#   x <- seq(0,720,by=5)
#   relaxation_model <- nlsLM(NPQ_values ~ c*(exp(-d*time_post_light_off)) + e, data=subset,
#                             start=list(c=3,d=0.005,e=0.5))
#   y <- predict(relaxation_model, list(time_post_light_off=x))
#   points(x,y, type="l")
# }
# 
# dev.off()

#initial slope NPQ----

initial_data <- filter(NPQdata, cumulative_time <= 40)
initial_data <- select(initial_data, plot_id,Plot,Repeat, NPQ_values, cumulative_time)

initial_slope_model <- function(dataset, plot) {
  
  subset <- subset(dataset, plot_id==plot)
  initial_model <- lm(NPQ_values ~ 0 + cumulative_time, data=subset) #the 0 tells it to fit through the origin
  initial_slope_fitting <- data.frame(plot=plot,
                                      gradient=coef(initial_model))
}

initial_slope_plot_ids <-as.character(unique(initial_data$plot_id))
names(initial_slope_plot_ids) <- initial_slope_plot_ids

initial_slope_fits <- lapply(initial_slope_plot_ids, function(x) initial_slope_model(initial_data,x))
initial_slope_fits <- rbindlist(initial_slope_fits)

initial_data <- filter(NPQdata, cumulative_time == 40)
initial_data <- select(initial_data, plot_id, Plot, Repeat)

#maximum amplitude NPQ----

max_amp <- filter(NPQdata, NPQ_names=="NPQ_Lss")
max_amp <- select(max_amp, plot_id, NPQ_values)
names(max_amp)[2] <- "max_amp"

#end point NPQ ----

end_point <- filter(NPQdata, NPQ_names == "NPQ_D8")
end_point <- select(end_point, plot_id, NPQ_values, `Fv/Fm_values`)
names(end_point)[2] <- "end_NPQ"
names(end_point)[3] <- "end_FvFm"

# phi PSII in the dark parameters----

names(NPQdata)[5] <- "FvFm_values"
# 
relaxation_data <- filter(NPQdata, cumulative_time >= 600)
# 
# test <- relaxation_data[1:9,]
# 
# PSII_model <- nlsLM(FvFm_values ~ f*(1-exp(-g*time_post_light_off)) + h,data=test, start = list(f=0.25,g=0.002,h=0.6))
# f_fit <- PSII_model$m$getAllPars()[1]
# g_fit <- PSII_model$m$getAllPars()[2]
# h_fit <- PSII_model$m$getAllPars()[3]
# 
# plot(test$time_post_light_off, test$FvFm_values)
# x <- seq(0,720,by=5)
# y <- predict(PSII_model, list(time_post_light_off=x))
# points(x,y, type="l")

#uses relaxation_data

PSII_dark_fitting <- function(dataset,plot,f=0.25,g=0.002,h=0.6){
  
  subset <- subset(dataset, plot_id == plot) #subsets for each unique plot/repeat
  PSII_model <- nlsLM(FvFm_values ~ f*(1-exp(-g*time_post_light_off)) + h,data=subset, 
                           start = list(f=f,g=g,h=h))
  f_fit <- PSII_model$m$getAllPars()[1]
  g_fit <- PSII_model$m$getAllPars()[2]
  h_fit <- PSII_model$m$getAllPars()[3]

  PSII_fitting <- data.frame(plot = plot,
                             f_fit = f_fit,
                             g_fit = g_fit,
                             h_fit = h_fit)
  
}

PSII_plot_ids <-as.character(unique(relaxation_data$plot_id))
names(PSII_plot_ids) <- PSII_plot_ids

all_PSII_fits <- lapply(PSII_plot_ids, function(x) PSII_dark_fitting(relaxation_data,x))
all_PSII_fits <- rbindlist(all_PSII_fits)


#get all parameters into one dataframe

NPQ_parms <- merge(all_induction_fits, all_relaxation_fits,by="plot")
NPQ_parms <- merge(NPQ_parms, all_PSII_fits, by = "plot")
NPQ_parms <- merge(NPQ_parms, initial_slope_fits, by = "plot")
names(NPQ_parms)[1] <- "plot_id"
NPQ_parms <- merge(NPQ_parms, max_amp, by="plot_id")
NPQ_parms <- merge(NPQ_parms, end_point, by="plot_id")
NPQ_parms <- merge(NPQ_parms, FvFm_values, by = "plot_id", all = TRUE)



#tidy data based on confint for induction parameters? ----

ggplot(NPQ_parms, aes(x=plot_id, y=a_confint)) + geom_point() 
ggplot(NPQ_parms, aes(x=factor(0), y=a_confint)) + geom_boxplot()
ggplot(NPQ_parms, aes(x=a_confint)) + geom_density() + geom_vline(aes(xintercept=2)) #suggest the cut off is 2 for a! 

ggplot(NPQ_parms, aes(x=plot_id, y=b_confint)) + geom_point() 
ggplot(NPQ_parms, aes(x=factor(0), y=b_confint)) + geom_boxplot() # suggest cut off is 0.01 for b!
ggplot(NPQ_parms, aes(x=b_confint)) + geom_density() + geom_vline(aes(xintercept=0.01)) #suggest the cut off is 0.01 for b!

# NPQ_parms <- filter(NPQ_parms, a_confint <= 2)
# NPQ_parms <- filter(NPQ_parms, b_confint <= 0.01)

#visualise! ----

ggplot(NPQ_parms, aes(x=e_fit, y=end_NPQ)) + geom_point() # should be y=x ish as e_fit is estimating the end point essentially
ggplot(NPQ_parms, aes(x=a_fit, y=max_amp)) + geom_point() # should be fairly positive correlation
ggplot(NPQ_parms, aes(x=b_fit, y=gradient)) + geom_point() # should be fairly positive correlation

#final dataset ----

NPQ_parms <- merge(initial_data, NPQ_parms, by="plot_id")
NPQ_parms <- select(NPQ_parms, -c(a_confint, b_confint, plot_id))

#raw_data <- dcast(melt(NPQ_parms, id.vars = c("Plot", "Repeat")), Plot~variable+Repeat)
#raw_data$a_fit_mean <- rowMeans(raw_data[,c('a_fit_1.00', 'a_fit_2.00', 'a_fit_3.00')], na.rm = T)
#raw_data$b_fit_mean <- rowMeans(raw_data[,c('b_fit_1.00', 'b_fit_2.00', 'b_fit_3.00')], na.rm = T)
#raw_data$c_fit_mean <- rowMeans(raw_data[,c('c_fit_1.00', 'c_fit_2.00', 'c_fit_3.00')], na.rm = T)
#raw_data$d_fit_mean <- rowMeans(raw_data[,c('d_fit_1.00', 'd_fit_2.00', 'd_fit_3.00')], na.rm = T)
#raw_data$e_fit_mean <- rowMeans(raw_data[,c('e_fit_1.00', 'e_fit_2.00', 'e_fit_3.00')], na.rm = T)
#raw_data$f_fit_mean <- rowMeans(raw_data[,c('f_fit_1.00', 'f_fit_2.00', 'f_fit_3.00')], na.rm = T)
#raw_data$g_fit_mean <- rowMeans(raw_data[,c('g_fit_1.00', 'g_fit_2.00', 'g_fit_3.00')], na.rm = T)
#raw_data$h_fit_mean <- rowMeans(raw_data[,c('h_fit_1.00', 'h_fit_2.00', 'h_fit_3.00')], na.rm = T)
#raw_data$gradient_mean <- rowMeans(raw_data[,c('gradient_1.00', 'gradient_2.00', 'gradient_3.00')], na.rm = T)
#raw_data$max_amp_mean <- rowMeans(raw_data[,c('max_amp_1.00', 'max_amp_2.00', 'max_amp_3.00')], na.rm = T)
#raw_data$end_NPQ_mean <- rowMeans(raw_data[,c('end_NPQ_1.00', 'end_NPQ_2.00', 'end_NPQ_3.00')], na.rm = T)
#raw_data$end_FvFm_mean <- rowMeans(raw_data[,c('end_FvFm_1.00', 'end_FvFm_2.00', 'end_FvFm_3.00')], na.rm = T)

#raw_data <- subset(raw_data, 
#               select = c("Plot", 
#                          "a_fit_1.00", "a_fit_2.00", "a_fit_3.00", "a_fit_mean",
#                          "b_fit_1.00", "b_fit_2.00", "b_fit_3.00", "b_fit_mean",
#                          "c_fit_1.00", "c_fit_2.00", "c_fit_3.00", "c_fit_mean",
#                          "d_fit_1.00", "d_fit_2.00", "d_fit_3.00", "d_fit_mean",
#                          "e_fit_1.00", "e_fit_2.00", "e_fit_3.00", "e_fit_mean",
#                          "f_fit_1.00", "f_fit_2.00", "f_fit_3.00", "f_fit_mean",
#                          "g_fit_1.00", "g_fit_2.00", "g_fit_3.00", "g_fit_mean",
#                          "h_fit_1.00", "h_fit_2.00", "h_fit_3.00", "h_fit_mean",
#                          "gradient_1.00", "gradient_2.00", "gradient_3.00", "gradient_mean",
#                          "max_amp_1.00", "max_amp_2.00", "max_amp_3.00", "max_amp_mean",
#                          "end_NPQ_1.00", "end_NPQ_2.00", "end_NPQ_3.00", "end_NPQ_mean",
#                          "end_FvFm_1.00", "end_FvFm_2.00", "end_FvFm_3.00", "end_FvFm_mean"))

#names(raw_data)[names(raw_data) == 'Plot'] <- 'PlotID'

#maize_ref <- read.csv("~/Google Drive/My Drive/Maize 2021 analyses/Maize_reference.csv",
#                      header = T,
#                      stringsAsFactors = F)
#maize_ref <- maize_ref[,1:6]
#maize_ref <- na.omit(maize_ref)
#maize_ref <- subset(maize_ref, select = c("Accession", "PlotID", "Replicate", "Column", "Row", "Measurement.Day"))

#raw_data <- merge(maize_ref, raw_data, by = "PlotID", all = TRUE)

#colnames(raw_data) <- c("Plot", "Accession", "Column", "Block", "Replication", "Date",
#                      "NPQ_ind_amp_1", "NPQ_ind_amp_2", "NPQ_ind_amp_3", "NPQ_ind_amp_mean", 
#                      "NPQ_ind_rate_1", "NPQ_ind_rate_2", "NPQ_ind_rate_3", "NPQ_ind_rate_mean",
#                      "NPQ_rel_amp_1", "NPQ_rel_amp_2", "NPQ_rel_amp_3", "NPQ_rel_amp_mean",
#                      "NPQ_rel_rate_1", "NPQ_rel_rate_2", "NPQ_rel_rate_3", "NPQ_rel_rate_mean",
#                      "NPQ_rel_res_1", "NPQ_rel_res_2", "NPQ_rel_res_3", "NPQ_rel_res_mean", 
#                      "phiPSII_ind_amp_1", "phiPSII_ind_amp_2", "phiPSII_ind_amp_3", "phiPSII_ind_amp_mean",
#                      "phiPSII_ind_rate_1", "phiPSII_ind_rate_2", "phiPSII_ind_rate_3", "phiPSII_ind_rate_mean",
#                      "phiPSII_ind_res_1", "phiPSII_ind_res_2", "phiPSII_ind_res_3", "phiPSII_ind_res_mean",
#                      "NPQ_ind_linear_1", "NPQ_ind_linear_1_2", "NPQ_ind_linear_1_3", "NPQ_ind_linear_1_mean",
#                      "maxNPQ_1", "maxNPQ_2", "maxNPQ_3", "maxNPQ_mean",
#                      "endNPQ_1", "endNPQ_2", "endNPQ_3", "endNPQ_mean",
#                      "endFvFm_1", "endFvFm_2", "endFvFm_3", "endFvFm_mean")

#raw_data <- raw_data[c(2,1:ncol(raw_data))]

#raw_data <- subset(raw_data, select = -c(Accession.1))

#write.csv(raw_data, 
#          "~/My Drive/Maize 2021 processed data/Raw data/2022_fluorescence_maize.csv", 
#          row.names = F)


colnames(NPQ_parms) <- c("Plot", "Repeat",
                      "NPQ_ind_amp", 
                      "NPQ_ind_rate",
                      "NPQ_rel_amp",
                      "NPQ_rel_rate",
                      "NPQ_rel_res", 
                      "phiPSII_ind_amp",
                      "phiPSII_ind_rate",
                      "phiPSII_ind_res",
                      "NPQ_ind_linear",
                      "maxNPQ",
                      "endNPQ",
                      "endFvFm",
                      "initialFvFm")

write.csv(NPQ_parms, 
          "~/My Drive/Maize 2022 analyses/Chlorophyll fluorescence/2022_fluorescence_maize.csv",
          row.names = F)


#add in genotype to dataframe ----

#genotype_data <- read_excel("~/OneDrive - University of Cambridge/MPhil/Phenotyping Campaign/FieldDesign.xlsx")
#names(genotype_data)[1] <- "Plot"

#NPQ_parms <- merge(NPQ_parms, genotype_data, by = "Plot")

#write_xlsx(NPQ_parms, "~/OneDrive - University of Cambridge/MPhil/GitLink/ExportedData/NPQparms.xlsx")
