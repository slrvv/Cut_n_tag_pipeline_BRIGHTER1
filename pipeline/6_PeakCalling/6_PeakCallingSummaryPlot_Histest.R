################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 6. Peak calling with SEACR                                                   #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: peak calling Histidine test summary plotting                        #
################################################################################

#-----------------------Paths and libraries------------------------------------#
library(ggplot2)
library(viridis)
library(dplyr)
library(ggpubr)
library(RColorBrewer)

#--------------------------functions-------------------------------------------#

peakNumberFig <- function(peakSEACR, title, facetExpression){
  fig7A = peakSEACR %>% ggplot(aes(x = IggType, y = peakN, fill = IggType)) +
    geom_boxplot() +
    geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
    facet_grid(facetExpession) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Number of Peaks") +
    xlab("")+
    ggtitle(title)
  
  return(fig7A)
}

peakRepFig <- function(peakSEACR, title, facetExpression){
  fig7C = peakSEACR %>% ggplot(aes(x = IggType, y = peakReprodRate, fill = IggType, label = round(peakReprodRate, 2))) +
    geom_bar(stat = "identity") +
    geom_text(vjust = 0.1) +
    facet_grid(facetExpression) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("% of Peaks Reproduced") +
    xlab("")+
    ggtitle(title)
  
  return(fig7C)
  
}

#--------------------------Peak calling assessment-----------------------------#

projPath <- "/project/ChromGroup_Seq_data/Celeste/2024_summer_complete"  #fill this variable with the same path as PROJECTROOT in 
# the cut_and_tag_pipeline.sh script

peakpathSEACR <- paste0(projPath, 
                        "/peakCalling/His_test/SEACR_His_test_peakcalling_summary.txt")
peakpathSEACRrelaxed <- paste0(projPath, 
                        "/peakCalling/His_test/SEACR_relaxed_His_test_peakcalling_summary.txt")
peakpathMACS2  <- paste0(projPath, 
                         "/peakCalling/His_test/His_test_summary_peak_calling_MACS2.txt")

peakSEACR <- read.table(peakpathSEACR, header=T)
peakSEACR$peakCaller <- rep("SEACR", times = nrow(peakSEACR))
peakSEACR$peakCaller <- paste(peakSEACR$peakCaller,peakSEACR$SEACR_parameter, sep = "_")
peakMACS2 <- read.table(peakpathMACS2, header=T) 
peakMACS2$peakCaller <- rep("MACS2", times = nrow(peakMACS2))
head(peakMACS2)
peakMACS2$peakCaller <- paste(peakMACS2$peakCaller,peakMACS2$MACS2_qvalue, sep = "_")

peakSEACR$HistoneIgg <- paste(peakSEACR$Histone,peakSEACR$IggType, sep = "_")
peakSEACR <- read.table(peakpathSEACRrelaxed, header=T) 

peakSEACRnorm <- peakSEACR[peakSEACR$SEACR_parameter == "norm",]
peakSEACRrelaxed <- read.table(peakpathSEACRrelaxed, header=T) 
peakSEACRrelaxed <- peakSEACR[peakSEACR$peakType %in% c("control", "control.rmDup"),]

peakSEACR_fig7A <- peakNumberFig(peakSEACR,"SEACR number of peaks HBD.His", "SEACR_parameter~peakType")

fig7A = peakSEACR %>% ggplot(aes(x = IggType, y = peakN, fill = IggType)) +
  geom_boxplot() +
  geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  facet_grid(SEACR_parameter~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") +
  xlab("")+
  ggtitle("SEACR number of peaks HBD.His")

fig7A

fig7A = peakSEACRrelaxed %>% ggplot(aes(x = IggType, y = peakN, fill = IggType)) +
  geom_boxplot() +
  geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  facet_grid(~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") +
  xlab("")+
  ggtitle("SEACR relaxed norm number of peaks HBD.His")

fig7A

fig7A = peakSEACRnorm %>% ggplot(aes(x = IggType, y = peakN, fill = IggType)) +
  geom_boxplot() +
  geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  facet_grid(~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") +
  xlab("")+
  ggtitle("SEACR norm number of peaks HBD.His")

fig7A


fig7C = peakSEACR %>% ggplot(aes(x = IggType, y = peakReprodRate, fill = IggType, label = round(peakReprodRate, 2))) +
  geom_bar(stat = "identity") +
  geom_text(vjust = 0.1) +
  facet_grid(SEACR_parameter+Replicate~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Peaks Reproduced") +
  xlab("")+
  ggtitle("SEACR replicate reproducibility HBD.His")

fig7C

fig7C = peakSEACRrelaxed %>% ggplot(aes(x = IggType, y = peakReprodRate, fill = IggType, label = round(peakReprodRate, 2))) +
  geom_bar(stat = "identity") +
  geom_text(vjust = 0.1) +
  facet_grid(Replicate~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Peaks Reproduced") +
  xlab("")+
  ggtitle("SEACR relaxed replicate reproducibility HBD.His")

fig7C


fig7A = peakMACS2 %>% ggplot(aes(x = IggType, y = peakN, fill = IggType)) +
  geom_boxplot() +
  geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  facet_grid(MACS2_qvalue~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") +
  xlab("")+
  ggtitle("MACS2 number of peaks HBD.His")

fig7A

fig7C = peakMACS2 %>% ggplot(aes(x = IggType, y = peakReprodRate, fill = IggType, label = round(peakReprodRate, 2))) +
  geom_bar(stat = "identity") +
  geom_text(vjust = 0.1) +
  facet_grid(MACS2_qvalue+Replicate~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Peaks Reproduced") +
  xlab("")+
  ggtitle("MACS2 replicate reproducibility HBD.His")

fig7C

#WKK
peakpathWKK <- paste0(projPath, 
                       "/peakCalling/His_test/SEACR_relaxed_WKKvHBD_peakcalling_summary.txt")
peakWKK <- read.table(peakpathWKK, header=T)
peakWKK <- peakWKK[peakWKK$peakType %in% c("control", "control.rmDup"),]

fig7A = peakWKK %>% ggplot(aes(x = IggType, y = peakN, fill = IggType)) +
  geom_boxplot() +
  geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  facet_grid(~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") +
  xlab("")+
  ggtitle("SEACR norm relaxed number of peaks WKK.His")

fig7A

fig7C = peakWKK %>% ggplot(aes(x = IggType, y = peakReprodRate, fill = IggType, label = round(peakReprodRate, 2))) +
  geom_bar(stat = "identity") +
  geom_text(vjust = 0.1) +
  facet_grid(Replicate~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Peaks Reproduced") +
  xlab("")+
  ggtitle("SEACR replicate reproducibility WKK.His")

fig7C

peakpathrpa70 <- paste0(projPath, 
                      "/peakCalling/His_test/SEACR_relaxed_rpa70_peakcalling_summary.txt")
peakrpa <- read.table(peakpathrpa70, header=T)
peakrpa <- peakrpa[peakrpa$peakType %in% c("control", "control.rmDup"),]

fig7A = peakrpa %>% ggplot(aes(x = Histone, y = peakN, fill = IggType)) +
  geom_boxplot() +
  geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  facet_grid(IggType~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") +
  xlab("")+
  ggtitle("SEACR norm relaxed number of RPA70 peaks")

fig7A

fig7C = peakrpa %>% ggplot(aes(x = Histone, y = peakReprodRate, fill = Histone, label = round(peakReprodRate, 2))) +
  geom_bar(stat = "identity") +
  geom_text(vjust = 0.1) +
  facet_grid(Replicate~peakType+ IggType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Peaks Reproduced") +
  xlab("")+
  ggtitle("SEACR replicate reproducibility rpa70")

fig7C

