################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 6. Peak calling with SEACR                                                   #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: peak calling summary plotting                                       #
################################################################################

#-----------------------Paths and libraries------------------------------------#
library(ggplot2)
library(viridis)
library(dplyr)
library(ggpubr)
library(RColorBrewer)

#--------------------------Peak calling assessment-----------------------------#

projectPath <- "/project/ChromGroup_Seq_data/Celeste/2024_summer_complete/" #fill this variable with the same path as PROJECTROOT in 
# the cut_and_tag_pipeline.sh script

peakpath <- paste0(projectPath, "/alignment/summary_peak_calling_norm.txt")
peakWidthPath <- paste0(projectPath, "/alignment/summary_peak_calling_norm_width_AmpR.txt")
fripsPath <- paste0(projectPath, "/alignment/summary_peak_calling_frips.txt")
peak <- read.table(peakpath, header=T)
peakWidth <- read.table(peakWidthPath, header=T)
frips <- read.table(fripsPath, header=T)
peak[is.na(peak)] <- 0.0

peak$Histone
split_results <- strsplit(peak$Histone, "[.]")
split_results
peak$cellLine <- sapply(split_results, `[`, 1)

split_results <- strsplit(frips$Histone, "[.]")
split_results
frips$cellLine <- sapply(split_results, `[`, 1)

fig7A = peak %>% ggplot(aes(x = Histone, y = peakN, fill = as.factor(cellLine))) +
  geom_boxplot() +
  geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  facet_wrap(~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
fig7A + ggtitle("Number of peaks SEACR norm relaxed")

fig7B = peakWidth %>% ggplot(aes(x = Histone, y = width, fill = Histone)) +
  geom_violin() +
  facet_grid(Replicate~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_y_continuous(trans = "log", breaks = c(400, 3000, 22000)) +
  theme_bw(base_size = 18) +
  ylab("Width of Peaks") +
  xlab("")
fig7B +ggtitle("Number of peaks SEACR norm relaxed")

fig7C = peak %>% ggplot(aes(x = Histone, y = peakReprodRate, fill =as.factor(cellLine), label = round(peakReprodRate, 2))) +
  geom_bar(stat = "identity") +
  geom_text(vjust = 0.1) +
  facet_grid(Replicate~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Peaks Reproduced") +
  xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

fig7C + ggtitle("Percentage of peaks overlapping replicate SEACR norm relaxed")
fig7D = frips %>% ggplot(aes(x = Histone, y = frip, fill = as.factor(cellLine), label = round(frip, 2))) +
  geom_boxplot() +
  geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Fragments in Peaks") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

fig7D
ggarrange(fig7A, fig7B, fig7C, fig7D, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")
