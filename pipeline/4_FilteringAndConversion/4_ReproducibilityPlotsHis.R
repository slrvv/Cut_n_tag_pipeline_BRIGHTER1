################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 4. Reproducibility assessment histidine                                      #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: create table for reproducibility assessment plot                    #
################################################################################

library(dplyr, quietly=TRUE)
library(corrplot)


projPath <- "/project/ChromGroup_Seq_data/Celeste/2024_summer_complete"
summaryPath <- "/project/ChromGroup_Seq_data/Celeste/2024_summer_complete/experiment_summary_align_formatted.csv"
sampletable <- read.table(summaryPath,
                          header = T, sep = ",")


cat("Reproducibility of the replicates\n")

sampleList <- sort(sampletable$SampleName)

sampleListHBD <- c("HBD.his_1","HBD.his_2","HBD.igg_1","HBD.igg_2","WKK.his_1","WKK.his_2","WKK.igg_1","WKK.igg_2")
sampleListHBD
sampleListHBD <- sampleListHBD[!grepl("WT*", sampleListHBD)]
sampleListHBD
reprod = c()
fragCount = NULL
for(hist in sampleListHBD){
  
  if(is.null(fragCount)){
    print(NULL)
    print(hist)
    fragCount = read.table(paste0(projPath, "/alignment/bed/", hist, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE) 
    colnames(fragCount) = c("chrom", "bin", hist)
    print(nrow(fragCount))
  }else{
    print(hist)
    fragCountTmp = read.table(paste0(projPath, "/alignment/bed/", hist, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE)
    colnames(fragCountTmp) = c("chrom", "bin", hist)
    fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
    print(head(fragCount))
  }
}
fragCofragCount$bin
length(fragCountHBD)
is.na(fragCount$HBD.flag_1)
fragCountHBD <- fragCount$bin

M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs")


corrplot(M, method = "color", outline = T, addgrid.col = "darkgray", order="hclust", addrect = 3, rect.col = "black", rect.lwd = 3,cl.pos = "b", tl.col = "indianred4", tl.cex = 1, cl.cex = 1, addCoef.col = "black", number.digits = 2, number.cex = 1, col = colorRampPalette(c("midnightblue","white","darkred"))(100))
write.table(M, paste0(projPath,"/alignment/rep_reproducibility_all_experiments.txt"),
            row.names = T)
