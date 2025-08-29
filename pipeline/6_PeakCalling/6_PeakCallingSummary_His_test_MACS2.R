################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 6. Peak calling for Histidine test                                           #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: MACS2 peak calling summary Histidine test                           #
################################################################################


#-------------------------Paths------------------------------------------------#
library(dplyr)
library(GenomicRanges, quiet = TRUE)
library(BRGenomics)

projPath <- "/project/ChromGroup_Seq_data/Celeste/2024_summer_complete"
summaryPath <- paste0(projPath, "/experiment_summary_peaks_same_cell.csv")


#------------------------Sequencing depth--------------------------------------#

sampletable <- read.table(summaryPath,
                          header = T, sep = ",")


sampleList <- sampletable$SampleName
sampleList
peakWidth = c()
peakN = c()
qvalues = c("0.01", "0.05", "0.1")
IggType = c("same_cell", "diff_cell", "merged_cell")
peakType = c("_peaks.broadPeak", "_peaks.narrowPeak", "_peaks.gappedPeak", 
             "_rmDup_peaks.broadPeak", "_rmDup_peaks.narrowPeak", "_rmDup_peaks.gappedPeak")
for(hist in sampleList){
  histInfo = unlist(strsplit(hist, "_", fixed = T))
  if(histInfo[1] != "IgG"){
    for(value in qvalues){
    for(igg in IggType){
    for(type in peakType){
      file <- paste0(projPath, "/peakCalling/His_test/", hist, "_", value, "_",
                     igg,type)
      print(file)
      if (file.info(file)$size != 0){
        peakInfo = read.table(file,
                              header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
        peakN = data.frame(peakN = nrow(peakInfo), peakType = type, 
                           IggType = igg,
                           MACS2_qvalue = value, 
                           Histone = histInfo[1], 
                           Replicate = paste(histInfo[-1], collapse = "_")) %>% rbind(peakN, .)
      } else {
        peakN = data.frame(peakN = 0, peakType = type, 
                           IggType = igg,
                           MACS2_qvalue = value, 
                           Histone = histInfo[1], 
                           Replicate = paste(histInfo[-1], collapse = "_")) %>% rbind(peakN, .)
      }
      
      
    }
    }
    }
  }
}



peakN %>% select(Histone, Replicate, peakType, IggType, MACS2_qvalue, peakN)
peakN
peakOverlap = c()

for(hist in sampleList){
  for(igg in IggType){
  for(value in qvalues){
  for(type in peakType){
    overlap.gr = GRanges()
    histInfo = unlist(strsplit(hist, "_", fixed = T))
    repL <- unique(peakN$Replicate[peakN$Histone == histInfo[1]])
    if(length(repL) > 1){
      for(rep in repL){
        file <- paste0(projPath,
                       "/peakCalling/His_test/",
                       histInfo[1],
                       "_",
                       rep,
                       "_",
                       value,
                       "_",
                       igg,
                       type)
        print(file)
        if (file.exists(file) & file.info(file)$size != 0){
          peakInfo = read.table(file, header = FALSE, fill = TRUE)
          peakInfo.gr = GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
          if(length(overlap.gr) >0){
            overlap.gr = subsetByOverlaps(overlap.gr, peakInfo.gr)
            peakReprod = length(overlap.gr)
          }else{
            overlap.gr = peakInfo.gr
            
          }
          
        } else {
          peakReprod <- NaN
        }
        
      } 
    } else {
      peakReprod <- NaN
    }
    peakOverlap = data.frame(peakReprod = peakReprod, 
                             Histone = histInfo[1],
                             Replicate = histInfo[-1],
                             IggType = igg,
                             MACS2_qvalue = value, 
                             peakType = type) %>% rbind(peakOverlap, .)
  }
  }
  }

}


peakOverlap
peakReprod = left_join(peakN, peakOverlap, 
                       by = c("Histone", 
                              "Replicate", 
                              "IggType", 
                              "MACS2_qvalue", 
                              "peakType")) %>% mutate(peakReprodRate = (peakReprod/peakN) * 100)
peakReprod
peakReprod %>% select(Histone, Replicate, IggType, MACS2_qvalue, peakType, peakN, peakReprodNum = peakReprod, peakReprodRate)


write.table(peakReprod, paste0(projPath, 
                               "/peakCalling/His_test/His_test_summary_peak_calling_MACS2.txt"),
            row.names = F)
