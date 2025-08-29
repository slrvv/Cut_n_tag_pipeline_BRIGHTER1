################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 6. Peak calling with SEACR                                                   #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: peak calling summary for Macs                                       #
################################################################################


#-------------------------Paths------------------------------------------------#
library(dplyr)
library(GenomicRanges, quiet = TRUE)
library(BRGenomics)
args <- commandArgs(trailingOnly=TRUE)
projPath <- args[1]
summaryPath <- args[2]


#------------------------Sequencing depth--------------------------------------#

sampletable <- read.table(summaryPath,
                          header = T, sep = ",")


sampleList <- sampletable$SampleName
sampleList
peakWidth = c()
peakN = c()
peakType = c("_peaks.broadPeak", "_peaks.narrowPeak", "_peaks.gappedPeak", 
            "_rmDup_peaks.broadPeak", "_rmDup_peaks.narrowPeak", "_rmDup_peaks.gappedPeak")
for(hist in sampleList){
  histInfo = unlist(strsplit(hist, "_", fixed = T))
  if(histInfo[1] != "IgG"){
    for(type in peakType){
      file <- paste0(projPath, "/peakCalling/MACS2/", hist, type)
      print(file)
      if (file.info(file)$size != 0){
        peakInfo = read.table(paste0(projPath, "/peakCalling/MACS2/", hist, type),
                              header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
        peakN = data.frame(peakN = nrow(peakInfo), peakType = type, 
                           Histone = histInfo[1], 
                           Replicate = paste(histInfo[-1], collapse = "_")) %>% rbind(peakN, .)
        peakWidth = data.frame(width = peakInfo$width, 
                               peakType = type, 
                               Histone = histInfo[1], 
                               Replicate = paste(histInfo[-1], collapse = "_"))  %>% rbind(peakWidth, .)
      } else {
        peakN = data.frame(peakN = 0, peakType = type, 
                           Histone = histInfo[1], 
                           Replicate = paste(histInfo[-1], collapse = "_")) %>% rbind(peakN, .)
        peakWidth = data.frame(width = 0, 
                               peakType = type, 
                               Histone = histInfo[1], 
                               Replicate = paste(histInfo[-1], collapse = "_"))  %>% rbind(peakWidth, .)
      }
      
      
    }
  }
}


write.table(peakWidth, paste0(projPath, 
                              "/alignment/summary_peak_calling_width_MACS2.txt"),
            row.names = F)

peakN %>% select(Histone, Replicate, peakType, peakN)

peakOverlap = c()

for(hist in sampleList){
  for(type in peakType){
    overlap.gr = GRanges()
    histInfo = unlist(strsplit(hist, "_", fixed = T))
    repL <- unique(peakN$Replicate[peakN$Histone == histInfo[1]])
    if(length(repL) > 1){
      for(rep in repL){
        file <- paste0(projPath,
                       "/peakCalling/MACS2/",
                       histInfo[1],
                       "_",
                       rep,
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
                             peakType = type) %>% rbind(peakOverlap, .)
  }
}


peakReprod = left_join(peakN, peakOverlap, by = c("Histone", "Replicate", "peakType")) %>% mutate(peakReprodRate = (peakReprod/peakN) * 100)
peakReprod
peakReprod %>% select(Histone, Replicate, peakType, peakN, peakReprodNum = peakReprod, peakReprodRate)


write.table(peakReprod, paste0(projPath, 
                               "/alignment/summary_peak_calling_MACS2.txt"),
            row.names = F)
