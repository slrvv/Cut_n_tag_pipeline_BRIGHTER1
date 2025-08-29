################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 6. Peak calling with SEACR                                                   #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: peak calling summary                                                #
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
peakType = c("control")
for(hist in sampleList){
  histInfo = unlist(strsplit(hist, "_", fixed = T))
  if(histInfo[1] != "IgG"){
    for(type in peakType){
      file <- paste0(projPath, "/peakCalling/SEACR/AmpR_calibrated/", hist, 
             "_seacr_norm_", type, ".peaks.AmpR.relaxed.bed")
      print(file)
      if (file.info(file)$size != 0){
        peakInfo = read.table(paste0(projPath, "/peakCalling/SEACR/AmpR_calibrated/", hist, 
                                     "_seacr_norm_", type, ".peaks.AmpR.relaxed.bed"),
                              header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
        print(nrow(peakInfo))
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
                               "/alignment/summary_peak_calling_norm_width_AmpR.txt"),
            row.names = F)

peakN %>% select(Histone, Replicate, peakType, peakN)

peakOverlap = c()

peakOverlap = c()
for(hist in unique(peakN$Histone)){
  print(hist)
  for(type in peakType){
    repL <- unique(peakN$Replicate[peakN$Histone == hist])
    print(length(repL))
    for(rep in repL){  
      if(length(repL) > 1){
        if (rep == "1"){
          file <- paste0(projPath,
                         "/peakCalling/SEACR/AmpR_calibrated/",
                         hist,
                         "_",
                         rep,
                         "_seacr_norm_",
                         type,
                         ".peaks.AmpR.relaxed.bed")
          repfile <- paste0(projPath,
                            "/peakCalling/SEACR/AmpR_calibrated/",
                            hist,
                            "_",
                            "2",
                            "_seacr_norm_",
                            type,
                            ".peaks.AmpR.relaxed.bed")
          
        } else {
          file <- paste0(projPath,
                         "/peakCalling/SEACR/AmpR_calibrated/",
                         hist,
                         "_",
                         rep,
                         "_seacr_norm_",
                         type,
                         ".peaks.AmpR.relaxed.bed")
          repfile <- paste0(projPath,
                            "/peakCalling/SEACR/AmpR_calibrated/",
                            hist,
                            "_",
                            "1",
                            "_seacr_norm_",
                            type,
                            ".peaks.AmpR.relaxed.bed")
          
          
        }
        
        if (file.exists(file) & file.info(file)$size != 0 & file.exists(repfile) & file.info(repfile)$size != 0){
          #load granges of the replicate in which we centered
          peakInfo = read.table(file, header = FALSE, fill = TRUE)
          peakInfo.gr = GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
          #load granges of the other replicate
          peakInforep = read.table(repfile, header = FALSE, fill = TRUE)
          peakInforep.gr = GRanges(peakInforep$V1, IRanges(start = peakInforep$V2, end = peakInforep$V3), strand = "*")
          
          #overlap the granges
          overlap.gr = subsetByOverlaps(peakInfo.gr, peakInforep.gr)
          peakReprod = length(overlap.gr)
          
          
        } else {
          peakReprod <- NaN
        }
      } else {
        peakReprod <- NaN
      }
      
      peakOverlap = data.frame(peakReprod = peakReprod, 
                               Histone = hist,
                               Replicate = rep, 
                               peakType = type) %>% rbind(peakOverlap, .)
    }
  }
}
peakReprod = left_join(peakN, peakOverlap, by = c("Histone", "Replicate", "peakType")) %>% mutate(peakReprodRate = (peakReprod/peakN) * 100)
peakReprod
peakReprod %>% select(Histone, Replicate, peakType, peakN, peakReprodNum = peakReprod, peakReprodRate)


write.table(peakReprod, paste0(projPath, 
                                 "/alignment/summary_peak_calling_norm_AmpR.txt"),
            row.names = F)
