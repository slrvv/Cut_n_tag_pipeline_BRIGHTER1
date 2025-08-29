################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 6. Peak calling with SEACR                                                   #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: peak calling summary FRips Calculation                              #
################################################################################


#-------------------------Paths------------------------------------------------#
library(dplyr, quietly = TRUE)
library(chromVAR, quietly = TRUE)
library(GenomicRanges, quietly = TRUE)
args <- commandArgs(trailingOnly=TRUE)
projPath <- args[1]
summaryPath <- args[2]



#------------------------Sequencing depth--------------------------------------#


sampletable <- read.table(summaryPath,
                          header = T, sep = ",")
print(summaryPath)
alignResult <- read.table(paste0(projPath, 
                                 "/alignment/summary_seq_depth_all_experiments.txt"), 
                          header=T, sep = " ")

alignResult$Replicate <- as.character(alignResult$Replicate)
sampleList <- sampletable$Sample


bamDir = paste0(projPath, "/alignment/bam")
inPeakData = c()
sampleList

for(hist in sampleList){
  ## overlap with bam file to get count for the not deduplicated peaks
  file <- paste0(projPath, "/peakCalling/SEACR/AmpR_calibrated/", 
         hist,"_seacr_norm_control.peaks.AmpR.relaxed.bed")
  print(hist)
  print(file)
  if (file.info(file)$size !=0 ){
    peakRes = read.table(file, header = FALSE, fill = TRUE)
    peak.gr = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*")
    bamFile = paste0(bamDir, "/", hist,".mapped.sorted.filtered.bam")
    histInfo <- unlist(strsplit(hist, "_", fixed = T))
    fragment_counts <- getCounts(bamFile, peak.gr, paired = TRUE, by_rg = FALSE, format = "bam")
    inPeakN = counts(fragment_counts)[,1] %>% sum
    
  } else { 
    histInfo <- unlist(strsplit(hist, "_", fixed = T))
    inPeakN <- 0
  }
  


  inPeakData = rbind(inPeakData, data.frame(inPeakN = inPeakN,
                                            Histone = histInfo[1],
                                            Replicate = paste(histInfo[-1], collapse = "_")))

}


frip = left_join(inPeakData, alignResult, by = c("Histone", "Replicate")) %>% 
  mutate(frip = inPeakN/MappedFragNum_hg38 * 100, )
frip %>% select(Histone, Replicate, SequencingDepth, MappedFragNum_hg38, AlignmentRate_hg38, FragInPeakNum = inPeakN, FRiPs = frip)
write.table(frip, paste0(projPath,
                               "/alignment/summary_peak_calling_frips_ampr.txt"),
            row.names = F)
