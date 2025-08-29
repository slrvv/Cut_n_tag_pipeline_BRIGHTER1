################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 6. Peak calling for Histidine test                                           #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: SEACR peak calling summary Histidine test                           #
################################################################################


#-------------------------Paths------------------------------------------------#
library(dplyr)
library(GenomicRanges, quiet = TRUE)
library(BRGenomics)

projPath <- "/project/ChromGroup_Seq_data/Celeste/2024_summer_complete"



#-------------------------Functions--------------------------------------------#

peak_number_comp <- function(sampleList, projPath, peakType, IggTyp, norm, relaxed){
  peakN = c()
  for(hist in sampleList){
    histInfo = unlist(strsplit(hist, "_", fixed = T))
    if(histInfo[1] != "IgG"){
      for (igg in IggType){
        for(type in peakType){
          suffix <-  ".peaks.stringent.bed"
          if(relaxed == TRUE) suffix <- ".peaks.relaxed.bed";
          if (norm == TRUE){
            file <- paste0(projPath, "/peakCalling/His_test/", hist, "_", igg,
                           "_seacr_norm_", type, suffix)
          } else {
            file <- paste0(projPath, "/peakCalling/His_test/", hist, "_", igg,
                           "_seacr_", type, suffix)
          }
          print(file)
          if (file.info(file)$size != 0){
            peakInfo = read.table(file,
                                  header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
            print(nrow(peakInfo))
            peakN = data.frame(peakN = nrow(peakInfo), IggType = igg, peakType = type, 
                               Histone = histInfo[1], 
                               Replicate = paste(histInfo[-1], collapse = "_")) %>% rbind(peakN, .)
          } else {
            peakN = data.frame(peakN = 0, peakType = type, 
                               Histone = histInfo[1], 
                               Replicate = paste(histInfo[-1], collapse = "_")) %>% rbind(peakN, .)
            
          }
          
        }
        
        
        
      }
    }
  }
  
  peakN %>% select(Histone, Replicate, peakType, IggType, peakN)
  return(peakN)
  
}



peak_rep <- function(sampleList, projPath, peakType, IggType, peakN, norm, relaxed){
  peakOverlap = c()
  
  for(hist in sampleList){
    for (igg in IggType){
      for(type in peakType){
        overlap.gr = GRanges()
        histInfo = unlist(strsplit(hist, "_", fixed = T))
        repL <- unique(peakN$Replicate[peakN$Histone == histInfo[1]])
        if(length(repL) > 1){
          for(rep in repL){
            suffix <-  ".peaks.stringent.bed"
            if(relaxed == TRUE) suffix <- ".peaks.relaxed.bed";
            if (norm == TRUE){
              file <- paste0(projPath,
                             "/peakCalling/His_test/",
                             histInfo[1],
                             "_",
                             rep,
                             "_",
                             igg,
                             "_seacr_norm_",
                             type,
                             suffix)
            } else {
              file <- paste0(projPath,
                             "/peakCalling/His_test/",
                             histInfo[1],
                             "_",
                             rep,
                             "_",
                             igg,
                             "_seacr_",
                             type,
                             suffix)
              
            }
            cat(paste0(file, "\n"))
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
                                 peakType = type) %>% rbind(peakOverlap, .)

      }
    }
  }
  
  peakReprod = left_join(peakN, peakOverlap, by = c("Histone", "Replicate", "peakType", "IggType")) %>% mutate(peakReprodRate = (peakReprod/peakN) * 100)
  
  
  peakReprod = peakReprod %>% select(Histone, Replicate, peakType, IggType, peakN, peakReprodNum = peakReprod, peakReprodRate)
  
  if (norm == TRUE){
    parameter <- rep("norm", times = nrow(peakReprod) )
    
    peakReprod$SEACR_parameter <- parameter
    
    return(peakReprod)
  } else {
    parameter <- rep("non", times = nrow(peakReprod) )
    
    peakReprod$SEACR_parameter <- parameter
    
    return(peakReprod)
  }
  
}



#------------------------Sequencing depth--------------------------------------#

## HBD.his 
summaryPath <- paste0(projPath, "/experiment_summary_peaks_same_cell.csv")
sampletable <- read.table(summaryPath,
                          header = T, sep = ",")
sampleList <- sampletable$SampleName
sampleList

peakType = c("control", "top0.01", "control.rmDup", "top0.01.rmDup")
IggType = c("same_cell", "diff_cell", "merged_cell")

#NON normalization
peakNnon <- peak_number_comp(sampleList, 
                             projPath, 
                             peakType,
                             IggTyp,
                             norm = FALSE, 
                             relaxed = FALSE)

peakReprodnon <- peak_rep(sampleList, 
                          projPath,
                          peakType, 
                          IggType,
                          peakNnon, 
                          norm=FALSE)

#with norm parameter
peakNnorm <- peak_number_comp(sampleList,
                              projPath,
                              peakType, 
                              IggTyp,
                              norm = TRUE, 
                              relaxed = FALSE)

peakReprodnorm <- peak_rep(sampleList, 
                           projPath,
                           peakType, 
                           IggType, 
                           peakNnorm, 
                           norm=TRUE, 
                           relaxed = FALSE)


peakReprodall <- rbind(peakReprodnon, peakReprodnorm)
peakReprodall

write.table(peakReprodall, paste0(projPath, 
                                  "/peakCalling/His_test/SEACR_His_test_peakcalling_summary.txt"),
            row.names = F)

#norm and relaxed

peakN <- peak_number_comp(sampleList,
                              projPath,
                              peakType, 
                              IggTyp,
                              norm = TRUE, 
                              relaxed = TRUE)

peakReprod <- peak_rep(sampleList, 
                           projPath,
                           peakType, 
                           IggType, 
                           peakN, 
                           norm=TRUE, 
                           relaxed = TRUE)



write.table(peakReprod, paste0(projPath, 
                               "/peakCalling/His_test/SEACR_relaxed_His_test_peakcalling_summary.txt"),
            row.names = F)


###WKK.his vs WKK.igg and HBD.igg

summaryPath <- paste0(projPath, "/experiment_summary_peaks_same_cell_WKK.csv")
sampletable <- read.table(summaryPath,
                          header = T, sep = ",")
sampleList <- sampletable$SampleName

sampleList
IggType = c("same_cell_WKK", "diff_cell_WKK")

peakN <- peak_number_comp(sampleList,
                          projPath,
                          peakType, 
                          IggTyp,
                          norm = TRUE, 
                          relaxed = TRUE)

peakReprod <- peak_rep(sampleList, 
                       projPath,
                       peakType, 
                       IggType, 
                       peakN, 
                       norm=TRUE, 
                       relaxed = TRUE)
peakReprod

peakN <- peak_number_comp(sampleList,
                          projPath,
                          peakType, 
                          IggTyp,
                          norm = TRUE, 
                          relaxed = FALSE)

peakReprod <- peak_rep(sampleList, 
                       projPath,
                       peakType, 
                       IggType, 
                       peakN, 
                       norm=TRUE, 
                       relaxed = FALSE)
peakReprod


write.table(peakReprod, paste0(projPath, 
                               "/peakCalling/His_test/SEACR_relaxed_WKKvHBD_peakcalling_summary.txt"),
            row.names = F)


## rpa_70 HBD vs WT

summaryPath <- paste0(projPath, "/rpa70_summary_peaks_same_cell.csv")
sampletable <- read.table(summaryPath,
                          header = T, sep = ",")
sampleList <- sampletable$SampleName

sampleList
IggType = c("same_cell", "diff_cell")

peakN <- peak_number_comp(sampleList,
                          projPath,
                          peakType, 
                          IggTyp,
                          norm = TRUE, 
                          relaxed = TRUE)

peakReprod <- peak_rep(sampleList, 
                       projPath,
                       peakType, 
                       IggType, 
                       peakN, 
                       norm=TRUE, 
                       relaxed = TRUE)
peakReprod

peakN <- peak_number_comp(sampleList,
                          projPath,
                          peakType, 
                          IggTyp,
                          norm = TRUE, 
                          relaxed = FALSE)

peakReprod <- peak_rep(sampleList, 
                       projPath,
                       peakType, 
                       IggType, 
                       peakN, 
                       norm=TRUE, 
                       relaxed = FALSE)
peakReprod
sampleList

write.table(peakReprod, paste0(projPath, 
                               "/peakCalling/His_test/SEACR_relaxed_rpa70_peakcalling_summary.txt"),
            row.names = F)
