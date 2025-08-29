################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 6. Merge replicates                                                          #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: merge peak files                                                     #
################################################################################

#------------------------------Paths-------------------------------------------#

projPath=$1
histName=$2
#-------------------------Peak calling-----------------------------------------#

bedtools intersect -a $projPath/peakCalling/SEACR/${histName}_1_seacr_norm_control.peaks.relaxed.bed \
-b $projPath/peakCalling/SEACR/${histName}_2_seacr_norm_control.peaks.relaxed.bed > \
$projPath/peakCalling/SEACR/${histName}_merged_seacr_norm_control.peaks.relaxed.bed

bedtools intersect -a $projPath/peakCalling/SEACR/${histName}_1_seacr_norm_control.rmDup.peaks.relaxed.bed \
-b $projPath/peakCalling/SEACR/${histName}_2_seacr_norm_control.rmDup.peaks.relaxed.bed > \
$projPath/peakCalling/SEACR/${histName}_merged_seacr_norm_control.rmDup.peaks.relaxed.bed

