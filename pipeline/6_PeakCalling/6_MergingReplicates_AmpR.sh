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

bedtools intersect -a $projPath/peakCalling/SEACR/AmpR_calibrated/${histName}_1_seacr_norm_control.peaks.AmpR.relaxed.bed \
-b $projPath/peakCalling/SEACR/AmpR_calibrated/${histName}_2_seacr_norm_control.peaks.AmpR.relaxed.bed > \
$projPath/peakCalling/SEACR/AmpR_calibrated/${histName}_merged_seacr_norm_control.peaks.AmpR.relaxed.bed

