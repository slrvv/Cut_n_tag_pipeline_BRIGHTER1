################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 6. Remove blacklist from peaks                                               #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: remove blacklist                                                    #
################################################################################

#------------------------------Paths-------------------------------------------#

projPath=$1
histName=$2
blacklist=$3
#-------------------------Peak calling-----------------------------------------#

bedtools intersect -v -a $projPath/peakCalling/SEACR/${histName}_merged_seacr_norm_control.peaks.relaxed.filtered.bed \
-b $projPath/GRCh38_unified_blacklist.bed > $projPath/peakCalling/SEACR/${histName}_merged_seacr_norm_control.peaks.relaxed.filtered.rmblck.bed


