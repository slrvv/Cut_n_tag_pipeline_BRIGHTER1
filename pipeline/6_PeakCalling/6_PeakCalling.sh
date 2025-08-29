################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 6. Peak calling with SEACR                                                   #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: peak calling                                                        #
################################################################################

#------------------------------Paths-------------------------------------------#

seacr="/project/ChromGroup/Serkan_Project/cut_and_tag_analysis/tools/SEACR/SEACR_1.3.sh"
projPath=$1
histControl=$2
histName=$3
dupRemove=$4
mkdir -p $projPath/peakCalling/SEACR


#-------------------------Peak calling-----------------------------------------#

echo "Peak calling on normalized experiments"


## produce peak files with norm parameter
 
bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph \
$projPath/alignment/bedgraph/${histControl}_bowtie2.fragments.normalized.bedgraph \
norm relaxed $projPath/peakCalling/SEACR/${histName}_seacr_norm_control.peaks
  
if [ "$dupRemove" = true ]; then
  echo "Peak calling on normalized and deduplicated experiments"
  bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.fragments.normalized.bedgraph \
  $projPath/alignment/bedgraph/${histControl}_bowtie2.rmDup.fragments.normalized.bedgraph \
  norm relaxed $projPath/peakCalling/SEACR/${histName}_seacr_norm_control.rmDup.peaks
    
fi