################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 6. Peak calling with SEACR Histidine test                                    #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: peak calling                                                        #
################################################################################

#------------------------------Paths-------------------------------------------#

seacr="/project/ChromGroup/Serkan_Project/cut_and_tag_analysis/tools/SEACR/SEACR_1.3.sh"
projPath=$1
histControl=$2
histName=$3
dupRemove=$4
norm=$5
summary=$6
relaxed=$7


#-------------------------Peak calling-----------------------------------------#

echo "Peak calling on normalized experiments"

# bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph \
# $projPath/alignment/bedgraph/${histControl}_bowtie2.fragments.normalized.bedgraph \
# non stringent $projPath/peakCalling/His_test/${histName}_${summary}_seacr_control.peaks
# 
# bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph \
# 0.01 non stringent $projPath/peakCalling/His_test/${histName}_${summary}_seacr_top0.01.peaks
# 
# if [ "$dupRemove" = true ]; then
# echo "Peak calling on normalized and deduplicated experiments"
# bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.fragments.normalized.bedgraph \
# $projPath/alignment/bedgraph/${histControl}_bowtie2.rmDup.fragments.normalized.bedgraph \
# non stringent $projPath/peakCalling/His_test/${histName}_${summary}_seacr_control.rmDup.peaks
# 
# bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.fragments.normalized.bedgraph \
# 0.01 non stringent $projPath/peakCalling/His_test/${histName}_${summary}_seacr_top0.01.rmDup.peaks
# fi
# 
# 
# ## produce peak files with norm parameter
if [ "$norm" = true ]; then
  bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph \
  $projPath/alignment/bedgraph/${histControl}_bowtie2.fragments.normalized.bedgraph \
  norm stringent $projPath/peakCalling/His_test/${histName}_${summary}_seacr_norm_control.peaks

  bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph \
  0.01 norm stringent $projPath/peakCalling/His_test/${histName}_${summary}_seacr_norm_top0.01.peaks

  if [ "$dupRemove" = true ]; then
  echo "Peak calling on normalized and deduplicated experiments"
  bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.fragments.normalized.bedgraph \
  $projPath/alignment/bedgraph/${histControl}_bowtie2.rmDup.fragments.normalized.bedgraph \
  norm stringent $projPath/peakCalling/His_test/${histName}_${summary}_seacr_norm_control.rmDup.peaks

  bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.fragments.normalized.bedgraph \
  0.01 norm stringent $projPath/peakCalling/His_test/${histName}_${summary}_seacr_norm_top0.01.rmDup.peaks
  fi
fi

if [ "$relaxed" = true ]; then
  bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph \
  $projPath/alignment/bedgraph/${histControl}_bowtie2.fragments.normalized.bedgraph \
  norm relaxed $projPath/peakCalling/His_test/${histName}_${summary}_seacr_norm_control.peaks

  bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph \
  0.01 norm relaxed $projPath/peakCalling/His_test/${histName}_${summary}_seacr_norm_top0.01.peaks

  if [ "$dupRemove" = true ]; then
  echo "Peak calling on normalized and deduplicated experiments"
  bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.fragments.normalized.bedgraph \
  $projPath/alignment/bedgraph/${histControl}_bowtie2.rmDup.fragments.normalized.bedgraph \
  norm relaxed $projPath/peakCalling/His_test/${histName}_${summary}_seacr_norm_control.rmDup.peaks

  bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.fragments.normalized.bedgraph \
  0.01 norm relaxed $projPath/peakCalling/His_test/${histName}_${summary}_seacr_norm_top0.01.rmDup.peaks
  fi
fi
