################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 6. Peak calling with MACS2                                                   #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: peak calling with macs2                                             #
################################################################################

#------------------------------Paths-------------------------------------------#

projPath=$1
histControl=$2
histName=$3
dupRemove=$4
qvalues=($5)

mkdir -p $projPath/peakCalling/MACS2

#---------------------------Peak Calling---------------------------------------#

outdir=$projPath/peakCalling/MACS2
bamdir=$projPath/alignment/bam

#g is 2.7e9 which is the parameter for human genome
for value in "${qvalues[@]}"
do
  echo $value
  echo "Calling narrow peaks with MACS2"
  macs2 callpeak -t $bamdir/${histName}.mapped.bam \
  	-c $bamdir/${histControl}.mapped.bam \
   	-f BAM -g hs -q $value \
  	-n ${histName}_$value \
  	--outdir $outdir 2> $outdir/${histName}.${value}.macs2.log

  echo "calling broad peaks with MACS2"
  macs2 callpeak -t $bamdir/${histName}.mapped.bam \
  	-c $bamdir/${histControl}.mapped.bam \
   	-f BAM -g hs -q $value --broad \
  	-n ${histName}_value \
  	--outdir $outdir 2> $outdir/${histName}.${value}.broad.macs2.log

  if [ "$4" = true ]; then
    echo "calling narrow peaks on deduplicated"
    macs2 callpeak -t $bamdir/${histName}.mapped.rmDup.bam \
  	-c $bamdir/${histControl}.mapped.rmDup.bam \
   	-f BAM -g hs -q $value \
  	-n ${histName}_${value}_rmDup \
  	--outdir $outdir 2> $outdir/${histName}.${value}.rmDup.macs2.log
    echo "calling broad peaks on deduplicated"
    macs2 callpeak -t $bamdir/${histName}.mapped.rmDup.bam \
  	-c $bamdir/${histControl}.mapped.rmDup.bam \
   	-f BAM -g hs -q $value --broad \
  	-n ${histName}_${value}_rmDup \
  	--outdir $outdir 2> $outdir/${histName}.${value}.rmDup.broad.macs2.log
  fi

done

