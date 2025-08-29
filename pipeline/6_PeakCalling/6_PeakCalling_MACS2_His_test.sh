################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 6. Peak calling with MACS2 Histidine test                                    #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: peak calling with macs2                                             #
################################################################################

#------------------------------Paths-------------------------------------------#

projPath=$1
histControl=$2
histName=$3
dupRemove=$4
qvalues=($5)
summary=$6



#---------------------------Peak Calling---------------------------------------#

outdir=$projPath/peakCalling/His_test
bamdir=$projPath/alignment/bam

#g is 2.7e9 which is the parameter for human genome
for value in "${qvalues[@]}"
do
echo $value
echo "Calling narrow peaks with MACS2"
macs2 callpeak -t $bamdir/${histName}.mapped.bam \
-c $bamdir/${histControl}.mapped.bam \
-f BAM -g hs -q $value \
-n ${histName}_${value}_${summary} \
--outdir $outdir 2> $outdir/${histName}.${value}.${summary}.macs2.log

echo "calling broad peaks with MACS2"
macs2 callpeak -t $bamdir/${histName}.mapped.bam \
-c $bamdir/${histControl}.mapped.bam \
-f BAM -g hs -q $value --broad \
-n ${histName}_${value}_${summary} \
--outdir $outdir 2> $outdir/${histName}.${value}.${summary}.broad.macs2.log

if [ "$4" = true ]; then
echo "calling narrow peaks on deduplicated"
macs2 callpeak -t $bamdir/${histName}.mapped.rmDup.bam \
-c $bamdir/${histControl}.mapped.rmDup.bam \
-f BAM -g hs -q $value \
-n ${histName}_${value}_${summary}_rmDup \
--outdir $outdir 2> $outdir/${histName}.${value}.${summary}.rmDup.macs2.log
echo "calling broad peaks on deduplicated"
macs2 callpeak -t $bamdir/${histName}.mapped.rmDup.bam \
-c $bamdir/${histControl}.mapped.rmDup.bam \
-f BAM -g hs -q $value --broad \
-n ${histName}_${value}_${summary}_rmDup \
--outdir $outdir 2> $outdir/${histName}.${value}.${summary}.rmDup.broad.macs2.log
fi

done

