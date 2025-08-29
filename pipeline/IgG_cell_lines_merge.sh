################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# IgG cell lines merge                                                         #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: merge all IgG for each cell lin                                     #
################################################################################

#--------------------------Project Paths---------------------------------------#
#project folder
PROJECTROOT=/project/ChromGroup_Seq_data/Celeste/2024_summer_complete
#rootfolder to the raw experimental data

PIPELINE=$PROJECTROOT/pipeline
#csv file with all of the experiments and a short id name
EXPSUMMARY=$PROJECTROOT/experiment_summary.csv

#-------------------------Module Paths-----------------------------------------#
SAMPATH=$PROJECTROOT/alignment/sam
BAMPATH=$PROJECTROOT/alignment/bam
BEDPATH=$PROJECTROOT/alignment/bed
BWPATH=$PROJECTROOT/alignment/bigwig
BEDGRAPHPATH=$PROJECTROOT/alignment/bedgraph
CHROMSIZES=$PROJECTROOT/pipeline/hg38.chrom.sizes
filename=merged.igg
DUPREMOVE=true
#---------------------------Make iGg merged------------------------------------#

#merge from bam
for rep in 1 2 ;
do
  samtools merge -f $PROJECTROOT/alignment/bam/${filename}_${rep}.mapped.bam \
  $PROJECTROOT/alignment/bam/HBD.igg_${rep}.mapped.bam \
  $PROJECTROOT/alignment/bam/WKK.igg_${rep}.mapped.bam \
  $PROJECTROOT/alignment/bam/WT.igg_${rep}.mapped.bam
  
  samtools sort -n $PROJECTROOT/alignment/bam/${filename}_${rep}.mapped.bam \
  -o $PROJECTROOT/alignment/bam/${filename}_${rep}.sorted.mapped.bam
  
  #file conversion and filtering
  
  ## Convert into bed file format
  echo "Convert to bed file formant"
  bedtools bamtobed -i $BAMPATH/${filename}_${rep}.sorted.mapped.bam -bedpe > \
  $BEDPATH/${filename}_${rep}_bowtie2.bed
  
  ## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
  echo "Clean the bed file"
  awk '$1==$4 && $6-$2 < 1000 {print $0}' $BEDPATH/${filename}_${rep}_bowtie2.bed > \
  $BEDPATH/${filename}_${rep}_bowtie2.clean.bed
  
  ## Only extract the fragment related columns
  echo "Extract fragment related columns"
  cut -f 1,2,6 $BEDPATH/${filename}_${rep}_bowtie2.clean.bed | \
  sort -k1,1 -k2,2n -k3,3n  > $BEDPATH/${filename}_${rep}_bowtie2.fragments.bed
  
  # Convert to bedgraph
  bedtools genomecov -bg -i $BEDPATH/${filename}_${rep}_bowtie2.fragments.bed \
  				-g $CHROMSIZES > $BEDGRAPHPATH/${filename}_${rep}_bowtie2.fragments.bedgraph
  				
  #rmDup
  
  samtools merge -f $PROJECTROOT/alignment/bam/${filename}_${rep}.mapped.rmDup.bam \
  $PROJECTROOT/alignment/bam/HBD.igg_${rep}.mapped.rmDup.bam \
  $PROJECTROOT/alignment/bam/WKK.igg_${rep}.mapped.rmDup.bam \
  $PROJECTROOT/alignment/bam/WT.igg_${rep}.mapped.rmDup.bam
  
  samtools sort -n $PROJECTROOT/alignment/bam/${filename}_${rep}.mapped.rmDup.bam \
  -o $PROJECTROOT/alignment/bam/${filename}_${rep}.sorted.mapped.rmDup.bam
  
  #file conversion and filtering
  
  ## Convert into bed file format
  echo "Convert to bed file formant"
  bedtools bamtobed -i $BAMPATH/${filename}_${rep}.sorted.mapped.rmDup.bam -bedpe > \
  $BEDPATH/${filename}_${rep}_bowtie2.rmDup.bed
  
  ## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
  echo "Clean the bed file"
  awk '$1==$4 && $6-$2 < 1000 {print $0}' $BEDPATH/${filename}_${rep}_bowtie2.rmDup.bed > \
  $BEDPATH/${filename}_${rep}_bowtie2.rmDup.clean.bed
  
  ## Only extract the fragment related columns
  echo "Extract fragment related columns"
  cut -f 1,2,6 $BEDPATH/${filename}_${rep}_bowtie2.rmDup.clean.bed | \
  sort -k1,1 -k2,2n -k3,3n  > $BEDPATH/${filename}_${rep}_bowtie2.rmDup.fragments.bed
  
  # Convert to bedgraph
  bedtools genomecov -bg -i $BEDPATH/${filename}_${rep}_bowtie2.rmDup.fragments.bed \
  				-g $CHROMSIZES > $BEDGRAPHPATH/${filename}_${rep}_bowtie2.rmDup.fragments.bedgraph
  
  echo "Calculate seq depth of spike in genome"

##calculate seq depth of spike in genome

seqDepth=1207123

echo "Rescale"
if [[ "$seqDepth" -gt "1" ]]; then

  mkdir -p $PROJECTROOT/alignment/bedgraph
  
  scale_factor=`echo "10000 / $seqDepth" | bc -l`
  echo "Scaling factor for $name is: $scale_factor!"
  bedtools genomecov -bg -scale $scale_factor -i \
  $PROJECTROOT/alignment/bed/${filename}_${rep}_bowtie2.fragments.bed -g $CHROMSIZES > \
  $PROJECTROOT/alignment/bedgraph/${filename}_${rep}_bowtie2.fragments.normalized.bedgraph

fi

if [ "$DUPREMOVE" = true ]; then
  echo "Rescale Dedup"
  if [[ "$seqDepth" -gt "1" ]]; then
  
    mkdir -p $PROJECTROOT/alignment/bedgraph
    
    scale_factor=`echo "10000 / $seqDepth" | bc -l`
    echo "Scaling factor for $name is: $scale_factor!"
    bedtools genomecov -bg -scale $scale_factor -i \
    $PROJECTROOT/alignment/bed/${filename}_${rep}_bowtie2.rmDup.fragments.bed -g $CHROMSIZES > \
    $PROJECTROOT/alignment/bedgraph/${filename}_${rep}_bowtie2.rmDup.fragments.normalized.bedgraph
  
  fi
fi

  bedGraphToBigWig $BEDGRAPHPATH/${filename}_${rep}_bowtie2.fragments.bedgraph \
  $CHROMSIZES $BWPATH/${filename}_${rep}_bowtie2.fragments.bw

done