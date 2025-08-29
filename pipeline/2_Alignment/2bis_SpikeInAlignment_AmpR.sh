################################################################################
#                                                                              #
# cut & tag analysis pipeline R-loop identification project                    #
# 2bis. Spike-In Alignment                                                     #
#                                                                              #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: Alignment to Spike-In genome for Spike-In calibration (only needed) #
# if the experiment was performed with a Spike-In)                             #
################################################################################



#-----------------------------Paths--------------------------------------------#

#root path of th whole project
projPath=$1
#fastqfiles for the replicates
rep1=$2
rep2=$3
#short name to id experiment
name=$4

#reference genome of the spike-in
spikeInRef=$5
#chromsize files for the reference genome of experiment (human)
chromSize=$6
DUPREMOVE=$7

#------------------------Spike-in calibration script---------------------------#
cores=8

echo "Align to ecoli"

echo $projPath/alignment/sam/${name}_bowtie2_spikeIn_AmpR.sam
## bowtie2-build path/to/Ecoli/fasta/Ecoli.fa /path/to/bowtie2Index/Ecoli
bowtie2  --local --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${spikeInRef} -1 $rep1 -2 $rep2 -S $projPath/alignment/sam/${name}_bowtie2_spikeIn_AmpR.sam &>  $projPath/alignment/sam/bowtie2_summary/${name}_bowtie2_spikeIn_AmpR.txt

echo "Calculate seq depth of spike in genome"
##calculate seq depth of spike in genome

samtools view -bS -F 0x04 $projPath/alignment/sam/${name}_bowtie2_spikeIn_AmpR.sam > $projPath/alignment/bam/AmpR_calibration/${name}_bowtie2_spikeIn_AmpR.bam
samtools sort $projPath/alignment/bam/AmpR_calibration/${name}_bowtie2_spikeIn_AmpR.bam -o $projPath/alignment/bam/AmpR_calibration/${name}_bowtie2_spikeIn_AmpR.sorted.bam
samtools index $projPath/alignment/bam/AmpR_calibration/${name}_bowtie2_spikeIn_AmpR.sorted.bam -o $projPath/alignment/bam/AmpR_calibration/${name}_bowtie2_spikeIn_AmpR.sorted.bam.bai
seqDepthDouble=`samtools view -c $projPath/alignment/bam/AmpR_calibration/${name}_bowtie2_spikeIn_AmpR.sorted.bam AMPR`
seqDepth=$((seqDepthDouble/2))
echo $seqDepth > $projPath/alignment/sam/bowtie2_summary/${name}_bowtie2_spikeIn_AmpR.seqDepth
total=`samtools view -c $projPath/alignment/bam/${name}.mapped.sorted.filtered.bam`
total=$(($total/2))
echo "Rescale"
if [[ "$seqDepth" -gt "1" ]]; then

  mkdir -p $projPath/alignment/bedgraph
  
  scale_factor=`echo "1000*$seqDepth/$total" | bc -l`
  echo "Scaling factor for $name is: $scale_factor!"
  echo $scale_factor > $projPath/alignment/sam/bowtie2_summary/${name}_bowtie2_spikeIn_AmpR.scaleFactor
  bedtools genomecov -bg -scale $scale_factor -i \
  $projPath/alignment/bed/${name}_bowtie2.fragments.filtered.bed -g $chromSize > \
  $projPath/alignment/bedgraph/${name}_bowtie2.fragments.filtered.normalized.AmpR.bedgraph

fi


