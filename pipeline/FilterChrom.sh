#-------------------------------------------------------------------------------
# Author:  Sara Lopez Ruiz de Vargas
# Email:  lopez_s@molgen.mpg.de
#
# Date:    2025-03-17
#
# Script Name: cut_and_tag_pipeline_tutorial.sh
#
# Original Pipeline: https://yezhengstat.github.io/CUTTag_tutorial/#I_Introduction
#
# Script Description: This is a script to filter out all mitochondrial and 
# non-standard chromosomes from the alignment and peaCalling files
#-------------------------------------------------------------------------------



#--------------------------Project Paths---------------------------------------#

#Reminder: In shell scripting variables are assigned using = without spaces
#Example: PROJECTROOT=/this/is/my/path will work.
#         PROJECTROOT = /this/is/my/path doesn't work.
# To access a variable we use '$',
# Example: echo $PROJECTROOT will print whatever is in that variable.
#          echo PROJECTROOT will throw an error or not work as expected.

#project folder
PROJECTROOT=/project/ChromGroup_Seq_data/Celeste/2024_summer_complete

#rootfolder to the raw experimental data (folder from SeqCore)

#Reminder: Raw data files are huge, there is no need to copy them into
#your own directory, this leads to IT related problems.
RAWROOT=/project/avitidata/pipelined/20240917_AV234501_A-PE75-Ad-FS-HO/Samples/V

#path to the folder with your pipeline
PIPELINE=$PROJECTROOT/pipeline

#csv file with all of the experiments but replicates are given in columns
#more suitable for the alignment step
EXPSUMMARYAL=$PROJECTROOT/experiment_summary_align_formatted.csv

#csv file for peaks
EXPSUMMARYPEAKS=$PROJECTROOT/experiment_summary_peaks.csv

EXPMERGE=$PROJECTROOT/experiment_summary_merge.csv

#-------------------------------Filtering--------------------------------------#

bam=$PROJECTROOT/alignment/bam
bed=$PROJECTROOT/alignment/bed
bedgraph=$PROJECTROOT/alignment/bedgraph
bigwig=$PROJECTROOT/alignment/bigwig
peakfiles=$PROJECTROOT/peakCalling/SEACR
chromSize=$PROJECTROOT/pipeline/hg38.chrom.sizes

##filtering bam files
# while read line ; do
#    set $line
#    IFS=$','; split=($line); unset IFS;
#    bamfile=$bam/${split[0]}.mapped.bam
#    bamfilesorted=$bam/${split[0]}.mapped.sorted.bam
#    bamfilefiltered=$bam/${split[0]}.mapped.sorted.filtered.bam
#    samtools sort $bamfile -o $bamfilesorted
#    samtools index $bamfilesorted
#    samtools idxstats $bamfilesorted | cut -f 1 | grep -Ev 'random|chrM|chrUn' | xargs samtools view -b $bamfilesorted > $bamfilefiltered
#    samtools index $bamfilefiltered
# done < <(tail -n +2 $EXPSUMMARYAL)

##filter the bed files

while read line ; do
   set $line
   IFS=$','; split=($line); unset IFS;
   input_bed=$bed/${split[0]}_bowtie2.fragments.bed
   out_bed=$bed/${split[0]}_bowtie2.fragments.filtered.bed
   awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/' $input_bed > $out_bed
done < <(tail -n +2 $EXPSUMMARYAL)

#filter bedgraph and convert to bigwig 
while read line ; do
   set $line
   IFS=$','; split=($line); unset IFS;
   input_bedgraph=$bedgraph/${split[0]}_bowtie2.fragments.normalized.bedgraph
   out_bedgraph=$bedgraph/${split[0]}_bowtie2.fragments.normalized.filtered.bedgraph
   awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/' $input_bedgraph > $out_bedgraph
   bedGraphToBigWig  $bedgraph/${split[0]}_bowtie2.fragments.normalized.filtered.bedgraph \
   $chromSize $bigwig/${split[0]}_bowtie2.fragments.normalized.filtered.bw
done < <(tail -n +2 $EXPSUMMARYAL)

##filter peakCalling
while read line ; do
   set $line
   IFS=$','; split=($line); unset IFS;
   input_peak=$peakfiles/${split[0]}_seacr_norm_control.peaks.relaxed.bed
   out_peak=$peakfiles/${split[0]}_seacr_norm_control.peaks.relaxed.filtered.bed
   awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/' $input_peak > $out_peak
done < <(tail -n +2 $EXPSUMMARYAL)

while read line ; do
   set $line
   IFS=$','; split=($line); unset IFS;
   input_peak=$peakfiles/${split[0]}_merged_seacr_norm_control.peaks.relaxed.bed
   out_peak=$peakfiles/${split[0]}_merged_seacr_norm_control.peaks.relaxed.filtered.bed
   awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/' $input_peak > $out_peak
done < <(tail -n +2 $EXPMERGE)


