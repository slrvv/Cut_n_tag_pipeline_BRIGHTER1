################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 4. Reproducibility assessment                                                #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: create deeptools heatmap for reproducibility assessment             #
################################################################################

#-----------------------------Paths--------------------------------------------#

projPath=/project/ChromGroup_Seq_data/Celeste/2024_summer_complete
experimentSummary=$projPath/experiment_summary_align_formatted.csv
namesummary=all_experiments_bam_summary
nameplot=all_experiments_corr_bam

#-------------------------Script-----------------------------------------------#

#populate bam files array
bamfiles=()
samples=()
while IFS=, read -r SampleName R1 R2
do
   bamfile=$projPath/alignment/bam/$SampleName.mapped.sorted.bam
   bamfilefiltered=$projPath/alignment/bam/$SampleName.mapped.sorted.filtered.bam
   bamfiles+=" $bamfilefiltered"
   samples+=" $SampleName"
done < <(tail -n +2 $experimentSummary )

echo $samples
multiBamSummary bins --bamfiles $bamfiles -l $samples --binSize 2000 \
  -o $projPath/alignment/bam/$namesummary.filtered.npz \
  --blackListFileName $projPath/GRCh38_unified_blacklist.bed \
  --outRawCounts $projPath/alignment/bam/$namesummary.filtered.rawCounts.tsv

  
plotCorrelation -in $projPath/alignment/bam/${namesummary}.filtered.npz -c pearson  \
--colorMap bwr --plotNumbers  --removeOutliers --plotHeight 18 --plotWidth  20 \
--plotFileFormat "pdf" -p heatmap -o $projPath/summary_figures_AmpR_calibrated/${nameplot}_2000_pearson_filtered.pdf \
--outFileCorMatrix $projPath/summary_figures_AmpR_calibrated/${nameplot}_2000_pearson_filtered.cormatrix.txt\

