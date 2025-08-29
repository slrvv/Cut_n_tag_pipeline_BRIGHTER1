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
nameplot=all_experiments_corr_heatmap
nameplotpca=all_experiments_pca
namesummarybw=all_experiments_bam_summarybw_spearman_2000
nameplotbw=all_experiments_corr_heatmapbw_spearman_2000
nameplotpcabw=all_experiments_pcabw_spearman_2000
#-------------------------Script-----------------------------------------------#

#populate bam files array
bamfiles=()
samples=()
while IFS=, read -r SampleName R1 R2
do
   bamfile=$projPath/alignment/bam/$SampleName.mapped.bam.sorted
   bamfilefiltered=$projPath/alignment/bam/$SampleName.mapped.bam.sorted.filtered
   samtools idxstats $bamfile | cut -f 1 | grep -Ev 'random|chrM|chrUn' | xargs samtools view -b $bamfile > $projPath/alignment/bam/$SampleName.mapped.bam.sorted.filtered
   samtools index $bamfilefiltered
   bamfiles+=" $bamfilefiltered"
   samples+=" $SampleName"
done < <(tail -n +2 $experimentSummary )

echo $samples
multiBamSummary bins --bamfiles $bamfiles -l $samples --binSize 2000 \
  -o $projPath/alignment/bam/$namesummary.filtered.npz --outRawCounts $projPath/alignment/bam/$namesummary.filtered.rawCounts.tsv
  
  
plotCorrelation -in $projPath/alignment/bam/${namesummary}.filtered.npz -c pearson  \
--colorMap RdYlBu --plotNumbers  --removeOutliers --plotHeight 18 --plotWidth  20 \
--plotFileFormat "pdf" -p heatmap -o $projPath/alignment/${nameplot}_2000_pearson_filtered_scatterplot.pdf \
--outFileCorMatrix $projPath/alignment/${nameplot}_2000_pearson_filtered_spearman.cormatrix.txt\

# bamfiles=()
# samples=()
# for SampleName in HBD.his_1 HBD.his_2 HBD.flag_1 HBD.flag_2 HBD.igg_1 HBD.igg_2
# do
#    bam=$projPath/alignment/bam/$SampleName.mapped.bam.sorted
#    bamfiles+=" $bam"
#    samples+=" $SampleName"
# done

# 
# echo $samples
# multiBamSummary bins --bamfiles $bamfiles -l $samples --binSize 2000 \
#   -o $projPath/alignment/bam/$namesummary.npz --outRawCounts $projPath/alignment/bam/$namesummary.rawCounts.tsv

#Figure out what these outliers are
# awk -F'\t' '{
#     for(i=4; i<=NF; i++) {
#         if($i > 5000) { print $0; next }
#     }
# }' $projPath/alignment/bam/$namesummary.rawCounts.tsv > $projPath/alignment/bam/$namesummary.rawCounts.outliers.tsv
#They seem to be mitochondrial DNA

##Filter out only the chromosomes be want from bam file

# bamfiles=()
# samples=()
# for SampleName in HBD.his_1 HBD.his_2 HBD.flag_1 HBD.flag_2 HBD.igg_1 HBD.igg_2
# do
#    bamfile=$projPath/alignment/bam/$SampleName.mapped.bam.sorted
#    bamfilefiltered=$projPath/alignment/bam/$SampleName.mapped.bam.sorted.filtered
#    samtools idxstats $bamfile | cut -f 1 | grep -Ev 'random|chrM|chrUn' | xargs samtools view -b $bamfile > $projPath/alignment/bam/$SampleName.mapped.bam.sorted.filtered
#    samtools index $bamfilefiltered
#    bamfiles+=" $bamfilefiltered"
#    samples+=" $SampleName"
#    
# 
# done
# 
# 
# echo $samples
# multiBamSummary bins --bamfiles $bamfiles -l $samples --binSize 2000 \
#   -o $projPath/alignment/bam/$namesummary.filtered.npz --outRawCounts $projPath/alignment/bam/$namesummary.filtered.rawCounts.tsv



# plotCorrelation -in $projPath/alignment/bam/${namesummary}.filtered.npz -c pearson  \
# --colorMap RdYlBu --plotNumbers  --removeOutliers --plotHeight 18 --plotWidth  20 \
# --plotFileFormat "pdf" -p scatterplot -o $projPath/alignment/${nameplot}_2000_pearson_filtered_scatterplot.pdf \
# --outFileCorMatrix $projPath/alignment/${nameplot}_2000_pearson_filtered_spearman.cormatrix.txt\

# plotPCA -in $projPath/alignment/bam/$namesummary.npz   \
#  --plotHeight 18 --plotWidth  20 \
#  --plotFileFormat "pdf"  -o $projPath/alignment/$nameplotpca.pdf


#populate bam files array
# bamfiles=()
# samples=()
# while IFS=, read -r SampleName R1 R2
# do
#    bam=$projPath/alignment/bigwig/${SampleName}_bowtie2.fragments.bw
#    bamfiles+=" $bam"
#    samples+=" $SampleName"
# done < <(tail -n +2 $experimentSummary )
# 
# echo $samples
# multiBigwigSummary bins -b $bamfiles -l $samples --binSize 5000 \
#   -o $projPath/alignment/bigwig/$namesummarybw.npz
# 
# plotCorrelation -in $projPath/alignment/bigwig/$namesummarybw.npz -c spearman  \
# --colorMap RdYlBu --plotNumbers --removeOutliers --plotHeight 18 --plotWidth  20 \
# --plotFileFormat "pdf" -p heatmap -o $projPath/alignment/$nameplotbw.pdf \
# --outFileCorMatrix $projPath/alignment/$nameplotbw.cormatrix.txt\
# 
# plotPCA -in $projPath/alignment/bigwig/$namesummarybw.npz   \
#  --plotHeight 18 --plotWidth  20 \
#  --plotFileFormat "pdf"  -o $projPath/alignment/$nameplotpcabw.pdf

Figure out what these outliers are
awk -F'\t' '{
    for(i=4; i<=NF; i++) {
        if($i > 5000) { print $0; next }
    }
}' all_experiments_bam_summary.filtered.rawCounts.tsv > all_experiments_bam_summary.filtered.rawCounts.outliers.tsv
