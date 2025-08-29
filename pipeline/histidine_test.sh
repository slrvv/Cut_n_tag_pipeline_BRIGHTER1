#-------------------------------------------------------------------------------
# Author:  Sara Lopez Ruiz de Vargas
# Email:  lopez_s@molgen.mpg.de
#
# Date:    2024-04-12
#
# Script Name: Histidine_test
#
# Original Pipeline: https://yezhengstat.github.io/CUTTag_tutorial/#I_Introduction
#
# Script Description: Tests using the histidine experiment to see which peak 
# caller performs best and test a few parameters
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

#path to the folder with your pipeline
PIPELINE=$PROJECTROOT/pipeline

#csv file with all of the experiments and a short id name
EXPSUMMARY=$PROJECTROOT/experiment_summary.csv

#csv file with all of the experiments but replicates are given in columns
#more suitable for the alignment step
EXPSUMMARYAL=$PROJECTROOT/experiment_summary_align_formatted.csv

#csv file for peaks
EXPSUMMARYPEAKS=$PROJECTROOT/experiment_summary_peaks


EXPMERGE=$PROJECTROOT/experiment_summary_merge.csv
DUPREMOVE=true #set to true if you want to have the duplicate removal files converted
IGGNEED=true #set to true if you want to remove IgG from bw files
MACS2=true
NORMPARAM=true #set to false if you only want to run SEACR on non stringent mode
RELAXED=true
#-------------------------Module Paths-----------------------------------------#
# No need to change anything from now on
#bash file that does the first step QC
#Peak calling modules
PEAKS=$PIPELINE/6_PeakCalling/6_PeakCalling_His_test.sh
PEAKSUMM=$PIPELINE/6_PeakCalling/6_PeakCallingSummary.R
PEAKSFRIP=$PIPELINE/6_PeakCalling/6_PeakCallingFrips.R
PEAKSMACS2=$PIPELINE/6_PeakCalling/6_PeakCalling_MACS2_His_test.sh
PEAKSFRIPMACS=$PIPELINE/6_PeakCalling/6_PeakCallingFrips_MACS2.R
PEAKSUMMACS=$PIPELINE/6_PeakCalling/6_PeakCallingSummary_MACS2.R


#---------------------Start Script---------------------------------------------#

#To make a decision on which peak caller to use and which parameters we focus on
# the Histidine experiment in the HBD cell line

# Things to test 
# 1. Should we merge IgG, use one for all or use the one corresponding to each cell
# line

# 2. SEACR parameters non vs. norm 

# 3. Possible q-values of MACS2

mkdir -p $PROJECTROOT/peakCalling/His_test

#for summary in merged_cell; do
# for summary in same_cell diff_cell merged_cell; do
#   echo ${EXPSUMMARYPEAKS}_${summary}.csv
#   while read line ; do
#       set $line
#       IFS=$','; split=($line); unset IFS;
#       echo "${split[0]}"
#       bash $PEAKS $PROJECTROOT ${split[1]} ${split[0]} $DUPREMOVE $NORMPARAM $summary $RELAXED
#   done < <(tail -n +2 ${EXPSUMMARYPEAKS}_${summary}.csv)
#   #
#   # 6_PeakCallingSummaryPlot.R for visualization
#   # if [ "$MACS2" = true ]; then
#   #    while read line ; do
#   #        set $line
#   #        IFS=$','; split=($line); unset IFS;
#   #        echo "${split[0]}"
#   #        bash $PEAKSMACS2 $PROJECTROOT ${split[1]} ${split[0]} $DUPREMOVE "0.01 0.05 0.1" $summary
#   #    done < <(tail -n +2 ${EXPSUMMARYPEAKS}_${summary}.csv)
#   # fi
# 
# done

## Test the same but with WKK vs HBD
# for summary in same_cell_WKK diff_cell_WKK; do
#   while read line ; do
#       set $line
#       IFS=$','; split=($line); unset IFS;
#       echo "${split[0]}"
#       bash $PEAKS $PROJECTROOT ${split[1]} ${split[0]} $DUPREMOVE $NORMPARAM $summary $RELAXED
#   done < <(tail -n +2 ${EXPSUMMARYPEAKS}_${summary}.csv)
# done

# EXPSUMMARYPEAKS=$PROJECTROOT/experiment_summary_peaks_WKK_ant.csv
# 
# while read line ; do
#     set $line
#     IFS=$','; split=($line); unset IFS;
#     echo "${split[0]}"
#     bash $PEAKS $PROJECTROOT ${split[1]} ${split[0]} $DUPREMOVE $NORMPARAM $summary $RELAXED
# done < <(tail -n +2 $EXPSUMMARYPEAKS )

EXPSUMMARYPEAKS=$PROJECTROOT/rpa70_summary_peaks
for summary in same_cell diff_cell; do
  while read line ; do
      set $line
      IFS=$','; split=($line); unset IFS;
      echo "${split[0]}"
      bash $PEAKS $PROJECTROOT ${split[1]} ${split[0]} $DUPREMOVE $NORMPARAM $summary $RELAXED
  done < <(tail -n +2 ${EXPSUMMARYPEAKS}_${summary}.csv)
done