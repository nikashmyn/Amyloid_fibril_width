#!/bin/sh

datadir=$1
cpus=$2
scriptsdir=$3 #path to the folder which includes all auxiliary scripts
outdir=$4

echo ${datadir}
echo ${cpus}

echo "Calculating h-bonds and generating visuals for the following files:"
echo $PWD
cd ${datadir}
#ls ${datadir}/*.pdb

#read in .pdb files to process
ls ${datadir}/*.pdb > ${outdir}/pdb_files.txt

#Reset the cmd list file
echo "" > ${outdir}/cmd_list.cmds

#Run script to generate h-bonds csvs
while read line  
do
  echo "python3 ${scriptsdir}/find_hbonds_04252021_debugged.py ${line} ${outdir}" >> ${outdir}/cmd_list.cmds;
  echo "" >> ${outdir}/cmd_list.cmds;
done < ${outdir}/pdb_files.txt

parallel --jobs ${cpus} < ${outdir}/cmd_list.cmds &> ${outdir}/run.macro_hbond.stdouterr


#read in successfully generated csvs 
#ls ${outdir}/*datatable.csv > ${outdir}/datatables_files.txt
#ls ${outdir}/*hbonds.csv > ${outdir}/hbonds_files.txt

#Run script to generate h-bonds plots
#while read line
#do
#  echo "python3 ${scriptsdir}/Generate_H-Bond_Plots.ipynb ${line} ${outdir}" >> ${outdir}/cmd_list.cmds;
#  echo "" >> ${outdir}/cmd_list.cmds;
#done < ${outdir}/pdb_files.txt

#parallel --jobs ${cpus} < ${outdir}/cmd_list.cmds &> ${outdir}/run.macro_hbond.stdouterr


#parallel --jobs ${cpus} bash ${scriptsdir}/CollectMultipleMetricsForMultipleBams.sh {}  ::: ${dstDir}/star/*.Aligned.sortedByCoord.out.bam &

