#!/bin/sh

datadir=$1
cpus=$2
scriptsdir=$3 #path to the folder which includes all auxiliary scripts
outdir=$4
distdir=$5
helixdir=$6

echo ${datadir}
echo ${cpus}

echo "Calculating h-bonds and generating visuals for the following files:"
echo $PWD

cd ${datadir}
ls ${datadir}/*.pdb

#clear old results from outdir
#i dont think we need to remove old results. It should save over. uncomment if im wrong. 
#rm ${outdir}/*

#read in .pdb files to process
ls ${datadir}/*.pdb > ${outdir}/pdb_files.txt

#Reset the cmd list file
echo "" > ${outdir}/cmd_list.cmds

#Generate Run script commands (cmds) for h-bonds csvs
while read line  
do
  echo "python3 ${scriptsdir}/find_hbonds_04252021_debugged.py ${line} ${outdir}" >> ${outdir}/cmd_list.cmds;
  echo "" >> ${outdir}/cmd_list.cmds;
done < ${outdir}/pdb_files.txt

parallel --jobs ${cpus} < ${outdir}/cmd_list.cmds &> ${outdir}/run.macro_hbond.stdouterr

#Reset the cmd list file
echo "" > ${outdir}/cmd_list_plots.cmds

#read in successfully generated csvs 
ls ${outdir}/*datatable.csv > ${outdir}/datatables_files.txt
ls ${outdir}/*hbonds.csv > ${outdir}/hbonds_files.txt
ls ${distdir}/*dist.csv > ${outdir}/dist_files.txt

#Generate run script commands (cmds) for h-bond plots
while read line1 <&1 && read line2 <&2 && read line3 <&3; 
do
  echo "python3 ${scriptsdir}/Generate_H-Bond_Plots.py ${line1} ${line2} ${line3} ${helixdir} ${outdir}" >> ${outdir}/cmd_list_plots.cmds;
  echo "" >> ${outdir}/cmd_list_plots.cmds;
done 1<${outdir}/datatables_files.txt 2<${outdir}/hbonds_files.txt 3<${outdir}/dist_files.txt

parallel --jobs ${cpus} < ${outdir}/cmd_list_plots.cmds &> ${outdir}/run.macro_hbond_plots.stdouterr


