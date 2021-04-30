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

ls ${datadir}/*.pdb > ${outdir}/pdb_files.txt

#Reset the cmd list file
echo "" > ${outdir}/cmd_list.cmds

#Run script to generate h-bonds csvs
while read line  
do
  echo "python3 ${scriptsdir}/find_hbonds_04252021_debugged.ipynb ${line} ${outdir}" >> ${outdir}/cmd_list.cmds;
  echo "" >> ${outdir}/cmd_list.cmds;
done < ${outdir}/pdb_files.txt

parallel --jobs ${cpus} < ${outdir}/cmd_list.cmds &> ${outdir}/run.macro_hbond.stdouterr




#for i in {1..22}
#do
#        echo "$HC_cmd -L chr${i} --out ${output_prefix}_${output_type}_HC.${i}.vcf" >> ${HC_output}
#        echo "$UG_cmd -L chr${i} --out ${output_prefix}_${output_type}_UG.${i}.vcf">> ${UG_output}
#        echo "$Discovery_cmd -L chr${i} --out ${output_prefix}_HC_Discovery.${i}.vcf" >> ${Discovery_output}
#done

#parallel --jobs ${cpus} < star.cmds &> star.cmds.stdouterr

#parallel --jobs ${cpus} bash ${scriptsdir}/CollectMultipleMetricsForMultipleBams.sh {}  ::: ${dstDir}/star/*.Aligned.sortedByCoord.out.bam &

