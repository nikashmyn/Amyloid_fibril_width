#!/bin/sh

no_h_pdb_dir=/Users/davidboyer/Dropbox/Tau_Project_EISENBERG_LAB/amyloid_width/python/AB-Fibril-Radius-vs-Hydrogen-Bond-Relationship/pdb_noh/
h_pdb_dir=/Users/davidboyer/Dropbox/Tau_Project_EISENBERG_LAB/amyloid_width/python/AB-Fibril-Radius-vs-Hydrogen-Bond-Relationship/pdb_h/
cpus=4
scripts_dir=/Users/davidboyer/Dropbox/Tau_Project_EISENBERG_LAB/amyloid_width/python/AB-Fibril-Radius-vs-Hydrogen-Bond-Relationship/ #path to the folder which includes all auxiliary scripts
output_graphs_dir=/Users/davidboyer/Dropbox/Tau_Project_EISENBERG_LAB/amyloid_width/python/AB-Fibril-Radius-vs-Hydrogen-Bond-Relationship/graphs/
helical_parameters=/Users/davidboyer/Dropbox/Tau_Project_EISENBERG_LAB/amyloid_width/width_pitch/helical_parameters.txt

echo "Calculating interstrand distances for the following files:"
echo `ls ${no_h_pdb_dir}/*.pdb`


#read in .pdb without hydrogens to process
cd ${no_h_pdb_dir}
ls *.pdb > pdb_files.txt

#Reset the cmd list file
echo "" > cmd_list.cmds

#Generate run script commands (cmds) for interstrand distance calculations

while read line
do
	echo "python3 ${scripts_dir}/find_interstrand_dist_rad.py ${line} ${no_h_pdb_dir}" >> ${no_h_pdb_dir}/cmd_list.cmds;
	echo "" >> ${no_h_pdb_dir}/cmd_list.cmds
done < ${no_h_pdb_dir}/pdb_files.txt

parallel --jobs ${cpus} < ${no_h_pdb_dir}/cmd_list.cmds &> ${no_h_pdb_dir}/run.macro_interstrand.stdouterr

echo "Done!"

echo "Calculating beta-sheet hydrogen bond lengths and off-axis tilts for the following files:"
echo `ls ${h_pdb_dir}/*.pdb` 

#read in .pdb files to process
cd ${h_pdb_dir}
ls *.pdb > pdb_files.txt

#Reset the cmd list file
echo "" > ${h_pdb_dir}/cmd_list.cmds

#Generate Run script commands (cmds) for h-bonds csvs
while read line  
do
  echo "python3 ${scripts_dir}/find_hbonds_04252021_debugged.py ${line} ${h_pdb_dir}" >> ${h_pdb_dir}/cmd_list.cmds;
  echo "" >> ${h_pdb_dir}/cmd_list.cmds;
done < ${h_pdb_dir}/pdb_files.txt

parallel --jobs ${cpus} < ${h_pdb_dir}/cmd_list.cmds &> ${h_pdb_dir}/run.macro_hbond.stdouterr

echo "Done!"
echo "I am now graphing your data :)"

#Reset the cmd list file
echo "" > ${output_graphs_dir}/cmd_list_plots.cmds

#read in successfully generated csvs 
ls ${h_pdb_dir}/*datatable.csv > ${output_graphs_dir}/datatables_files.txt
ls ${h_pdb_dir}/*hbonds.csv > ${output_graphs_dir}/hbonds_files.txt
ls ${no_h_pdb_dir}/*dist.csv > ${output_graphs_dir}/dist_files.txt

#Generate run script commands (cmds) for h-bond plots
while read line1 <&1 && read line2 <&2 && read line3 <&3; 
do
  echo "python3 ${scripts_dir}/Generate_H-Bond_Plots.py ${line1} ${line2} ${line3} ${helical_parameters} ${output_graphs_dir}" >> ${output_graphs_dir}/cmd_list_plots.cmds;
  echo "" >> ${output_graphs_dir}/cmd_list_plots.cmds;
done 1<${output_graphs_dir}/datatables_files.txt 2<${output_graphs_dir}/hbonds_files.txt 3<${output_graphs_dir}/dist_files.txt

parallel --jobs ${cpus} < ${output_graphs_dir}/cmd_list_plots.cmds &> ${output_graphs_dir}/run.macro_hbond_plots.stdouterr

echo "Done!"
