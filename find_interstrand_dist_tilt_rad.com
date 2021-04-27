#!/bin/csh -f
set pdbindir =  "../../../fixed_pdbs/"
set dataoutdir = "data"
set workdir = `pwd`
rm -rf $dataoutdir
mkdir $dataoutdir
rm find_interstrand_dist_tilt_rad.log

cd $pdbindir
ls -1 *.pdb > $workdir/pdblist
cd $workdir
set list = `cat pdblist`
foreach pdbfile ($list )
echo 'preparing pdb file' $pdbfile

cp $pdbindir/$pdbfile .

python find_interstrand_dist_tilt_rad.py <<eof | tee -a find_interstrand_dist_tilt_rad.log
$pdbfile
eof

mv *.csv $dataoutdir
rm $pdbfile
end
echo "Bye" `whoami`
exit
