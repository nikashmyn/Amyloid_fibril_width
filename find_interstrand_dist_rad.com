#!/bin/csh -f
set pdbindir =  "pdb_noh/"
set dataoutdir = "pdb_noh/"
set workdir = `pwd`
mkdir $dataoutdir

cd $pdbindir
ls -1 *.pdb > $workdir/pdblist
cd $workdir
set list = `cat pdblist`
foreach pdbfile ($list )
echo 'preparing pdb file' $pdbfile

cp $pdbindir/$pdbfile .

python find_interstrand_dist_tilt_rad.py <<eof | tee -a $dataoutdir/find_interstrand_dist_tilt_rad.log
$pdbfile
eof

mv *.csv $dataoutdir
rm $pdbfile
end
echo "Bye" `whoami`
exit
