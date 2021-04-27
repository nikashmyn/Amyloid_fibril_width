#!/bin/csh -f
set pdbindir =  "../../../fixed_pdbs_h/"
set dataoutdir = "data"
set workdir = `pwd`
rm -rf $dataoutdir
mkdir $dataoutdir
rm find_hbonds.log

cd $pdbindir
ls -1 *.pdb > $workdir/pdblist
cd $workdir
set list = `cat pdblist`
foreach pdbfile ($list )
echo 'preparing pdb file' $pdbfile

cp $pdbindir/$pdbfile .

python find_hbonds_04252021.py <<eof | tee -a find_hbonds.log
$pdbfile
eof

mv *.csv $dataoutdir
rm $pdbfile
end
echo "Bye" `whoami`
exit
