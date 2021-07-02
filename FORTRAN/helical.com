#!/bin/csh -f

#Program written and freely distributed by Michael Sawaya

# Generate twisted filament by applying nops copies

# input file with asymmetric unit
set pdbin = RealSpaceRefine_1/6lni_origin_5layers_real_space_refined.pdb

# output file with fibril coordinates
set pdbout = 6lni_origin_5layers_refined.pdb

# total number of increments of the symmetry operator to apply to the asymmetric unit
set nops = 9  

# number of degrees to rotate per step
set rotdegincrement = 179.45       

# number of angstroms to translate along z in angstroms
set translatez = 2.4       

set opcount = 0   

rm -rf temp*.pdb
rm -rf in_rechain*.pdb

# grep only bottom layer chain A
 grep ' A ' $pdbin > unit.pdb

loopofgreatrewards:
 set outfiledigit = `echo $opcount |awk '{printf("%.4d \n" , $1);}'`
 set pdbfileoutm = tempm-$outfiledigit.pdb


set m1 = `echo "c( $opcount * $rotdegincrement * 3.1418 / 180)"  | bc -l`
set m2 = `echo "s(-1 * $opcount * $rotdegincrement * 3.1418 / 180)"  | bc -l`
set m3 = 0
set m4 = `echo "s( $opcount * $rotdegincrement * 3.1418 / 180)"  | bc -l`
set m5 = `echo "c( $opcount * $rotdegincrement * 3.1418 / 180)"  | bc -l`
set m6 = 0
set m7 = 0
set m8 = 0
set m9 = 1
set transx = 0
set transy = 0
set transz = `echo "$opcount*$translatez" | bc -l`

pdbset xyzin unit.pdb xyzout $pdbfileoutm  <<eof-1
symgen NCS
transform $m1 $m2 $m3    $m4 $m5 $m6    $m7 $m8 $m9   $transx  $transy $transz
eof-1


echo $opcount
@ opcount = ($opcount + 1)

if ( $opcount ==  $nops + 1 ) goto exitstageleft
goto loopofgreatrewards

exitstageleft:

grep CRYST1 $pdbin > in.pdb
grep -h ATOM tempm--0010.pdb >> in.pdb
grep -h ATOM tempm--0009.pdb >> in.pdb
grep -h ATOM tempm--0008.pdb >> in.pdb
grep -h ATOM tempm--0007.pdb >> in.pdb
grep -h ATOM tempm--0006.pdb >> in.pdb
grep -h ATOM tempm--0005.pdb >> in.pdb
grep -h ATOM tempm--0004.pdb >> in.pdb
grep -h ATOM tempm--0003.pdb >> in.pdb
grep -h ATOM tempm--0002.pdb >> in.pdb
grep -h ATOM tempm--0001.pdb >> in.pdb
grep -h ATOM tempm-0000.pdb >> in.pdb
grep -h ATOM tempm-0001.pdb >> in.pdb
grep -h ATOM tempm-0002.pdb >> in.pdb
grep -h ATOM tempm-0003.pdb >> in.pdb
grep -h ATOM tempm-0004.pdb >> in.pdb
grep -h ATOM tempm-0005.pdb >> in.pdb
grep -h ATOM tempm-0006.pdb >> in.pdb
grep -h ATOM tempm-0007.pdb >> in.pdb
grep -h ATOM tempm-0008.pdb >> in.pdb
grep -h ATOM tempm-0009.pdb >> in.pdb
grep -h ATOM tempm-0010.pdb >> in.pdb
grep -h ATOM tempm-0011.pdb >> in.pdb
grep -h ATOM tempm-0012.pdb >> in.pdb
grep -h ATOM tempm-0013.pdb >> in.pdb
grep -h ATOM tempm-0014.pdb >> in.pdb
grep -h ATOM tempm-0015.pdb >> in.pdb
grep -h ATOM tempm-0016.pdb >> in.pdb


#automatically increment chain IDs
beyond26chainsdiscont <<eof
in.pdb
eof

cat in_rechain.pdb > $pdbout
rm in_rechain.pdb 

rm temp*.pdb
rm in.pdb 
echo 'Auf wiedersehen' `whoami`
exit
#
