#!/bin/csh -f

#Script and distrecord2phenix_withsidechains_planarity3.f program written and freely distributed by Michael Sawaya

set pdbin = 6lni_origin_5layers.pdb
set hbplusout = `echo $pdbin | sed "s/pdb/hb2/"`
set pymolout = `echo $pdbin | sed "s/pdb/pml/"`
set cgoout = `echo $pdbin | sed "s/pdb/cgo/"`
set phenixout = `echo $pdbin | sed "s/pdb/edits/"`


rm -rf $hbplusout
rm -rf $pymolout 
rm -rf $phenixout
rm -rf $cgoout
/joule2/programs/hbplus/hbplus3.0_linux/hbplus -h 3.0 $pdbin

distrecord2phenix_withsidechains_planarity3 <<eof
$hbplusout
$phenixout
$pymolout
$cgoout
$pdbin
eof

echo ' good bye'
whoami

