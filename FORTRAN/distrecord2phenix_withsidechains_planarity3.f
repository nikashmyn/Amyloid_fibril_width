c Convert HBLUS records to external distance restraints for use in phenix
c sawaya 23 Apr 2018
c
c 2019-06-24 In version 3 I added conditions to restrict planarity restraints to inter-chain H-bonds between residue i and i+/-1. 
c That is normal hydrogen bonding pattern for parallel in-register sheets.
c It appears to have eliminated the weird H-bonds we saw with the previous version.
c
	character hbplusfile*90,header*90,resn1*3,resn2*3
        character editsfile*90, pmlfile*90,cgofile*90
        character pdbfile*90
        character ato1*3,ato2*3,cha1*1,cha2*1,DA*2
        integer resid1,resid2,noatom,chi,ch
        real x0(9000),y0(9000),z0(9000)
        character pdbline(9000)*80
        character atopdb(9000)*3,chpdb(9000)*1
        integer residpdb(9000),iresidpdb(9000)
        print *, 'type name of input hbplus file '
        read(*,'(a90)')hbplusfile
	open (1,file=hbplusfile,status='old',form='formatted')

        print *, 'type name of output distance restraints'
        read(*,'(a90)')editsfile
	open (2,file=editsfile,status='new',form='formatted')

        print *, 'type name of output hbonds for pymol display'
        read(*,'(a90)')pmlfile
	open (3,file=pmlfile,status='new',form='formatted')

        print *, 'type name of output cgo file for pymol display'
        read(*,'(a90)')cgofile
	open (4,file=cgofile,status='new',form='formatted')

        print *, 'type name of input pdb file'
        read(*,'(a90)')pdbfile
        open (5,file=pdbfile,status='old',form='formatted')
c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        noatom=0
c      remove header
c      remove header
c write the load pdb file record to pymol file
       write (3,10)pdbfile
 10    format('load ',a90,', protearella')
       write (3,13)cgofile
 13    format('@',a90)
c write the import cmd record to cgo file
       write (4,14)
 14    format('from pymol.cgo import *')
       write (4,11)
 11    format('from pymol import cmd')
       write (4,12)
 12    format('obj = [ BEGIN, TRIANGLES, COLOR, 0.5, 0.5, 0.5, \')
c
c      Read pdb coordinates
        n=1
304     read(5,303,end=306)pdbline(n)
303     format(a80)
        if(pdbline(n)(1:4).ne.'ATOM')goto 304
        if(pdbline(n)(18:20).eq.'HOH')goto 304
        chpdb(n)=pdbline(n)(22:22)
        atopdb(n)=pdbline(n)(14:16)
        read (pdbline(n)(23:26),*)residpdb(n)
        iresidpdb(n)=residpdb(n)
c       Convert character string to real (x)
        read (pdbline(n)(31:38),*)x0(n)
        read (pdbline(n)(39:46),*)y0(n)
        read (pdbline(n)(47:54),*)z0(n)
        n=n+1
        goto 304
306     n=n-1
        print*, n, 'Atoms read from ',pdbfile
        if(n.gt.9000)then
         print*,'Redimension for >9000 atoms'
         go to 999
        endif
c
c      Skip over hbplus header
       do 6, i=1,8
        read (1,7)header
  7     format(a90)
  6    continue
c
  4     read (1,5,end=99)cha1,resid1,resn1,ato1,cha2,resid2,resn2,ato2,
     1  DA
  5     format(a1,i4,x,a3,x,a3,x,a1,i4,x,a3,x,a3,6x,a2)
c       n   atom  resd res      DA  || num        DHA   H-A  angle D-A-AA Bond
c       s   type  num  typ     dist DA aas  dist angle  dist       angle   num
c  A0210-ALA N   A0240-TYR O   2.88 MM  30  4.90 103.9  2.47 123.5 143.0     1
c  A0240-TYR N   A0210-ALA O   3.15 MM  30  4.90 145.9  2.27 118.0 128.2     2
c  A0216-SER N   A0214-VAL O   2.83 MM   2  5.66 143.1  1.97 101.7  93.0     3
c  etc.
c NOTE to MIKE, you might want to comment out the next line for general purposes
c       if((resid1.lt.186).and.(resid2.lt.186))go to 4
c
c       Option to skip over sidechain atoms
        if((DA.ne.'MM').and.(DA.ne.'SS'))goto 4

        noatom=noatom+1

c now write the distance record to .editsfile
        write (2,101)
101     format('refinement.geometry_restraints.edits {')
        write (2,102)
102     format('  bond {')
        write (2,103)
103     format('    action = *add')
        write (2,104)ato1,cha1,resn1,resid1
104     format('    atom_selection_1 = name ',a4,' and chain ',a1,
     1  ' and resname ',a4,' and resseq ',i4)
        write (2,105)ato2,cha2,resn2,resid2
105     format('    atom_selection_2 = name ',a4,' and chain ',a1,
     1  ' and resname ',a4,' and resseq ',i4)
        write (2,106)
106     format('    distance_ideal = 2.80')
        write (2,107)
107     format('    sigma = 0.200')
        write (2,108)
108     format('  }')
        write (2,109)
109     format('}')

c refinement.geometry_restraints.edits {
c   bond {
c     action = *add
c     atom_selection_1 = name  S   and chain B and resname PMS and resseq  201
c     atom_selection_2 = name  OG  and chain A and resname SER and resseq  225
c     distance_ideal = 2.80
c     sigma = 0.200
c   }
c }

c now write the planarity restraint record to .edits file
       if((ato1.eq.'N  ').and.(ato2.eq.'O  '))then
       if(cha1.eq.cha2)goto 4
       if(abs(resid1-resid2).ne.1.0)goto 4
        write (2,201)
201     format('refinement.geometry_restraints.edits {')
        write (2,202)
202     format(' planarity {')
        write (2,203)
203     format('  action = *add')
        write (2,204)cha1,resn1,resid1,
     1               cha1,resn2,resid2,
     1               cha1,resn2,resid2,
     1               cha2,resn1,resid1,
     1               cha2,resn2,resid2,
     1               cha2,resn2,resid2 
204     format('   atom_selection = name N and chain ',a1,
     1  ' and resname ',a4,' and resseq ',i4,
     1  ' or name C and chain ',a1,' and resname ',a4,' and resseq ',i4,
     1  ' or name O and chain ',a1,' and resname ',a4,' and resseq ',i4,
c
     1  ' or name N and chain ',a1,' and resname ',a4,' and resseq ',i4,
     1  ' or name C and chain ',a1,' and resname ',a4,' and resseq ',i4,
     1  ' or name O and chain ',a1,' and resname ',a4,' and resseq ',i4)
        write (2,207)
207     format('    sigma = 0.010')
        write (2,208)
208     format('  }')
        write (2,209)
209     format('}')
c      Write coordinates to illustrate cgo plane in .cgo file
        do 308, i=1,n
         if((chpdb(i).eq.cha1).and.
     1      (residpdb(i).eq.resid1).and.
     1      (atopdb(i).eq.'N  '))then
          x1=x0(i)
          y1=y0(i)
          z1=z0(i)
         endif
         if((chpdb(i).eq.cha1).and.
     1      (residpdb(i).eq.resid2).and.
     1      (atopdb(i).eq.'C  '))then
          x2=x0(i)
          y2=y0(i)
          z2=z0(i)
         endif
         if((chpdb(i).eq.cha1).and.
     1      (residpdb(i).eq.resid2).and.
     1      (atopdb(i).eq.'O  '))then
          x3=x0(i)
          y3=y0(i)
          z3=z0(i)
         endif
         if((chpdb(i).eq.cha2).and.
     1      (residpdb(i).eq.resid1).and.
     1      (atopdb(i).eq.'N  '))then
          x4=x0(i)
          y4=y0(i)
          z4=z0(i)
         endif
         if((chpdb(i).eq.cha2).and.
     1      (residpdb(i).eq.resid2).and.
     1      (atopdb(i).eq.'C  '))then
          x5=x0(i)
          y5=y0(i)
          z5=z0(i)
         endif
         if((chpdb(i).eq.cha2).and.
     1      (residpdb(i).eq.resid2).and.
     1      (atopdb(i).eq.'O  '))then
          x6=x0(i)
          y6=y0(i)
          z6=z0(i)
         endif
308     continue
        write(4,309)x1,y1,z1,x2,y2,z1,x3,y3,z3
        write(4,309)x1,y1,z1,x2,y2,z1,x6,y6,z6
        write(4,309)x1,y1,z1,x5,y5,z5,x6,y6,z6
        write(4,309)x1,y1,z1,x4,y4,z4,x6,y6,z6
309     format('VERTEX,',f9.1,',',f9.1,',',f9.1,',' 
     1         'VERTEX,',f9.1,',',f9.1,',',f9.1,','
     1         'VERTEX,',f9.1,',',f9.1,',',f9.1,',\')
       endif

c refinement.geometry_restraints.edits {
c   planarity {
c      action = *add delete change
c      atom_selection = name  N  and chain A and resname  VAL and resseq   55 
c                    or name  C  and chain A and resname  THR and resseq   54 
c                    or name  O  and chain A and resname  THR and resseq   54 
c
c                    or name  N  and chain C and resname  VAL and resseq   55 
c                    or name  C  and chain C and resname  THR and resseq   54 
c                    or name  O  and chain C and resname  THR and resseq   54
c      sigma = 0.01
c   }
c }
c now write the distance record to pymol file
       write (3,9)cha1,resid1,ato1,cha2,resid2,ato2
  9    format('distance (chain ',a1,' and resid ',i4,' and name ',a3, 
     1   '), (chain ',a1,' and resid ',i4,' and name ',a3,')')
        go to 4
 99     print*, noatom, ' distance records written'
c     Final touch on the .cgo file
       write(4,312)
312    format('   END, \')
       write(4,313)
313    format('   END ]')
       write(4,314)'cmd.load_cgo(obj, ', '''planes''' ,')'
314    format(a18,a8,a2)
999	stop
	end
