c Renumber chains when you have more than 26.
c Use upper and lower case chain IDs.
c Advance chain ID each if resid n+1 < resid 1.
c Mikey 2018-05-17
c
c version 3, eliminate segid in right side of pdb file
c Put TER and END records at their proper locations in output
c Mikey 2020-08-20
c 
c version 4, carry over the alternate conformation ID in column 17.
c Mikey 2020-08-20
c
c version 5, propagate the occupancy faithfully.
c
	character pdbinn*80,chain*1,atnam*80
        character chainout*2, chainin(99999)*2
        character chainid(62)*2
        character response*1
        integer residprev
        integer neworder(62)
        integer zcomp(62)
        integer nochain(99999)
        character pdbout*80
        character nama(99999)*3,resn(99999)*3,elem(99999)*14
        character altid(99999)*1
        real occ(99999),bf(99999)
        integer noatomin
        integer resida(99999)
        real xa(99999),ya(99999),za(99999)
        real dist(99999)
        chainid(1)=' A'
        chainid(2)=' B'
        chainid(3)=' C'
        chainid(4)=' D'
        chainid(5)=' E'
        chainid(6)=' F'
        chainid(7)=' G'
        chainid(8)=' H'
        chainid(9)=' I'
        chainid(10)=' J'
        chainid(11)=' K'
        chainid(12)=' L'
        chainid(13)=' M'
        chainid(14)=' N'
        chainid(15)=' O'
        chainid(16)=' P'
        chainid(17)=' Q'
        chainid(18)=' R'
        chainid(19)=' S'
        chainid(20)=' T'
        chainid(21)=' U'
        chainid(22)=' V'
        chainid(23)=' W'
        chainid(24)=' X'
        chainid(25)=' Y'
        chainid(26)=' Z'
        chainid(27)=' 0'
        chainid(28)=' 1'
        chainid(29)=' 2'
        chainid(30)=' 3'
        chainid(31)=' 4'
        chainid(32)=' 5'
        chainid(33)=' 6'
        chainid(34)=' 7'
        chainid(35)=' 8'
        chainid(36)=' 9'
        chainid(37)=' a'
        chainid(38)=' b'
        chainid(39)=' c'
        chainid(40)=' d'
        chainid(41)=' e'
        chainid(42)=' f'
        chainid(43)=' g'
        chainid(44)=' h'
        chainid(45)=' i'
        chainid(46)=' j'
        chainid(47)=' k'
        chainid(48)=' l'
        chainid(49)=' m'
        chainid(50)=' n'
        chainid(51)=' o'
        chainid(52)=' p'
        chainid(53)=' q'
        chainid(54)=' r'
        chainid(55)=' s'
        chainid(56)=' t'
        chainid(57)=' u'
        chainid(58)=' v'
        chainid(59)=' w'
        chainid(60)=' x'
        chainid(61)=' y'
        chainid(62)=' z'
        print *,'REMARK type name of input pdb file'
        read(*,'(a80)')pdbinn
	open (1,file=pdbinn,status='old',form='formatted')
        print*,'You said ', pdbinn
        do 532, i=1,80
         if(pdbinn(i:i+3).eq.'.pdb')then
          nlastchar=i-1
          go to 533
         endif
532     continue
        print*, 'The coordinate file must end with .pdb'
        go to 999

533     pdbout=pdbinn(1:nlastchar)//'_rechain.pdb'
        print*, 'output pdbfile is called ',pdbout
        open (2,file=pdbout,status='new',form='formatted')
c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        noatomin=0

  4     read (1,2,end=91)atnam
  2     format(a80)
c       print*, noatomin, atnam
        if(atnam(1:4).ne.'ATOM')then
         if(atnam(1:3).eq.'TER')goto 4
         if(atnam(1:3).eq.'END')goto 4
         write(2,2)atnam
         goto 4
        endif
        if(atnam(18:20).eq.'HOH')goto 4
        noatomin=noatomin+1
        if(noatomin.ge.99990)then
         print*, 'Redimension for more than 99999 atoms'
         go to 999
        endif
        read (atnam(14:16),*)nama(noatomin)
        altid(noatomin)=atnam(17:17)
        read (atnam(18:20),*)resn(noatomin)
        chainin(noatomin)=atnam(21:22)
        read (atnam(23:26),*)resida(noatomin)
        read (atnam(31:38),*)xa(noatomin)
        read (atnam(39:46),*)ya(noatomin)
        read (atnam(47:54),*)za(noatomin)
        read (atnam(55:60),*)occ(noatomin)
        read (atnam(61:66),*)bf(noatomin)
        read (atnam(67:80),14)elem(noatomin)
14      format(a14)
        go to 4
 91     print*,'Finished reading coordinates'
        print*,'REMARK ',noatomin,' atoms in', pdbinn
c
c
c       Here we go
       noatomout=0
c
c      First atom
       print*, 'Starting with chain A'
       nchain=1
       nout=1
       chainout=chainid(nchain)
c      Make segid = chainIDchainIDchainIDchainID
c      elem(1)(7:7)=chainout(2:2)
c      elem(1)(8:8)=chainout(2:2)
c      elem(1)(9:9)=chainout(2:2)
c      elem(1)(10:10)=chainout(2:2)
c      David Boyer says that segid messes up structural alignments in pymol.
c      Replace segid with blanks.
       elem(i)(7:10)='    '
       write(2,96)nout,nama(1),altid(1),resn(1),chainout,
     1     resida(1),xa(1),ya(1),za(1),occ(1),bf(1),elem(1)
c      write(6,96)nout,nama(1),altid(1),resn(1),chainout,
c    1     resida(1),xa(1),ya(1),za(1),occ(1),bf(1),elem(1)
       residprev=resida(1)
c 
c      All other atoms
         do 86, i=2,noatomin
          if(resida(i).ge.residprev)goto 85
          if(resida(i).lt.residprev+1)then
           print*, 'We have come to the following line'
           write(6,96)i,nama(i),altid(i),resn(i),chainin(i),
     1    resida(i),xa(i),ya(i),za(i),occ(i),bf(i),elem(i)
           print*,'I will advance chain ID whether you like it or not!'
           write(2,95)'TER'
 95        format(a3)
c84        read(*,'(a1)')response
c          if(response.eq.'n')goto 85
c          if(response.eq.'y')then 
            nchain=nchain+1
            chainout=chainid(nchain)
            print*, 'Starting chain',chainout
c           goto 85
c          endif
c          print*, 'I did not understand. Type y or n'
c          goto 84
          endif
c         Make segid = chainIDchainIDchainIDchainID
c85       elem(i)(7:7)=chainout(2:2)
c         elem(i)(8:8)=chainout(2:2)
c         elem(i)(9:9)=chainout(2:2)
c         elem(i)(10:10)=chainout(2:2)
c         David Boyer says that segid messes up structural alignments in pymol.
c         Replace segid with blanks.
 85       elem(i)(7:10)='    '
          write(2,96)i,nama(i),altid(i),resn(i),chainout,
     1    resida(i),xa(i),ya(i),za(i),occ(i),bf(i),elem(i)
 96       format('ATOM',i7,2x,a3,a1,a3,a2,i4,4x,3f8.3,2f6.2,a14)
          residprev=resida(i)
 86      continue
         write(2,97)'END'
 97      format(a3)
         print*,'Congratulations! you have',nchain,'chains.'
999	stop 
	end
