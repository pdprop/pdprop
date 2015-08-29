C  plateLongFAD.f   vers. 0.60     FAC aug.22 2015
      SAVE
C  Computes the FAD data from outputs of plateLongSurfFlaw.f programs
C  Compile:  gfortran  -g -w -fbounds-check plateLongFAD.f  -o plateLongFAD
C  Usage:   plateLongFAD    >outputFile

C   The inputfile is a random access (direct) recl=72 file.
C   This program also reads items from the    pdprop.env   file.
C   Items include  B=   Plate thickness in mm
C                  W=   Plate width in mm
C                  Kmat = Fracture stress intensity
C                  PmEOL= End of Life Pm (membrane stress)
C                  PbEOL= End of Life Pb (bending stress)
C                  PinJoint= 0 or 1 for pin jointed structure
C  
C   FAD boundaries are read from file "limitsFAD"  in which are contained
C   the data for  FAD1, FAD2a and FAD2b boundaries.

C  The program is made availble to help students develop advances in crack
C  propagation software.

C  Copyright (C) 2013  Al Conle
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published by
C the Free Software Foundation; either version 2 of the license, or (at
C your option) any later version.
C  This  file is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTA-
C BILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
C License for more details.
C  You should have received a copy of the GNU General PUblic License along
C with this program; if not, write to the Free Software Foundation, Inc.,
C 59 Temple Place -Suite 330, Boston, MA 02111-1307, USA. Try also their
C web site: http://www.gnu.org/copyleft/gpl.html

C vers 0.6  Fix the indexing for plot interval loops   Aug 22 2015
C           (loops start near statement no. 1050 )
C vers 0.5  Try fix when records <100000.  maxintervals = smaller Aug.19 2015
C           Look for lbruteforce in code
C vers. 0.40  Fix major arithm. error in compute for Sigref and SigrefEOL
C             (forgot to include "term2" )                     Aug. 28 2013

C   Expected input format:
C  The numbers were printed by plateLongSurfFlaw.f  statement:
C        nrecord=nrecord+1
C        write(60,rec=nrecord)nrev,totdam90,nblk,nact,
C     &             lobj90,xMm90,xMb90, xMkm90,xMkb90,xfw,
C     &             stsMembrane,stsBending

C  Output file for plotting etc is:


C  Process:
C    1. find the number of the last record in the input file.
C       This is written in the FIRST record, first element.
C       (or one could use an "INQUIRE" function in Fortran.)
C    2. Divide the history into  "maxintervals= 1000"  sections
C       If history contains less than 100,000 points just go to 
C       bruteforce mode and check every peak for FAD.
C    3. Step through each record and fetch the Kmax, PmMax and PbMax etc
C       values. At end of each interval take the max values of variables
C       in that interval and compute the FAD variables for that interval.
C       Also compute the Extreme Lr using PmEOL and PbEOL values.
C    4. Using these resulting max value tables for the intervals, go back
C       thru the suspect intervals and compute exact numbers 

C    (it is not clear it this last point (4.) is needed)
C    


C      Keep same variable names as in plateLongSurfFlaw.f for record reads

      character*300  inp300,jnp300      ! used to read in lines as chars.
      character*1    inpone(300),jnpone(300)
      character*10   inpten(30)
      equivalence (inpone(1), inpten(1), inp300)
      equivalence (jnpone(1),jnp300)

      character*30 firstfield, ctail
      character*80 probType
      integer*4 iargc, nargc
      character*80 argv

      character*30 ctypefad   ! reads in which type of FAD data is.

      real*4 totdam90,lobj90,
     &             xMm90,xMb90,
     &             xMkm90,xMkb90,xfw,
     &             stsMembrane,stsBending
      integer*4  nrev,nblk,nact,nrecord

C     Storage for each interval. The extra 1 in 1001 is for remainders
      real*4  amaxS(1001),cmaxS(1001),
     &    xlobj90maxS(1001),xlobj90minS(1001),
     &    stsMmaxS(1001),stsMminS(1001),stsBmaxS(1001),stsBminS(1001),
     &    xKr90maxS(1001), SrS(1001),xLrS(1001),
     &    yKrEOL90S(1001), SrEOLS(1001),xLrEOLS(1001)
      logical logiFAD1(1001),logiFAD2a(1001),logiFAD2b(1001),
     &        logiEOL(1001)
      integer nrevS(1001)

      real   xKmat, PmEOL, PbEOL, Syield, Sult, Sflow, Emod

      logical logiRunFAD1, logiRunFAD2a, logRunFAD2b
      integer ipinjoint ! =0 if not a pin joint structure.
      logical lbruteforce

C     Storage for the FAD  diagram points
      real*4  xLrfad1(100),yKrfad1(100)
      real*4  xLrfad2a(100),yKrfad2a(100)
      real*4  xLrfad2b(100),yKrfad2b(100)

      maxintervals=1000  !  the dimensions of the interval storage.
      maxfads=100        !  dimension of FAD point storage.
      pi=3.1415927




      write(0,180)
      write(6,180)
  180 format("# plateLongFAD.f  vers. 0.6 starts...")
  190 continue
      write(6,191)
      write(0,191)
  191 format("#Opening pdprop.env file...")
C     Initilize the things to be read to checkable items
C     to make sure they have been entered.
      probType= " "
      Bthick= -1.
      Width=  -1.

      xKmat = -1 ! Fracture stress intensity
      ipinjoint= -1  ! switch for pin or non pin joint structure
      PmEOL=  -1 ! End of Life Pm (membrane stress)
      PbEOL=  -1 ! End of Life Pb (bending stress)


C---------------------------------pdprop.env file reads----------------------
C -----------   Open and read in the pdprop.env  file
C     In this file all lines should begin with a #  or are blank lines
      open(unit=10,file="pdprop.env")
      ninput=0
  800 continue
c     Loop back to here for next input line.

      read(10,"(a300)",end=380)inp300
      ninput=ninput+1
Cdebug      write(0,*)" read input line ",ninput

C     Check for blank line
      if(inp300.eq." ")then
C        write(6,"(a1)")" "
        go to 800
      endif

      if(inpone(1).ne."#")then
C        We may have a data line, or someone screwed up and put the # later in line.
C        See if 1st char is a # later in line
         do 805 i=1,300
           if(inpone(i).eq." ") go to 805
C           No? 1st non-blank found, if # its not a data line
            loc=i
           if(inpone(i).eq."#") then
C            Shift the whole mess over to begining of field
             inew=0
             do 803 j=loc,300
              inew=inew+1
              inpone(inew)=inpone(j)
  803        continue
Cdebug        write(0,*)"Shifted a comment: ",(inpone(j),j=1,inew)
C            Ok, now its a nice comment. Go play with it.
             go to 820
          else
C         The first non blank is not a #
          write(0,*)"# Skipped garbage in *.env file: line no.= ",ninput
          write(0,*)"# Text is= ",inp300
          write(6,*)"# Skipped garbage in *.env file: line no.= ",ninput
          write(6,*)"# Text is= ",inp300
          go to 800
          endif
  805    continue


      endif
  820 continue
C     First char was a #
C     pdprop.env file input line has a # in 1st col.  
C     See if it has a keyword.
      if(inpone(1).ne."#")then
C       something is bad in program
        write(0,*)" ERROR 820, sorry prog. messed up call ? admin"
        stop
      endif

C     Ok, its a nice comment.  Figure out if its a special tag.
      read(inp300,*)firstfield


      if(firstfield .eq."#TYPE=" .or.
     &   firstfield .eq."#Type=" .or.
     &   firstfield .eq."#type=" )then
         read(inp300,*) firstfield, probType
         write(6,833)probType
  833    format("#TYPE= ",a80)
         go to 800
      endif

      if(firstfield .eq."#B=" .or.
     &   firstfield .eq."#b=" )then
         read(inp300,*) firstfield, Bthick
         write(6,844)Bthick
  844    format("#B= ",E14.7," # Thickness, mm.")
         go to 800
      endif

      if(firstfield .eq."#W=" .or.
     &   firstfield .eq."#w=" )then
         read(inp300,*) firstfield, Width
         write(6,845)Width
  845    format("#W= ",E14.7," # width, mm.")
         go to 800
      endif

C      if(firstfield .eq."#ri=" .or.
C     &   firstfield .eq."#Ri" .or.
C     &   firstfield .eq."#RI" )then
C         read(inp300,*) firstfield, rinternal
C         write(6,846)rinternal
C  846    format("#ri= ",E14.7," # interal pipe diam, mm.")
C         go to 800
C      endif

C      Kmat = -1 ! material  Fracture stress intensity
      if(firstfield .eq."#KMAT=" .or.
     &   firstfield .eq."#Kmat=" .or.
     &   firstfield .eq."#kmat=" )then
         read(inp300,*) firstfield, xKmat
         write(6,852)xKmat
  852    format("#Kmat= ",e14.7)
         go to 800
      endif

C      PinJoint = -1 ! Structure is either pin jointed=1 or not=0
      if(firstfield .eq."#PinJoint=" .or.
     &   firstfield .eq."#Pinjoint=" .or.
     &   firstfield .eq."#PINJOINT=" .or.
     &   firstfield .eq."#pinjoint=" .or.
     &   firstfield .eq."#pinJoint=" )then
         read(inp300,*) firstfield, ipinjoint
         write(6,853)ipinjoint
  853    format("#PinJoint= ",i2)
         go to 800
      endif

C      PmEOL= -1 ! End of Life Pm (membrane stress)
      if(firstfield .eq."#PMEOL=" .or.
     &   firstfield .eq."#PmEOL=" .or.
     &   firstfield .eq."#Pmeol=" .or.
     &   firstfield .eq."#pmeol=" )then
         read(inp300,*) firstfield, PmEOL
         write(6,854)PmEOL
  854    format("#PmEOL= ",e14.7)
         go to 800
      endif

C      PbEOL= -1 ! End of Life Pb (bending stress)
      if(firstfield .eq."#PBEOL=" .or.
     &   firstfield .eq."#PbEOL=" .or.
     &   firstfield .eq."#Pbeol=" .or.
     &   firstfield .eq."#pbeol=" )then
         read(inp300,*) firstfield, PbEOL
         write(6,855)PbEOL
  855    format("#PbEOL= ",e14.7)
         go to 800
      endif

  350 continue
C     None of the above?  Then it must be a plain old
C     comment line.   Write it out too.
C       Count backwards and see where the last char is
        do 360 i=1,300
          j=300-(i-1)
          if(inpone(j).ne." ")then
C           found last char
            lastloc=j
            go to 362
          endif
  360   continue

  362   continue
      write(6,"(300a1)")(inpone(i),i=1,lastloc)
C     Go read another line
      go to 800

C     End of pdprop.env file reached
  380 continue
C     Check if critical items have been read in.
      close(unit=10)
      istop=0
      if(probType .eq. " ")then
        write(0,*)"#ERROR: #Type= not found."
        write(6,*)"#ERROR: #Type= not found."
        istop=1
      endif
      if(Bthick .eq. -1.0)then
        write(0,*)"ERROR:  #B=  not found"
        write(6,*)"ERROR:  #B=  not found"
        istop=1
      endif
      if(Width .eq. -1.0)then
        write(0,*)"ERROR:  #W=  not found"
        write(6,*)"ERROR:  #W=  not found"
        istop=1
      endif
C      if(rinternal .eq. -1.0)then
C        write(0,*)"ERROR:  #ri=  not found"
C        write(6,*)"ERROR:  #ri=  not found"
C        istop=1
C      endif
      if(Kmat .eq. -1)then
        write(0,*)"ERROR: #Kmat= not found"
        write(6,*)"ERROR: #Kmat= not found"
        istop=1
      endif
      if(ipinjoint .eq. -1)then
        write(0,*)"ERROR: #PinJoint= not found"
        write(6,*)"ERROR: #PinJoint= not found"
        istop=1
      endif
      if(PmEOL .eq. -1)then
        write(0,*)"ERROR: #PmEOL= not found"
        write(6,*)"ERROR: #PmEOL= not found"
        istop=1
      endif
      if(PbEOL .eq. -1)then
        write(0,*)"ERROR: #PbEOL= not found"
        write(6,*)"ERROR: #PbEOL= not found"
        istop=1
      endif

      if(istop.eq. 1)then
         write(0,*)"# Stopping..."
         write(6,*)"# Stopping..."
         stop
      endif



C-------------- Read in the FAD  file-----------------------------
C    In addition to FADs  we need to get Sy, Su, Emod 

      write(0,701)
      write(6,701)
  701 format("# #Opening fads.table  ..."/  )
      open(unit=10,file="fads.table")

      ninput=0
      nfad1=0
      nfad2a=0
      nfad2b=0
      Syield= -999.  ! use to check if specified in input file.
      Sult= -999.
      Emod= -999.

C     Loop back to here for next input line.
  700 continue

      read(10,"(a300)",end=750)inp300
      ninput=ninput+1
Cdebug      write(0,*)" read input line ",ninput

C     Check for blank line
      if(inp300.eq." ")then
C        write(6,"(a1)")" "
        go to 700
      endif

      if(inpone(1).ne."#")then
C        We may have a data line, or someone screwed up and put the # later in line.
C        See if 1st char is a # later in line
         do 705 i=1,300
           if(inpone(i).eq." ") go to 705
C           No? 1st non-blank found, if # its not a data line
            loc=i
           if(inpone(i).eq."#") then
C            Shift the whole mess over to begining of field
             inew=0
             do 703 j=loc,300
              inew=inew+1
              inpone(inew)=inpone(j)
  703        continue
Cdebug        write(0,*)"Shifted a comment: ",(inpone(j),j=1,inew)
C            Ok, now its a nice comment. Go play with it.
             go to 730
           else
C            The first non blank is not a #
             go to 710
           endif
  705      continue

  710      continue
C        1st non blank is not a #.  Check if its a number.
         if(inpone(loc).ne."+" .and.
     &      inpone(loc).ne."-" .and.
     &      inpone(loc).ne."." .and.
     &      inpone(loc).ne."0" .and.
     &      inpone(loc).ne."1" .and.
     &      inpone(loc).ne."2" .and.
     &      inpone(loc).ne."3" .and.
     &      inpone(loc).ne."4" .and.
     &      inpone(loc).ne."5" .and.
     &      inpone(loc).ne."6" .and.
     &      inpone(loc).ne."7" .and.
     &      inpone(loc).ne."8" .and.
     &      inpone(loc).ne."9"  )then
C           It must be a letter, not a number. Time to bomb out.
            write(0,714)ninput
            write(6,714)ninput
  714       format("# ERROR fads.table: input line no. ",I5,
     &      " not a # and not a number. Edit your FADs table file. ")
            stop
         endif
C        Ok, its a number.
C        Each data line in a fads.table file should have 3 elements. e.g.:
C                  0.0      0.707  #FAD1
C                        or
C                  1.02345  0.300  #FAD2a 
C                        or
C                  1.02345  0.300  #FAD2b 

C         Where the numbers are  Lr  and Kr  points on the FAD boundary.

C         Thus read in the two numbers and the string and then save in
C         the correct storage according to the string.
        read(inp300,*)xvalue,yvalue,ctypefad
        if(ctypefad .eq. "#FAD1" )then
C         type FAD 1   data
          nfad1=nfad1+1
          xLrfad1(nfad1)=xvalue
          yKrfad1(nfad1)=yvalue
          write(6,715)xvalue,yvalue
  715     format("#FAD ",f7.5,1x,f8.4," #FAD1")
          goto 700   !get next point
        endif

        if(ctypefad .eq. "#FAD2a" )then
          nfad2a=nfad2a+1
          if(nfad2a .gt. maxfads)then
             write(0,716)
             write(6,716)
  716        format("#Error: plateLongFAD: too many FAD2a data points."/
     &            "# you need to edit your FAD table file, or "
     &            " recompile the program plateLongFAD.f   Stopping...")
             stop
          endif
          xLrfad2a(nfad2a)=xvalue
          yKrfad2a(nfad2a)=yvalue
          write(6,718)xvalue,yvalue
  718     format("#FAD ",f7.5,1x,f8.4," #FAD2a")
          goto 700   !get next point
        endif

        if(ctypefad .eq. "#FAD2b" )then
          nfad2b=nfad2b+1
          if(nfad2b .gt. maxfads)then
             write(0,721)
             write(6,721)
  721        format("#Error: plateLongFAD: too many FAD2b data points."/
     &            "# you need to edit your FAD table file, or "
     &            " recompile the program plateLongFAD.f   Stopping...")
           stop
          endif
          xLrfad2b(nfad2b)=xvalue
          yKrfad2b(nfad2b)=yvalue
          write(6,722)xvalue,yvalue
  722     format("#FAD ",f7.5,1x,f8.4," #FAD2b" )
          goto 700   !get next point
        endif
      endif  !end of the:     if first char is not a "#" statement

C     We have a comment line begining with #
  730 continue
         read(inp300,*)firstfield

         if(firstfield .eq. "#Sy=" .or.
     &      firstfield .eq. "#SY=" )then
           read(inp300,*)firstfield,Syield
           write(6,732)Syield
  732      format("#fads.table #Found Sy= ",e14.7," (assume MPa!)" )
         endif
         if(firstfield .eq. "#Su=" .or.
     &      firstfield .eq. "#su=" .or.
     &      firstfield .eq. "#Sult=" .or.
     &      firstfield .eq. "#SULT=" )then
           read(inp300,*)firstfield,Sult
           write(6,733)Sult
  733      format("#fads.table #Found Su= ",e14.7," (assume MPa!)" )
         endif
         if(firstfield .eq. "#E=" .or.
     &      firstfield .eq. "#e=" .or.
     &      firstfield .eq. "#EMOD=" .or.
     &      firstfield .eq. "#Emod=" .or.
     &      firstfield .eq. "#emod=" )then
           read(inp300,*)firstfield,Emod
           write(6,734)Emod
  734      format("#fads.table #Found E= ",e14.7," (assume MPa!)" )
         endif
       goto 700

  750 continue  ! end of file comes here
C        All data is in.
      if(Syield .eq. -999.)then
         write(0,752)
  752    format("#Error: #Sy=   not found in  fads.table file",
     &          " Stopping...")
         stop
      endif
      if(Sult .eq. -999.)then
         write(0,754)
  754    format("#Error: #Su=   not found in  fads.table file",
     &          " Stopping...")
         stop
      endif
      if(Emod .eq. -999.)then
         write(0,756)
  756    format("#Error: #E=   not found in  fads.table file",
     &          " Stopping...")
         stop
      endif
C     compute some items for FAD check
      Sflow=(Syield+Sult)/2.0
      if(Sflow .gt. (1.2*Syield) ) Sflow=1.2*Syield
      write(6,760)Sflow
  760 format("#Sflow= ",f7.1," =  (Syield+Sult)/2.0")

      close(unit=10)  !   close fads.table  file



C------------------  scan the fadInput.rand file for max values
      open(unit=60, file="fadInput.rand", access="direct",
     &     recl=48, status="old")
      irec=1
      read(60,rec=irec)nrev,totdam90,nblk,nact,
     &             lobj90,xMm90,xMb90,
     &             xMkm90,xMkb90,xfw,
     &             stsMembrane,stsBending
C     In this first record  nrev   should actually be the number of
C     records written in the file, including this first one.
C     Also, as a check,  totdam90 should be 0.0  and nrev should have
C     been saved as a -ve number (no. of records)
      if(totdam90 .ne. 0.0 .or. totdam00 .ne. 0.0 .or.
     &   nrev .ge. 0)then
C         something is wrong with this random access file
          write(0,1010)nrev,totdam90
          write(6,1010)nrev,totdam90
 1010     format("# Error: 1st rec. of file fadInput.rand  is wrong:"/
     &           "# nrev = ",i10," (should be -ve)"/
     &           "#    a = ",e14.7," should be = 0.0"/
     &           "# Stopping now...")
      endif
      nrev=-nrev  
      maxrecords=nrev

C     If we get to here the 1st rec is ok.  Go read the last rec. and
C     determine what the total reversals of the test was:
      read(60,rec=maxrecords)nrevmax,totdam90,nblk,nact,
     &             lobj90,xMm90,xMb90,
     &             xMkm90,xMkb90,xfw,
     &             stsMembrane,stsBending
      write(0,1012)maxrecords,nrevmax,nrevmax,nblk,nact
      write(6,1012)maxrecords,nrevmax,nrevmax,nblk,nact
 1012 format("#MAXRECORDS= ",i10," In last rec. Reversal= ",i10/
     &       "#MAXREVERSALS= ",i10/
     &       "#MAXBLOCKS= ",i10/
     &       "#NACT= ",i10)

C     Divide the data into about 1000 intervals for quick plotting
      nitems= maxrecords-1
      nIntervals= maxIntervals 
      if( nitems .le. nIntervals)then
          nIntervals= 1
C         There are not that many recs items in this history.
          write(0,1015)maxrecords
          write(6,1015)maxrecords
 1015     format("# Warning! plateLongFAD: Max Recs ",i10," < 1000  ",
     &           "Will try to continue  FAD...")
      endif

C     We are dividing the history into nIntervals 
C     At this point nIntervals is either 1 or 1000
      nrecsPerInt= nitems / nIntervals
C     Due to roundoff, there may be extra recs. past last interval
      nremainder= nitems - ( nIntervals * nrecsPerInt )
      if(nremainder .lt. 0)nremainder=0 !  only 1 interval exists.
      write(0,1020)nitems, nIntervals,nrecsPerInt,nremainder
      write(6,1020)nitems, nIntervals,nrecsPerInt,nremainder
 1020 format("#plateLongFAD: For plotting: "/
     &      "#   No. Recs with reversals: nitems=      ",i10/
     &      "#   No. of plot Intervals:   nIntervals=  ",i10/
     &      "#   No. Recs per Interval:   nrecsPerInt= ",i10/
     &      "#   No. of leftover items:   nremainder=  ",i10)
      
      interval=0
 1050 continue  !-------------------- Begin of loop for each interval
      interval=interval+1
      if(interval .eq. 1)then  !initilize 1st interval limits
         jrecStart=2
         jrecEnd= jrecStart+nrecsPerInt
         goto 1060
      endif

      if(interval .gt. nIntervals)then 
C       All the full intervals have been done.  Check if remaider?
        if(nremainder .eq. 0)goto 3900 !all done, exit interval loop
C       There are one or more remainders
        jrecStart=jrecEnd+1
        jrecEnd=maxrecords
        if(jrecStart .gt.maxrecords)then !we are on remainders 2nd time.
          write(0,*)"#Done with remainders"
          goto 3900  !finished scans
        endif
        write(0,*)"#scanning remainders, jrecStart,jrecEnd= ",
     &     jrecStart,jrecEnd
C       It is possible, with a remainder=1 that jrecStart=jrecEnd
        goto 1060
      endif

C     Its just a regular interval
      jrecStart=jrecEnd+1
      jrecEnd= jrecStart+nrecsPerInt-1


 1060 continue  ! limits are set start interval or remainder loop
C        Scan this interval of records
C        Read in the first rec of the interval to set the max mins
         read(60,rec=jrecStart)nrev,totdam90,nblk,nact,
     &             lobj90,xMm90,xMb90,
     &             xMkm90,xMkb90,xfw,
     &             stsMembrane,stsBending
             xlobj90max=lobj90 !stress intensities
             xlobj90min=lobj90

             stsMmax=stsMembrane !Stresses
             stsMmin=stsMembrane
             stsBmax=stsBending
             stsBmin=stsBending
      
        if(jrecStart .eq. jrecEnd)then ! happens when remainder=1
           jrec=jrecEnd
           write(0,*)"#plotLongFAD: jrecStart=jrecEnd= ",jrecStart
           write(6,*)"#plotLongFAD: jrecStart=jrecEnd= ",jrecStart
           read(60,rec=jrec)nrev,totdam90,nblk,nact,
     &             lobj90,xMm90,xMb90,
     &             xMkm90,xMkb90,xfw,
     &             stsMembrane,stsBending
           xlobj90max=lobj90
           xlobj90min=lobj90
           stsMmax=stsMembrane
           stsMmin=stsMembrane
           stsBmax=stsBending
           stsBmin=stsBending
           goto 1901 !skip the interval scan
        endif 

C       Ok, its a regular set of interval points or a set of remainders
        do 1900 jrec=jrecStart+1,jrecEnd
           read(60,rec=jrec)nrev,totdam90,nblk,nact,
     &             lobj90,xMm90,xMb90,
     &             xMkm90,xMkb90,xfw,
     &             stsMembrane,stsBending

           if(xlobj90max .lt. lobj90)xlobj90max=lobj90
           if(xlobj90min .gt. lobj90)xlobj90min=lobj90

           if(stsMmax .lt. stsMembrane)stsMmax=stsMembrane
           if(stsMmin .gt. stsMembrane)stsMmin=stsMembrane

           if(stsBmax .lt. stsBending)stsBmax=stsBending
           if(stsBmin .gt. stsBending)stsBmin=stsBending

 1900    continue  !done scans in this interval
 1901    continue  !come here if remainder =1 and no loop required


         amax=totdam90 
         nrevIntMax=nrev

C        Now for this interval, given the various maxima
C        compute  Kr, Lr, and the KrEOL, and LrEOL

         if(interval .gt. maxIntervals+1)then
           write(0,*)"#Interval Storage Error: ",
     &     "interval, jrecStart,jrecEnd,jrec= ",
     &     interval, jrecStart,jrecEnd,jrec
         endif
         amaxS(interval)=amax
         nrevS(interval)=nrevIntMax
         xlobj90maxS(interval)= xlobj90max
         xlobj90minS(interval)= xlobj90min

         stsMmaxS(interval)= stsMmax
         stsMminS(interval)= stsMmin
         stsBmaxS(interval)= stsBmax
         stsBminS(interval)= stsBmin

C        Compute Kr  (no correction for secondary stresses)
         xKr90maxS(interval)= xlobj90max/xKmat

C        compute Sigma_ref:  (depends on if structure is pin jointed)
C        See BS7910 Annex P.3.3 page 240
C        (same as plateSurf flaw except for   a2prime   calc. )
           a2prime=(amax/Bthick)
           term2= 3.0 *( (1.0-a2prime)**2 )

         if(ipinjoint .eq. 0) then
C          No, not pin jointed
           Sigref=stsBmax + sqrt( stsBmax**2 + (9.0*stsMmax**2) *
     &            (1.0-a2prime)**2 )
           Sigref=Sigref/term2
         else
           Sigref=stsBmax +3.0*stsMmax*a2prime +
     &         sqrt (  (stsBmax + 3.0*stsMmax*a2prime)**2  +
     &                 (9.0*stsMmax**2) * (1.0-a2prime)**2
     &              )
           Sigref=Sigref/term2
         endif
         SrS(interval)= Sigref/Sflow
         xLrS(interval)= Sigref/Syield


C        Compute End of Life (severe storm near end of life) Kr Sr Lr
C        EOL stresses are in PmEOL and PbEOL
C        We will use the amax, cmax, and the various Mm, Mkm  etc factors
C        from the last record  read in this interval.
         rootPiA=sqrt(pi*totdam90)
         xK90=xfw*(xMkm90*xMm90*PmEOL + xMkb90*xMb90*PbEOL )
         xK90=xK90*rootPiA


         yKrEOL90S(interval)=xK90/xKmat

C        compute Sigma_refEOL:  (depends on if structure is pin jointed)
C        a2prime is still same
         if(ipinjoint .eq. 0) then
C          No, not pin jointed
           SigrefEOL=PbEOLstsBmax + sqrt( PbEOL**2 + (9.0*PmEOL**2) *
     &            (1.0-a2prime)**2 )
           SigrefEOL=SigrefEOL/term2
         else
           SigrefEOL=PbEOL +3.0*PmEOL*a2prime +
     &         sqrt (  (PbEOL + 3.0*PmEOL*a2prime)**2  +
     &                 (9.0*PmEOL**2) * (1.0-a2prime)**2
     &              )
           SigrefEOL=SigrefEOL/term2
         endif
         SrEOLS(interval)= SigrefEOL/Sflow
         xLrEOLS(interval)= SigrefEOL/Syield

C        nblk and nact are simply the last values of the interval


C        write out one long row for everything in this interval
         i=interval
         write(6,3810)amaxS(i),nrevS(i),nblk,nact
     &     ,xlobj90maxS(i),xlobj90minS(i)
     &     ,stsMmaxS(i),stsMminS(i),stsBmaxS(i),stsBminS(i)
     &     ,xkr90maxS(i), SrS(i),xLrS(i)
     &     ,yKrEOL90S(i), SrEOLS(i),xLrEOLS(i)
Cdebug     &     ,jrecStart,jrecEnd,i
     &     ,xMm90,xMb90
     &     ,xMkm90,xMkb90
     &     ,xfw
 3810    format("#FADints ",e10.3,3(1x,i10)
     &      ,2(f6.0,1x)
     &      ,4(f6.1,1x)
     &      ,3(f6.3,1x)
     &      ,3(f6.3,1x)
Cdebug     &      ,1x,i7,1x,i7,1x,i7
     &      ,2(f6.4,1x)
     &      ,2(f6.4,1x)
     &      ,f6.4
     &   )

C     Now check for exceeding  FADs

C     All done with this interval. go to next one
      goto 1050 !go back up for next interval

 3900 continue  !exit from  interval loop


C     Detailed interval inspection 
 5000 continue
C     It is not clear if this is necessary

 9000 continue
      close(unit=60)
      stop
      end

