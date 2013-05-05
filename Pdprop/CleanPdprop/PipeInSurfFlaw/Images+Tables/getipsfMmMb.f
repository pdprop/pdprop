C======================================================================
      SUBROUTINE getMmMb00a90(ivec,iret,aoc,aob,xMm00,xMb00,xMm90,xMb90)
      SAVE
C     ivec= function instruction  =1=interpolate         (input)
C                                 =2 turn on debug
C     iret= return value = 0 if no errors                (return)
C     aoc=  a/c  crack depth / half crack surface length (input)
C     aob=  a/B  crack depth / Plate thickness           (input)
C     xMm00=  multiplication factor for membrane stress  (return)
C     xMb00=  multiplication factor for bending stress   (return)
C     xMm90=  multiplication factor for membrane stress  (return)
C     xMb90=  multiplication factor for bending stress   (return)

C     s/r to interpolate for Mkm and Mkb for: both 00 deg and 90 deg
C      (where 0 deg is along surface of plate, and 90 deg is deepest pt. of crack)
C     As long as the dimensions of the storage matrices are the same we
C     will save some computation time in calculating the location of the solutions
C     in the matrix because aoc and aob  are the same for both 00deg and 90deg.
C     i.e. the four pts in the matrix will have the same j,i  co-ordinates.

C    The input file for 90deg looks like this:      ( 00deg is similar)
C# getMmMb vers. 0.4 starts...
C#ROWS=   51
C#COLS=   11
C#theta=  0.1570796E+01   #radians =  90.00 deg.
C#iaoc iaob a/c a/B         Mm          Mb
C  1  1 0.000  0.00  0.1130000E+01  0.1130000E+01
C  1  2 0.000  0.10  0.1170396E+01  0.1034045E+01
C  1  3 0.000  0.20  0.1307138E+01  0.1016954E+01
C  1  4 0.000  0.30  0.1586888E+01  0.1084638E+01
C  1  5 0.000  0.40  0.2087415E+01  0.1252449E+01
C  1  6 0.000  0.50  0.2917596E+01  0.1539032E+01
C  1  7 0.000  0.60  0.4217416E+01  0.1965316E+01
C  1  8 0.000  0.70  0.6157967E+01  0.2558635E+01
C  1  9 0.000  0.80  0.8941448E+01  0.3361984E+01
C  1 10 0.000  0.90  0.1280117E+02  0.4448406E+01
C  1 11 0.000  1.00  0.1800154E+02  0.5940509E+01
C
C  2  1 0.040  0.00  0.1122352E+01  0.1122352E+01
C  2  2 0.040  0.10  0.1154351E+01  0.1018274E+01
C etc...

C The data was read into common by   readMmMb00.f   and readMmMb90.f


C     Matrices for storing Mm, Mb
      real*4 mb00(51,11),mm00(51,11),mb90(51,11),mm90(51,11)
      integer*4 nmmdatarow,maxmmdatarow, nmmdatacol, maxmmdatacol
C     Limits and intervals of input data used for quicker interpolation
      real*4 startaoc, endaoc, deltaaoc
      real*4 startaob, endaob, deltaaob
      logical debugMm,lactivateMmMb
      common/MMDATA/   mb00,mm00,mb90,mm90,
     &       nmmdatarow,maxmmdatarow, nmmdatacol, maxmmdatacol,
     &      startaoc,endaoc,deltaaoc, startaob,endaob,deltaaob,debugMm,
     &      lactivateMmMb

      logical debug
      debug= .false.
      if(ivec.eq.2)debug= .true.   !set to true if testing code
      if(debugMm)debug= .true.

C---------------------Interpolate ------------------
C
C     Use 2D interpolation: http://en.wikipedia.org/wiki/Bilinear_interpolation
C     Given the aoc and aob inputs, we need to establish the four known
C     matrix points that surround this  x,y pair.   It may also be that we
C     are directly on one of the four points.

C     The matrix row and cols are equal sized increments.  Thus we can "jump"
C     to the relative co-ords without  step wise searching.

C      real*4 startaoc, endaoc, deltaaoc
C      real*4 startaob, endaob, deltaaob
      i=ifix((aoc-startaoc)/deltaaoc)+1
      j=ifix((aob-startaob)/deltaaob)+1
      if(i.lt.1 .or. j.lt.1 .or. j.gt.maxmmdatacol .or.
     &   i.gt. maxmmdatarow) then
C        Out of bounds error
         write(0,1010)startaoc,aoc,endaoc,
     &                startaob,aob,endaob,i,j
 1010    format("#ERROR: S/R getMmMb00a90 a/c or a/b out of bounds:"/
     &    "#      start a/c=",f6.3," a/c= ",e14.7," end a/c= ",f6.3/
     &    "#      start a/b=",f6.3," a/b= ",e14.7," end a/b= ",f6.3/
     &    "#      i,j=  ",i5,1x,i5/
     &    "# Stopping...")
          iret=1010    ! set the error flag to format sta. no.
          return
      endif

C     j,i  is the corner (see internet address above)  of x1,y1
      ip1=i+1
      jp1=j+1
C     We need to use x1, x2, y1, y2 to calculate Mm and Mb
C     x,y  are aob,aoc
      x=aob
      y=aoc

C     Things appear to be within bounds of matrices, continue
      x1=startaob+float(j-1)*deltaaob
      x2=x1+deltaaob
      y1=startaoc + float(i-1)*deltaaoc
      y2=y1+deltaaoc
      if(debug)then
        write(0,1015)x1,y1,  x2,y1, x,y, x1,y2, x2,y2
        write(6,1015)x1,y1,  x2,y1, x,y, x1,y2, x2,y2
 1015   format("#TargetCo-ords:"/
     &        "#x1,y1       x2,y1  :",2f7.4,14x,2f7.4/
     &        "#       x,y         :",15x,2f7.4/
     &        "#x1,y2       x2,y2  :",2f7.4,14x,2f7.4/
     &  )
        write(0,1016)i,j,  i,jp1,  ip1,j, ip1,jp1
        write(6,1016)i,j,  i,jp1,  ip1,j, ip1,jp1
 1016   format("# i,j       i,j+1    : "i3,",",i3, 10x, i3,",",i3/
     &         "# i+1,j     i+1,j+1  : "i3,",",i3, 10x, i3,",",i3/
     &  )
      endif

      Q11=mm90(i,j)
      Q12=mm90(ip1,j)
      Q21=mm90(i,jp1)
      Q22=mm90(ip1,jp1)

C     Ok, lets find  R1 on the line between q11 and q21
C     Original formula is:
C      R1= ( (x2-x)/(x2-x1))*Q11 + ((x-x1)/(x2-x1))*Q21
C      R2= ( (x2-x)/(x2-x1))*Q12 + ((x-x1)/(x2-x1))*Q22
C     Since we use same for Mm and Mb compute once
      x2mx=x2-x
      x2mx1 = x2-x1
      xmx1= x-x1
      frac1= x2mx/x2mx1
      frac2= xmx1/x2mx1
      
C      R1=(x2mx/x2mx1)*Q11 + (xmx1/x2mx1)*Q21
C      R2=(x2mx/x2mx1)*Q12 + (xmx1/x2mx1)*Q22
      R1 = frac1*Q11 + frac2*Q21
      R2 = frac1*Q12 + frac2*Q22
      
      y2my1=y2-y1
      frac1y= (y2-y)/y2my1 
      frac2y= (y-y1)/y2my1
C      xmm = ( (y2-y)/y2my1 )*R1 + ( (y-y1)/y2my1 )*R2
      xMm90 = frac1y*R1 + frac2y*R2
      if(debug)then
        write(0,1020)aoc,aob,i,j, Q11,R1,Q21, xMm90, Q12,R2,Q22
        write(6,1020)aoc,aob,i,j, Q11,R1,Q21, xMm90, Q12,R2,Q22
 1020   format("#Solution for a/c= ",f6.3," a/b= ",f6.3," (i,j)=",2i3/
     &        "# Q11,  R1,  Q21 :  ",3(f7.4,1x)/
     &        "#      Mm90=     :  ",10x,f7.4/
     &        "# Q12,  R2,  Q22 :  ",3(f7.4,1x)/
     &        )
      endif

C     Now do the  Mb  interp.
      Q11=mb90(i,j)
      Q12=mb90(ip1,j)
      Q21=mb90(i,jp1)
      Q22=mb90(ip1,jp1)

      R1 = frac1*Q11 + frac2*Q21
      R2 = frac1*Q12 + frac2*Q22
      xMb90 = frac1y*R1 + frac2y*R2
      if(debug)then
        write(0,1022)aoc,aob,i,j, Q11,R1,Q21, xMb90, Q12,R2,Q22
        write(6,1022)aoc,aob,i,j, Q11,R1,Q21, xMb90, Q12,R2,Q22
 1022   format("#Solution for a/c= ",f6.3," a/b= ",f6.3," (i,j)=",2i3/
     &        "# Q11,  R1,  Q21 :  ",3(f7.4,1x)/
     &        "#      Mb90=     :  ",10x,f7.4/
     &        "# Q12,  R2,  Q22 :  ",3(f7.4,1x)/
     &        )
      endif

C----------Now do the  Mb00  and Mm00   interpolations
C          Note!!! this only works if the matrices have same dimensions !!
C                  and same aoc and aob  for each entry !!
      Q11=mm00(i,j)
      Q12=mm00(ip1,j)
      Q21=mm00(i,jp1)
      Q22=mm00(ip1,jp1)

      R1 = frac1*Q11 + frac2*Q21
      R2 = frac1*Q12 + frac2*Q22
      
      xMm00 = frac1y*R1 + frac2y*R2
      if(debug)then
        write(0,1024)aoc,aob,i,j, Q11,R1,Q21, xMm00, Q12,R2,Q22
        write(6,1024)aoc,aob,i,j, Q11,R1,Q21, xMm00, Q12,R2,Q22
 1024   format("#Solution for a/c= ",f6.3," a/b= ",f6.3," (i,j)=",2i3/
     &        "# Q11,  R1,  Q21 :  ",3(f7.4,1x)/
     &        "#      Mm00=     :  ",10x,f7.4/
     &        "# Q12,  R2,  Q22 :  ",3(f7.4,1x)/
     &        )
      endif

C     Now do the  Mb  interp.
      Q11=mb00(i,j)
      Q12=mb00(ip1,j)
      Q21=mb00(i,jp1)
      Q22=mb00(ip1,jp1)

      R1 = frac1*Q11 + frac2*Q21
      R2 = frac1*Q12 + frac2*Q22
      xMb00 = frac1y*R1 + frac2y*R2
      if(debug)then
        write(0,1026)aoc,aob,i,j, Q11,R1,Q21, xMb00, Q12,R2,Q22
        write(6,1026)aoc,aob,i,j, Q11,R1,Q21, xMb00, Q12,R2,Q22
 1026   format("#Solution for a/c= ",f6.3," a/b= ",f6.3," (i,j)=",2i3/
     &        "# Q11,  R1,  Q21 :  ",3(f7.4,1x)/
     &        "#      Mb00=     :  ",10x,f7.4/
     &        "# Q12,  R2,  Q22 :  ",3(f7.4,1x)/
     &        )
      endif



      iret=0
      return
      end

C==============================================================================
      SUBROUTINE readMmMb00(ivec,iret)
      SAVE
C     ivec= function instruction =0=read in,  =2=debug
C     iret= return value = 0 if no errors
C     aoc=  a/c  in mm. of crack depth / half crack surface length (input)
C     aob=  a/B  in mm. of crack depth / Plate thickness (input)
C     xMm=  multiplication factor for membrane stress (return)
C     xMb=  multiplication factor for bending stress (return)

C     s/r to read in the Mkm and Mkb for:
C------00 deg   00 deg   00 deg   00 deg   00 deg   00 deg   00 deg   00 deg   00 deg
C      (where 00 deg is along surface of plate, and 90 deg is deepest pt. of crack)

C    The input file (t_MmMb_Surflaw_00) looks like this:
C# getMmMb vers. 0.4 starts...
C#ROWS=   51
C#COLS=   11
C#theta=  0.0000000E+00   #radians =   0.00 deg.
C#iaoc iaob a/c a/B         Mm          Mb
C  1  1 0.000  0.00  0.0000000E+00  0.0000000E+00
C  1  2 0.000  0.10  0.0000000E+00  0.0000000E+00
C  1  3 0.000  0.20  0.0000000E+00  0.0000000E+00
C  1  4 0.000  0.30  0.0000000E+00  0.0000000E+00
C  1  5 0.000  0.40  0.0000000E+00  0.0000000E+00
C  1  6 0.000  0.50  0.0000000E+00  0.0000000E+00
C  1  7 0.000  0.60  0.0000000E+00  0.0000000E+00
C  1  8 0.000  0.70  0.0000000E+00  0.0000000E+00
C  1  9 0.000  0.80  0.0000000E+00  0.0000000E+00
C  1 10 0.000  0.90  0.0000000E+00  0.0000000E+00
C  1 11 0.000  1.00  0.0000000E+00  0.0000000E+00
C
C  2  1 0.040  0.00  0.2469174E+00  0.2469174E+00
C  2  2 0.040  0.10  0.2547652E+00  0.2459911E+00
C  2  3 0.040  0.20  0.2797244E+00  0.2604570E+00
C  2  4 0.040  0.30  0.3261514E+00  0.2924534E+00
C  2  5 0.040  0.40  0.4016670E+00  0.3463333E+00
C  2  6 0.040  0.50  0.5176973E+00  0.4285498E+00
C  2  7 0.040  0.60  0.6902308E+00  0.5476015E+00
C  2  8 0.040  0.70  0.9407914E+00  0.7139853E+00
C  2  9 0.040  0.80  0.1297628E+01  0.9401052E+00
C  2 10 0.040  0.90  0.1797119E+01  0.1240084E+01
C  2 11 0.040  1.00  0.2485397E+01  0.1629426E+01
C
C  3  1 0.080  0.00  0.3454381E+00  0.3454381E+00
C etc...

      character*300  inp300,jnp300,Cwebpage
      character*1    inpone(300),jnpone(300)
      character*10   inpten(30)
      equivalence (inpone(1), inpten(1), inp300)
      equivalence (jnpone(1),jnp300)

C     Matrices for storing Mm, Mb
      real*4 mb00(51,11),mm00(51,11),mb90(51,11),mm90(51,11)
      integer*4 nmmdatarow,maxmmdatarow, nmmdatacol, maxmmdatacol
C     Limits and intervals of input data used for quicker interpolation
      real*4 startaoc, endaoc, deltaaoc
      real*4 startaob, endaob, deltaaob
      logical debugMm,lactivateMmMb
      common/MMDATA/   mb00,mm00,mb90,mm90,
     &       nmmdatarow,maxmmdatarow, nmmdatacol, maxmmdatacol,
     &      startaoc,endaoc,deltaaoc, startaob,endaob,deltaaob,debugMm,
     &      lactivateMmMb


      ndatalines=0
      ninput=0


C------------------------ Read Input file --------------------------------
C     The matrix max limit values have been set at begin of mainline
      write(0,402)
      write(6,402)
  402 format("# Opening Mm_Mb_00deg. table, file= t_MmMb_Surflaw_00")
      open(unit=10,file= "t_MmMb_Surflaw_00",status="old")

  400 continue     !Loop back to here for next input line.
      read(10,"(a300)",end=450)inp300
      ninput=ninput+1

      if(inp300.eq." ")then    ! Check for blank line
C        write(6,"(a1)")" "
        go to 400
      endif

      if(inpone(1).ne."#")then
C        We may have a data line, or someone screwed up and put the # later in line.
C        See if 1st char is a # later in line
         do 405 i=1,300
           if(inpone(i).eq." ") go to 405
C           No? 1st non-blank found, if # its not a data line
            loc=i
           if(inpone(i).eq."#") then
C            Shift the whole mess over to begining of field
             inew=0
             do 403 j=loc,300
              inew=inew+1
              inpone(inew)=inpone(j)
  403        continue
Cdebug        write(0,*)"Shifted a comment: ",(inpone(j),j=1,inew)
C            Ok, now its a nice comment. Go play with it.
             go to 430
           else
C            The first non blank is not a #
              go to 410
           endif
  405      continue
           go to 400  !assume garbage in line.

  410      continue
C        1st non blank is not a #  Check if its a number.
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
            write(0,415)ninput
            write(6,415)ninput
  415       format(" ERROR MmMbfile00: input line no. ",I5," 1st Char.",
     &      " not a # and not a number. Fix data in file MmMbfile00")
            stop
         endif
C      real*4 mb00(51,11),mm00(51,11),mb90(51,11),mm90(51,11)
C      integer*4 nmmdatarow,maxmmdatarow, nmmdatacol, maxmmdatacol
C        Ok, its a number
         read(inp300,*) iaoc,iaob,aoc,aob, xm,xb
         ndatalines=ndatalines+1
         if(iaoc .gt. maxmmdatarow .or. iaob .gt. maxmmdatacol)then
           write(0,416)iaoc,iaob
           write(6,416)iaoc,iaob
  416      format(" Error: MmMbfile00: too many data points:",
     &            "rows: iaoc= ",i5," cols: iaob= "i5/
     &            " recompile pipeIntSurfFlaw.f or reduce file data")
           stop
         endif
C        Data is within size of matrix. Place into matrices
         mm00(iaoc,iaob)= xm
         mb00(iaoc,iaob)= xb
         write(6,420)iaoc,iaob,aoc,aob,xm,xb
  420    format("#MmMbfile00 ",i4,i4,2(1x,f6.3),2(1x,e14.7))
         if(iaoc.eq.1 .and. iaob.eq.1)then
            startaoc=aoc
            startaob=aob
            write(0,424)startaoc,startaob
            write(6,424)startaoc,startaob
  424       format("#MmMbfile00 #STARTAOC= ",e14.7/
     &             "#MmMbfile00 #STARTAOB= ",e14.7)
         endif
         if(iaoc.eq.maxmmdatarow .and. iaob .eq. maxmmdatacol)then
            endaoc=aoc
            endaob=aob
            write(0,425)endaoc,endaob
            write(6,425)endaoc,endaob
  425       format("#MmMbfile00 #ENDAOC= ",e14.7/
     &             "#MmMbfile00 #ENDAOB= ",e14.7)
         endif
         go to 400
      endif


  430    continue
C        All comment lines are just written out, or deleted
         write(6,431)(inpone(i2),i2=1,80)
  431    format("#MmMbfile00 ",80a1)
         go to 400

  450    continue
C        All data is in. Close file
         close(unit=10)
         if(iaoc .ne. maxmmdatarow .or. iaob .ne. maxmmdatacol)then
           write(0,453)iaoc,maxmmdatarow,iaob,maxmmdatacol
           write(6,453)iaoc,maxmmdatarow,iaob,maxmmdatacol
 453       format("#WARNING! WARNING! WARNING!  Matrix not full"/
     &            "#iaoc= ",i4,"maxmmdatarow= ",i4/
     &            "#iaob= ",i4,"maxmmdatacol= ",i4/
     &            "# RESETTING maxmmdatarow and maxmmdatacol"/)
           maxmmdatarow=iaoc
           maxmmdatacol=iaob
           endaoc=aoc
           endaob=aob
           write(0,425)endaoc,endaob
           write(6,425)endaoc,endaob
         endif

         deltaaoc=(endaoc-startaoc)/(maxmmdatarow-1)
         deltaaob=(endaob-startaob)/(maxmmdatacol-1)
         write(0,451)deltaaoc,deltaaob
         write(6,451)deltaaoc,deltaaob
  451    format("#MmMbfile00 #DELTAAOC= ",e14.7/
     &          "#MmMbfile00 #DELTAAOB= ",e14.7)
         write(6,452)
         write(0,452)
  452    format("#MmMbfile00: #MmMb00 data input completed. ")
         if(ndatalines .eq. 0)then   ! ? no data lines ?
           write(0,455)
           write(6,455)
  455      format("ERROR: No data lines in file t_MmMb_Surflaw_00a"/
     &            " Stopping")
           iret=455
           return
         endif

      iret=0  !no errors
      return
      end

C====================================================================================
      SUBROUTINE readMmMb90(ivec,iret)
      SAVE
C     ivec= function instruction =0=read in
C     iret= return value = 0 if no errors
C     aoc=  a/c  in mm. of crack depth / half crack surface length 
C     aob=  a/B  in mm. of crack depth / Plate thickness 

C     s/r to read in the Mkm and Mkb for:
C------90 deg   90 deg   90 deg   90 deg   90 deg   90 deg   90 deg   90 deg   90 deg
C      (where 0 deg is along surface of plate, and 90 deg is deepest pt. of crack)

C    The input file looks like this:
C# getMmMb vers. 0.4 starts...
C#ROWS=   51
C#COLS=   11
C#theta=  0.1570796E+01   #radians =  90.00 deg.
C#iaoc iaob a/c a/B         Mm          Mb
C  1  1 0.000  0.00  0.1130000E+01  0.1130000E+01
C  1  2 0.000  0.10  0.1170396E+01  0.1034045E+01
C  1  3 0.000  0.20  0.1307138E+01  0.1016954E+01
C  1  4 0.000  0.30  0.1586888E+01  0.1084638E+01
C  1  5 0.000  0.40  0.2087415E+01  0.1252449E+01
C  1  6 0.000  0.50  0.2917596E+01  0.1539032E+01
C  1  7 0.000  0.60  0.4217416E+01  0.1965316E+01
C  1  8 0.000  0.70  0.6157967E+01  0.2558635E+01
C  1  9 0.000  0.80  0.8941448E+01  0.3361984E+01
C  1 10 0.000  0.90  0.1280117E+02  0.4448406E+01
C  1 11 0.000  1.00  0.1800154E+02  0.5940509E+01
C
C  2  1 0.040  0.00  0.1122352E+01  0.1122352E+01
C  2  2 0.040  0.10  0.1154351E+01  0.1018274E+01
C etc...

      character*300  inp300,jnp300,Cwebpage
      character*1    inpone(300),jnpone(300)
      character*10   inpten(30)
      equivalence (inpone(1), inpten(1), inp300)
      equivalence (jnpone(1),jnp300)

C     Matrices for storing Mm, Mb
      real*4 mb00(51,11),mm00(51,11),mb90(51,11),mm90(51,11)
      integer*4 nmmdatarow,maxmmdatarow, nmmdatacol, maxmmdatacol
C     Limits and intervals of input data used for quicker interpolation
      real*4 startaoc, endaoc, deltaaoc
      real*4 startaob, endaob, deltaaob
      logical debugMm,lactivateMmMb
      common/MMDATA/   mb00,mm00,mb90,mm90,
     &       nmmdatarow,maxmmdatarow, nmmdatacol, maxmmdatacol,
     &      startaoc,endaoc,deltaaoc, startaob,endaob,deltaaob,debugMm,
     &      lactivateMmMb


      ndatalines=0
      ninput=0


C------------------------ Read Input file --------------------------------
C     The matrix max limit values have been set at begin of mainline
      write(0,402)
      write(6,402)
  402 format("# Opening Mm_Mb_90deg. table, file= t_MmMb_Surflaw_90")
      open(unit=10,file= "t_MmMb_Surflaw_90",status="old")

  400 continue     !Loop back to here for next input line.
      read(10,"(a300)",end=450)inp300
      ninput=ninput+1

      if(inp300.eq." ")then    ! Check for blank line
C        write(6,"(a1)")" "
        go to 400
      endif

      if(inpone(1).ne."#")then
C        We may have a data line, or someone screwed up and put the # later in line.
C        See if 1st char is a # later in line
         do 405 i=1,300
           if(inpone(i).eq." ") go to 405
C           No? 1st non-blank found, if # its not a data line
            loc=i
           if(inpone(i).eq."#") then
C            Shift the whole mess over to begining of field
             inew=0
             do 403 j=loc,300
              inew=inew+1
              inpone(inew)=inpone(j)
  403        continue
CdebugMm        write(0,*)"Shifted a comment: ",(inpone(j),j=1,inew)
C            Ok, now its a nice comment. Go play with it.
             go to 430
           else
C            The first non blank is not a #
              go to 410
           endif
  405      continue
           go to 400  !assume garbage in line.

  410      continue
C        1st non blank is not a #  Check if its a number.
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
            write(0,415)ninput
            write(6,415)ninput
  415       format(" ERROR readMmMb90: input line no. ",I5," 1st Char.",
     &      " not a # and not a number. Fix data in file mkfile00")
            stop
         endif
C      real*4 mb00(51,11),mm00(51,11),mb90(51,11),mm90(51,11)
C      integer*4 nmmdatarow,maxmmdatarow, nmmdatacol, maxmmdatacol
C        Ok, its a number
         read(inp300,*) iaoc,iaob,aoc,aob, xm,xb
         ndatalines=ndatalines+1
         if(iaoc .gt. maxmmdatarow .or. iaob .gt. maxmmdatacol)then
           write(0,416)iaoc,iaob
           write(6,416)niaoc,iaob
  416      format(" Error: readMmMb90: too many data points:",
     &            "rows: iaoc= ",i5," cols: iaob= "i5/
     &            " recompile pipeIntSurfFlaw.f or reduce file data")
           stop
         endif
C        Data is within size of matrix. Place into matrices
         mm90(iaoc,iaob)= xm
         mb90(iaoc,iaob)= xb
         write(6,420)iaoc,iaob,aoc,aob,xm,xb
  420    format("#MmMbfile90 ",i4,i4,2(1x,f6.3),2(1x,f7.4))
         if(iaoc.eq.1 .and. iaob.eq.1)then
            startaoc=aoc
            startaob=aob
            write(0,424)startaoc,startaob
            write(6,424)startaoc,startaob
  424       format("#MmMbfile90 #STARTAOC= ",e14.7/
     &             "#MmMbfile90 #STARTAOB= ",e14.7)
         endif
         if(iaoc.eq.maxmmdatarow .and. iaob .eq. maxmmdatacol)then
            endaoc=aoc
            endaob=aob
            write(0,425)endaoc,endaob
            write(6,425)endaoc,endaob
  425       format("#MmMbfile90 #ENDAOC= ",e14.7/
     &             "#MmMbfile90 #ENDAOB= ",e14.7)
         endif
         go to 400
      endif

  430    continue
C        All comment lines are just written out, or deleted
         write(6,431)(inpone(i2),i2=1,80)
  431    format("#MmMbfile90 ",80a1)
         go to 400

  450    continue
C        All data is in. Close file
         close(unit=10)
        if(iaoc .ne. maxmmdatarow .or. iaob .ne. maxmmdatacol)then
           write(0,453)iaoc,maxmmdatarow,iaob,maxmmdatacol
           write(6,453)iaoc,maxmmdatarow,iaob,maxmmdatacol
 453       format("#WARNING! WARNING! WARNING!  Matrix not full"/
     &            "#iaoc= ",i4,"maxmmdatarow= ",i4/
     &            "#iaob= ",i4,"maxmmdatacol= ",i4/
     &            "# RESETTING maxmmdatarow and maxmmdatacol"/)
           maxmmdatarow=iaoc
           maxmmdatacol=iaob
           endaoc=aoc
           endaob=aob
           write(0,425)endaoc,endaob
           write(6,425)endaoc,endaob
         endif

         deltaaoc=(endaoc-startaoc)/(maxmmdatarow-1)
         deltaaob=(endaob-startaob)/(maxmmdatacol-1)
         write(0,451)deltaaoc,deltaaob
         write(6,451)deltaaoc,deltaaob
  451    format("#MmMbfile90 #DELTAAOC= ",e14.7/
     &          "#MmMbfile90 #DELTAAOB= ",e14.7)
         write(6,452)
         write(0,452)
  452    format("#MmMbfile90: #MmMb90 data input completed. ")
         if(ndatalines .eq. 0)then   ! ? no data lines ?
           write(0,455)
           write(6,455)
  455      format("ERROR: No data lines in file t_MmMb_Surflaw_90a"/
     &            "Stopping")
           iret=455
           return
         endif

      iret=0  !no errors
      return
      end

C=====================================================================
