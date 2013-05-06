c  pdrain.f  v2.1   prog to rainflow count a column of data. Apr.12 2013
      SAVE
c Compile:  gfortran -g -Wuninitialized  pdrain.f  -o pdrain
c Usage:     pdrain ichan  <infile >outfile
C                   ichan = which column of data to process
c
c  Input file is read free form, output is SAE standard rainflow file.
C  Data points are read into memory, (larger file requires a recompile).
C  Points are scanned once for max, min then rainflowed starting at abs. max.
C  
C---------------------------------------------------------------------------
C  Copyright (C) 2013 A.Conle
C  This program is free software; you can redistribute it and/or
C  modify it under the terms of the GNU General Public License as
C  published by the Free Software Foundation; either version 2 of the
C  license, or (at your option) any later version.
C  This program is distributed in the hope tha it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General PUblic License for more details.
C  You should have received a copy of the GNU General PUblic License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA
C  Try also their web site: http://www.gnu.org/copyleft/gpl.html
C---------------------------------------------------------------------------

C  vers. 0.9 forked from hilo2.f v2.1 and pdprop.f v0.2
C  hilo2.f is an opensource program from FatigueDes+Eval. Comm.
C  pdprop.f is an opensource prog. from A.Conle

C Input file head example:
C  Stage, Time
C  Point, X, Y, Z, dX, dY, dZ
C  Stage 1, -1.495000
C  1000 -1897.499 -651.101 459.407 -0.1715 0.0422 -0.1059
C  1001 -477.977 -449.566 264.794 -0.0505 0.2622 -0.1382
C  1002 -1647.488 -398.552 256.436 -0.0166 0.4729 -0.0119
C  #other comment
C  1003 1000.750 -424.300 265.610 -0.0650 0.4323 0.0243
C  1004 1919.852 -409.564 467.080 -0.0231 0.0216 -0.2006
C  ....etc
C
C  Thus skip lines that start with any non-numeric character in 1st field.
C
      logical debug

      real*4 xvalue(2000000)   ! storage for data
      real*4 mcount(64,64)  ! rainflow count matrix

      character*200  jnp200
      character*5    jnp5(40)
      character*1    jnp1(200)
      equivalence    (jnp200, jnp5(1), jnp1(1))

      CHARACTER*80 INP80, argv, field1
      CHARACTER*5  INP5(16)
      CHARACTER*1  INP1(80)
      EQUIVALENCE (INP80,INP5(1),INP1(1))

C     Right now only one of these will be used:
      REAL*4     xhigh(100), xlow(100), xdata(100)
      integer*4  nhigh(100), nlow(100)

C     These are the push-down list members:
      real*4 tliml(20000), climl(20000)
      real*4 ldo,lobj,ldn,ldnn,xmean
    
      INTEGER UNO,UNI
      uni=5  ! input 
      UNO=6  ! output

C
      debug= .false.
C      debug=.true.
      nbins=64
      nlistmax=100
      maxpoints= 2000000
      maxpd=20000
C     if you change the dimensions of xhigh, xlow,etc above, change above too.

C     Test if -Wuninitilized  warning works in gfortran:
C      i=ifix(x)
C     Works!

      write(0,*) "#pdrain vers. 0.9 Starts. Usage: ",
     &" pdrain ichan <in  >out"

C
      nargc = iargc()
      if(nargc .ne. 1)then
           write(0,*) "#Error: Incorrect args list. Should be: ",
     &     " pdrain ichan <in  >out"
           stop
      endif
      jvect=1
      call getarg(jvect,argv)
      ncol=0
      read(argv,*,err=15)ncol
      write(i0,10)ncol
      write(UNO,10)ncol
   10 format("# Rainflow count will be on Data column= ",i4)
      go to 16
C
   15 write(0,*)"#Error: Could not read 1st arg: column no., ichan"
      write(0,*) "# Usage should be:  pdrain ichan <in  >out"
      stop
C
   16 continue
      NL=0          !counter for line no. in file
      ndatalines=0

C     Find the first data line, and initilize max min registers
  110 continue
      read(UNI,111,end=112)jnp200
  111 FORMAT(A200)
      goto 113

  112 write(0,*)"#Error: End of File before 1st data line read." 
      stop

  113 continue
      NL=NL+1
      if(jnp200.eq." ")go to  110

c     Find the first non blank character.
      n=0
      do  115 i=1,200
           if(jnp1(i).eq." ")go to 115
c          it is a non blank char
           n=i
           go to 120
  115 continue

  120    continue
c        If 1st char is a number or + or - or .  then assume it is a point number
c        otherwise it is a character word and of no interest
         if(jnp1(n).eq."." .or.
     &      jnp1(n).eq."-" .or.
     &      jnp1(n).eq."+" .or.
     &      jnp1(n).eq."0" .or.
     &      jnp1(n).eq."1" .or.
     &      jnp1(n).eq."2" .or.
     &      jnp1(n).eq."3" .or.
     &      jnp1(n).eq."4" .or.
     &      jnp1(n).eq."5" .or.
     &      jnp1(n).eq."6" .or.
     &      jnp1(n).eq."7" .or.
     &      jnp1(n).eq."8" .or.
     &      jnp1(n).eq."9" )then
C           yes it is a number

C     Allow nlistmax columns max.  This is the first data line.  
C          We only want data column  "ncol"
        ndatalines=ndatalines+1

C       Repeatedly read the same line, one exta arg each time
C       when we run out or args we know how many there are
        read(jnp200,*,err=130,end=130)(xdata(i),i=1,ncol)
        write(0,108)ncol
        write(6,108)ncol
  108   format("# First data line found. Read ",i4," columns.")
        go to 135
C
C       If we reach this point there is a comment before column "ncol"
  130   write(0, 132)jnp200
  132   format('#Error: in 1st data line: '/a200)
        write(0,*)"#At file line number= ",NL
        write(0,*)"# Unable to read ",ncol," columns of data. Stopping"
        stop
C
  135   continue
C       We now have read the 1st data line. It had at least ncol columns.
C       For all the rest of the data lines assume ncol  is valid.
C       Initialize the highs and lows registers:
        do  137 i=1,ncol
          xhigh(i)=xdata(i)
          nhigh(i)=ndatalines
          xlow(i)=xdata(i)
          nlow(i)=ndatalines
  137   continue
        xvalue(ndatalines)= xdata(ncol)  !save for rainflow count
        go to  140
      endif
c     We have not found the first data line yet. Keep looking
      go to  110




C     First data line is in, Look for rest of data.
  140 continue
c        Fetch a new line. It may be comment or data.
         read(UNI,111,end=190)jnp200
         NL=NL+1
Cdebug         write(0,*)"Line ",NL,jnp200
         if(jnp200.eq." ")go to 140

c        Find the first non blank character.
         n=0
         do  145 i=1,200
           if(jnp1(i).eq." ")go to 145
c          it is a non blank char
           n=i
           go to 146
  145    continue

  146    continue
c        If 1st char is a number or + or -, then assume it is a point number
c        otherwise it is a character word and of no interest
         if(jnp1(n).eq."." .or.
     &      jnp1(n).eq."-" .or.
     &      jnp1(n).eq."+" .or.
     &      jnp1(n).eq."0" .or.
     &      jnp1(n).eq."1" .or.
     &      jnp1(n).eq."2" .or.
     &      jnp1(n).eq."3" .or.
     &      jnp1(n).eq."4" .or.
     &      jnp1(n).eq."5" .or.
     &      jnp1(n).eq."6" .or.
     &      jnp1(n).eq."7" .or.
     &      jnp1(n).eq."8" .or.
     &      jnp1(n).eq."9" )then
C           yes it is a number
C           Read in the whole line of ncol numbers
            read(jnp200,*,err= 170)(xdata(i),i=1,ncol)
            ndatalines=ndatalines+1
            if(ndatalines .gt. maxpoints)then
              write(0,*)"#Error: Too many data points= ",ndatalines
              write(0,*)"# You need to re-compile with bigger storage",
     &                  " for ""xvalue()"" in pdrain.f. Stopping."
              stop
            endif
            xvalue(ndatalines)= xdata(ncol)

            do  160 i=1,ncol
             if( xhigh(i).lt.xdata(i) )then
                xhigh(i)=xdata(i)
                nhigh(i)=ndatalines
             endif
             if( xlow(i).gt.xdata(i) )then
                xlow(i)=xdata(i)
                nlow(i)=ndatalines
             endif
  160       continue
         endif
c        if we get to here it is not a data line, or we are done
c        with this data line. Go fetch the next line
         go to  140



  170   continue
        write(0, 172)NL,jnp200
        write(6, 172)NL,jnp200
  172   format("#ERROR: non-number in file line no=",i6,
     &  " : "/a200/"#***Assuming line is a comment...")
C       Keep going
        go to  140


C      Done with all data input.
  190 continue
       IER=0
       write(0,*)"#Done max,min hunt. Scanned ",ndatalines,
     &           " data lines."
       write(0,*)"#Total lines in file= ",NL
       write(0,191)xhigh(ncol),nhigh(ncol)
       write(0,192)xlow(ncol),nlow(ncol)

       write(6,191)xhigh(ncol),nhigh(ncol)
  191  format("#HIGH= ",e14.7," inLine= ",i9)
       write(6,192)xlow(ncol),nlow(ncol)
  192  format("#LOW=  ",e14.7," inLine= ",i9)

       iupdown=+1   ! assume up half-cycle or reversal
       absmax=abs(xhigh(ncol) )
       nstart=nhigh(ncol)
       if(absmax.lt. (abs(xlow(ncol))) )then
C         the -ve side is bigger than the +ve
          nstart=nlow(ncol)
          iupdown=-1
       endif
       write(0,193)nstart,xvalue(nstart)
       write(6,193)nstart,xvalue(nstart)
  193  format("# Counting will start at DataPt no.= "/
     &        "#NSTART= ",i9/
     &        "#STARTVALUE= ",e14.7)
          
C      Check for silly user input
       if((xhigh(ncol) .eq. 0.0) .and. (xlow(ncol) .eq. 0.0 ))then
         write(0,195)
         write(6,195)
  195    format("# Error: Your Max. = Min. = 0.0  !"/
     &   "# There is no history to cycle count.  Stopping.")
         stop
       endif
C      If the load history is only a single non-zero point, things should
C      still work out ok, since we would go from 0,0  to the max, 
C      back to the 0,0 origin and then back to the max. Thus creating
C      a single large cycle.
       if(ndatalines .eq. 0)then
         write(0,196)ndatalines
         write(6,196)ndatalines
  196    format("# Error: No Data?  ndatalines= ",i9,"  Stopping.")
       endif

       if(uno.ne.6)close(unit=uno)


  200 continue
C  pdprop.f   vers. 0.2  Notched Spec Crack Prop.  FAC apr.22 2012
C  Push-Down List crack propagation program.
C  Explanation of primary variables used=
C    tliml, climl = 1 dimensional arrays of vectors containing
C       the push-ddown list loads.  TLIML is the tensile, CLIML is compr.
C  nptt,nptc =  push-down list pointers.
C  lorg   The load origin of a given 1/2 cycle
C  lobj        "           "      end point        "    "
C  dld   the change in "    "  during a 1/2 cycle
C  maxpd = maximum no of points possible in PD list
C  maxpoints = max no. of pts in Load history

C      Initilize some variables.
C      Origin and zero
      nrev=0  ! count half cycles
      nptc=0
      nptt=0
      do 20 i=1,maxpd
        tliml(i)=0.
        climl(i)=0.
   20 continue
      do 25 i=1,nbins
      do 25 j=1,nbins
      mcount(i,j)=0.
   25 continue

C     Rainflow counts are saved in matrix mcount()
C     vmax and vmin are the outer boundary limits of the matrix.
C     They are slightly (0.5 bin) above and below xhigh() and xlow().
C     For a given counted closed "loop" the maximum tip is saved
C     in column jmax,  while the minimum tip is saved in row, imin.
C     jmax and imin are computed by the distance between the value
C     and the vmin matrix limit.  
C     When closed loops are very small, less than 1 bin in range, 
C     jmax and imin will be equal.  These are saved but not printed
C     out in final table. If you wish to activate them you will need
C     to increase the number of matrix bins, or make some assumption about
C     their size, and change the output code.
      deltabin= (xhigh(ncol)-xlow(ncol))/(nbins-1)
      vmax=xhigh(ncol)+deltabin/2.0
      vmin=xlow(ncol) -deltabin/2.0
      write(0,26)vmax,vmin,nbins
   26 format("#Matrix Max= ",e14.7," Min= ",e14.7," nbins= ",i3)

      npoint=nstart ! pointer to the xvalue(npoint) being examined.


c
c---------------------------------------------------------------------


C  This push-down list rainflow counter is based on material memory models
C  such as described in: Conle, A., T.R.Oxland, T.H.Topper, "Computer-Based 
C  Prediction of Cyclic Deformation and Fatigue Behavior," Low Cycle Fatigue 
C  ASTM STP 942, 1988, pp.1218-1236. 
C  and in multi dimension models:
C  Chu, C.-C., "A Three-Dimensional Model of Anisotropic Hardening in Metals
C  and Its Application to the Analysis of Sheet Metal Formability," 
C  J.of Mech. Phys. of Solids, Vol.32, 1984, pp.197-212.

C  The plan is to:
C  Get values from xvalue().     We will start at the largest excursion
C  defined by nstart  i.e.:   xvalue(nstart)   We will then look at each
C  point that follows until we get to the last one xvalue(ndatalines)
C  at which time we will return back to xvalue(1)  and then read successive
C  points up to and including xvalue(nstart).  The last point at nstart 
C  will clear out the push-down list and thus close all hysteresis loops.


      lobj=0.

C     Return to this point after each half cycle is completed.
 3000 continue
      ldo=lobj
      nrev=nrev+1
      if(nrev .eq. 1)then
        iphase=1
C       This is the very first half cycle. Located at nstart.  It should
C       be the largest excursion from 0,0 origin  and it is definitely a
C       reversal point,  so we don't have to look ahead to see if the next
C       point is even "bigger".  Note that it may be in either tension or
C       compression.
C       The direction indicator "iupdown" was set previously. It should
C       still be correct for 1st ramp.
        lobj=xvalue(nstart)
        go to 3050
      endif

C     Ok, we have had a reversal. Now we are hunting for the next half cycle's
C     new reversal point.
      iupdown=-iupdown

 3012 continue     !  Loop  to find next reversal value.--------
      npoint=npoint+1
      if(npoint .gt. ndatalines)then
        write(0,*)"#Wrap npoint around to xvalue(1)..."
        iphase=2
        if(debug)write(6,*)"#-----------iphase=2 -------------"
        npoint=1
      endif

      if(iphase .eq. 2 .and. npoint .eq. nstart)then
C        Finish the ramp to the largest excursion xvalue(nstart)
        lobj=xvalue(nstart)
        go to 3050
      endif

      if(iphase .eq. 2 .and. npoint .gt. nstart)then
C       we have hit the end of history,
        goto 5000
      endif
C     Ok, we have the next potential data point
      trialpoint1= xvalue(npoint)

C     Now check on the next to next
      npointplus1=npoint+1
      if(npointplus1 .gt. ndatalines)then
        write(0,*)"#Wrap npointplus1 around to xvalue(1)..."
        npointplus1=1
      endif
      trialpoint2= xvalue(npointplus1)



      if(iupdown .eq. -1)go to 3014

C     We are trying to go in tensile direction. Figure out where the next
C     reversal point is. ----------------------------------------------------
C     iupdown = +1
      if(trialpoint1 .eq. ldo)go to 3012  !new point is same as old, skip it.

      if(trialpoint1 .lt. ldo)then
C       New point is not in tensile direction (this should not actually happen)
        write(0,3015)npoint,ldo,iupdown
 3015   format("#Error: logic: near statement 3015"/
     &         "# dataline No.= ",i9/
     &         "# Old Reversal = ",e14.7,", New direction, iupdown= ",
     &         i2,"  Stopping.")
         stop
       endif

C      Ok, the trial point is above the ldo old reversal.  Whats next?
C      It potentially is a reversal peak. Depends on the following point value
       if(trialpoint2 .ge. trialpoint1)then
C          The 2nd point is above the first, so ignore the first
C          and keep looking
           go to 3012
       endif
C      If we got to here then
C      trialpoint2 indicates that trialpoint 1 is a valid reversal point,
C      so go ahead and run it through the counter stuff.
C      The new target is trialpoint1
       lobj=trialpoint1
       go to 3050

C     Going Down:  iupdown = -1  and we are hunting for the next compressive
C     peak. ---------------------------------------------------------------
 3014 continue
      if(trialpoint1 .eq. ldo)then
C        Next point is same as where we are now,  skip it
         go to 3014
      endif

      if(trialpoint1 .gt. ldo)then
C        New point is not in compressive direction, this should not happen
         write(0,3015)npoint,ldo,iupdown
 3016    format("#Error: logic: near statement 3016"/
     &          "# dataline No.= ",i9/
     &          "# Old Reversal = ",e14.7,", New direction, iupdown= ",
     &          i2,"  Stopping.")
          stop
       endif

C      Ok, trial pt is below ldo (old reversal), it is in the right direction.
C      It potentially is a reversal peak. Depends on the following point value
       if(trialpoint2 .le. trialpoint1 )then
C         2nd new point is equal to or below the 1st new point, thus trialpoint1
C         is not a reversal.  Keep looking.
          go to 3012
       endif

C      If we get to here then
C      trialpoint2 indicates that trialpoint1 is our desired new compressive peak,
C      go ahead and run it through the rainflow logic.
       lobj=trialpoint1
       go to 3050


 3050 continue
      if(debug)write(6,*)"#lobj= ",lobj,"  ldo= ",ldo
      dld=abs(lobj-ldo)
C ----------------------------------------------------------------

C     Start Cycling
C     Going up or down ?
  500 continue
      if(lobj .lt. ldo) go to 2100

C************ Going UP,  Tensile Direction **************************
C     Are we on the monotonic curve? (on very 1st half cycle for sure)
 1100 if(nptt .eq. 0) go to 1500
c     No. Has a exceedence  occurred?
 1110 if(lobj .gt. tliml(nptt)) go to 1120

c     Is this Same as previous amplitude?
 1130 if(lobj .eq. tliml(nptt) .and. nptt .ne. 1) go to 1135
      go to 1350

c     Is a loop being closed without the monotonic?
 1120 if(nptt .ne. 1) go to 1400

c     Is closure in connection  with the monotonic?
 1150 if(nptc .eq. 2) go to 1170

c     Check for impossible
 1160 if(nptc .eq. 1) go to 1165
      idump=1160
      go to 9999

C     Same load level as previous level.
 1135 continue
C 1135 dam=SMITH( (lobj-ldo),tlsts(nptt),clsts(nptc),tlstr(nptt),
C     &            clstr(nptc),totdam,nrev)
C     Count this closed half cycle. The other side has already been counted
      jmax= ifix((lobj-vmin)/deltabin)+1
      imin= ifix((ldo-vmin)/deltabin)+1
      mcount(imin,jmax)=mcount(imin,jmax)+0.5
      if(debug)write(6,*)"#1135 same ld, jmax,imin= ",jmax,imin
      nptc=nptc-1
c     The closed loop is erased and the point of rev. is already there
      go to 3000


C     An unmatched 1/2 cycle is being forgotten and a return to the 
c     monotonic curve is occurring. Count the unmatched 1/2 cycle.
 1165 continue
C 1165 dam=SMITH ( (tliml(1)-climl(1)), tlsts(1), clsts(1),
C     &         tlstr(1),clstr(1),totdam,nrev)
      jmax= ifix((tliml(1)-vmin)/deltabin)+1
      imin= ifix((climl(1)-vmin)/deltabin)+1
      mcount(imin,jmax)=mcount(imin,jmax)+0.5 !this should not happen
C     because we started at the largest excursion.
      if(debug)write(6,*)"#1165 same ld, jmax,imin= ",jmax,imin
      write(0,*)"#Error: unmatched 1/2 cycle,  npoint= ",npoint
      write(6,*)"#Error: unmatched 1/2 cycle,  npoint= ",npoint
      nptc=0
      nptt=0
      go to 1500

C     A loop is being closed and a return to the monotonic 
c     curve is occurring.  Damage is the same as the other half of
c     the cycle being closed.
 1170 continue
C 1170 totdam=totdam+SMITH( (tliml(1)-climl(2)),tlsts(1),clsts(2),
C     &    tlstr(1),clstr(2), totdam,nrev)
      jmax= ifix((tliml(1)-vmin)/deltabin)+1
      imin= ifix((climl(2)-vmin)/deltabin)+1
      mcount(imin,jmax)=mcount(imin,jmax)+0.5
      if(debug)write(6,*)"#1170 : jmax,imin= ",jmax,imin
      write(0,*)"#Error: return to monotonic,  npoint= ",npoint
      write(6,*)"#Error: return to monotonic,  npoint= ",npoint
c     Eliminate the closed loop.
      nptc=0
      nptt=0
      go to 1500

C     A new entry in the P.D. list is being made.
 1350 continue
C 1350 dld=abs(lobj-climl(nptc))
C      dam=SMITH( dld,sobj,so,eobj,eo,totdam,nrev)
      jmax= ifix((lobj-vmin)/deltabin)+1
      imin= ifix((climl(nptc)-vmin)/deltabin)+1
      mcount(imin,jmax)=mcount(imin,jmax)+0.5
      if(debug)write(6,*)"#1350 , jmax,imin= ",jmax,imin
c     Enter the rev in the P.D. list
      nptt=nptt+1
      tliml(nptt)=lobj
      go to 3000

C     Deformation is occurring on the monotonic curve again.
 1500 continue
C 1500 dld=abs(lobj)
c     Subtract damage of previous use of monotonic curve (=0 in cyc )
      go to 600

C     A loop is being closed, count the remaining half of the
c     loop and then eliminate the loop.
 1400 continue
C 1400 totdam=totdam+SMITH( (tliml(nptt)-climl(nptc)), tlsts(nptt),
C     &      clsts(nptc),tlstr(nptt),clstr(nptc),totdam,nrev)
      jmax= ifix((tliml(nptt)-vmin)/deltabin)+1
      imin= ifix((climl(nptc)-vmin)/deltabin)+1
      mcount(imin,jmax)=mcount(imin,jmax)+0.5
      if(debug)write(6,*)"#1400 : jmax,imin= ",jmax,imin
      nptc=nptc-1
C     Subtract the old half cycle damage of the Re-incurred
c     stress-strain path.
C     totdam=totdam-tldam(nptt)
C     Use the same   jmax   as the one above
      imin= ifix((climl(nptc)-vmin)/deltabin)+1
      mcount(imin,jmax)=mcount(imin,jmax)-0.5    !subtract !
      if(debug)write(6,*)"#1400 subtract: jmax,imin= ",jmax,imin
      nptt=nptt-1
c     Reset the local stress-strain origin to the old origin.
      ldo=climl(nptc)
      go to 1100

C*****  Going Down, Compressive Direction ****************************

C     Are we on the monotonic curve?
 2100 if(nptc .eq. 0) go to 2500

c     No. Has a exceedence  occurred?
 2110 if(lobj .lt. climl(nptc)) go to 2120

c     Is this Same as previous amplitude?
 2130 if(lobj .eq. climl(nptc) .and. nptc .ne. 1) go to 2135
      go to 2350

c     Is a loop being closed without the monotonic?
 2120 if(nptc .ne. 1) go to 2400

c     Is closure in connection  with the monotonic?
 2150 if(nptt .eq. 2) go to 2170

c     Check for impossible
 2160 if(nptt .eq. 1) go to 2165
      idump=2160
      go to 9999

C     Same load level as previous level.
 2135 continue
C 2135 dam=SMITH( (ldo-lobj),tlsts(nptt),clsts(nptc),
C     &      tlstr(nptt),clstr(nptc),totdam,nrev)
      jmax= ifix((ldo-vmin)/deltabin)+1
      imin= ifix((lobj-vmin)/deltabin)+1
      mcount(imin,jmax)=mcount(imin,jmax)+0.5
      if(debug)write(6,*)"#2135 : jmax,imin= ",jmax,imin
c     The closed loop is erased and the point of rev. is already there
      nptt=nptt-1
      go to 3000

C     An unmatched 1/2 cycle is being forgotten and a return to the
c     monotonic curve is occurring. Count the unmatched 1/2 cycle.
 2165 continue
C 2165 dam=SMITH( (tliml(1)-climl(1)),tlsts(1),clsts(1),tlstr(1),
C     &        clstr(1),totdam,nrev)
      jmax= ifix((tliml(1)-vmin)/deltabin)+1
      imin= ifix((climl(1)-vmin)/deltabin)+1
      mcount(imin,jmax)=mcount(imin,jmax)+0.5  ! this should not happen
      if(debug)write(6,*)"#2165 : jmax,imin= ",jmax,imin
C       because we started with the max excursion
      write(0,*)"#Error: unmatched 1/2 cycle,  npoint= ",npoint
      write(6,*)"#Error: unmatched 1/2 cycle,  npoint= ",npoint
      nptc=0
      nptt=0
      go to 2500


C     A loop is being closed and a return to the monotonic
c     curve is occurring.  Damage is the same as the other half of
c     the cycle being closed.
 2170 continue
C 2170 totdam=totdam+SMITH( (tliml(2)-climl(1)),tlsts(2),
C     &       clsts(1),tlstr(2),clstr(1),totdam,nrev)
      jmax= ifix((tliml(2)-vmin)/deltabin)+1
      imin= ifix((climl(1)-vmin)/deltabin)+1
      mcount(imin,jmax)=mcount(imin,jmax)+0.5
      if(debug)write(6,*)"#2170 : jmax,imin= ",jmax,imin
      nptc=0
      nptt=0
      go to 2500

C     A new entry in the P.D. list is being made.
 2350 continue
C 2350 dld=abs(lobj-tliml(nptt))
C      dam=SMITH(dld,so,sobj,eo,eobj,totdam,nrev)
      jmax= ifix((tliml(nptt)-vmin)/deltabin)+1
      imin= ifix((lobj-vmin)/deltabin)+1
      mcount(imin,jmax)=mcount(imin,jmax)+0.5
      if(debug)write(6,*)"#2350 : jmax,imin= ",jmax,imin
      nptc=nptc+1
      climl(nptc)=lobj
      go to 3000

C     Deformation is occurring on the monotonic curve again.
 2500 continue
C 2500 dld=abs(lobj)
      write(0,*)"#Error: we are back on monotonic,  npoint= ",npoint
      write(6,*)"#Error: we are back on monotonic,  npoint= ",npoint
c     Subtract damage of previous use of monotonic curve (=0 in cyc )
C      totdam=totdam - tldam(1)


C     Add the present damage. Both Tens and Comp comes here.
C     Note that the monotonic does not get counted.
  600 continue
C  600 dsx=abs(sobj)
C      dex=abs(eobj)
C     dsx and dex are amplitudes
C      dsamp=dsx
C      deamp=dex
C      dam= SMITH(dld,dsamp,-dsamp,deamp,-deamp,totdam,nrev)
C      jmax= ifix((lobj-vmin)/deltabin)+1
C      imin= ifix((0.0-vmin)/deltabin)+1
C      mcount(imin,jmax)=mcount(imin,jmax)+0.5  ! this should probably not be used

c     Set the appropriate push-down list arrays.
      climl(1)=-abs(lobj)
      tliml(1)=abs(lobj)
      nptt=1
      nptc=1
      go to 3000

C     A loop is being closed, count the remaining half of the
c     loop and then eliminate the loop.
 2400 continue
C 2400 totdam=totdam+SMITH( (tliml(nptt)-climl(nptc)), tlsts(nptt),
C     &       clsts(nptc),tlstr(nptt),clstr(nptc),totdam,nrev)
      jmax= ifix((tliml(nptt)-vmin)/deltabin)+1
      imin= ifix((climl(nptc)-vmin)/deltabin)+1
      mcount(imin,jmax)=mcount(imin,jmax)+0.5
      if(debug)write(6,*)"#2400 : jmax,imin= ",jmax,imin
      nptt=nptt-1
C     Subtract the old half cycle damage of the Re-incurred
c     stress-strain path.
C       totdam=totdam-cldam90(nptc)
      jmax= ifix((tliml(nptt)-vmin)/deltabin)+1
C     Use the same   imin    from above
      mcount(imin,jmax)=mcount(imin,jmax)-0.5     !   subtract  !!
      if(debug)write(6,*)"#2400 + subtract: jmax,imin= ",jmax,imin
      nptc=nptc-1
c     Reset the local stress-strain origin.
      ldo=tliml(nptt)
      go to 2100

 5000 continue
C     End of data input, empty the left-overs in push-down list

C     This code is not necessary,  these items have already been counted.
C 5010 continue
CC     Take one out from tliml() and one from climl() until  we
CC     get back to the monotonic of tliml(1) or climl(1).
CC     nptc and nptt are the pointers to the entries in tliml and climl
C      if(nptc.eq.0 .or. nptt .eq. 0) go to 5020 ! done
CC       not done yet
C        jmax= ifix((tliml(nptt)-vmin)/deltabin)+1
C        imin= ifix((climl(nptc)-vmin)/deltabin)+1
C        mcount(imin,jmax)=mcount(imin,jmax)+1
C        nptt=nptt-1
C        nptc=nptc-1
C        go to 5010

C     End of data.  Put out the accumulated Rainflow Matrix.
 5020 continue
C     The reference point is at the center of the lowest bin:
      xminbin= vmin+deltabin/2.0
      do 5050 i=1,nbins
      do 5045 j=1,nbins
      xcount=mcount(i,j)
      if(xcount .eq. 0.0)go to 5045 ! nothing found
      top=((j-1)*deltabin)+xminbin
      bot=((i-1)*deltabin)+xminbin
      vrange=top-bot
      vmean=(top+bot)/2.0
      if(i.eq.j)go to 5045
      write(6,5040)vrange,vmean,xcount,top,bot
 5040 format(2x,g10.3,2x,g10.3,2x, f10.1 ,2x,g10.3,2x,g10.3)
 5045 continue
 5050 continue
C     write out end of rainflow  numbers
      xcount=0.
      vrange=0.
      vmean=0.
      top=0.
      bot=0.
      write(6,5040)vrange,vmean,xcount,top,bot
      stop



 9999 write(6,100) idump
  100 format("# *** Error: near statement label no. ",i5/)
      stop
      end


