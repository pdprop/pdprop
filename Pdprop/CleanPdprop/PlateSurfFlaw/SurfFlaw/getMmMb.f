C getMmMb.f   vers 0.4  Surface flaw Membrane    July 27 2012
      SAVE
C Create the table of Mb values for a surface flaw as per BS7910 pg.186
C "Mm"  is the membrane (axial) stress intensity multiplier
C "Mb"  is then calculated from Mm. It is the Plate in bending multiplier.
C compile:   gfortran -g -w  getMmMb.f  -o getMmMb
C Usage:    getMmMb theta >tableMmSurfFlaw
C                 theta=  angle in degrees

C  Copyright (C) 2012  F.A.Conle
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

C vers. 0.4 Correction for 2nd G1  calculation. Jul27-12
C vers. 0.3 Puts out both Mm  and Mb  Jul27-12
C vers. 0.2 Change output to tabulate traditional a/c rather than a/2c Jul27-12
C There are two variables that define the matrix of Mm
C   1. ao2c :  a/2c    where a= crack depth,  2c= crack width at surface
C   2. aob  :  a/B     where a is crack depth, B= plate thickness

C matrices are created by this program for the theta= angle from surface. e.g.:
C   theta=  0    where the crack meets the surface
C   theta=  pi/2    90 deg. from surface (where crack is deepest)

C Note that the output file is not in table form, but rather one Mm value per line
C If you want the table form activate the write statment at lines that
C   begin with "Ctable" 


      real*4 Mm(51,11), Mb(51,11)
      character*80 argv,fname
      integer*4  iargc,argc

      nrows=51
      ncols=11

      pi=3.1415927
      write(6,10)nrows,ncols
   10 format("# getMmMb vers. 0.4 starts..."/
     &   "#ROWS= ",i4/"#COLS= ",i4)
      argc= iargc()
      if(argc .ne. 1)then
        write(0,*)"#Error: Angle argument needed"
        write(6,*)"#Error: Angle argument needed:"
        write(0,*)"#e.g.:  getMmMb  45  > outfile"
        write(0,*)"#Stopping now..."
        stop
      endif
      jvect=1
      call getarg(jvect,argv)
      read(argv,*,err=30)angleDeg
      if(angleDeg .lt.0.0 .or. angleDeg .gt. 180)then
        write(0,*)"Error: Angle must  be between 0 and 180 deg."
        write(0,*)"Error: Angle must  be between 0 and 180 deg."
        stop
      endif
      go to 35

   30 continue
      write(0,*)"#Error: reading angle:",argv
      write(6,*)"#Error: reading angle:",argv
      stop


   35 continue
C     convert to radians:
      theta= (angleDeg/360.)*pi*2.0
      sintheta= sin(theta)
      costheta= cos(theta)
      sin2theta= sintheta*sintheta
      cos2theta= costheta*costheta

      write(0,20)theta,angleDeg
      write(6,20)theta,angleDeg
   20 format("#theta= ", e14.7,"   #radians = ",f6.2," deg.")

      write(6,25)
C   25 format("# iao2c  iaob   a/2c   a/B     Mm   phi")
   25 format("#iaoc iaob a/c a/B         Mm          Mb")

C     Loop for the a/2c variable.
      ao2cinterval=0.020
      ao2c=0.0001
C     ao2c gets incremented at end of do loop
      do 2000 iaoc=1,51    ! ----------------------loop  a/c----------------
      aoc= ao2c*2.0
      coa=1.0/aoc

C     Loop for the a/B  variable.
C     Distance between a/B points is:
      aobinterval= 0.1
      aob=0.0
      do 1000 iaob=1,11    ! ----------------------loop  a/b----------------
C      aob will be incremented at end of loop

C     Ok, generate all the parameters

C ------------- Case A:    0 .ge. a/c .le. 1.0  ---------------------------------
C               note that ao2c = aoc/2
      if(ao2c .ge. 0. .and. ao2c .le. 0.5)then
        xM1= 1.13-0.09*aoc
        xM2= (0.89 / (0.2+aoc) )-0.54
        xM3= 0.50 -1.0 / (0.65+aoc) + 14* (1.0-aoc)**24

        g= 1.0 + (0.1+0.35*(aob**2) ) *(1-sintheta )**2
        ftheta= ( aoc*aoc * cos2theta + sin2theta)**(0.25) ! NOTE: =0 when aoc=0
        phi= sqrt( 1.0 + 1.464* (aoc**(1.65)) )

        xMm= ( (xM1 + (xM2*aob*aob) + xM3*aob**4) * g *ftheta )/phi
        Mm(iaoc,iaob)= xMm

C       Ok, we have Mm,   now get Mb
        q= 0.2 +aoc + 0.6*aob
        H1= 1.0 - 0.34*aob - 0.11*aoc*aob
        G1= -1.22 - 0.12*aoc
        G2= 0.55 - 1.05*(aoc**(0.75)) + 0.47*(aoc**(1.5))
        H2= 1.0 + G1*aob + G2*(aob**2)

        H= H1+ (H2-H1)* ((sintheta)**(q) )
        xMb= H * xMm
        Mb(iaoc,iaob)= xMb

CpointPlot
Cold        write(6,400)iao2c,iaob,ao2c,aob,xMm,phi

        write(6,400)iaoc,iaob,aoc,aob,xMm,xMb
C        write(6,400)iaoc,iaob,aoc,aob,xMm,xMb,xM1,xM2,xM3,g,
C     &              ftheta,phi,q,H1,G1,G2,H2,H
 400   format(i3,i3,f7.4,1x,f5.2,1x, e14.7,1x,e14.7)
C 400   format(i3,i3,f6.3,1x,f5.2,1x, 14(e10.3,1x))
        if(iaoc.eq.1 )then
          write(6,60)iaoc,iaob,aoc,aob,xMm,xMb,xM1,xM2,xM3,g,ftheta,phi
  60   format("#",i3,i3,f7.4,1x,f5.2,1x, e14.7,1x,e14.7, 6(1x,e14.7))
        endif
        go to 950
      endif

C ------------- Case B:   1 .gt. a/c .le. 2.0  ---------------------------------
C      if(ao2c .gt. 0.50 .and. ao2c .le. 1.0) then
      if(ao2c .gt. 0.50 ) then  ! technically it must be <= 1.0, but ignore
C        if(ao2c .gt. 1.0)then
C          ao2c=1.
C          aoc=2.
C          coa=1.0/aoc
C        endif
        xM1= sqrt(coa) * (1 + 0.04*coa)
        xM2= 0.2*(coa**4)
        xM3= -0.11 * (coa**4)

        g= 1.0 + (0.1 + 0.35*coa*(aob**2) ) * ( (1.0-sintheta)**2)
        ftheta= ( ( coa**2)*sin2theta + cos2theta)**(0.25)
        phi= sqrt( 1.0 + 1.464* (coa**(1.65)) )

        xMm= ( (xM1 + (xM2*aob*aob) + xM3*aob**4) * g *ftheta )/phi
        Mm(iaoc,iaob)= xMm

C       Ok, we have Mm,   now get Mb
        q= 0.2 +coa + 0.6*aob
        H1= 1.0 - (0.04 + 0.41*coa) * aob +
     &         (0.55-1.93*(coa**(0.75)) + 1.38*(coa**(1.5)) )*(aob**2)
        G1= -2.11 + 0.77*coa
        G2= 0.55 - 0.72*(coa**(0.75)) + 0.14*(coa**(1.5))
        H2= 1.0 + G1*aob + G2*(aob**2)

        H= H1+ (H2-H1)* ((sintheta)**(q) )
        xMb= H * xMm
        Mb(iaoc,iaob)= xMb

CpointPlot
        write(6,400)iaoc,iaob,aoc,aob,xMm,xMb
C        write(6,400)iaoc,iaob,aoc,aob,xMm,xMb,xM1,xM2,xM3,g,
C     &              ftheta,phi,q,H1,G1,G2,H2,H
        go to 950
      endif

C     We should not ever end up here
  450  write(0,*)"# Error near Sta.No. 450 in prog getMb...Stopping."
       stop

C     End of loop for ao2c
  950 continue
      aob=aob+aobinterval
 1000 continue

C     One row is done, write it out for the table
Ctable      write(6,470)aoc,(Mb(iaoc,j),j=1,11)
Ctable  470 format(f4.2,2x,11(e13.6,1x))

C     Put out a blank line for better gnuplot in non-table output
      write(6,*)
      ao2c=ao2c+ao2cinterval
 2000 continue


      write(0,*)"# Finished theta = ",theta," Radians"

      stop
      end
