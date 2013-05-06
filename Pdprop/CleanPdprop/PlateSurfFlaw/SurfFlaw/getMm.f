C getMm.f   vers 0.2  Surface flaw Membrane    July 22 2012
      SAVE
C Create the table of Mm values for a surface flaw as per BS7910 pg.186
C "Mm"  is the membrane (axial) stress intensity multiplier
C compile:   gfortran -g -w  getMm.f  -o getMm
C Usage:    getMm theta >tableMmSurfFlaw
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


      real*4 Mm(51,11)
      character*80 argv,fname
      integer*4  iargc,argc

      nrows=51
      ncols=11

      pi=3.1415927
      write(6,10)nrows,ncols
   10 format("# getMm vers. 0.2 starts..."/
     &   "#ROWS= ",i4/"#COLS= ",i4)
      argc= iargc()
      if(argc .ne. 1)then
        write(0,*)"#Error: Angle argument needed"
        write(6,*)"#Error: Angle argument needed:"
        write(0,*)"#e.g.:  getMn  45  > outfile"
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
   25 format("# iaoc  iaob   a/c   a/B     Mm   phi")

C     Loop for the a/2c variable.
      ao2cinterval=0.020
      ao2c=0.0
C     ao2c gets incremented at end of do loop
      do 2000 iao2c=1,51
      aoc= ao2c*2.0
      coa=1.0/aoc

C     Loop for the a/B  variable.
C     Distance between a/B points is:
      aobinterval= 0.1
      aob=0.0
      do 1000 iaob=1,11
C      aob will be incremented at end of loop

C     Ok, generate all the parameters
      if(ao2c .ge. 0. .and. ao2c .le. 0.5)then
        xM1= 1.13-0.09*aoc
        xM2= (0.89 / (0.2+aoc) )-0.54
        xM3= 0.50 -1.0 / (0.65+aoc) + 14* (1.0-aoc)**24

        g= 1.0 + (0.1+0.35*(aob**2) ) *(1-sintheta )**2
        ftheta= ( aoc*aoc * cos2theta + sin2theta)**(0.25)
        phi= sqrt( 1.0 + 1.464* (aoc**(1.65)) )

        xMm= ( (xM1 + (xM2*aob*aob) + xM3*aob**4) * g *ftheta )/phi
        Mm(iao2c,iaob)= xMm

CpointPlot
C        write(6,400)iao2c,iaob,ao2c,aob,xMm,phi
        write(6,400)iao2c,iaob,aoc,aob,xMm,phi
 400   format(i3,i3,f6.3,1x,f5.2,1x, e14.7,1x,e14.7)
        go to 950
      endif

      if(ao2c .gt. 0.50 .and. ao2c .le. 1.0) then
        xM1= sqrt(coa) * (1 + 0.04*coa)
        xM2= 0.2*(coa**4)
        xM3= -0.11 * (coa**4)

        g= 1.0 + (0.1 + 0.35*coa*(aob**2) ) * ( (1.0-sintheta)**2)
        ftheta= ( ( coa**2)*sin2theta + cos2theta)**(0.25)
        phi= sqrt( 1.0 + 1.464* (coa**(1.65)) )

        xMm= ( (xM1 + (xM2*aob*aob) + xM3*aob**4) * g *ftheta )/phi
        Mm(iao2c,iaob)= xMm

CpointPlot
C        write(6,400)iao2c,iaob,ao2c,aob,xMm,phi
        write(6,400)iao2c,iaob,aoc,aob,xMm,phi
        go to 950
      endif

C     We should not ever end up here
  450  write(0,*)"# Error near Sta.No. 450 in prog getMm...Stopping."
       stop

C     End of loop for ao2c
  950 continue
      aob=aob+aobinterval
 1000 continue

C     One row is done, write it out for the table
Ctable      write(6,470)ao2c,(Mm(iao2c,j),j=1,11)
Ctable  470 format(f4.2,2x,11(e13.6,1x))
C     Put out a blank line for better gnuplot
      write(6,*)
      ao2c=ao2c+ao2cinterval
 2000 continue


      write(0,*)"# Finished theta = ",theta," Radians"

      stop
      end
