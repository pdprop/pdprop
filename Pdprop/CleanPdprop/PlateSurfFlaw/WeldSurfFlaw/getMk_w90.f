C getMk_w90.f   vers 0.1  Surface flaw near Weld Aug 10 2012
      SAVE
C Create the table of Mkm and Mkb values for a surface flaw as per BS7910 pg.224
C "Mkm"  is the membrane (axial) stress intensity multiplier
C "Mkb"  is the Plate in bending multiplier.
C  90  refers to the Deepest location on flaw
C compile:   gfortran -g -w -fbounds-check getMk_w90.f  -o getMk_w90
C Usage:    getMk_w90 LoB > tableMk_w90
C                     LoB is L/B ratio

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

C There are three variables that define the matrix of Mkm and Mkb
C   1. aoc :  a/c    where a= crack depth,  2c= crack width at surface
C   2. aob :  a/B     where a is crack depth, B= plate thickness
C   3. LoB :  L/B    where L is weld joint width, B = flawed plate thickness

C Note that the output file is not in table form, but rather one Mkm 
C      and one Mkb value per line
C If you want the table form activate the write statment at lines that
C   begin with "Ctable" 


      real*4 Mkm(200,18), Mkb(200,18)
      real*4 LoB
      character*80 argv,fname
      integer*4  iargc,argc

      nrows=200
      ncols=18

      pi=3.1415927
      write(6,10)nrows,ncols
   10 format("# getMk_w90 vers. 0.1 starts..."/
     &   "#ROWS= ",i4/"#COLS= ",i4)
      argc= iargc()
      if(argc .ne. 1)then
        write(0,*)"#Error: L/B argument needed"
        write(6,*)"#Error: L/B argument needed:"
        write(0,*)"#e.g.:  getMk_w90  2.5  > outfile"
        write(0,*)"#Stopping now..."
        stop
      endif
      jvect=1
      call getarg(jvect,argv)
      read(argv,*,err=30)LoB
      if(LoB .lt.0.0 )then
        write(0,*)"Error: LoB must  be between >= 0.5 "
        write(6,*)"Error: LoB must  be between >= 0.5 "
        stop
      endif
      if(LoB .gt.2.75)then
        write(0,15)LoB
        write(6,15)LoB
   15   format("#Setting your L/B= ",f6.3,"  to max allowed = 2.75")
      endif
      write(6,16)LoB
   16 format("#LoB= ",f6.3)
      go to 35

   30 continue
      write(0,*)"#Error: reading arg. L/B:",argv
      write(6,*)"#Error: reading arg. L/B:",argv
      stop


   35 continue

      write(6,25)
   25 format("#iaoc iaob a/c a/B      L/B     Mkm          Mkb")
C     LoB is constant for each run of this program

C     Loop for the a/c variable.
      aocinterval=0.050
      aoc=0.1
C     aoc gets incremented at end of do loop
      do 2000 iaoc=1,18

C     Loop for the a/B  variable.
C     Distance between a/B points is:
      aobinterval= 0.005
C       small interval since coming out of a notch
      aob=0.005
      do 1000 iaob=1,200
C      aob will be incremented at end of loop

C     Ok, generate all the parameters

      if(aoc .ge. 0.1 .and. aoc .le. 1.0 .and.
     &   aob .ge. 0.005 .and. aob .lt. 1.0)then

C ------------- Case A: Axial stress    ---------------------------------
C       Get Mkm  membrane stress Factor
        g1= -1.0343*(aoc**2) -0.15657*(aoc) +1.3409

        g2= 1.3218*( aoc**(-0.61153) )

        g3= -0.87238*aoc +1.2788

        g4= -0.46190*(aoc**3) -0.67090*(aoc**2) -0.37571*aoc +4.6511

        fm1= 0.43358*(  aob**(g1 + (g2*aob)**g3 ) )  +
     &      0.93163* exp( aob**(-0.050966) )    + g4


        fm2= -0.21521*( (1-aob)**(176.4199) )   +
     &       2.8141* ( aob**(-0.10740*aob) )

        g5= -0.015647*(LoB**3) +0.090889*(LoB**2) -0.17180*LoB -0.24587

        g6= -0.20136*(LoB**2) +0.93311*LoB -0.41496

        g7= 0.20188*(LoB**2) -0.97857*LoB +0.068225

        g8= -0.027338*(LoB**2) +0.12551*LoB -11.218

        fm3= 0.33994*(aob**(g5)) +1.9493*( aob**(0.23003) ) +
     &      (g6*(aob**2) +g7*aob +g8)

        xMkm= fm1 +fm2 +fm3
        if( xMkm .lt.1.0) xMkm=1.0
        Mkm(iaob,iaoc)= xMkm


C  -------------- CASE B   Mkb -----------------------------------------
C       Ok, we have Mm,   now get Mb

        g1= -0.014992*(aoc**2) -0.021401*aoc -0.23851

        g2= 0.61775*(aoc**(-1.0278) )

        g3= 0.00013242*aoc -1.4744

        g4= -0.28783*( aoc**3 ) +0.58706*(aoc**2) -0.37198*aoc -0.89887

        f1= 0.065916*( aob**(g1+ (g2*aob)**(g3) ) )   +
     &      0.52086*(exp(aob**(-0.10364) ) )  +g4

        g5= -17.195*( aob**2 ) +12.468*aob -0.51662

        f2= -0.021950*( (1.0 -aob)**(2.8086) )  +
     &      0.021403*( aob**(g5) )

        g6= -0.059798*(LoB**3) +0.38091*(LoB**2) 
     &      -0.8022020*LoB +0.31906

        g7= -0.35848*(LoB**2) +1.3975*LoB -1.7535

        g8= 0.31288*(LoB**2) -1.3599*LoB +1.6611

        g9= -0.001470*(LoB**2) -0.0025074*Lob -0.0089846

        f3= 0.23344*( aob**(g6) ) -0.14827*( aob**(-0.20077) ) +
     &      +( g7*(aob**2) +g8*aob +g9 )

        xMkb= f1+f2+f3
        if(xMkb .lt. 1.0)xMkb=1.0
        Mkb(iaob,iaoc)=xMkb


CpointPlot
        write(6,*)iaoc,iaob,aoc,aob,LoB,xMkm,xMkb
  400   format(i3,1x,i3, f6.3,1x,f6.3,1x,f5.2, 2(1x,f6.3) )
C        write(6,400)iaoc,iaob,aoc,aob,LoB,xMkm,xMkb,
C     &              fm1,fm2,fm3,f1,f2,f3
C 400   format(i3,1x,i3, f6.3,1x,f6.3,1x,f5.2, 8(1x,f6.3) )
        go to 950
      endif

C     We should not ever end up here
  450  write(0,*)"# Error near Sta.No. 450 in prog getMb...Stopping."
       stop

C     End of loop for ao2c
  950 continue
      aob=aob+aobinterval
      if(aob .ge. 1.0) go to 1900
 1000 continue

C     One row is done, write it out for the table
Ctable      write(6,470)aoc,(Mb(iaoc,j),j=1,11)
Ctable  470 format(f4.2,2x,11(e13.6,1x))

C     Put out a blank line for better gnuplot in non-table output
 1900 continue
      write(6,*)
      aoc=aoc+aocinterval
      if(aoc .gt. 1.0)go to 9000
 2000 continue

 9000 continue
      write(0,*)"# Finished L/B = ",LoB
      stop
      end
