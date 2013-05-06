C getMk_w00.f   vers 0.1  Surface flaw near Weld Aug 10 2012
      SAVE
C Create the table of Mkm and Mkb values for a surface flaw as per BS7910 pg.224
C "Mkm"  is the membrane (axial) stress intensity multiplier
C "Mkb"  is the Plate in bending multiplier.
C  00 is for flaw edge at surface, 90  refers to the Deepest location on flaw
C compile:   gfortran -g -w -fbounds-check getMk_w00.f  -o getMk_w00
C Usage:    getMk_w00 LoB > tableMk_w00
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
      real*4 LoB,LoBsq,LoBcube
      character*80 argv,fname
      integer*4  iargc,argc

      nrows=200
      ncols=18

      pi=3.1415927
      write(6,10)nrows,ncols
   10 format("# getMk_w00 vers. 0.1 starts..."/
     &   "#ROWS= ",i4/"#COLS= ",i4)
      argc= iargc()
      if(argc .ne. 1)then
        write(0,*)"#Error: L/B argument needed"
        write(6,*)"#Error: L/B argument needed:"
        write(0,*)"#e.g.:  getMk_w00  1.5  > outfile"
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
        LoBsq= LoB**2
        LoBcube=LoB**3

C     Loop for the a/c variable.
      aocinterval=0.050
      aoc=0.1
C     aoc gets incremented at end of do loop
      do 2000 iaoc=1,18
        coa= 1.0/aoc
        coasq= coa**2

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

        g1= 0.0078157*coasq -0.070664*coa +1.8508

        g2= -0.000054546*LoBsq +0.00013651*LoB -0.00047844

        g3= 0.00049192*LoBsq -0.0013595*LoB +0.011400

        g4= 0.0071654*LoBsq -0.033399 -0.25064

        g5= -0.018640*coasq +0.24311*coa -1.7644

        g6= -0.0016713*LoBsq +0.0090620*LoB -0.016479

        g7= -0.0031615*LoBsq -0.010944*LoB +0.13967

        g8= -0.045206*LoBcube +0.32380*LoBsq -0.68935*LoB +1.4954

        fm1= g1*(aob**(g2*coasq +g3*coa +g4) ) +
     &       g5*(  (1.0-aob)**( g6*coasq +g7*coa +g8 ) )

        g9= -0.25473*(aoc**2) +0.40928*aoc +0.0021892

        g10= 37.423*(aoc**2) -15.741*aoc +64.903

        fm2= (-0.28639*(aoc**2) +0.35411*aoc +1.6430 )*(aob**(g9)) +
     &        0.27449*( (1-aob)**(g10) )

        g11= 0.10553*LoBcube +0.59895*LoBsq -1.0942*LoB -1.2650

        g12= 0.043891*LoBcube -0.24898*LoBsq +0.44732*LoB +0.60136

        g13= -0.011411*(aoc**2) +0.004369*aoc +0.51732

        fm3= g11*(aob**(0.75429)) +g12*exp( aob**(g13) )


        xMkm= fm1 *fm2 *fm3
C       BS7910 2005 does not call for this, however put it in anyway:
        if( xMkm .lt.1.0) xMkm=1.0
        Mkm(iaob,iaoc)= xMkm


C  -------------- CASE B   Mkb -----------------------------------------
C       Ok, we have Mkm,   now get Mkb

        g1= 0.0023232*coasq -0.00037156*coa +4.5985

        g2= -0.000044010*LoBsq +0.00014425*LoB -0.00086706

        g3= 0.00039951*LoBsq -0.0013715*LoB +0.014251

        g4= 0.0046169*LoBsq -0.017917*LoB -0.16335

        g5= -0.018524*coasq +0.27810*coa -5.4253

        g6= -0.00037981*LoBsq +0.0025078*LoB +0.00014693

        g7= -0.0038508*LoBsq +0.0023212*LoB -0.026862

        g8= -0.011911*LoBcube +0.082625*LoBsq -0.16086*LoB +1.2302

        g9= 0.27798*(aob**3) -1.2144*(aob**2) -2.4860*aob +0.099981

        f1= g1*( aob**(g2*coasq +g3*coa +g4) ) +
     &      g5*(  (1.0-aob)**(g6*coasq +g7*coa +g8) ) +g9

        g10= -0.25922*(aoc**2) +0.39566*aoc +0.011759

        g11= 6.5964*(aoc**2) +55.787*aoc +37.053

        f2= ( -0.35006*(aoc**2) +0.40768*aoc +1.7053)*(aob**(g10)) +
     &        0.24988*( (1.0-aob)**(g11) )

        g12= -0.14895*LoBcube +0.815*LoBsq -1.4795*LoB -0.89808

        g13= 0.055459*LoBcube -0.30180*LoBsq +0.54154*LoB +0.53433

        g14= -0.01343*(aoc**2) +0.0066702*aoc +0.75939

        f3= g12*( aob**(0.94761) ) +g13*exp( aob**(g14) )


        xMkb= f1*f2*f3
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
