C getfw.f   vers 0.4  Surface flaw Finite Width Corr Factor  Oct 28 2012
      SUBROUTINE fwSurfFlaw(ivec,iret,cow,aob,xfw)
      SAVE
C  ivec=0=create Matrix,  ivec=1=interpolate, ivec=2=debug ON
C Create the table of fw values for a surface flaw as per BS7910 pg.185
C    and if ivec=1  interpolate for xfw  and return to mainline
C c/W  is 1/2 surf crack length / plate Width
C a/B  is   crack depth/ thickness
C compile:   gfortran -g -w  getfw.f  -o getfw
C Usage:    getfw  >  t_fw_SurfFlaw

C For plot output file :
C       grep #plotfw t_fw_SurfFlaw | delete1arg |delete1arg >fwData4Plot
C   then in gnuplot do:
C       set grid
C       set term x11 enhanced font "arial,15"
C       set xlabel "c/W"
C       set ylabel "a/B"
C       set zlabel "fw" rotate by 90
C       set title "Surface Flaw Finite Width Correction"
C       set title "Plate Surface Flaw Finite Width Correction"
C       splot "fwData4Plot" u 3:4:5
C For Matrix viewing:
C       grep -v #plotfw t_fw_SurfFlaw | delete1arg > t_fw_SurfFlawMatrix

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

C There are two variables that define the matrix of fw
C   1. cow :   c/W     where W= plate width,  2c= crack width at surface
C   2. aob  :  a/B     where a is crack depth, B= plate thickness

C A matrix is created by this program, also the plot file information.


      integer*4 ivec,iret

C     Matrix for storing  fw
      real*4 fw(26,26)
      integer*4 nfwdatarow,maxfwdatarow, nfwdatacol, maxfwdatacol
C     Limits and intervals of input data used for quicker interpolation
      real*4 startcow, endcow, deltacow
      real*4 startaob3, endaob3, deltaaob3
      common/FWDATA/ fw
     &      nfwdatarow,maxfwdatarow, nfwdatacol, maxfwdatacol,
     &      startcow,endcow,deltacow, startaob3,endaob3,deltaaob3

      logical debug
      debug= .false.
      if(ivec.eq.2)debug= .true.   !set to true if testing code
      if(ivec.eq.1 .or. ivec.eq.2)go to 3000  !go and interpolate

C     Check on BS standard limits for aob and cow

C     ivec=0  Initilize the matrix and other variables
      nrows=26
      ncols=26
      maxfwdatarow=nrows  !save in common block
      maxfwdatacol=ncols

      startaob3=0.0
      endaob3=1.0
      deltaaob3= (endaob3-startaob3)/(float(ncols-1))  ! (1.0-0.)/25 = 0.04

      startcow=0
      endcow=0.40
      deltacow= 0.40/(float(nrows-1))  ! = 0.016

      pi=3.1415927
      write(6,10)nrows,ncols
   10 format("#fw # getfw vers. 0.4 starts..."/
     &   "#fw #ROWS= ",i4/"#fw #COLS= ",i4/
     &   "##fw a/b runs between 0.0  and 1.0, c/W between 0.0 and 0.4")

      write(6,25)
C   25 format("#fw #iao2c  iaob   a/2c   a/B     Mm   phi")
   25 format("#fw #icow c/W                                     a/B ")
      write(6,26)startaob3,endaob3
   26 format("#fw         ",f5.3,196x,f5.3)

C     Distance between a/B points is:
C     If there are 26 locations, there are 26-1 intervals between them
      write(6,30)startaob3,endaob3,deltaaob3
   30 format("#fw #startaob3,endaob3,deltaaob3 = ",3(f6.4,1x))

C     Loop for the c/W  variable.
      write(6,32)deltacow
   32 format("#fw #deltacow= ",f6.4)
      cow=startcow
C     cow gets incremented at end of do loop
      do 2000 icow=1,nrows    ! ----------------------loop  cow----------------

C     Loop for the a/B  variable.
      aob=startaob3
      do 1000 iaob=1,ncols    ! ----------------------loop  a/B----------------
C      aob will be incremented at end of loop

C     Ok, generate all the parameters

      if(cow .ge. 0. .and. cow .le. 1.0)then
        temp=(pi*cow)*sqrt(aob)
        ysecant=1.0/cos(temp)   ! note that sec(a)=  1.0/cos(a)
        xfw=sqrt(ysecant)
        fw(icow,iaob)= xfw
        write(6,50)icow,iaob,cow,aob,xfw
   50   format("#fw #plotfw ",i4,1x,i4,1x,f6.4,1x,f6.4,2x,e13.6)

        if(icow .eq. 1 .and. iaob .eq. 1)then
          write(6,54)startcow,startaob3
   54     format("#fw #STARTCOW= ",e14.7/"#fw #STARTAOB= ",e14.7)
        endif
      if(icow.eq.maxfwdatarow .and. iaob.eq.maxdatacol)then
         write(0,2050)startcow,endcow,deltacow
         write(6,2050)startcow,endcow,deltacow
 2050    format("#fw #STARTCOW= ",e14.7/"#fw #ENDCOW= ",e14.7/
     &          "#fw #DELTACOW= ",e14.7)
         write(0,2052)startaob3,endaob3,deltaaob3
         write(6,2052)startaob3,endaob3,deltaaob3
 2052    format("#fw #STARTAOB= ",e14.7/"#fw #ENDAOB= ",e14.7/
     &          "#fw #DELTAAOB= ",e14.7)
      endif

        go to 950
      endif

C     We should not ever end up here
  450  write(0,*)"# Error near Sta.No. 450 in prog getfw...Stopping."
       stop

C     End of loop for aob
  950 continue
      aob=aob+deltaaob3
 1000 continue

C     One row is done, write it out for the table
      write(6,470)icow,cow,(fw(icow,j),j=1,ncols)
  470 format("#fw",i3,1x,f4.3,2x,50(f6.3,1x))  ! not all are used, hopefully

      cow=cow+deltacow
      write(6,475)         !write out a blank line for cow plot seperation
  475 format("#fw #plotfw")
 2000 continue

      write(0,*)"#fw  Finished creating fw Matrix "
      iret=0
      return



C================== Interpolation ===================================
 3000 continue
      if(debug)then
        write(0,3010)
 3010   format(" Enter c/W  and  a/B : ",$)
        read(5,*,end=9000)cow,aob
      endif
C     if not debug,  cow and  aob  are passed into s/r

      if(cow .le.0. .and. aob .le.0.)then
        write(0,3020)cow,aob
 3020   format("#fw #Error: fwSurfFlaw: c/W,a/B = ",e14.7,",",e14.7/
     &         "#fw #Stopping...")
        stop
      endif

      i=ifix((cow-startcow)/deltacow)+1
      j=ifix((aob-startaob3)/deltaaob3)+1
      if(i.lt.1 .or. j.lt.1 .or. j.gt.maxfwdatacol .or.
     &   i.gt. maxfwdatarow) then
C        Out of bounds error
         write(0,3025)startcow,cow,endcow,
     &                startaob3,aob3,endaob3,i,j
 3025    format("#ERROR: S/R getfw c/W or a/B out of bounds:"/
     &    "#      start c/W=",f6.3," c/W= ",e14.7," end c/W= ",f6.3/
     &    "#      start a/B=",f6.3," a/B= ",e14.7," end a/B= ",f6.3/
     &    "#      i,j=  ",i5,1x,i5/
     &    "# Stopping...")
          iret=3025    ! set the error flag to format sta. no.
          return
      endif

C     j,i  is the corner (see internet address above)  of x1,y1
      ip1=i+1
      jp1=j+1
C     We need to use x1, x2, y1, y2 to calculate Mm and Mb
C     x,y  are aob,cow
      x=aob
      y=cow

C     Things appear to be within bounds of matrices, continue
      x1=startaob3+float(j-1)*deltaaob3
      x2=x1+deltaaob3
      y1=startcow + float(i-1)*deltacow
      y2=y1+deltacow
      if(debug)then
        write(0,3035)x1,y1,  x2,y1, x,y, x1,y2, x2,y2
        write(6,3035)x1,y1,  x2,y1, x,y, x1,y2, x2,y2
 3035   format("#fw #TargetCo-ords:"/
     &        "#fw #x1,y1       x2,y1  :",2f7.4,14x,2f7.4/
     &        "#fw #       x,y         :",15x,2f7.4/
     &        "#fw #x1,y2       x2,y2  :",2f7.4,14x,2f7.4/
     &  )
        write(0,3036)i,j,  i,jp1,  ip1,j, ip1,jp1
        write(6,3036)i,j,  i,jp1,  ip1,j, ip1,jp1
 3036   format("#fw # i,j       i,j+1    : "i3,",",i3, 10x, i3,",",i3/
     &         "#fw # i+1,j     i+1,j+1  : "i3,",",i3, 10x, i3,",",i3/
     &  )
      endif

      Q11=fw(i,j)
      Q12=fw(ip1,j)
      Q21=fw(i,jp1)
      Q22=fw(ip1,jp1)

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
      xfw = frac1y*R1 + frac2y*R2
      if(debug)then
        write(0,3040)cow,aob,i,j, Q11,R1,Q21, xfw, Q12,R2,Q22
        write(6,3040)cow,aob,i,j, Q11,R1,Q21, xfw, Q12,R2,Q22
 3040   format("#fw #Solution for c/W= ",f6.3," a/B= ",f6.3,
     &        " (i,j)=",2i3/
     &        "#fw # Q11,  R1,  Q21 :  ",3(f7.4,1x)/
     &        "#fw #       fw =     :  ",10x,f7.4/
     &        "#fw # Q12,  R2,  Q22 :  ",3(f7.4,1x)/
     &        )
      endif

      iret=0
      return
      go to 3000

 9000 stop
      end
