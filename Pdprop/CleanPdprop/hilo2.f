c  hilo2.f  v2.1   prog to read lines with columns of numbers  rev. May2008
c              and find highs and lows.
C Compile:  gfortran  -g -w -fbounds-check hilo2.f  -o hilo2
c Usage:     hilo2   <infile >outfile
c
c  Input file is read free form, output is list of unique first column elements
C  
C---------------------------------------------------------------------------
C  Copyright (C) 2000 SAE Fatigue Design and Evaluation Committee
C  This program is free software; you can redistribute it and/or
C  modify it under the terms of the GNU General Public License as
C  published by the Free Software Foundation; either version 2 of the
C  license, or (at your option) any later version.
C
C  This program is distributed in the hope tha it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General PUblic License for more details.
C
C  You should have received a copy of the GNU General PUblic License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA
C  Try also their web site: http://www.gnu.org/copyleft/gpl.html
C---------------------------------------------------------------------------

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
C  Thus skip lines that start with "Stage,"  "Point," 
C  or any non-numeric character in 1st field.
C
COct5.08  ver2.1 : fixed ndatalines bug.  Took out excess write statem.

      character*200  jnp200
      character*5    jnp5(40)
      character*1    jnp1(200)
      equivalence    (jnp200, jnp5(1), jnp1(1))

      CHARACTER*80 INP80, argv, field1
      CHARACTER*5  INP5(16)
      CHARACTER*1  INP1(80)
      EQUIVALENCE (INP80,INP5(1),INP1(1))

      REAL xhigh(50), xlow(50), xdata(50)

      INTEGER UNO,UNI

C
      uni=5
      UNO=6
      nlistmax=50
C     if you change the dimensions of xhigh, xlow,etc above, change this too

      write(0,*) "hilo2 Starts. Usage: ",
     &" hilo2  <in  >out"
C
      if(iargc().ne. 0)then
           write(0,*) " Error, extra args, should be: ",
     &     " hilo2  <in  >out"
           stop
      endif
C
C
      NL=0
      ndatalines=0

 2010    continue
 1010    FORMAT(A200)
         read(UNI,1010,end=1190)jnp200
         NL=NL+1
         if(jnp200.eq." ")go to 2010

c        Find the first non blank character.
         n=0
         do 2015 i=1,200
           if(jnp1(i).eq." ")go to 2015
c          it is a non blank char
           n=i
           go to 2020
 2015    continue

 2020    continue
c        If 1st char is a number or + or -, then assume it is a point number
c        otherwise it is a character word and of no interest
         if(jnp1(n).eq."-" .or.
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

C     Allow nlistmax columns max.  This is the first data line.  Figure
C          out how many columns.
        ndatalines=ndatalines+1
        ncols=0
        do 2025 j=1,nlistmax+1
          ncols=ncols+1
          if(ncols.gt.nlistmax)then
            write(0,*)"Error: Too many columns (",
     &     nlistmax," max). Recompile hilo2.f  Stopping now."
            stop
          endif
C         Repeatedly read the same line, one exta arg each time
C         when we run out or args we know how many there are
          read(jnp200,*,err=2030,end=2035)(xdata(i),i=1,ncols)
          write(0,2008)ncols
 2008     format('#The 1st data line has ',I2,' columns of numbers')
 2025   CONTINUE
C       If we reach this point
C
 2030   write(0,2032)jnp200
 2032   format('ERROR: non-number in 1st data line: ',A200)
        write(0,*)"At file line number= ",NL
        stop
C
 2035   continue
        ncols=ncols-1
C       We now have read the 1st data line. It had ncol columns.
C       For all the rest of the data lines assume ncol  is valid.
C       Initialize the highs and lows registers:
        do 2037 i=1,ncols
          xhigh(i)=xdata(i)
          xlow(i)=xdata(i)
 2037   continue
        go to 2040
      endif
c     We have not found the data line yet. Keep looking
      go to 2010

C     Look for rest of data.
 2040 continue
c        Fetch a new line. It may be comment.
         read(UNI,1010,end=1190)jnp200
         NL=NL+1
Cdebug         write(0,*)"Line ",NL,jnp200
         if(jnp200.eq." ")go to 2040

c        Find the first non blank character.
         n=0
         do 2045 i=1,200
           if(jnp1(i).eq." ")go to 2045
c          it is a non blank char
           n=i
           go to 2046
 2045    continue

 2046    continue
c        If 1st char is a number or + or -, then assume it is a point number
c        otherwise it is a character word and of no interest
         if(jnp1(n).eq."-" .or.
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
C           Read in the whole line of ncols numbers
            read(jnp200,*,err=2070)(xdata(i),i=1,ncols)
Cdebug            write(0,*)(xdata(i),i=1,ncols)
            ndatalines=ndatalines+1
            do 2060 i=1,ncols
             if(xhigh(i).lt.xdata(i))xhigh(i)=xdata(i)
             if(xlow(i).gt.xdata(i))xlow(i)=xdata(i)
 2060       continue
         endif
c        if we get to here it is not a data line, or we are done
c        with this data line. Go fetch the next line
         go to 2040



 2070    write(0,*)"#ERROR: non-number? in data line no=",NL,
     &  " : ",jnp200
        write(0,*)"#***Assume line is a comment"

        write(6,2072)NL,jnp200
 2072   format("#ERROR: non-number in data line no=",i6,
     &  " : ",a200/"#***Assume line is a comment")
C       Keep going
        go to 2040



 1190  Continue
           IER=0
           write(uno,*)(xhigh(i),i=1,ncols)
           write(uno,*)(xlow(i),i=1,ncols)
           write(0,*)"#Done. Scanned ",ndatalines," data lines."
           write(0,*)"#Total lines= ",NL
           write(0,*)"#Highs:",(xhigh(i),i=1,ncols)
           write(0,*)"#Lows :",(xlow(i),i=1,ncols)
   
           if(uno.ne.6)close(unit=uno)
       stop  
 9000 Continue
          write(0,*)"Encountered EOF during 1st line read" 
      stop
      end   

C==============================================================
      SUBROUTINE WR(UN,INP1)
C S/R TO WRITE A LINE BUT CUT THE TRAILING BLANKS
      CHARACTER*1 INP1(200)
      INTEGER UN
      LONG=200
C FIND LAST NON BLANK CHAR
      DO 10 I=LONG,1,-1
        IF(INP1(I).EQ.' ')GO TO 10
C       NO? THIS IS IT
        N=I
        GO TO 20
   10 CONTINUE
C COMPLETELY BLANK LINE
      WRITE(UN,11)
   11 FORMAT(' ')
      RETURN
   20 WRITE(UN,24)(INP1(J),J=1,N)
   24 FORMAT(200A1)
      RETURN
      END

