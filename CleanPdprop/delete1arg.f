C delete1arg.f  Strip off the first field of each line in a table
C Compile:  gfortran  -g -w -fbounds-check delete1arg.f  -o delete1arg
C Usage:    delete1arg  <input  >output
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

C  Program is intended as a post processer to the saefcalc1.f program.
C  saefcalc1.f puts out special tagged lines in its output stream, such
C  as #xcalc1,  #xcalc2,   #plotloops   etc.
C  Using "grep" one can strip out these tagged lines, and then remove the
C  tags using this program.

      CHARACTER*1   INP1(80)
      CHARACTER*5   INP5(16),JNP5(16)
      CHARACTER*10  INP10(8),JNP10(56)

      CHARACTER*80  JNP80(7),INP80
      EQUIVALENCE  (INP1(1),INP5(1),INP10(1),INP80),
     &             (JNP5(1),JNP10(1),JNP80(1))

      character*300  inp300,jnp300
      character*1    inpone(300),jnpone(300)
      character*10   inpten(30)
      equivalence (inpone(1), inpten(1), inp300)
      equivalence (jnpone(1),jnp300)


      ninput=0

C---------------------------- Read in from stdin-------------
  800 continue
c     Loop back to here for next input line.


      read(5,"(a300)",end=9000)inp300
      ninput=ninput+1
Cdebug      write(0,*)" read input line ",ninput
   
C     Check for all blank line
      if(inp300.eq." ")then
C        write(6,"(a1)")" "
        go to 800
      endif

      if(inpone(1).eq." ")then
C        We may have a line, but the 1st char is later in line.
C        Find the first non-blank
         do 805 i=1,300
           if(inpone(i).eq." ") go to 805
C           No? 1st non-blank found,
            loc=i
C            Shift the whole mess over to begining of field
             inew=0
             do 803 j=loc,300
              inew=inew+1
              inpone(inew)=inpone(j)
  803        continue
Cdebug             write(0,*)"Shifted all chars left for blanks: ",
Cdebug               (inpone(j),j=1,inew)
C            Ok, now its a nice line with something in 1st col. Go play with it.
             go to 810
  805    continue
      endif

      
  810    continue
C        1st non blank is now in col1.  
C        Find the end of the first field
         lastlocfield1=0
         do 850 i=1,300
          if(inpone(i).ne. " ")go to 850
C         ok, found end of 1st field
          lastlocfield1=i-1
C         So where does the next non-blank field start?
          locfield2=0
          do 840 j=i,300
           if(inpone(j).eq." ")go to 840
C          Found 2nd field
           locfield2=j
           go to 860
  840     continue
C         We should only end up here if there is no 2nd field (rest is blank)
          go to 8000

  850    continue
C        we would end up here if all 300 chars are non-blank, thus only 1 field total
         go to 8000

  860    continue
C        Now shift the rest over to the 1st col
         inew=0
         do 863 j=locfield2,300
          inew=inew+1
          inpone(inew)=inpone(j)
  863    continue



  950 continue
C     Write out the comment line
C       Count backwards and see where the last char is
        do 960 i=1,300
          j=300-(i-1)
          if(inpone(j).ne." ")then
C           found last char
            lastloc=j
            go to 962
          endif
  960   continue

  962   continue
      write(6,"(300a1)")(inpone(i),i=1,lastloc)
C     Go read another line
      go to 800





 8000 continue
C     After 1st field removal, line is blank
      write(6,*)
      go to 800


 9000 continue
      stop
      end

