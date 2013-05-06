C  replaceline.f   used to read a line from replaceline.env  and
C     depending on first string, replace the line(s) in input file
      SAVE
C Compile   gfortran  -g -w -fbounds-check replaceline.f  -o replaceline
C  Usage:   replaceline  <infile  >outputFile
C           ( substitute line is in file: replaceline.env )
C        

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

      character*50   string1,string2
      character*300  inp300,jnp300
      character*1    inpone(300),jnpone(300)
      character*10   inpten(30)
      equivalence (inpone(1), inpten(1), inp300)
      equivalence (jnpone(1),jnp300)

      maxchars=300

      open(unit=10,file="replaceline.env")
  100 continue
      read(10,"(a300)",end=8000) inp300
      if(inp300 .eq. " ")then
       write(0,105)
  105  format("#Error: replaceline.env: Replacement line is empty  !!!")
       stop
      endif
      read(inp300,*)string1
      close(unit=10)


      ninput=0
      nfound=0
      ioutdevice=6    !stdo
  800 continue
C     Loop back here for next input line
      read(5,"(a300)",end=1000)jnp300
      ninput=ninput+1
C     Now get first item from line
      if(jnp300 .eq." ")then   ! in case of blank line
        write(6,*)
        go to 800
      endif
      read(jnp300,*)string2
      if(string2 .eq. string1)then
C        Yes, we have found a match line. Replace it
         nfound=nfound+1
         call subCutBlanks(ioutdevice,inpone,maxchars)
         go to 800
      endif
C     No, line does not match,  write out the line as-is
      call subCutBlanks(ioutdevice,jnpone,maxchars)
      goto 800


 1000 continue
C     File has been completed.  Check if anything happened.
      write(0,1010)ninput, nfound
 1010 format("#replaceline: read in: ",i9," lines, replaced: ",i9,
     &       " lines.")
      if(nfound .ne. 1)then
        write(0,1020)nfound
 1020   format(/"#Warning!: replacelines: found ",i9,
     &         " lines to replace"/)
      endif
      stop

 8000 continue
      write(0,8010)
 8010 format(/"#Error: replacelines: No input lines in file= ",
     &       "replacelines.env"/)
      stop

      end


C====================================================================
      SUBROUTINE subCutBlanks(idevice,inp1,idimen)
      SAVE
C  subCutBlanks.f   s/r to remove the leading and trailing blanks for output
C  Usage:    call subCutBlanks(idevice,inp1,idimen)
C                              idevice= output device no. in fortran
C                              inp1=    char storage vector
C                              idimen=  length of incoming char string
C                                       (or how many to look at ?)

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

      integer*4  idevice,idim, ncharmax/400/
      integer*4  ifirst,ilast
      character*1 inp1(400)  ! if you change this, change also ncharmax above

C     Check for programming error:
      if(idimen .gt. ncharmax)then
       write(0,100)idimen,ncharmax
  100  format("#Error: s/r subCutBlanks(): dimension of char varible= ",
     & i10," bigger than max allowed =",i5," stopping...")
       stop
      endif
      if(idimen .lt. 1 )then
       write(0,100)idimen
  101  format("#Error: s/r subCutBlanks(): dimension of char varible= ",
     & " is less than 1,  stopping...")
       stop
      endif
         
C     Find the first non-blank character
      do 200 i=1,ncharmax
         if(inp1(i) .eq. " ")goto 200
         ifirst=1
         goto 210
  200 continue
C     if we end up here there are no chars within inp(1) and inp(idimen)
C     Just write out one blank for this line
      write(idevice,205)
  205 format(" ")
      return

  210 continue  !  now hunt for the last char by starting at idimen and 
C                  working backwards.
      do 300 i=idimen,1,-1
         if(inp1(i) .eq. " ")goto 300
         ilast=i  !  No?, then we have found an non blank
         go to 310
  300 continue
C     If we end up here, the whole line was blanks, write out one blank
      write(idevice,205)
      return

  310 continue  ! ifirst and ilast are set, write out the stuff
      write(idevice,315)(inp1(i),i=ifirst,ilast)
  315 format(400a1)
      return
      end




