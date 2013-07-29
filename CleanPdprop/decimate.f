C decimate.f    remove lines from a file
      SAVE
C  Compile:  gfortran  -g -w -fbounds-check decimate.f  -o decimate
C  Usage eg.:   decimate  99  <infile >outfile
C                                            saves every 100th line
C  Max no. of chars per line is 300
C  Comment lines are not counted in decimation process.

C  Copyright (C) 2012  Al Conle
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


      integer idimen/300/   ! the dimension of char vector inpone()
      character*300  inp300 ! used to read in lines as chars.
      character*1    inpone(300)
      character*10   inpten(30)
      equivalence (inpone(1), inpten(1), inp300)
      

      character*80 argv
      integer*4 iargc,nargc

      write(6,165)
      write(0,165)
  165 format("# decimate.f vers. 1.0"/
     & "#Usage e.g.: decimate 99 <infile  >outfile"/
     & "#        saves every 100 th  line.")

      nargc = iargc()
C     Note that in HP fortran the arg numbers must be changed by -1
C     and that iargc probably includes the "decimate" as an arg.
      if( nargc .ne. 1)then
        write(0,*)" decimate:  usage ERROR"
        write(0,*)"# Usage e.g.: decimate 99  <infile  >outfile"
        write(6,*)" decimate:  usage ERROR"
        write(6,*)"# Usage e.g.: decimate 99  <infile  >outfile"
        stop
      endif


C       The first arg is the no. of points removed between saved point
        jvect=1
        call getarg(jvect,argv)
        read(argv,*,err= 178)nremove
        write(6,*)"#nremove= ",nremove
        if(nremove .gt. 1)go to 180
C     Bad arguments in command line
  178 write(0,*)"# ERROR: bad removal argument=",argv
      write(0,*)"# Usage e.g.: decimate 99  <infile  >outfile"
      stop

  180 continue
      write(0,185)nremove
      write(6,185)nremove
  185 format("#decimate by removing ",i7," lines, for each 1 saved.")



  200 continue
C     Get the first data line and keep it. Then delete "nremove" data lines
      
      ninput=0
      ndelete=0
      nsavedata=0
  700 continue
c     Loop back to here for next input line.

      read(5,"(a300)",end=750)inp300
      ninput=ninput+1
Cdebug      write(0,*)" read input line ",ninput

C     Check for blank line
      if(inp300.eq." ")then
C        write(6,"(a1)")" "
        go to 700
      endif

      if(inpone(1).ne."#")then
C        We may have a data line, or someone screwed up and put the # later in line.
C        See if 1st char is a # later in line
         do 705 i=1,300
           if(inpone(i).eq." ") go to 705
C           No? 1st non-blank found, if # its not a data line
            loc=i
           if(inpone(i).eq."#") then
C            Shift the whole mess over to begining of field
             inew=0
             do 703 j=loc,300
              inew=inew+1
              inpone(inew)=inpone(j)
  703        continue
Cdebug        write(0,*)"Shifted a comment: ",(inpone(j),j=1,inew)
C            Ok, now its a nice comment. Go play with it.
             go to 730
           else
C            The first non blank is not a #
             go to 710
           endif
  705      continue
  710      continue
C        1st non blank is not a #.  Check if its a number.
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
            write(0,715)ninput
            write(6,715)ninput
  715       format("# ERROR : input line no. ",I5,
     &      " not a # and not a number. Fix input file. ")
            call subCutBlanks(0,inpone,idimen)
            stop
         endif

C        Ok, its a number
         ndata=ndata+1

         if(nsavedata .eq.0) then  ! its the first data line, save
           nsavedata=nsavedata+1
           call subCutBlanks(6,inpone,idimen)
           go to 700
          endif

C       Potential delete point
        ndelete=ndelete+1
        if(ndelete .gt. nremove)then  ! enough deleted, save one
          nsavedata=nsavedata+1
          call subCutBlanks(6,inpone,idimen)
          ndelete=0
          go to 700
         endif
C        ok,  decimate or skip this point
         go to 700
         
      endif
  730    continue
C        Simply write out the comment line
         call subCutBlanks(6,inpone,idimen)
         go to 700

  750    continue   ! end of file comes here
C      All data is in. 
       write(0,755)ninput,ndata,nsavedata
  755  format("#Done. Total input lines including comments= ",i8/
     &  "#of which ",i8," were data lines.  Lines saved= ",i8)
       stop
       end



C==============================================================

      SUBROUTINE subCutBlanks(idevice,inp1,idimen)
      SAVE
C  subCutBlanks.f  s/r to cut the leading and trailing blanks from output line
C     Used to write out long char. vector lines, usually comments.
C  Usage:    call subCutBlanks(idevice,inp1,idimen)
C                              idevice= output device no. in fortran
C                              inp400=  char storage vector
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


