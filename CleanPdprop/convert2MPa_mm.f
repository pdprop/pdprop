C   convert2MPa_mm.f   vers. 1.4   dec 10 2013   FAC
      SAVE
C     Convert points from a digital dadn vs DeltaK  plot into 
C     da/dn (mm) vs. Delta_K table (mpa_mm)
C     Used mostly to do a unit conversion  table->convert->table
C     Compile:  gfortran -g -w -fbounds-check convert2MPa_mm.f -o convert2MPa_mm
C     Usage:    convert2MPa_mm  <inputfile   >dadnFileName.dadn
C
C  Copyright (C) 2013 Fatigue Des. and Eval. Comm. and F.A.Conle
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

C vers 1.4 Add output for #deltaKunits= mpa_mm and #dadnUnits= mm Dec10 2013
C vers 1.3 forks  mkparis2.f  to convert2Mpa_mm.  Uses seperate definition
C          of the units for deltaK and dadn.
C vers 1.2 increases no. of data pts to  1250
C vers 1.1  Set output format to be   N/(mm**(3/2) ) which is same as MPa*sqrt(mm)

C  Input file format:
C    ______ first column in file
C   |
C   v
C   #NAME= G40.11_steel
C   #DELTAKUNITS= ksi_in   # can be: ksi_in,  mpa_m,  mpa_mm
C   #DADNUNITS= inch       # can be: inch,  meter, mm
C   #  It is better not to show deltaKc  or deltaK threshold values, since
C   #  they are ignored in the simulations anyway.
C   8.64953 1.10478e-07
C   11.0953 4.18774e-07
C   11.7759 4.04328e-08
C    etc...
C                      # File "haddadG40.11dadn_no-Lo.txt"  is an example file
C---------------------------------------------------------------------------

C e.g.:   Output file format is a table, preceeded by comments:
C   #Name= G40.21-50A                    
C   #Got Original #units=                     ksi_in    
C    #Got Original #dKth=    18.000000    
C    #Got Original #dKc=    100.00000    
C    #Got Original  #A=   4.99999997E-09
C    #Got Original  #m=    1.9000000    
C   #Conversions:
C   # 1 ksi*sqrt(inch) = 1.0989*MPa*sqrt(m)
C   # 1 MPa*sqrt(m)  =   31.6228*N/(mm**(3/2) )
C   # 1 MPa*sqrt(mm) =         1 N/(mm**(3/2) )
C   
C   #All inputs converted to    MPa*sqrt(mm)   and  mm/cycle
C   #Note that this is same as  N/(mm**(3/2))  and  mm/cycle
C   
C   #deltaKunits= mpa_mm
C   #dadnunits=  mm
C   #   mm/cycle    MPa*sqrt(mm)
C   #     dk            da/dn  
C   0.6255053E+03  0.3081921E-04
C   0.3475030E+04  0.8013158E-03
C etc....
C---------------------------------------------------------------------------

      CHARACTER*1   INP1(80)
      CHARACTER*5   INP5(16),JNP5(16)
      CHARACTER*10  INP10(8),JNP10(56)

      CHARACTER*80  JNP80(7),INP80, INP80temp
      EQUIVALENCE  (INP1(1),INP5(1),INP10(1),INP80),
     &             (JNP5(1),JNP10(1),JNP80(1))

      character*300  inp300,jnp300
      character*1    inpone(300),jnpone(300)
      character*10   inpten(30)
      equivalence (inpone(1), inpten(1), inp300)
      equivalence (jnpone(1),jnp300)

      character*10 dadnunits,deltaKunits
      character*30 name30, firstfield, ctail
      character*30 cfiletype, Cdatatype

 
      integer maxpoints,npts
      real*4  dadn0, dadn1(1250), dadn2(1250)

      real*4  dk0,   dk1(1250),   dk2(1250)

      real*4  A,m,dKth,dKc

      maxpoints=1250

      write(6,701)
      write(0,701)
  701 format("# convert2MPa_mm  vers. 1.4 starts..."/
     &   "# Program convert to table of mpa_mm, from other units table")

      XMPAS=6.894759
C     Set some values to check if user forgets inputs
      deltaKunits=" "
      dadnunits=" "
      name30= " "
      dKth=-999.
      dKc=-999.
C      A=-999.   not used here.
C      m=-999.   not used here.

C---------------------------- Read in Paris parameters from stdin-------------
      ninput=0
      ndata=0
  800 continue
c     Loop back to here for next input line.

C     Input lines, in this version, should be #Comment lines or blank
C        It will be hard to distinguish the real comment from the junk comment.
C        Leave it up to the user to edit the comment section.

      read(5,"(a300)",end=980)inp300
      ninput=ninput+1
Cdebug      write(0,*)" read input line ",ninput

C     Check for blank line
      if(inp300.eq." ")then
C        write(6,"(a1)")" "
        go to 800
      endif

      if(inpone(1).ne."#")then
C        We may have a data line, or someone screwed up and put the # later in line.
C        See if 1st char is a # later in line
         do 805 i=1,300
           if(inpone(i).eq." ") go to 805
C           No? 1st non-blank found, if # its not a data line
            loc=i
           if(inpone(i).eq."#") then
C            Shift the whole mess over to begining of field
             inew=0
             do 803 j=loc,300
              inew=inew+1
              inpone(inew)=inpone(j)
  803        continue
Cdebug             write(0,*)"Shifted a comment: ",(inpone(j),j=1,inew)
C            Ok, now its a nice comment. Go play with it.
             go to 880
           else
C            The 1st non blank is not a #
             go to 810
           endif
  805     continue

  810    continue
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
C            It must be a letter, not a number. Time to bomb out.
             write(0,815)ninput
             write(6,815)ninput
  815        format(" ERROR convert2MPa_mm: input line no. ",I5,
     &       " not a # and not a number. Edit input file")
             stop
           endif
C          Ok, it should be like this one:
C          11.0953 4.18774e-07
         read(inp300,*)dk0,dadn0
         ndata=ndata+1
         if(ndata.gt.maxpoints)then
           write(0,816)
           write(6,816)
  816      format(" Error too many data points:",i5,
     &            " recompile convert2MPa_mm.f")
           stop
         endif

         if(deltaKunits .eq. " " .or. dadnunits .eq. " ")then
           write(0,820)
           write(6,820)
  820      format("#ERROR: convert2MPa_mm: Got a data point but ",
     &           "#deltaKunits= or #dadnunits= not defined yet."/
     &           "#Edit your dadn vs deltaK table file...Stopping")
           stop
         endif

         if(deltaKunits .eq. "ksi_in")then
C           # 1 ksi*sqrt(inch) = 1.0989*MPa*sqrt(m)
C           # 1 MPa*sqrt(m)  =   31.6228*N/(mm**(3/2) )
C           # 1 MPa*sqrt(mm) =         1 N/(mm**(3/2) )
            dk1(ndata)= dk0
            dk2(ndata)= 1.0989 * 31.6228 * dk0
         endif
         if(deltaKunits .eq. "mpa_m")then
C           deltaK data  is in MPA*sqrt(m)
            dk1(ndata)= dk0
            dk2(ndata)=   31.6228*dk0
         endif
         if(deltaKunits .eq. "mpa_mm")then
C           deltaK data  is in MPA*sqrt(mm)
            dk1(ndata)= dk0
            dk2(ndata)= dk0   ! no change
         endif

         if(dadnunits .eq. "inch")then
            dadn1(ndata)=dadn0
            dadn2(ndata)= 25.4*dadn0
         endif
         if(dadnunits .eq. "meter")then
            dadn1(ndata) = dadn0
            dadn2(ndata) = 1000. * dadn0
         endif
         if(dadnunits .eq. "mm")then
            dadn1(ndata) = dadn0
            dadn2(ndata) = dadn0  ! no change
         endif

C        End of data line processing, write it out later, when all are in.
         go to 800
      endif

C     Input line has a # in first col, check it just in case
  880 continue

      if(inpone(1).ne."#")then
C       We should not be here.  something is bad in program
       write(0,*)" ERROR 880, sorry prog. messed up call ? admin"
        stop
      endif

C     Ok, its a nice comment.  Figure out if its a special tag.
      read(inp300,*)firstfield

      if(firstfield .eq."#NAME=" .or.
     &   firstfield .eq."#Name=" .or.
     &   firstfield .eq."#name=" )then
         read(inp300,*) firstfield, name30
         write(6,935)name30
  935    format("#Name= ",a30)
         go to 950
      endif

      if(firstfield .eq."#DELTAKUNITS=" .or.
     &   firstfield .eq."#deltakunits=" .or.
     &   firstfield .eq."#deltaKunits=" .or.
     &   firstfield .eq."#DELTAKunits=" .or.
     &   firstfield .eq."#DeltaKUnits=" )then
         read(inp300,*) firstfield, deltaKunits
         if(deltaKunits.eq. "ksi_in" .or. deltaKunits .eq. "mpa_m" .or.
     &      deltaKunits .eq. "mpa_mm")then
           write(6,940)deltaKunits
  940      format("#Got Original #deltaKunits= ",a30)
           go to 950
         else
           write(0,942)deltaKunits
           write(6,942)deltaKunits
  942      format("#Error:convert2MPa_mm: unknown #deltaKunits= ",
     &            "descriptor:",a30/"#Should be: ksi_in  or  ",
     &            "mpa_m  or mpa_mm"/
     &            "# Edit your  dadn table  file. Stopping...")
           stop
         endif
         goto 950
      endif

      if(firstfield .eq."#DADNUNITS=" .or.
     &   firstfield .eq."#dadnunits=" .or.
     &   firstfield .eq."#dadnUnits=" .or.
     &   firstfield .eq."#DadnUnits=" .or.
     &   firstfield .eq."#DADNunits=" )then
         read(inp300,*) firstfield, dadnunits
         if(dadnunits.eq. "inch" .or. dadnunits .eq. "meter" .or.
     &      dadnunits .eq. "mm" .or. dadnunits .eq. "millimeter")then
           write(6,945)dadnunits
  945      format("#Got Original #dadnunits= ",a30)
           go to 950
         else
           write(0,947)dadnunits
           write(6,947)dadnunits
  947      format("#Error:convert2MPa_mm: unknown #dadnunits= ",
     &            "descriptor:",a30/"#Should be: inch  or  ",
     &            "meter  or mm"/
     &            "# Edit your  dadn table  file. Stopping...")
           stop
         endif
         goto 950
      endif



C     We do not really need dKth and dKc  in convert2MPa_mm.f  It is
C     a table conversion program only
C      if(firstfield .eq."#dKth=" .or.
C     &   firstfield .eq."#dkth=" .or.
C     &   firstfield .eq."#DKTH=" .or.
C     &   firstfield .eq."#DKth=" )then
C         read(inp300,*) firstfield, dKth
C         write(6,*)"#Got Original #dKth= ",dKth
C         go to 950
C      endif
C      if(firstfield .eq."#dKc=" .or.
C     &   firstfield .eq."#dkc=" .or.
C     &   firstfield .eq."#DKC=" .or.
C     &   firstfield .eq."#DKc=" )then
C         read(inp300,*) firstfield, dKc
C         write(6,*)"#Got Original #dKc= ",dKc
C         go to 950
C      endif

C     If it is  none of the above, then assume it is a plain ol' comment
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
  950 continue
      go to 800



C     All input lines have been read.---------------------------------------
  980 continue
C     See if some of the critical values are missing:
C      if(dKth.eq.-999.)then
C        write(0,*)" ERROR #dKth= line missing"
C        write(6,*)" ERROR #dKth= line missing"
C        stop
C      endif
C      if(dKc.eq.-999.)then
C        write(0,*)" ERROR #dKc= line missing"
C        write(6,*)" ERROR #dKc= line missing"
C        stop
C      endif
      if(name30.eq." ")then
        write(0,*)" ERROR #Name= line missing"
        write(6,*)" ERROR #Name= line missing"
        stop
      endif
      if(deltaKunits.eq." ")then
        write(0,*)" ERROR #DELTAKUNITS= line missing"
        write(6,*)" ERROR #DELTAKUNITS= line missing"
        stop
      endif
      if(dadnunits.eq." ")then
        write(0,*)" ERROR #DADNUNITS= line missing"
        write(6,*)" ERROR #DADNUNITS= line missing"
        stop
      endif

C     Create the table.
      write(6,1006)
 1006 format("#Conversions:"/
     & "# 1 ksi*sqrt(inch) = 1.0989*MPa*sqrt(m)"/
     & "# 1 MPa*sqrt(m)  =   31.6228*N/(mm**(3/2) )"/
     & "# 1 MPa*sqrt(mm) =         1 N/(mm**(3/2) )"
     & //
     & "#All inputs converted to    MPa*sqrt(mm)   and  mm/cycle"/
     & "#Note that this is same as  N/(mm**(3/2))  and  mm/cycle"
     & /"#deltaKunits=  mpa_mm"/"#dadnUnits= mm"
     & /)


      if(ndata.eq.0)then
        write(0,2030)
 2030   format("#Error: convert2MPa_mm: No data lines in input file.")
      endif

      write(6,1008)
 1008 format("#   MPa*sqrt(mm)  mm/Cycle")
      do 2012 i=1,ndata
        write(6,2016)dk2(i),dadn2(i)
 2016   format(2(e14.7,1x))
 2012 continue

      stop
      end

