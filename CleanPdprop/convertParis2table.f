C   convertParis2table.f   vers. 1.2   mar 07 2013   FAC
      SAVE
C     make Paris Eq. into da/dn vs. Delta_K  table
C     Compile:  gfortran -g -w -fbounds-check convertParis2table.f -o convertParis2table
C     Usage:    convertParis2table  <inputfile.paris   >dadnFileName.dadn
C
C  Copyright (C) 2012 Fatigue Des. and Eval. Comm. and F.A.Conle
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

C vers 1.2  is a fork from mkparis.f  mainly for a name change but also
C           due to splitup of #UNITS=  to #DELTAKUNITS= and #DADNUNITS=.  march2013
C vers 1.1  Set output format to be   N/(mm**(3/2) ) which is same as MPa*sqrt(mm)

C  Input file format:
C    ______ first column in file
C   |
C   v
C   #NAME= G40.21-50A  
C   #deltaKunits= ksi_in     # specify:   mpa_m   or   mpa_mm   or   ksi_in
C   #dadnunits=    inch      # specify:   meter   or    mm      or   inch
C   #dKth=  18.0       # threshold stress intensity value
C   #dKc=  100.0       # Upper cut-off or Paris equation, or KIc
C   #A= 5.0E-09        # Paris co-efficient
C   #m= 1.90           # Paris exponent
C                      # Use file "g40.21-50A.paris"   as an example file
C---------------------------------------------------------------------------
C Output file format
C   #
C   # convertParis2table  vers. 1.2 starts...
C   #Name= G40.21-50A                    
C   #Got Original #DELTAKUNITS=                     ksi_in    
C   #Got Original #DADNUNITS=                       inch
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
C   #DELTAKUNITS= mpa_mm
C   #DADNUNITS= mm
C   #dKth=   626.
C   #dKc=    3475.
C   #   MPa*sqrt(mm)   mm/cycle
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

      character*10 deltaKunits, dadnunits
      character*30 name30, firstfield, ctail
      character*30 cfiletype, Cdatatype

 
      integer maxpoints,npts
      real*4  dadn0

      real*4  dk0

      real*4  A,m,dKth,dKc

      maxpoints=250

      write(6,701)
      write(0,701)
  701 format("# convertParis2table  vers. 1.2 starts...")

      XMPAS=6.894759
C     Set some values to check if user forgets inputs
      deltaKunits=" "
      dadnunits=" "
      name30= " "
      dKth=-999.
      dKc=-999.
      A=-999.
      m=-999.

C---------------------------- Read in Paris parameters from stdin-------------
      ninput=0
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
  805    continue

  810    continue
C        1st non blank is not a #.  
C        It must be a letter, not a number. Time to bomb out.
         write(0,815)ninput
         write(6,815)ninput
  815    format(" ERROR convertParis2table: input line no. ",I5,
     &   " not a # and not a number. Edit fitted input file")
         stop
      endif
C     end of check for 1st char non-#


C     Input line has a # in first col, check it just in case
  880 continue

      if(inpone(1).ne."#")then
C       something is bad in program
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

      if(firstfield .eq."#dKth=" .or.
     &   firstfield .eq."#dkth=" .or.
     &   firstfield .eq."#DKTH=" .or.
     &   firstfield .eq."#DKth=" )then
         read(inp300,*) firstfield, dKth
         write(6,*)"#Got Original #dKth= ",dKth
         go to 950
      endif
      if(firstfield .eq."#dKc=" .or.
     &   firstfield .eq."#dkc=" .or.
     &   firstfield .eq."#DKC=" .or.
     &   firstfield .eq."#DKc=" )then
         read(inp300,*) firstfield, dKc
         write(6,*)"#Got Original #dKc= ",dKc
         go to 950
      endif
      if(firstfield .eq."#A=" .or.
     &   firstfield .eq."#a=" )then
         read(inp300,*) firstfield, A
         write(0,*)"#Got Original  #A= ",A
         write(6,*)"#Got Original  #A= ",A
         go to 950
      endif
      if(firstfield .eq."#m=" .or.
     &   firstfield .eq."#M=" )then
         read(inp300,*) firstfield, m
         write(0,*)"#Got Original  #m= ",m
         write(6,*)"#Got Original  #m= ",m
         go to 950
      endif

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

  980 continue
C     All input lines have been read.
C     See if some of the critical values are missing:
      if(A.eq.-999.)then
        write(0,*)" ERROR #A= line missing"
        write(6,*)" ERROR #A= line missing"
        stop
      endif
      if(m.eq.-999.)then
        write(0,*)" ERROR #m= line missing"
        write(6,*)" ERROR #m= line missing"
        stop
      endif
      if(dKth.eq.-999.)then
        write(0,*)" ERROR #dKth= line missing"
        write(6,*)" ERROR #dKth= line missing"
        stop
      endif
      if(dKc.eq.-999.)then
        write(0,*)" ERROR #dKc= line missing"
        write(6,*)" ERROR #dKc= line missing"
        stop
      endif
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
         write(6,1002)
 1002    format("#DELTAKUNITS= mpa_mm")
         write(6,1001)
 1001    format("#DADNUNITS= mm")

C     Compute each end of the straight line Paris eq. in the 
C     original units.
C     Lowest possible (?) point first
C     Assume that Kth is close to zero,  note that it cannot be Zero
C       in the Paris equation.
      dk0=dKth
      if(dKth .eq. 0.0)then
        dk0= 0.1
        write(0,1010)dk0
        write(6,1010)dk0
 1010   format("#Warning: your value of dKth= 0. has been reset to ",
     &          e14.7," to create the 1st Table point")
       endif

      dadn0= A*(dk0**m )
      write(0,*)" #Got: dk0= ",dk0," dadn0= ",dadn0

C     Highest possible point next:   KIc
      if(dKc .le. dKth)then
        write(0,1020)dKc,dKth
        write(6,1020)dKc,dKth
 1020   format("#Error: convertParis2table: value of dKc <= dKth,  ",
     &         "dKc=",e14.7,"  dKth= ",e14.7)
        stop
      endif
      dk1=dKc
      dadn1= A*(dKc**(m) )
      write(0,*)" #Got: dk1= ",dk1," dadn1= ",dadn1


C     Create the table.  Only the 2 end points are in the table.
C     Convert units if necessary.
      write(6,1006)
 1006 format("#Conversions:"/
     & "# 1 ksi*sqrt(inch) = 1.0989*MPa*sqrt(m)"/
     & "# 1 MPa*sqrt(m)  =   31.6228*N/(mm**(3/2) )"/
     & "# 1 MPa*sqrt(mm) =         1 N/(mm**(3/2) )"
     & //
     & "#All inputs converted to    MPa*sqrt(mm)   and  mm/cycle"/
     & "#Note that this is same as  N/(mm**(3/2))  and  mm/cycle"
     & /)
      if(deltaKunits .eq. "ksi_in")then
C        Original units in ksi-in, new will be mpa-mm
         dKthnew= dKth*1.0989 * 31.6228
         write(6,1003)dKthnew
 1003    format("#dKth= ",f6.0)
         dKcnew= dKc*1.0989*31.6228
         write(6,1004)dKcnew
 1004    format("#dKc= ",f8.0)
         write(6,1007)
 1007    format("# MPa*sqrt(mm)  mm/cycle")
         dk0= dk0 * 1.0989 * 31.6228
         dk1= dk1 * 1.0989 * 31.6228
      endif

      if(deltaKunits .eq. "mpa_m")then
         write(6,1002)
         dKthnew= dKth*31.6228
         write(6,1003)dKthnew
         dKcnew= dKc*31.6228
         write(6,1004)dKcnew
         write(6,1007)
         dk0= dk0 * 31.6228
         dk1= dk1 * 31.6228
      endif

      if(deltaKunits .eq. "mpa_mm")then
C        Original units are same as final
         write(6,1002)
         write(6,1003)dKth
         write(6,1004)dKc
         write(6,1007)
      endif



      if(dadnunits .eq. "inch")then
         dadn0= dadn0*25.4
         dadn1= dadn1*25.4
      endif

      if(dadnunits .eq. "meter")then
         dadn0= dadn0*1000.
         dadn1= dadn1*1000.
      endif

C      if(dadnunits .eq. "mm") do nothing

C     Write out the 2 point table
      write(6,*)dk0,dadn0
      write(6,*)dk1,dadn1

      stop
      end

