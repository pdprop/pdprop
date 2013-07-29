C getFADs.f  ver. 0.5  read in the material tensile curve file and create the
C                      FAD file of Kr vs Lr 
      SAVE
C Linux Compile:  gfortran -g -w -fbounds-check   getFADs.f  -o getFADs

C Usage: getFADs  <matl_tensile_file  >transform_FADs_file
C---------------------------------------------------------------------------
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

C vers0.5 Adds code for elimination of upper yield and estimating
C         Luder's plateau size. (noise window = 5 ksi)  Jun22 2013

C Assumptions and other items of interest:

C The material input file is assumed to be a cleaned up tensile test file 
C  i.e.: the computer will assume that a line could be drawn through the 
C  successive input data points and that this line would make a reasonable 
C  curve, reasonalble enough to allow one to compute values between the points. 
C NOTE!: For  BS7910 FADs required Tags are:
C                #FileType= strain_stress
C                #DataType= engineering    (or "true" ) 
C                #Sy=   
C                #Su=   
C                #E= 
C                #Stress_units= mpa   (or ksi)
C                #Strain_units= strain      (no other options at this time)
C  will appear in the input stream.  If not the program should error out.

C  An example file:
C   #
C   #FileType= strain_stress
C   #DataType= engineering     #Can be "engineering" or "true"
C   #NAME= ASTM-A36
C   #NAME= HotRolled
C   #NAME= Steel
C   #Stress_units= mpa
C   #Strain_units= strain
C   #Sy= 224. mpa = 32.5 ksi
C   #Su= 414. mpa = 60.0 ksi
C   #eu= 0.14  #strain at Su
C   #E= 190000 mpa  27500 ksi
C   #FractureStrain=  1.19
C   #FractureStress=  952 mpa = 138 ksi
C   #monotonic_K= 780 mpa, 113 ksi
C   #monotonic_n= 0.258
C   #BHN= 0. not reported
C   #%RA= 69.7 %
C   
C   #From initial MattosStress-Strain plot
C   0  0
C   0.00146  285
C   0.00206  220
C   0.00274  225
C    ....etc.


C  If you are changing this routine please consider that sooner or later
C  someone will send a different type of file by mistake.


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

      character*10 name1, name2, stressunits
      character*10 strainunits 
      character*30 names30(10), firstfield, ctail
      character*30 cfiletype, Cdatatype
      character*80 argv, fname

      integer*4 iargc, argc

C     Save the digital curves here:
      real xStrain(250),xStress(250)  ! only used for read
      real EngStress(250),EngStrain(250)
      real TrueStress(250),TrueStrain(250)
      integer ndata, nxyFAD1, nxyFAD2a, nxyFAD2b, iFADtype
      real xLrmax,yKrmax
      real xLr,yKr2a, xKr2b
      

      idim=250    !  = max dimension of digital curve stuff

      write(0,100)
      write(6,100)
  100 format("#getFADs.f   vers 1.5 starts...")

      ndata=0
      nxyFAD1=0
      nxyFAD2a=0
      nxyFAD2b=0
      iFADtype=0   ! could be  1, 2, -2,  or 3

      XMPAS=6.894759
C     Set some default values in case user forgets
      stressunits="none"
      strainunits="none"
      EMOD=0.
      Cdatatype="none"
      cfiletype="none"
      Sult=0.
      Syield=0.
      FractureStress=0.
      FractureStrain=0.
C---------------------------- Read in monotonic tensile file from stdin-------------
  800 continue
c     Loop back to here for next input line.

C     Input lines may either be 1. data lines,  2. #Comment lines or 3. blank
C        It will be hard to distinguish the real comment from the junk comment.
C        Leave it up to the user to edit the comment section.
C     Special comment lines begin with "#string= value"   tags.  They are
C        special.

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
C           The first non blank is not a #
             go to 810
           endif
  805    continue

      
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
C           It must be a letter, not a number. Time to bomb out.
            write(0,815)ninput
            write(6,815)ninput
  815       format(" ERROR: getFADs: input line no. ",I5,
     &      " not a # and not a number. Edit fitted input file")
            stop
         endif


C        Ok, it should be like this one:
C        0.0125          180   279.   .0       0.0030        30100.    #specimen comment
C        # Total Strain   2Nf  Stress  Mean   Plastic Strain   Initial
C        #    Amp               Amp   Stress      Amp        Elastic Mod.
C        Note that the trailing comment field "#specimen comment" may or may not be there.

C        We should be able to read the first 5 fields as real numbers.  Since 2Nf may have
C        more than 7 signif. digits, it needs to be read in as a double precision.  It may
C        be an integer.  The value can be changed later to whatever the local format is.

C      Since this is all "Fitted" data expressing a curve, we only need the first 3 items.

         ndata=ndata+1
         if(ndata.gt.idim)then
           write(0,816)
           write(6,816)
  816      format(" Error too many data points:",i5,
     &            " recompile getFADs.f")
           stop
         endif
C        Place into temporary storage
         read(inp300,*)xStrain(ndata), xStress(ndata)

C          If we have not crashed by here, the data from the line has been read in. 
C          Now figure out if there is comment at the end of line.
C          In getFADs we don't really need to do this, but the code is here for
C          some possible future use. ? :)
C          Brute force it. Hunt for a #

           do 820 i=1,300
             if(inpone(i).eq."#")then
C               we found it in col i
Cdebug                write(0,*)"found # trailer in data line"
                loc=i
                go to 822
             endif
  820      continue
C          If we got to here, the trailing field is empty
           go to 855

  822      continue
C          Now count backwards and see where the last char is 
           do 830 i=1,300
             j=300-(i-1)
             if(inpone(j).ne." ")then
C              found last char
               lastloc=j
               go to 832
             endif
  830      continue

  832      continue
C          It is possible that the #field is an html tag, so it may be real long
C          Move the stuff to begin of a new field for decoding
           j=0
           do 840 i=loc,lastloc
             j=j+1
             jnpone(j)=inpone(i)
  840      continue
Cdebug           write(0,*)"Trailer # on data: ",(jnpone(i),i=1,j)

C          Check to see if any special tags are in this field
C          If the field is blank, then what?
           read(jnp300,*)ctail
C           if(ctail.eq."#Runout" .or.
C     &        ctail.eq."#RUNOUT" .or.
C     &        ctail.eq."#runout" )then
C            The data point was a runout
C            In the local company's format, runouts are -ve nos.
C            Change it to whatever you like !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C             longlife=-longlife
C           endif
C          #knifeEdge and #outsideGage  are also known tags
C           but we dont do anything with them right now.
C           So just write out, and say which point.
Cdebug           write(6,850)ndata,(jnpone(i),i=1,j)
  850      format(" Pt No.", I5, " : ", 300a1)

C          Now save the data in the local program fields
  855    continue


C       End of data line processing, write it out later, when all are in.
        go to 800
       endif

C     Input line has a # in first col, check it just in case
  880 continue

      if(inpone(1).ne."#")then
C       something is bad in program
       write(0,*)" ERROR 880, sorry prog. messed up call ? admin"
        stop
      endif
C     Ok, its a nice comment.  Figure out if its a special tag.
      read(inp300,*)firstfield
      if(firstfield .eq."#FileType=" .or.
     &   firstfield .eq."#filetype=" .or.
     &   firstfield .eq."#FILETYPE=" )then
         read(inp300,*) firstfield, cfiletype
         if(cfiletype.ne."strain_stress")then
C          We have the wrong kind of sae file here folks
           write(0,*)" WARNING: wrong type of SAE std file. ",
     &               " Not #FileType= strain_stress"
         endif
C        filetype is ok.  Write it out as a comment
         go to 950
      endif

      if(firstfield .eq."#DataType=" .or.
     &   firstfield .eq."#DATATYPE=" .or.
     &   firstfield .eq."#datatype=" )then
         read(inp300,*) firstfield, Cdatatype
C        correct caps input:
         if(Cdatatype .eq. "TRUE" .or.
     &      Cdatatype .eq. "True")then
            Cdatatype = "true"
         endif
         if(Cdatatype .eq. "ENGINEERING" .or.
     &      Cdatatype .eq. "Engineering")then
            Cdatatype = "engineering"
         endif
         
         if(Cdatatype .ne. "engineering" .and.
     &      Cdatatype .ne. "true")then
C          We have the wrong kind of sae file here folks
           write(0,*)" ERROR, wrong type of data file. ",
     &               " Not #DataType= engineering or true"
           stop
         endif
C        DataType is ok.  Write it out as a comment
         go to 950
      endif

      if(firstfield .eq."#NAME=" .or.
     &   firstfield .eq."#Name=" .or.
     &   firstfield .eq."#name=" )then
         numnames=numnames+1
         read(inp300,*) firstfield, names30(numnames)
         go to 950
      endif

      if(firstfield .eq."#UNITS=" .or.
     &   firstfield .eq."#units=" .or.
     &   firstfield .eq."#Units=" .or.
     &   firstfield .eq."#Stress_units=" .or.
     &   firstfield .eq."#STRESS_UNITS=" .or.
     &   firstfield .eq."#Stress_Units=" .or.
     &   firstfield .eq."#stress_units=" )then
         read(inp300,*) firstfield, stressunits
         go to 950
      endif

      if(firstfield .eq."#STRAIN_UNITS=" .or.
     &   firstfield .eq."#Strain_units=" .or.
     &   firstfield .eq."#STRAIN_units=" .or.
     &   firstfield .eq."#Strain_Units=" .or.
     &   firstfield .eq."#strain_units=" )then
         read(inp300,*) firstfield, strainunits
         go to 950
      endif

      if(firstfield .eq."#Su=" .or.
     &   firstfield .eq."#SU=" )then
         read(inp300,*) firstfield, Sult
         go to 950
      endif
       if(firstfield .eq."#Sy=" .or.
     &    firstfield .eq."#SY=" )then
          read(inp300,*) firstfield, Syield
          go to 950
      endif
       if(firstfield .eq."#E=" .or. firstfield .eq."#e=" .or.
     &    firstfield.eq."#EMOD=" .or. firstfield .eq."#Emod=".or.
     &    firstfield .eq."#emod=" .or.
     &    firstfield .eq."#MODULUS=" .or.
     &    firstfield .eq."#Modulus=" .or.
     &    firstfield .eq."#modulus="  )then
          read(inp300,*) firstfield, EMOD
          write(6,*)"# Found EMOD =", EMOD
          go to 950
      endif
  
      if(firstfield .eq."#FractureStress=" .or.
     &   firstfield .eq."#fracturestress=" .or.
     &   firstfield .eq."#FRACTURESTRESS=" .or.
     &   firstfield .eq."#Fracturestress=" )then
         read(inp300,*) firstfield, FractureStress
         go to 950
      endif

      if(firstfield .eq."#FractureStrain=" .or.
     &   firstfield .eq."#Fracturestrain=" .or.
     &   firstfield .eq."#fracturestrain=" .or.
     &   firstfield .eq."#FRACTURESTRAIN=" )then
         read(inp300,*) firstfield, FractureStrain
         go to 950
      endif

C       Look for things we should skip over in the output
        write(INP80temp,939)
  939   format("#Here is the bottom part of the ",
     &         "html  graph/calc wrapper-----------")
        if(INP300.eq. INP80temp)go to 800

        if(INP300.eq.
     &"#</textarea>"
     &  )go to 800

        if(INP300.eq.
     &"#</pre></DL></FORM></body></html>"
     &  )go to 800



CC     It was NOT a special comment. Save it. If it was a # blank line, discard it.
CC     Count backwards  to find last non-blank char of string
C      do 940 i=1,300
C        j=300-(i-1)
C        if(inpone(j).ne." ") then
CC          found last char
C           lastloc=j
C           go to 943
C        endif
C  940 continue
C
CC     We also need to save the no of chars for each comm line
CC     If it was a #blank forget it
C  943 continue
C      if(lastloc.eq.1)go to 800
CC     Look for some other dumb lines and discard them
C      read(inp300,*) dummy
C      if(dummy.eq."#saeinput")go to 800
C      if(dummy.eq."#---------")go to 800
C
C      numcomment=numcomment+1
C      numcomchars(numcomment)=lastloc
CC     save the first 80 chars
C      read(inp300,"(a80)")comment80(numcomment)
C      go to 800



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

C-------------------------------------------------read  is done-------------
  980 continue
C     All input lines have been read.  The comments were
C     put out along the way.  

C     See if some of the critical values are missing:
      if(EMOD.eq.0.)then
        write(0,*)" ERROR #EMOD= missing value"
        write(6,*)" ERROR #EMOD= missing value"
        stop
      endif
      if(Cfiletype.eq."none")then
        write(0,*)" ERROR #FileType= missing value"
        write(6,*)" ERROR #FileType= missing value"
        stop
      endif
      if(Cdatatype.eq."none")then
        write(0,*)" ERROR #DataType= missing value"
        write(6,*)" ERROR #DataType= missing value"
        stop
      endif
      if(stressunits.eq."none")then
        write(0,*)" ERROR #stressunits= missing value"
        write(6,*)" ERROR #stressunits= missing value"
        stop
      endif
      if(strainunits.eq."none")then
        write(0,*)" ERROR #strainunits= missing value"
        write(6,*)" ERROR #strainunits= missing value"
        stop
      endif
      if(Syield.eq.0.)then
        write(0,*)" ERROR: #Sy= missing value"
        write(6,*)" ERROR: #Sy= missing value"
        stop
      endif
      if(Sult.eq.0.)then
        write(0,*)" ERROR: #Su= missing value"
        write(6,*)" ERROR: #Su= missing value"
        stop
      endif

      write(0,*)"#CHECK THESE Inputs !:"
      write(0,*)"#Stress_units=",stressunits," Strain_units=",
     &          strainunits
      write(0,*)"#EMOD=",EMOD, " Sult=",Sult," Syield=",Syield
      write(6,*)"#CHECK THESE Inputs:"
      write(6,*)"#Stress_units=",stressunits," Strain_units=",
     &          strainunits
      write(6,*)"#EMOD=",EMOD, " Sult=",Sult," Syield=",Syield

C     Ok,  the critical Tags have all been defined. Lets move the
C     temp data into the correct storage
      if(stressunits .eq."mpa" .or.
     &   stressunits .eq."MPA" .or.
     &   stressunits .eq."MPa" .or.
     &   stressunits .eq."Mpa" )then
         stressunits="MPa"
         do 982 i=1,ndata
           if(Cdatatype .eq. "engineering")then
              EngStress(i)=xStress(i)
              EngStrain(i)=xStrain(i)
              TrueStrain(i)=log( 1.0+EngStrain(i) )
              TrueStress(i)=EngStress(i)*( 1+EngStrain(i) )
           endif
           if(Cdatatype .eq. "true")then
              TrueStrain(i)=xStrain(i)
              TrueStress(i)=xStress(i)
C             And probably not used, but:
              EngStrain(i)=  exp(xStrain(i)) -1.0
              EngStress(i)= TrueStress(i)/(1.0+EngStrain(i) )
           endif
  982    continue
         go to 1010
      endif
      
      if(stressunits .eq."ksi" .or.
     &   stressunits .eq."KSI" .or.
     &   stressunits .eq."Ksi" )then
C        We really do not need to change units?? What if loads
C        are in mpa and material in ksi?? hm, Ok make everything mpa
         EMOD=EMOD*XMPAS
         Sult=Sult*XMPAS
         Syield=Syield*XMPAS
         do 983 i=1,ndata
           if(Cdatatype .eq. "engineering")then
              EngStress(i)=xStress(i)*XMPAS  ! also convert to mpa
              EngStrain(i)=xStrain(i)
              TrueStrain(i)=log( 1.0+EngStrain(i) )
              TrueStress(i)=EngStress(i)*( 1+EngStrain(i) )
           endif
           if(Cdatatype .eq. "true")then
              TrueStrain(i)=xStrain(i)
              TrueStress(i)=xStress(i)*XMPAS
C             And probably not used, but:
              EngStrain(i)=  exp(xStrain(i)) -1.0
              EngStress(i)= TrueStress(i)/(1.0+EngStrain(i) )
           endif
  983    continue
         write(6,*)"#Tensile file ksi -> MPa."
         write(0,*)"#Tensile file ksi -> MPa."
         stressunits="MPa"
         go to 1010
      endif
      
      if(stressunits .eq."psi" .or.
     &   stressunits .eq."PSI" .or.
     &   stressunits .eq."Psi" )then
         xtemp =XMPAS/1000.
         EMOD=EMOD*xtemp
         Sult=Sult*xtemp
         Syield=Syield*xtemp
         FractureStress=FractureStress*xtemp
         do 985 i=1,ndata
           if(Cdatatype .eq. "engineering")then
              EngStress(i)=xStress(i)*xtemp  ! also convert to mpa
              EngStrain(i)=xStrain(i)
              TrueStrain(i)=log( 1.0+EngStrain(i) )
              TrueStress(i)=EngStress(i)*( 1+EngStrain(i) )
           endif
           if(Cdatatype .eq. "true")then
              TrueStrain(i)=xStrain(i)
              TrueStress(i)=xStress(i)*xtemp
C             And probably not used, but:
              EngStrain(i)=  exp(xStrain(i)) -1.0
              EngStress(i)= TrueStress(i)/(1.0+EngStrain(i) )
           endif
  985    continue
         write(6,*)"#Material file: psi -> MPa."
         write(0,*)"#Material file: psi -> MPa."
         stressunits="MPa"
      endif
 1010 continue

      if(strainunits .eq."Microstrain" .or.
     &   strainunits .eq."MICROSTRAIN" .or.
     &   strainunits .eq."MicroStrain" )then
         FractureStrain=FractureStrain/1.0E6
         do 1020 i=1,ndata
           TrueStrain(i)=TrueStrain(i)/1.0E+6
C          and for consistency:
           EngStrain(i)=EngStrain(i)/1.0E+6
           xStrain(i)=xStrain(i)/1.0E+6
 1020   continue
         write(6,*)"#Material file Microstrain -> strain"
         write(0,*)"#Material file Microstrain -> strain"
         strainunits="strain"
      endif

      if(strainunits .eq."Percent" .or.
     &   strainunits .eq."PERCENT" .or.
     &   strainunits .eq."percent" )then
         FractureStrain=FractureStrain/1.0E6
         do 1021 i=1,ndata
           TrueStrain(i)=TrueStrain(i)/1.0E+6
C          and for consistency:
           EngStrain(i)=EngStrain(i)/1.0E+6
           xStrain(i)=xStrain(i)/1.0E+6
 1021   continue
         write(6,*)"#Material file Microstrain -> strain"
         write(0,*)"#Material file Microstrain -> strain"
         strainunits="strain"
      endif



C     Now sort the stress,strain,  values  ------------------
C     (this is silly. We should not need to do this. -but if
C      someone has dumped experimental data with noise, then strains
C      may not be monotonically increasing.)  Thus:
      do 1090 i=1,ndata
        ifind=0
        smallstrain= +999.

        do 1080 j=i,ndata
         if(smallstrain.lt.TrueStrain(j)) go to 1080
C        Nope, its smaller
         smallstrain=TrueStrain(j)
         ifind=j
 1080   continue

Cdebug       write(6,*)" Sorting ",ifind, " into ",i
C      smallstrain is now the smallest & is at location j
C      Save it in a temp location
       tempEngStrain=EngStrain(ifind)
       tempEngStress=EngStress(ifind)
       tempTrueStrain=TrueStrain(ifind)
       tempTrueStress=TrueStress(ifind)
       tempxStrain=xStrain(ifind)
       tempxStress=xStress(ifind)

C      Now move the stuff pointed to by the big loop
C      into "bigstrain's" location
       EngStrain(ifind)  =EngStrain(i)
       EngStress(ifind)  =EngStress(i)
       TrueStrain(ifind)  =TrueStrain(i)
       TrueStress(ifind)  =TrueStress(i)
       xStrain(ifind)  =xStrain(i)
       xStress(ifind)  =xStress(i)

C      Put the temp stuff into the "i" location
       EngStrain(i)=tempEngStrain
       EngStress(i)=tempEngStress
       TrueStrain(i)=tempTrueStrain
       TrueStress(i)=tempTrueStress
       xStrain(i)=tempxStrain
       xStress(i)=tempxStress
 1090 continue

C     Check for an upper yield in the tensile curve:
C     Assume that an Upper yield would occur before the  0.2% offset
C     or Lower yield = Syield
C     Some of the initial points in the curve may be on the "elastic"
C     modulus line .   (before we get to the Upper yield.
C     Thus any points before the Syield that are below  Syield can
C     be ignored as valid.  If we can assume that Syield is at the
C     0.2% offset strain, then we can stop the hunt at 
C           TrueStrain() > Syield/Emod +0.002
      strainAtYield=  Syield/Emod +0.002
      write(6,1104)strainAtYield, Syield
 1104 format("# Checking for Upper Yld: Strain at yield= ",f8.5,
     &       " Syield= ",f8.1)
      isUpper=0 ! will be set to 1  if there is an upper yield detected
      do 1110 i=1,ndata
Cdebug        write(6,*)"#checking: ",TrueStrain(i),strainAtYield
        if(TrueStrain(i) .gt. strainAtYield )goto 1115 ! done
C         Yes, its less than the yield strain, check for an Sy upper
          if(TrueStress(i) .gt. Syield)then
C           Yes this is an "upper" or "rolling" yield
C           ("rolling" upper yields occur in some  HSLAs)
            isUpper=1
            write(0,1108)i,TrueStrain(i),TrueStress(i),
     &                   TrueStrain(i),Syield
            write(6,1108)i,TrueStrain(i),TrueStress(i),
     &                   TrueStrain(i),Syield
 1108       format("#Warning!: Upper yield detected and changed:"/
     &      "#Old point: ",i4,1x,f8.5,1x,f5.0," mpa"/
     &      "#Changed to:       ",f8.5,1x,f5.0," mpa"/)
            TrueStress(i)=Syield  !  Set it back to Lower yield
          endif
C       No its just a lower stress strain point
 1110 continue

      if(isUpper .eq. 0)then
        write(0,1112)
        write(6,1112)
 1112   format("# No Upper yield points detected before 0.2% Syield")
        write(6,1113)
 1113   format("#FADinput #Lueders= no")
      endif

 1115 continue
      if(isUpper .eq. 1)then
C       We had an upper yield detected, so lets use BS7910 to estimate
C       length of the Luder's region:
        strainLuders= 0.0375*(1.0-Syield/1000.)
        write(0,1118)strainLuders
        write(6,1118)strainLuders
 1118   format("#BS7910 Lueders Strain:  eL = ",f8.5)
        write(6,1119)
 1119   format("#FADinput #Lueders= yes")
        slopeN= 0.30*(1.0-Syield/Sult)  !used down below in FAD2a
        write(6,1121)slopeN
 1121   format("#FADinput #slopeN= ",e14.7)
      endif

C     Check for a Lueders region (data may not have an upper yield)
C     HSLA  Lrmax  = 1.15 typical.  
C     Thus if Lrmax is less than about 1.2  we probably do not need
C     to check for a Lueders.  HSLA  are very flat.  Lueders are 
C     (usually) followed by a strain hardening region
      xLrmax= (Syield+Sult)/(2.0*Syield)
      if(xLrmax .le. 1.20)then
        write(0,1130)xLrmax
        write(6,1130)xLrmax
 1130   format("#Lrmax= ",f5.2," < 1.20 indicates an",
     &         "  HSLA type material. No Lueders check needed.")
        goto 1170
      endif

C     BS7910 estimates the length of the Luders plateau with the formula
      if(isUpper .eq. 0)then
         write(0,1135)xLrmax
         write(6,1135)xLrmax
 1135    format("#Lrmax= ",f5.2," > 1.2 indicates a STRAIN HARDENING",/
     &         "# material.  No upper Syield was indicated, but curve"/
     &         "# could be Flat at the yield. Checking for Lueders...:")
         estimateLuders=  0.0375*(1.0 - (Syield/1000.) )
         write(0,1138)estimateLuders
         write(6,1138)estimateLuders
 1138    format("#BS7910 Lueders length estimate epsL= ",
     &          f8.5," strain")

C        Lets start at the strain at 0.2% offset yield
C        strainAtYield was computed above
         Snoise= 5.0*6.895  !  5ksi =  34.225 mpa  Detect Window size
         write(0,1139)Snoise
         write(6,1139)Snoise
 1139    format("#Lueders stress noise window  Snoise= ",f8.2," mpa")
         endLuders= strainAtYield
         sizeLuders= 0.
         do 1150 i=1,ndata
           if(TrueStress(i) .gt. (Syield + Snoise) )then 
              write(0,1141)i,TrueStrain(i),TrueStress(i)
              write(6,1141)i,TrueStrain(i),TrueStress(i)
 1141         format("# Lueders window exceeded at pt ",i4,
     &         "  trueStrain= ",f8.5," trueStress= ",f8.0," mpa")
C             We have risen outside the noise window.
C             Assume that we are now in the strain hardening region
C             Thus the previous SS point is the last in the Luders
              endLuders=TrueStrain(i-1)
              sizeLuders= endLuders-strainAtYield
              write(0,1142)strainAtYield,Syield,
     &                     endLuders,TrueStress(i-1),
     &                     sizeLuders,estimateLuders
              write(6,1142)strainAtYield,Syield,
     &                     endLuders,TrueStress(i-1),
     &                     sizeLuders,estimateLuders
 1142         format("# Check for Lueders plateau length:"/
     &        "# start Lueders= ",f8.5,1x,f8.0/
     &        "# end   Lueders= ",f8.5,1x,f8.0/
     &        "# length= ",f8.5,", compares with BS7910= ",f8.5)
C             Ok,  if the length is bigger than 0.005 we will assume
C             That it is a Lueders   material.
              if(sizeLuders .ge. 0.005 )then
                 write(0,1144)sizeLuders
                 write(6,1144)sizeLuders
 1144            format("# Lueders length= ",f8.5," is bigger than ",
     &                  "0.005 strain"/
     &                  "# conclusion: "/
     &                  "#FADinput #Lueders= yes")
                 slopeN= 0.30*(1.0-Syield/Sult)  !used down below in FAD2a
                 write(6,1145)slopeN
 1145            format("#FADinput #slopeN= ",e14.7)
                 isUpper=1  ! indicates Lueders
C                 strainLuders=estimateLuders ! BS7910 approx.
                 strainLuders=sizeLuders     !  observed
                 write(0,1146)strainLuders
 1146            format("# We will use Lueders length epsL=",f8.5)
                 slopeN= 0.30*(1.0-Syield/Sult)  !used down below in FAD2a
                 write(6,1121)slopeN
                 goto 1170
              endif
           endif

 1150      continue
C          if we get to here, no Lueders plateau was found
           write(0,1152)sizeLuders
           write(6,1152)sizeLuders
 1152      format("# No Lueders larger than 0.005 was found"/
     &       "# largest flat region was= ",f8.5," strain")
      endif ! If was for "if(isUpper .eq. 0) "


C     All is well with input data (hopefully). 
 1170 continue

      write(6,1172)EMOD
 1172 format("#FADinput #E= ",f8.0," mpa")
      write(6,1173)Syield
 1173 format("#FADinput #Sy= ",f5.0," mpa")
      write(6,1174)Sult
 1174 format("#FADinput #Su= ",f6.0," mpa")

      write(6,1176)
 1176 format("#FADinput TrueStrain TrueStress_mpa")

      do 1180 i=1,ndata
      write(6,1178) TrueStrain(i), TrueStress(i)
 1178 format("#FADinput ",f8.5,1x,f6.1)
 1180 continue
      

C analysis  for the FADs   ------------------------------------------------ 

C     FAD 1   (boundaries are not material dependent)
C            but in the actual analysis we will need  Kmat and Syield
      write(6,1202)
 1202 format("#FAD1= 1"/"#Krmax= 0.707   #FAD1"/
     &       "#Srmax= 0.8     #FAD1"/
     &       "0.0 0.707       #FAD1"/
     &       "0.8 0.707       #FAD1"/
     &       "0.8 0.0         #FAD1"/)


C     FAD 2a   Assume only yield and ultimate are known
      nxyFAD2a= 25
C      write(6,1205)nxyFAD2a
C 1205 format(/"#FAD2a= ",i8,"       #FAD2a")
      xLrmax=(Syield+Sult)/(2.0*Syield)
      xinc=xLrmax/float(nxyFAD2a-1)
      xLr=0.
      yKr2a=1.0
      write(6,*)xLr,yKr2a,"       #FAD2a"

      if(isUpper .eq. 0) then ! Not a Lueders material----------------
        do 1210  i=2,24
          xLr=xLr+xinc
          yKr2a=(1.0-0.14*(xLr*xLr) ) * ( 0.3 + 
     &          0.7*( exp(-0.65*( xLr**6 ) ) ) )
          write(6,*)xLr,yKr2a,"       #FAD2a"
          goto 1210
 1210   continue
C       After last point drop the line down to zero Kr:
        write(6,*)xLr," 0.0           #FAD2a"
      endif

      if(isUpper .eq. 1)then  ! Lueders material---------------------
C        Compute up to but not including xLr=1.0
        do 1215 i=2,24
          xLr=xLr+xinc
          if(xLr .ge. 1.0)then
            isave=i
            goto 1216
          endif
          yKr2a=(1.0-0.14*(xLr*xLr) ) * ( 0.3 + 
     &          0.7*( exp(-0.65*( xLr**6 ) ) ) )
          write(6,*)xLr,yKr2a,"       #FAD2a"
 1215   continue

 1216   continue
C       Now compute for Lr=1.0
C       The top of the drop-off is still defined by the eq. above
        zLr=1.0  !  use zLr just for visibility
        yKr2a= (1.0-0.14*(zLr*zLr) ) * ( 0.3 +
     &          0.7*( exp(-0.65*( zLr**6 ) ) ) )
          write(6,*)zLr,yKr2a,"       #FAD2a"
C           BS7910 says to cut off at Lr=1.0, or extend a bit if
C           the material has Lueders
C           Set above:  slopeN= 0.30*(1.0-Syield/Sult)
        yKrAtLrOne= (1.0 + (Emod*strainLuders/Syield) + 
     &             (1.0/(2.0*(1 + (Emod*strainLuders/Syield) )) )
     &            )**(-0.5)
C       This is the bottom of the drop-off line:
        write(6,*)"1.0",yKrAtLrOne,"       #FAD2a Lueders change here"

C       Now continue on for Lr > 1   with incs to xLr etc
        do 1220 i=isave,24
           if(xLr .eq. 1.0)goto 1220  !  this point was done above
           xLr=xLr+xinc
C          A Lueder's material, then at xLr >= 1 has special equation:
           yKr2a= yKrAtLrOne*xLr**( (sopeN-1.0)/(2.0*slopeN) )
           write(6,*)xLr,yKr2a,"       #FAD2a"
 1220   continue
C       Put in the last point at Lrmax
        write(6,*)xLr," 0.0           #FAD2a"
      endif    !end of if for isUpper .eq. 1


C--------------------------------------------------------------------
C     FAD2b   Use the monotonic strain stress points as input
      nxyFAD2b= ndata
      write(6,1235)nxyFAD2b
 1235 format(/"#FAD2b= ",i8,"      #FAD2b")
      xLrmax=(Syield+Sult)/(2.0*Syield)
      xLr=0.
      yKr2b=1.0
      write(6,*)xLr,yKr2b,"       #FAD2b"
C     Assume that the first point of the tensile data i (0,0)
      istart=2
      if(TrueStress(1).ne.0.0)istart=1

      do 1240 i=istart,ndata
        sref=TrueStress(i)
        eref=TrueStrain(i)
        xLr=sref/Syield
        if(xLr .gt. xLrmax)goto 1250 !  end of curves
        yKr2b= 1.0/sqrt( ( (Emod*eref)/(xLr*Syield) )  +
     &             ( (xLr**3 *Syield )/( 2.0*Emod*eref) )  )
        write(6,*)xLr,yKr2b,"       #FAD2b"
 1240 continue
C     After last point drop the line down to zero Kr:
      write(6,*)xLr," 0.0           #FAD2b"
      goto 9000

 1250 continue    !we have exceeded Lrmax
      xLr= xLrmax
      yKr2b= 1.0/sqrt( ( (Emod*eref)/(xLr*Syield) )  +
     &             ( (xLr**3 *Syield )/( 2.0*Emod*eref) )  )
      write(6,*)xLr,yKr2b,"       #FAD2b"
C     Now drop the curve down to Kr = 0.
      write(6,*)xLr," 0.0           #FAD2b"

 9000 continue
      stop
      end


