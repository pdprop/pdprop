C saefcalc2.f  Ver. 1.8 Prog to calculate life using sae std. fitted (updated Nov 1 2011)
      SAVE
C              fatigue data file and a sae standard rainflow file as input. (fac Oct25/04)
C Linux Compile:  g77 -g saefcalc2.f  -o saefcalc2
C              :  gfortran -g -w -fbounds-check saefcalc2.f  -o saefcalc2
C freeBSD      :  g95 -g  saefcalc2.f -o saefcalc2
C Sun Compile:  f77 -g -Bstatic saefcalc2.f  -o saefcalc2
C Usage: 
C saefcalc2  matl_fitted_file multfactor <rainflow_file  >outfile
C                 ^             ^------- is a multiplying factor to scale the rainflow file
C                 |______________________is the digital fitted fatigue curve file
C
C In the rainflow file the data lines are assumed to be:
C   Srange  Smean  N    Smax   Smin
C        Where Smax, Smin are Elastic Stresses in MPa
C   and    N_ is number of cycles. The "triples" sets can be repeated as
C   long as the argument list will allow. Only N, Smax, and Smin are actually
C   used in this program. Srange and Smean are only provided for visual assessment.

C---------------------------------------------------------------------------
C  Copyright (C) 2004 SAE Fatigue Design and Evaluation Committee
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

C             Ver.1.8 Added a "SAVE" to each SUBROUTINE to preserve data
C                     between calls  (a new gfortran "feature" ) Feb 2012
C             Ver.1.7 Fixed nbrinell format in sta            Nov 2011
C             Ver.1.6 fixes exact match interpolation bug. See Cbugfix   May2005
C saefcalc2.f Ver.1.5 fixes format size problem for large local stresses
C             and a Solaris f77 inability to read a number into a char varible
C             using a * format type. See "Cbugfix  Oct24/04
C saefcalc2.f  Ver. 1.4 fixes %dam SWT bug: See "Cbugfix     May31/04"
C Derived from GPL program:
C saefcalc1.f  Ver. 1.3 Program to calculate life using sae std. fatigue data file
C Version 1.3 fixes bug in SWT damage calc s/r (FAC May 22/04)
C Version 1.1 has additions to allow for plotting of history and damage.
C Search for "#plothist "  in code below. (FAC Feb.12/02)

C Assumptions and other items of interest:
C  Interpolation, on curves such as stress-life or strain-life is done in
C   log-log space.  Generally anything that plots "better" on log-log
C   co-ords is also interpolated there.

C The material input file is assumed to be a "fitted" file i.e.: the computer
C  will assume that a line could be drawn through the successive input data
C  points and that this line would make a reasonable curve, reasonalble 
C  enough to allow one to compute values between the points. 
C  It is expected that the tag 
C                #DataType= fitted
C  will appear in the input stream.  If not the program should error out.
C  If you are changing this routine please consider that sooner or later
C  someone will send a different type of file by mistake.
C  It is also expected that the tag:
C                #FileType= strain_life
C  will appear in the input stream.  But, in the intrests of possibly 
C  allowing other types of data files (e.g.: plastic materials), I'm not
C  certain this rule will always be enforced.

C Other important identifiers expected in the input stream along with
C defaults assumed and others allowed:
C
C   Tag        Example      Default  Others_possible......................
C  #name=  SAE1010          Unknown  anything that starts with a letter

C  #E=     29000.           None(error out)

C  #Stress_Units= MPA       MPA      ksi,   psi

C  #Strain_Units= strain    strain   microstrain

C  #Life_Units= Reversals   Reversals  Cycles

C  #Su=  119.9              none
C  #Sy=   92.0              none
C  #BHN= 243                none

C Future feature?:  Allow material file tag
C #Sort= life (default)   or  strain,  stress
C and sort the digital curves accordingly.  Life sort should have biggest
C first, Stress & strain sorts should have smallest first.
C This would accomodate stress life curves that have shorter lives at low
C stress levels. (i.e.: Cup towards left, or Cup right etc)
C Cup or Cap type data is not allowed.

C The first part of the program was adapted from saedigcurve.f.  It reads
C in the information from the material file, and decodes the input arguments.



C-----------------------------------------------------------------------------
C Eg. of SAE standard form fatigue data file:

C   ______ first column in file
C  |
C  v
C
C  # SAE Exchange File Format.   
C  # data collected from ASTM E606 axial fatigue test data.
C  # Note: The data below is not a real file
C   
C  #DataType= fitted
C  
C  #NAME= SAE1045
C  #NAME= SAE350X
C  #NAME= SAE050X
C  #Ford= 34   
C  #Stress_Units= KSI
C  #Strain_Units= strain
C  #Life_Units= reversals
C  #E=  30000.
C  #Su= 89.
C  #Sy= 50.
C  #%RA= 85.
C  #BHN= 325
C  #WebPage= http://fde.uwaterloo.ca/Fde/Materials/Steels/ASTM-A588C/g40.21-50A_non_os.html
C   
C  # Total Strain   2Nf  Stress  Mean   Plastic Strain   Initial
C  #    Amp               Amp   Stress      Amp        Elastic Mod.
C  0.0125          180   279.   .0       0.0030        30100.    #Fitted_point
C  0.0095          490   253.   .0       0.0011        29400.    #Fitted_point
C  0.0090          950   229.   .0       0.0007        29800.    #Fitted_point
C  0.0075         2260   220.   .0       0.0002        30050.    #Fitted_point
C  0.0050        38000.  149.   .0       0.0           29900.    #Fitted_point
C  0.0040       770000   119.   .0       0             30700.    #Fitted_point
C  #
C-----------------------------------------------------------------------------
C SAMPLE RAINFLOW INPUT FILE:
C  # grpltorweld.rain  Used for Weld Challeng 2
C  #Data format is:
C  #   Range      Mean         Cycles    Max        Min
C  #BEGIN DATA
C      560.        40.0             1    320.       -240.
C      452.        40.0             1    266.       -186.
C      434.        31.0             1    248.       -186.
C      397.        31.0             1    230.       -168.
C      361.        31.0             1    212.       -150.
C      361.        12.9             1    194.       -168.
C      307.        40.0             2    194.       -114.
C      307.        21.9             2    175.       -132.
C      289.        85.2             2    230.       -59.4
C-----------------------------------------------------------------------------


      CHARACTER*1   INP1(80)
      CHARACTER*5   INP5(16),JNP5(16)
      CHARACTER*10  INP10(8),JNP10(56)
      character*11  clifeout(10)

      CHARACTER*80  JNP80(7),INP80, INP80temp 
      EQUIVALENCE  (INP1(1),INP5(1),INP10(1),INP80),
     &             (JNP5(1),JNP10(1),JNP80(1))

      character*300  inp300,jnp300,Cwebpage
      character*1    inpone(300),jnpone(300)
      character*10   inpten(30)
      equivalence (inpone(1), inpten(1), inp300)
      equivalence (jnpone(1),jnp300)

      character*10 name1, name2, stressunits
      character*10 strainunits, lifeunits, sorttype
      character*30 names30(10), firstfield, ctail
      character*30 cfiletype, Cdatatype
      character*80 argv, fname

      integer uno
      integer*4 iargc, argc

C     Save the arg line values here:
      integer nargsets
      real Smaxin(250), Sminin(250), Cyclesin(250)
C     And each input will create one of these:
      real Srangein(250), StressReps(250)
      real SwatReps(250), StrainReps(250), MorReps(250),GoodReps(250),
     &     Sigmaxin(250),Sigminin(250),Epsmaxin(250),Epsminin(250),
     &    DeltaSigin(250),DeltaEpsin(250)
      real pcswatdam(250),pcstraindam(250),pcstressdam(250),
     & pcgooddam(250), pcmordam(250)

      real Morlife


C     Save the digital curves here:
C     Check dimensions in the various subroutines too    !!!!!!!!!!!!

      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
     &     SigmaxStrainAmp(250)
      real tstrain(250),tstress(250)
C          tstrain & tstress are used only for loop plotting and passing
c          from s/r getloop
      integer ndata

      Common/Material/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
     &                 ndata,FractureStress,FractureStrain



Cdidntwork to s/rs      PARAMETER ( IDIM=250 )
      idim=250
C         = max dimension of digital curve stuff


      XMPAS=6.894759
C     Set some default values in case user forgets
      stressunits="MPA"
      strainunits="strain"
      lifeunits="reversals"
      EMOD=0.
      Cdatatype=" "
      cfiletype=" "
      Sult=0.
      Syield=0.
      percentRA=0.
      Cwebpage=" "
      FractureStress=0.
      FractureStrain=0.
      ifordnumber=0
      nbrinell=0
      sorttype="life"
      Xmagfactor= 1.0
      Xmeanadd= 0.

C  Check the number of args. Must be in sets of three
C
C  Note !!!!!!!!!!!!!    for HP fortran the arg numbers must be changed by -1
C      and the iragc probably includes the "saehbook" as an arg.

       argc = iargc()
C       In saefcalc1.f rainflow data was in args. NOT in saefcalc2.f
C       thus comment this stuff out:
C        argc=argc-1
C        nargsets= argc/3
C       Make sure there are triples:
C        If ( (nargsets *3).ne.argc .or.
C     &        argc .eq. 0) then
C         WRITE(0,*)" saefcalc1:  usage ERROR"
C         WRITE(6,*)" saefcalc1:  usage ERROR"
C         write(0,*)"Usage e.g.: saefcalc1 matl_filename 250 -250 1 ",
C     &             ">OutputFile"
C         write(0,*)" Where [ Smax Smin N ] can be multiple sets."

         if(argc.ne. 2)then
           write(0,*)"#Error: saefcalc2 : wrong no. of arguments"
           write(6,*)"#Error: saefcalc2 : wrong no. of arguments"
           write(0,*)"#Usage e.g.: saefcalc2 matl_file  multFactor",
     &               " <rainfile >outfile"
         STOP
        ENDIF

C       The first arg is material file name
        jvect=1
        call getarg(jvect,fname)
        write(6,*)"#Opening material file= ",fname," as unit 10"
        open(unit=10,file=fname)

C       2nd arg should be the multiplier factor
        jvect=2
        call getarg(jvect,INP80)
        read(INP80,*,err=600)xmultf
        write(6,*)"#Found History multiply factor: ",xmultf
        write(0,*)"#Found History multiply factor: ",xmultf
        go to 610

  600   write(0,*)"#Error: Could not read multFactor:",INP80
        write(0,*)"#Usage e.g.: saefcalc2 matl_file  multFactor",
     &               " <rainfile >outfile"
        stop

  610   continue

C      Code for saefcalc1.f and saefcalc3.f
CC     It appears that there are nargsets of triples.
C      do 750 i=1,nargsets
C        jvect=jvect+1
C        call getarg(jvect,argv)
C        read (argv,*,err=760)Smaxin(i)
C        jvect=jvect+1
C        call getarg(jvect,argv)
C        read (argv,*,err=760)Sminin(i)
C        jvect=jvect+1
C        call getarg(jvect,argv)
C        read (argv,*,err=760)Cyclesin(i)
CC        Ensure that Cycles is +ve
C         if(Cyclesin(i).le.0)then
C          write(6,741)i, Cyclesin(i)
C          write(0,741)i, Cyclesin(i)
C  741     format(" ERROR: saefcalc2: Cyclesin(",i4,") bad value =",
C     &           e14.7)
C         endif

C      This section replaces the bit directly above. Here read stdin
       i=0
       do 750 j=1,idim
C         Read in the rainflow file. Expect comments with #xxxx text
         read(5,742)inp300
  742    format(a300)
C        Check for blank line. 
         if(inp300.eq." ")then
           write(6,*)
           go to 750
         endif
Cbugfix  Oct24/05 : Was: read(inp300,*)INP80  (Solaris bombs out on data line)
         read(inp300,"(a80)")INP80

         if(INP1(1).eq."#")then
C          comment line, write and skip
           uno=6
           call WR(uno,inp300)
           go to 750
         endif
           
C        Must be a data line
         i=i+1
         read(inp300,*,end=751,err=760)rangein, xmeanin,
     &            Cyclesin(i),Smaxin(i),Sminin(i)
         if(Smaxin(i).eq.0. .and. Sminin(i).eq.0. .and.
     &      Cyclesin(i).eq.0. )then
C           Data line is all zeros, old method of exit
            i=i-1
            go to 751
         endif
C       Make sure max is max
        if(Smaxin(i).lt.Sminin(i))then
          xtemp=Smaxin(i)
          Smaxin(i)=Sminin(i)
          Sminin(i)=xtemp
        endif

C       See if someone put in equal max & min
        Srangein(i)=Smaxin(i)-Sminin(i)
        if(Srangein(i).eq.0)then
          write(0,745)j
          write(6,745)j
  745     format(" Error: in rainflow file, line= ",i4," Smax=Smin")
          write(0,*)"#",inp300
          stop
        endif
C       Apply the multiplication factor from command line
        Smaxin(i)=Smaxin(i)*xmultf
        Sminin(i)=Sminin(i)*xmultf
        Srangein(i)=Smaxin(i)-Sminin(i)

        write(6,*)"#Input x MultF : Smax=",Smaxin(i)," Smin=",Sminin(i),
     &            "Srange=",Srangein(i)," Cycles=",Cyclesin(i)
        if(i.eq.1)then
          overallsmax=Smaxin(i)
          overallsmin=Sminin(i)
        else
          if(overallsmax.lt.Smaxin(i))overallsmax=Smaxin(i)
          if(overallsmin.gt.Sminin(i))overallsmin=Sminin(i)
        endif
  750 continue

C     End of rainflow input file
  751 continue
      nargsets=i
      go to 770

  760 write(0,*)" saefcalc2: Bad Smax,Smin,N. =", inp300
      write(6,*)" saefcalc2: Bad Smax,Smin,N. =", inp300
      stop

  770 continue
C     All the comand line triplets are in. Now is a good time to
C     sort them by range size.
      if(nargsets.eq.1)go to 791
      do 790 i=1,nargsets
        ifind=0
        bigrange=0.

        do 780 j=i,nargsets
         if(bigrange.gt.Srangein(j)) go to 780
C        Nope, its bigger
         bigrange=Srangein(j)
         ifind=j
  780  continue

C      bigrange is now the biggest & is at location j
C      Save it in a temp location
       tsmax=Smaxin(ifind)
       tsmin=Sminin(ifind)
       tcyc= Cyclesin(ifind)
       trange=Srangein(ifind)

C      Now move the stuff pointed to by the big loop
C      into "bigrange's" location
       Smaxin(ifind)  =Smaxin(i)
       Sminin(ifind)  =Sminin(i)
       Cyclesin(ifind)=Cyclesin(i)
       Srangein(ifind)=Srangein(i)

C      Put the temp stuff into the "i" location
       Smaxin(i)=tsmax
       Sminin(i)=tsmin
       Srangein(i)=trange
       Cyclesin(i)=tcyc
  790 continue
 
  791 continue
      write(6,*)"#Inputs after sorting:"
      write(6,*)"#  Smax_in   Smin_in   SRange_in   Cycles_in"
      do 795 i=1,nargsets
        write(6,793)Smaxin(i),Sminin(i),Srangein(i),
     &     Cyclesin(i)
  793   format("#Inputs= ",3(1x,f6.1),1x,f11.1)
  795 continue
       

C     Initilize some counters and flags for reading stuff out of material file
      ifitted=0
C            =1 for each "#DataType= fitted" tag
      ninput=0
C           = Counts input lines
      ndata=0
C          =  Counts input data lines
      numnames=0
C             = counts no. #NAMES
      ifordnumber= 0
C             = contains Ford file number ID
      numcomment=0
C             = counts no. # comment lines
      

C---------------------------- Read in Fitted material file from stdin-------------
  800 continue
c     Loop back to here for next input line.

C     Input lines may either be 1. data lines,  2. #Comment lines or 3. blank
C        It will be hard to distinguish the real comment from the junk comment.
C        Leave it up to the user to edit the comment section.

      read(10,"(a300)",end=980)inp300
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
  815       format(" ERROR saefcalc2: input line no. ",I5,
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
     &            " recompile saefcalc2.f")
           stop
         endif
         read(inp300,*)StrainAmp(ndata), Lifecycles(ndata), 
     &                 StressAmp(ndata)
Cdebug         write(6,*)StrainAmp(ndata), Longlife(ndata), 
Cdebug      &            StressAmp(ndata)

C          If we have not crashed by here, the data from the line has been read in. 
C          Now figure out if there is comment at the end of line.
C          In saefcalc2 we don't really need to do this, but the code is here for
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
         if(cfiletype.ne."strain_life")then
C          We have the wrong kind of sae file here folks
           write(0,*)" WARNING: wrong type of SAE std file. ",
     &               " Not #FileType= strain_life"
         endif
C        filetype is ok.  Write it out as a comment
         go to 950
      endif

      if(firstfield .eq."#DataType=" .or.
     &   firstfield .eq."#DATATYPE=" .or.
     &   firstfield .eq."#datatype=" )then
         read(inp300,*) firstfield, Cdatatype
         if(Cdatatype.ne."fitted")then
C          We have the wrong kind of sae file here folks
           write(0,*)" ERROR, wrong type of SAE std file. ",
     &               " Not #DataType= fitted"
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

      if(firstfield .eq."#LIFEUNITS=" .or.
     &   firstfield .eq."#Life_units=" .or.
     &   firstfield .eq."#LIFE_UNITS=" .or.
     &   firstfield .eq."#Life_Units=" .or.
     &   firstfield .eq."#life_units=" )then
         read(inp300,*) firstfield, lifeunits
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
          go to 950
      endif

       if(firstfield .eq."#%RA=" .or.
     &    firstfield .eq."#%Ra=" )then 
          read(inp300,*) firstfield, percentRA
          go to 950 
      endif 
  
      if(firstfield .eq."#WebPage=" .or.
     &   firstfield .eq."#Webpage=" .or.
     &   firstfield .eq."#WEBPAGE=" .or.
     &   firstfield .eq."#webpage=" )then
         read(inp300,*) firstfield, Cwebpage
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


      if(firstfield .eq."#Ford=" .or.
     &   firstfield .eq."#FORD=" .or.
     &   firstfield .eq."#ford=" )then
         read(inp300,*) firstfield, ifordnumber
         go to 950
      endif


      if(firstfield .eq."#BHN=" .or.
     &   firstfield .eq."#bhn=" .or.
     &   firstfield .eq."#Bhn=" .or.
     &   firstfield .eq."#HBN=" .or.
     &   firstfield .eq."#HB=" )then
         read(inp300,*) firstfield, xnbrinell
         nbrinell=IFIX(xnbrinell)
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

  980 continue
C     All input lines have been read.  The comments were
C     put out along the way.  Its time to dump out the data
C     in whatever format the local machine wants it.  This bit
C     is probably site specific.

C     See if some of the critical values are missing:
      if(EMOD.eq.0.)then
        write(0,*)" ERROR #EMOD= missing value"
        write(6,*)" ERROR #EMOD= missing value"
        stop
      endif
      if(Cdatatype.eq." ")then
        write(0,*)" ERROR #DataType= missing value"
        write(6,*)" ERROR #DataType= missing value"
        stop
      endif
      if(Sult.eq.0.)then
        write(0,*)" WARNING: #Sult= missing value"
        write(6,*)" WARNING: #Sult= missing value"
        write(6,*)" Cannot compute Goodman Damage."
        stop
      endif
      write(6,*)"#CHECK THESE Assumptions:"
      write(6,*)"#Stress_units=",stressunits," Strain_units=",
     &          strainunits," Life_units=",lifeunits
      write(6,*)"#EMOD=",EMOD, " Sult=",Sult," Syield=",Syield,
     &          " %RA=",percentRA," Fracture_Stress=",FractureStress,
     &          " Fracture_Strain=",FractureStrain
      write(6,*)"#Ford=",ifordnumber," BHN=",nbrinell

      
      if(stressunits .eq."mpa" .or.
     &   stressunits .eq."MPA" .or.
     &   stressunits .eq."MPa" .or.
     &   stressunits .eq."Mpa" )then
         stressunits="MPa"
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
         FractureStress=FractureStress*XMPAS
         do 983 i=1,ndata
           StressAmp(i)=StressAmp(i)*XMPAS
  983    continue
         write(6,*)"#Material file ksi -> MPa."
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
           StressAmp(i)=StressAmp(i)*xtemp
  985    continue
         write(6,*)"#Material file: psi -> MPa."
         stressunits="MPa"
      endif
 1010 continue

      if(strainunits .eq."Microstrain" .or.
     &   strainunits .eq."MICROSTRAIN" .or.
     &   strainunits .eq."MicroStrain" )then
         FractureStrain=FractureStrain/1.0E6
         do 1020 i=1,ndata
           StrainAmp(i)=StrainAmp(i)/1.0E+6
 1020   continue
         write(6,*)"#Material file Microstrain -> strain"
         strainunits="strain"
      endif

      if(lifeunits .eq."reversals" .or.
     &   lifeunits .eq."REVERSALS" .or.
     &   lifeunits .eq."Reversals" )then
         do 1030 i=1,ndata
           Lifecycles(i)=Lifecycles(i)/2.
 1030   continue
         write(6,*)"#Material file: Reversals -> Cycles"
         lifeunits="cycles"
      endif



C     Now sort the stress,strain, life values  ------------------
      do 1090 i=1,ndata
        ifind=0
        biglife=0.

        do 1080 j=i,ndata
         if(biglife.gt.Lifecycles(j)) go to 1080
C        Nope, its bigger
         biglife=Lifecycles(j)
         ifind=j
 1080  continue

Cdebug       write(6,*)" Sorting ",ifind, " into ",i
C      biglife is now the biggest & is at location j
C      Save it in a temp location
       tempstrain=StrainAmp(ifind)
       tempstress=StressAmp(ifind)
       templife=Lifecycles(ifind)

C      Now move the stuff pointed to by the big loop
C      into "biglife's" location
       StrainAmp(ifind)  =StrainAmp(i)
       StressAmp(ifind)  =StressAmp(i)
       Lifecycles(ifind)=Lifecycles(i)

C      Put the temp stuff into the "i" location
       StrainAmp(i)=tempstrain
       StressAmp(i)=tempstress
       Lifecycles(i)=templife
 1090 continue
      do 1092 i=1,ndata
C        Look for FractureStress & Strain  value at N=0.5
         if(Lifecycles(i) .eq. 0.5)then
C          We have a good approx of  FractureStress & Strain
           FractureStress=StressAmp(i)
           write(6,*)"#Took SigmaF_primed value =",
     &                FractureStress," at N=0.5"
           FractureStrain=StrainAmp(i)
           write(6,*)"#Took Frac.Strain value =",
     &                FractureStrain," at N=0.5"
         endif
Cdebug      write(6,*)StrainAmp(i),StressAmp(i),Lifecycles(i)
 1092 continue
 
      if(FractureStress.eq. 0.)then
C       If no N=0.5 occured we need to extrapolate for FractureStress
        xslope=
     &   ( alog10(StressAmp(ndata-1 )) -alog10(StressAmp(ndata )) ) /
     &   ( alog10(Lifecycles(ndata-1)) -alog10(Lifecycles(ndata)) )
        xlogSigfp= alog10( StressAmp(ndata-1)) -
     &    xslope*( alog10(Lifecycles(ndata-1)) -alog10(0.5) )
        FractureStress=10.0**xlogSigfp
        write(6,*)"# Extrapolated to FractureStress =",
     &             FractureStress
        xslope=
     &   ( alog10(StrainAmp(ndata-1)) -alog10(StrainAmp(ndata)) ) /
     &   ( alog10(Lifecycles(ndata-1)) -alog10(Lifecycles(ndata)) )
        xlogS= alog10(StrainAmp(ndata-1)) -
     &    xslope*( alog10(Lifecycles(ndata-1)) -alog10(0.5) )
        FractureStrain=10.0**xlogS
        write(6,*)"# Extrapolated to FractureStrain =",
     &             FractureStrain
      endif


C     All is well with input data. Fill in the other "columns"
C     for the material file matrix (some folks might 
C     call this a spreadsheet).
Cbugfix  Oct24/04: change to #xcalc1 tag :
      write(6,1096)
 1096 format("#xcalc1 Strain_Amp     Cycles  Stress_Amp",
     &  "  Elas_Str_Amp ",
     &  " Plas_Str_Amp  Smax*Str_Amp  Snominal_Amp")
      do 1100 i=1,ndata
C       Put out in reverse order to storage
        j=ndata+1-i
        ElstrainAmp(j)=StressAmp(j)/EMOD
        PlstrainAmp(j)=StrainAmp(j)-ElstrainAmp(j)
        SigmaxStrainAmp(j)=StressAmp(j)*StrainAmp(j)
        SelasticAmp(j)=SQRT( StressAmp(j)*StrainAmp(j)*EMOD)
      write(6,1098) StrainAmp(j),Lifecycles(j),
     &   StressAmp(j),ElstrainAmp(j),PlstrainAmp(j),
     &   SigmaxStrainAmp(j),SelasticAmp(j)
Cbugfix  Oct24/04: change to #xcalc1 tag :
 1098 format("#xcalc1 ",f8.5,1x,f11.1,1x,f6.1,
     &       1x,f8.5,1x,f8.5,1x,E14.7,1x,f7.1)
 1100 continue
      



C analysis ------------------------------------------------ 

C We are not going to bother with other assumptions, but will deal
C only with hanging the loops on the largest Tensile 1/2 cycle.
C If for some reason the loop lies outside of the bounds of the 
C biggest cycle (it should not if rainlow counted), then hang it
C on the skeleton curve again.  In some cases this can cause a substantial
C diff. in the damage due to altered mean stress, so maybe send a warning 
C to the user.

C     Init total damage counters
      totswatdam=0.
      totstraindam=0.
      totstressdam=0.
      totgooddam=0.
      totmordam=0.

C     Loop for each user Smax,Smin,N triplet
      do 1700 iloop=1,nargsets

C       This is same for all loops:
        Samp=(Smaxin(iloop)-Sminin(iloop))/2.
        call getLoad2StressStrain(Samp,Sigamp,Epsamp,iexit)
        write(6,*)"getL2S-S:",Samp,Sigamp,Epsamp,iloop
        DeltaSigin(iloop)=Sigamp*2
        DeltaEpsin(iloop)=Epsamp*2
C       Check for overload beyond all data
        if(iexit.eq.1)then
          write(6,1240)Smaxin(iloop),Sminin(iloop)
 1240     format(" ERROR: Overload beyond all property",
     &    " values, Smaxin,Sminin=",2f6.1)
          go to 8000
        endif

C      The 1st has the biggest range due to previous sort, 
C      so it must hang on the skeleton curve.
       if(iloop .eq. 1) go to 1300

C      No, then it may fit inside the biggest which is no.1
C      Check for fit:
       if(Smaxin(iloop).le.Smaxin(1) .and.
     &    Sminin(iloop).ge.Sminin(1) )then
C         Yes it fits
C         The biggest loop was previously "placed", we already know
C         the two triples: Smaxin(1),Sigmaxin(1),Epsmaxin(1)  
C            and           Sminin(1),Sigminin(1),Epsminin(1)
C         Now find this smaller loops stuff

C         Now place the upper tip. Its attached to the biggest 1/2 cycle (1)
	  DiffS=Smaxin(iloop)-Sminin(1)
C         Get the DS & DE to get us there
          Samp=DiffS/2.
          call getLoad2StressStrain(Samp,Sigamp,Epsamp,iexit)
          write(6,*)"getL2S-S:",Samp,Sigamp,Epsamp,iloop
          if(iexit.eq.1)then
            write(6,1240)Smaxin(iloop),Sminin(iloop)
            go to 8000
          endif
          Sigmaxin(iloop)=Sigminin(1)+(Sigamp*2.)
          Epsmaxin(iloop)=Epsminin(1)+(Epsamp*2.)
C         Now place the lower tip: Upper tip less deltas
          Sigminin(iloop)=Sigmaxin(iloop)-DeltaSigin(iloop)
          Epsminin(iloop)=Epsmaxin(iloop)-DeltaEpsin(iloop)
          write(6,*)"NonSkeleton tips:",Sigmaxin(iloop),
     &      Epsmaxin(iloop),Sigminin(iloop),Epsminin(iloop)
     &     ,"  Nom:",Smaxin(iloop),Sminin(iloop)," loop=",iloop
C         Loop is placed. Go get the damage
          go to 1500
        endif

C       No, loop does not fit inside. Stick it on the skeleton curve.
 1300   continue
        Saverage=(Smaxin(iloop)+Sminin(iloop))/2.0
        if(Saverage .ge. 0.)then
C         The loop is more on the tensile side. Attach to tensile
C         skeleton curve.
          Samp= Smaxin(iloop)
          call getLoad2StressStrain(Samp,Sigamp,
     &             Epsamp,iexit)
          write(6,*)"getL2S-S:",Samp,Sigamp,Epsamp,iloop
          if(iexit.eq.1)then
            write(6,1240)Smaxin(iloop),Sminin(iloop)
            go to 8000
          endif
          Sigmaxin(iloop)=Sigamp
          Epsmaxin(iloop)=Epsamp
          Sigminin(iloop)=Sigmaxin(iloop)-DeltaSigin(iloop)
          Epsminin(iloop)=Epsmaxin(iloop)-DeltaEpsin(iloop)
          write(6,*)"TensSkeleton:",Sigmaxin(iloop),
     &      Epsmaxin(iloop),Sigminin(iloop),Epsminin(iloop)
     &     ,"Nom:",Smaxin(iloop),Sminin(iloop)
C         Loop is placed, go get damaged
          go to 1500
        else
C         Loop is more into compression. Thus minimum tip is attaced
          Samp=ABS(Sminin(iloop))
          call getLoad2StressStrain(Samp,Sigamp,
     &             Epsamp,iexit)
          if(iexit.eq.1)then
            write(6,1240)Smaxin(iloop),Sminin(iloop)
            go to 8000
          endif
          Sigminin(iloop)=-Sigamp
          Epsminin(iloop)=-Epsamp
C         Use deltas to get to upper tip
          Sigmaxin(iloop)=Sigminin(iloop)+DeltaSigin(iloop)
          Epsmaxin(iloop)=Epsminin(iloop)+DeltaEpsin(iloop)
          write(6,*)"CompSkeleton:",Sigmaxin(iloop),
     &      Epsmaxin(iloop),Sigminin(iloop),Epsminin(iloop)
     &     ," Nom:",Smaxin(iloop),Sminin(iloop)
C         Loop is placed, go get damaged
          go to 1500
        endif

 1500   continue
C       Damage assessment
C       We know where each loop's tip is. Sigmaxin Epsmaxin DeltaSigin DeltaEpsin
C                                         Sigminin Epsminin Deltasigin DeltaEpsin
C       We also know that the biggest loop is in (1), and could thus deploy
C       some overload compensation code- although we may already have that as
C       the input curve.

C       As usual, we have the problem of what is infinity?  i.e. if the level
C       is below the last user value of the material file, what value should we
C       pass back at inf?   Lets try the -ve of the last material curve point.
C       Thus a negative returned life means it was below the last data point on
C       the material curve.  And the absolute value is the value of the last
C       point on the material curve.

C       It is also possible to have flat spots in the life curves.  If we
C       hit a flat spot, the lowest life value will be used.  What about 
C      

C       Compute the Smith/Watson/Topper product Smax * Strain_Amp
C       Cyclesin() was checked for +ve
        xamp=DeltaEpsin(iloop)/2.
        xN=Cyclesin(iloop)

        Swatproduct=Sigmaxin(iloop) * xamp
        call getSwat2Life(Swatproduct,Swatlife,iexit)
        write(6,*)"getSwatl:",Swatproduct,Swatlife
C       That was per loop, now how many loops in this input triplet?
C       Infinite life will be a -ve number
        if(Swatlife .gt.0)then
           SwatReps(iloop)=Swatlife/xN
           totswatdam=totswatdam+(xN/Swatlife)
        else 
           SwatReps(iloop)=Swatlife
        endif

        call getStrain2Life(xamp,Strainlife,iexit)
        write(6,*)"getStrl:",xamp,Strainlife
        if(Strainlife .gt.0.)then
           StrainReps(iloop)=Strainlife/xN
           totstraindam=totstraindam+(xN/Strainlife)
        else 
Cbugfix    May31/04:
           StrainReps(iloop)=Strainlife
        endif

        Smean= (Sigmaxin(iloop)+Sigminin(iloop) )/2.
        Sa=DeltaSigin(iloop)/2.0
        call getStress2Life(Sa,Stresslife,iexit)
        write(6,*)"getStsl:",Sa,Stresslife," Sts"
        if(Stresslife .gt. 0.)then
           StressReps(iloop)=Stresslife/xN
           totstressdam=totstressdam+(xN/Stresslife)
        else 
Cbugfix    May31/04:
           StressReps(iloop)=Stresslife
        endif

C       Compute equiv stress for Goodman
        SeqGoodman=Sa * (Sult/(Sult-Smean))
        call getStress2Life(SeqGoodman,Goodlife,iexit)
        write(6,*)"getStsl:",SeqGoodman,Goodlife," Goo"
        if(Goodlife .gt. 0.)then
           GoodReps(iloop)=Goodlife/xN
           totgooddam=totgooddam+(xN/Goodlife)
        else 
Cbugfix    May31/04:
           GoodReps(iloop)=Goodlife
        endif

C       Compute Morrow equiv.
        SeqMorrow=Sa *(FractureStress/(FractureStress-Smean))
        call getStress2Life(SeqMorrow,Morlife,iexit)
        write(6,*)"getStsl:",SeqMorrow,Morlife," Mor"
        if(Morlife .gt. 0.)then
           MorReps(iloop)=Morlife/xN
           totmordam=totmordam+(xN/Morlife)
        else 
Cbugfix    May31/04:
           MorReps(iloop)=Morlife
        endif



C     End of loop selction do-loop
 1700 continue

C     --------------------------------------------------------------
C     Calc the % damage contribution of each input level
      write(6,*)"Debug totdams:",totswatdam,totstraindam,
     &  totstressdam,totgooddam,totmordam

      do 1800 iloop=1,nargsets


       if(SwatReps(iloop).gt.0)then
         pcswatdam(iloop)= (1.0/SwatReps(iloop))/totswatdam
     &   *100.
       else
         pcswatdam(iloop)=0.
       endif

       if(StrainReps(iloop).gt.0)then
         pcstraindam(iloop)= (1.0/StrainReps(iloop))/totstraindam
     &   *100.
       else
         pcstraindam(iloop)=0.
       endif

       if(StressReps(iloop).gt.0)then
         pcstressdam(iloop)= (1.0/StressReps(iloop))/totstressdam
     &   *100.
       else
         pcstressdam(iloop)=0.
       endif

       if(GoodReps(iloop).gt.0)then
         pcgooddam(iloop)= (1.0/GoodReps(iloop))/totgooddam
     &   *100.
       else
         pcgooddam(iloop)=0.
       endif

       if(MorReps(iloop).gt.0)then
         pcmordam(iloop)= (1.0/MorReps(iloop))/totmordam
     &   *100.
       else
         pcmordam(iloop)=0.
       endif

      write(6,*)"Debug %dams:",
     & pcstraindam(iloop),pcswatdam(iloop),
     & pcstressdam(iloop),pcgooddam(iloop),pcmordam(iloop)
 1800 continue

C     All done with calcs.


C    ---------------------------------------------- Output----------
C     Write out the title line
C     Get the first 10 chars from each of first two names
      name1=" "
      name2=" "
      if(numnames.ge.2)then
        read(names30(1),*)name1
        read(names30(2),*)name2
      endif
      if(numnames.eq.1)then
        read(names30(1),*)name1
      endif
      write(6,1490),name1,name2, nbrinell,ifordnumber
 1490 format("# ",A10," ",A10," BHN= ",I5," Fn= ",I4)

      write(6,2010)
 2010 format("#xcalc2 Loop   Smax    Smin         N  Sigmax ",
     & " Sigmin Delta Epsmax Epsmin DeltaEps ",
     & " %Eps %SWaT %Sts %Morr %Goodm ")
      do 2020 iloop=1,nargsets
      write(6,2014)iloop, Smaxin(iloop),Sminin(iloop),Cyclesin(iloop),
     & Sigmaxin(iloop),Sigminin(iloop),DeltaSigin(iloop),
     & Epsmaxin(iloop),Epsminin(iloop),DeltaEpsin(iloop),
     & pcstraindam(iloop),pcswatdam(iloop),pcstressdam(iloop),
     & pcmordam(iloop),pcgooddam(iloop)
Cbugfix  Oct24/04 Widen local stress format fields to allow -1999. :
 2014 format("#xcalc2 ",i4, 1x,f7.1,1x,f7.1,f11.1,
     & 1x,f6.0,1x,f6.0,1x,f5.0, 3(1x,f7.5), 5(1x,f5.1) )

 2020 continue
 
      if(totstraindam.eq.0.)then
        clifeout(1)=" Infinity  "
      else
       write(clifeout(1),2025)1.0/totstraindam
 2025  format(f11.1)
      endif

      if(totswatdam.eq.0.)then
        clifeout(2)=" Infinity  "
      else
       write(clifeout(2),2025)1.0/totswatdam
      endif

      if(totstraindam.eq.0.)then
        clifeout(3)=" Infinity  "
      else
       write(clifeout(3),2025)1.0/totstressdam
      endif

      if(totstraindam.eq.0.)then
        clifeout(4)=" Infinity  "
      else
       write(clifeout(4),2025)1.0/totmordam
      endif

      if(totstraindam.eq.0.)then
        clifeout(5)=" Infinity  "
      else
       write(clifeout(5),2025)1.0/totgooddam
      endif

      write(6,2028)
 2028 format("#xcalc3  StrainLife_Reps SWaT_Life_Reps ",
     &   "StressLife_Reps   Morrow_Reps   Goodman_Reps ",
     &   "(Reps= Repetions)" )
      write(6,2030)(clifeout(i),i=1,5)
 2030 format("#xcalc3 ",10(a11,5x))


C-------------- Plotting ---------------------------------------------

      do 2500 iloop=1,nargsets

C     If it is the 1st loop, we should also do the skeleton curve
      if(iloop.eq.1)then
C       Find the points for the skeleton curve
        Saverage=(Sigmaxin(iloop)+Sigminin(iloop))/2.
        if(Saverage.ge.0.)then
C         Its on the tensile side
          Samp=Sigmaxin(iloop)
          call getloop(Samp,tstrain,tstress,nvalues,iexit)
C         By this time iexit should always be 0
C         The vectors tstrain & tstress have the points along the 1/2 cycle
          write(6,2110)
 2110     format("#plotloops"/
     &       "#plotloops #Tensile Skeleton curve")
c         The first point is always(?)  0,0
          do 2150 i=1,nvalues
           write(6,2115)tstrain(i),tstress(i)
Cbugfix    Oct 24/04: increase size of stress format field:
 2115      format("#plotloops ",f8.5,1x,f8.1)
 2150      continue
          go to 2200
        else
C         Its on the Comp. side
          Samp=-Sigminin(iloop)
          call getloop(Samp,tstrain,tstress,nvalues,iexit) 
C         The vectors tstrain & tstress have the points along the 1/2 cycle 
          write(6,2120)
 2120     format("#plotloops"/
     &          "#plotloops #Compression Skeleton curve") 
c         The first point is always(?)  0,0 
          do 2160 i=1,nvalues 
           write(6,2115)-tstrain(i),-tstress(i) 
 2160      continue 
          go to 2200 
        endif
      endif

C     Its a loop
 2200 continue
C     The loops have all been previously located in space, thus
c     all we need is the "fill" between the tips
      Samp=(Sigmaxin(iloop)-Sigminin(iloop) )/2.
      call getloop(Samp,tstrain,tstress,nvalues,iexit)
      write(6,2205)
 2205 format("#plotloops")
      

C     Start at the compressive tip and go up:
      write(6,2207)iloop
 2207 format("#plotloops #Tensile side of loop ",i4)
C     tstrain & tstress are in amplitudes, thus:
      estart=Epsminin(iloop)
      sstart=Sigminin(iloop)
      do 2050 i=1,nvalues
        x=estart+tstrain(i)*2.
        s=sstart+tstress(i)*2.
        write(6,2115)x,s
 2050 continue

C     Start at the tensile tip and go down: 
      write(6,2257)iloop 
 2257 format("#plotloops #Compressive side of loop ",i4)
C     tstrain & tstress are in amplitudes, thus:
      estart=Epsmaxin(iloop)
      sstart=Sigmaxin(iloop)
      do 2060 i=1,nvalues  
        x=estart-tstrain(i)*2.   
        s=sstart-tstress(i)*2.
        write(6,2115)x,s 
 2060 continue 
 
C     end of loop plotting output
 2500 continue



C---------- History/Damage Plotting-------------------------------

C  In this section, we will put out a "box" for each cycle set.
C  The box is defined by Smaxin() on top, Sminin() on bottom, and
C  its X axis co-ords are defined by Cyclesin() = addcycles2 - addcycles1
C  Where addcycles1 is the cumulative cycle count as one loops through the
C  sorted (largest ranges first) cycle sets.
C  An additional output is the percent damage for each set. From experience
C  we should use the max of pcswatdam() and pcmordam(), but both will be plotted.
C  The cumulative cycle boxes will be plotted on semi-log co-ords,
C  with the X axis minimum of 10**-1 or 0.1 cycles.

C     Initilize the cumulative cycle counters
      addcycles2=0.

      do 2600 iloop=1,nargsets
        addcycles1=addcycles2
        addcycles2=addcycles1+Cyclesin(iloop)
C       The 1st time is special:
        if(iloop .eq.1)addcycles1=0.1

        write(6,2548)
 2548   format("#plothist")
        write(6,2550)addcycles1,Smaxin(iloop)
 2550   format("#plothist ",f11.1,1x,f11.2)
        write(6,2550)addcycles2,Smaxin(iloop)
        write(6,2550)addcycles2,Sminin(iloop)
        write(6,2550)addcycles1,Sminin(iloop)
        write(6,2550)addcycles1,Smaxin(iloop)

C       The x axis for damage will be the center of the cycle box
        x=10**( (ALOG10(addcycles1) + ALOG10(addcycles2) )/2.)
        write(6,2560)x,pcswatdam(iloop),pcmordam(iloop)
 2560   format("#plotdam ",f11.1, 1x,f6.1, 1x,f6.1)

 2600 continue


      go to 9000


C
 8000 continue
      write(6,*)"Failure in <=1 cycle."
      write(0,*)"Failure in <=1 cycle."

 9000 CONTINUE
C      write(6,9005)
C 9005 format(
C     &"#Here is the bottom part of the html  "
C     & ,"graph/calc wrapper-----------"/
C     &  "#</textarea>"/
C     &  "#</pre></DL></FORM></body></html>"
C     & )
      STOP
      END


C---------------------------------------------------------
      SUBROUTINE getStress2Life(Stress,Cycles,iexit)
      SAVE
C     Given Stress Ampl., interpolate Cycles

      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
     &     SigmaxStrainAmp(250)
      integer ndata
      real FractureStress,FractureStrain

      Common/Material/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
     &                 ndata,FractureStress,FractureStrain
      iexit=0

C       We are assuming that the data has been sorted to largest life,
C       and thus smallest stress, strain etc, first.
C       If the stress is below the first stress point, then it is below
C       the fatigue limit (the first life value in the table).

C      check if data is below the fatigue limit
      if(Stress .lt. StressAmp(1))then
C       Yes, its below. Set it equal to - fat_limit
        Cycles= -Lifecycles(1)
        return
      endif
      
C     Check if data is above FractureStress
      if(Stress .ge.FractureStress)then
C       The specimen fails in first 1/2 cycle. Damage is >= 1.
        Cycles=0.5
        iexit=1
        return
      endif

C     Ok, the input appears to lie between the two extremes. Lets
C     hunt & peck.  Basically we are going to hunt our way up the stress
C     ladder.  When our point is bigger than the last rung, and smaller than
C     the next rung, we can interpolate.

      do 1000 idat=1,ndata
        if(Stress - StressAmp(idat) ) 100,200,300
C         In-elegant but fast:          -ve  0  +ve
C         When the result is -ve, then the curve point is just above
C         the point we want, and we can interpolate.


  100   continue
C       Our stress is less than the point on the curve, we have arrived
C       Interpolate.  Since this is Stress vs Life, use log-log interpl.
        xslope=
     &     (ALOG10(StressAmp(idat)) -ALOG10(StressAmp(idat-1))  )
     &    /(ALOG10(Lifecycles(idat))-ALOG10(Lifecycles(idat-1)) )
        Clog=  (ALOG10(Stress)-ALOG10(StressAmp(idat-1)) ) /xslope
     &           +ALOG10(Lifecycles(idat-1))
        Cycles= 10.0**Clog
        return

  200   continue
C       The stress is exactly equal to the curve point(idat)
C       There could be a "flat" spot in the stress-life curve, and
C       thus the next point would have the same stress value. Check
Cbugfix   May2005 : im didnt increment
        im=idat
  201   im=im+1
        if(Stress .eq. StressAmp(im))then
C         yes, another is equal. is there more?
          go to 201
        endif
C       No, now the (im) point has bigger stress, use the last
C       one that was equal as a match
Cbugfix   May2005 : Should be im-1
        Cycles=Lifecycles(im-1)
        return

  300   continue
C       Ok, Point(idat) on the curve is still smaller, keep going
        go to 1000

 1000   continue
        write(0,*)" ERROR:saefcalc2:getStress2Life.endofdoloop"
        iexit=1
C       We should not be here.  Assume fracture strain
        Cycles=0.5
        return
        end





C---------------------------------------------------------
      SUBROUTINE getStrain2Life(Strain,Cycles,iexit)
      SAVE
C     Given Strain Ampl., interpolate Cycles

      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
     &     SigmaxStrainAmp(250)
      integer ndata
      real FractureStress,FractureStrain

      Common/Material/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
     &                 ndata,FractureStress,FractureStrain
      iexit=0

C       We are assuming that the data has been sorted to largest life,
C       and thus smallest stress, strain etc, first.
C       If the strain is below the first strain point, then it is below
C       the fatigue limit (the first life value in the table).

C      check if data is below the fatigue limit
      if(Strain .lt. StrainAmp(1))then
C       Yes, its below. Set it equal to - fat_limit
        Cycles= -Lifecycles(1)
        return
      endif
      
C     Check if data is above FractureStrain
      if(Strain .ge.FractureStrain)then
C       The specimen fails in first 1/2 cycle. Damage is >= 1.
        Cycles=0.5
        iexit=1
        return
      endif

C     Ok, the input appears to lie between the two extremes. Lets
C     hunt & peck.  Basically we are going to hunt our way up the stress
C     ladder.  When our point is bigger than the last rung, and smaller than
C     the next rung, we can interpolate.

      do 1000 idat=1,ndata
        if(Strain - StrainAmp(idat) ) 100,200,300
C         In-elegant but fast:          -ve  0  +ve
C         When the result is -ve, then the curve point is just above
C         the point we want, and we can interpolate.


  100   continue
C       Our stress is less than the point on the curve, we have arrived
C       Interpolate.  Since this is Strain vs Life, use log-log interpl.
        xslope=
     &     (ALOG10(StrainAmp(idat)) -ALOG10(StrainAmp(idat-1))  )
     &    /(ALOG10(Lifecycles(idat))-ALOG10(Lifecycles(idat-1)) )
        Clog=  (ALOG10(Strain)-ALOG10(StrainAmp(idat-1)) )/xslope
     &          +ALOG10(Lifecycles(idat-1))
        Cycles= 10.0**Clog
        return

  200   continue
C       The Strain is exactly equal to the curve point(idat)
C       There could be a "flat" spot in the Strain-life curve, and
C       thus the next point would have the same Strain value. Check
Cbugfix   May2005 : im didnt increment
        im=idat
  201   im=im+1
        if(Strain .eq. StrainAmp(im))then
C         yes, another is equal. is there more?
          go to 201
        endif
C       No, now the (im) point has bigger Strain, use the last
C       one that was equal as a match
Cbugfix   May2005 : Should be im-1
        Cycles=Lifecycles(im-1)
        return

  300   continue
C       Ok, Point(idat) on the curve is still smaller, keep going

 1000   continue
        write(0,*)" ERROR:saefcalc2:getStrain2Life.endofdoloop"
        iexit=1
C       We should not be here.  Assume fracture strain
        Cycles=0.5
        return
        end




C---------------------------------------------------------
      SUBROUTINE getStress2Strain(Stress,Strain,iexit)
      SAVE

C     Given Stress Ampl., interpolate Strain

      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
     &     SigmaxStrainAmp(250)
      integer ndata
      real FractureStress,FractureStrain

      Common/Material/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
     &                 ndata,FractureStress,FractureStrain
      iexit=0

C       We are assuming that the data has been sorted to largest life,
C       and thus smallest stress, strain etc, first.
C       If the stress is below the first stress point, then it is below
C       the fatigue limit (the first life value in the table).

      if(Stress .lt. StressAmp(1))then
C       Yes, its below. Assume straight line from (0,0)
        Strain=(Stress/StressAmp(1))*StrainAmp(1)
        return
      endif
      
C     Check if data is above Sigfp
      if(Stress .ge. Fracturestress)then
C       The specimen fails in first 1/2 cycle. Damage is >= 1.
        Strain=Fracturestrain
        iexit=1
        return
      endif

C     Ok, the input appears to lie between the two extremes. Lets
C     hunt & peck.  Basically we are going to hunt our way up the stress
C     ladder.  When our point is bigger than the last rung, and smaller than
C     the next rung, we can interpolate.

      do 1000 idat=1,ndata
        if(Stress - StressAmp(idat) ) 100,200,300
C         In-elegant but fast:          -ve  0  +ve
C         When the result is -ve, then the curve point is just above
C         the point we want, and we can interpolate.


  100   continue
C       Our stress is less than the point on the curve, we have arrived
C       Interpolate.  Try log-log interpl.  If we ever get -ve stresses
C       or strains we will definitely have to go to linear interpl.
        xslope=
     &     (ALOG10(StressAmp(idat)) -ALOG10(StressAmp(idat-1)))
     &    /(ALOG10(StrainAmp(idat))-ALOG10(StrainAmp(idat-1)) )
        Clog=  (ALOG10(Stress)-ALOG10(StressAmp(idat-1)) )/xslope
     &          +ALOG10(StrainAmp(idat-1))
        Strain= 10.0**Clog
        return

  200   continue
C       The stress is exactly equal to the curve point(idat)
C       There could be a "flat" spot in the stress-strain curve, and
C       thus the next point would have the same stress value. Check
Cbugfix   May2005 : im didnt increment
        im=idat
  201   im=im+1
        if(Stress .eq. StressAmp(im))then
C         yes, another is equal. is there more?
          go to 201
        endif
C       No, now the (im) point has bigger stress, use the last
C       one that was equal as a match
Cbugfix   May2005 : Should be im-1
        Strain=StrainAmp(im-1)
        return

  300   continue
C       Ok, Point(idat) on the curve is still smaller, keep going
        go to 1000

 1000   continue

        write(0,*)" ERROR:saefcalc2:getStress2Strain.endofdoloop"
        iexit=1
C       We should not be here.  Assume fracture strain
        Strain=Fracturestrain
        return
        end


C---------------------------------------------------------
      SUBROUTINE getLoad2StressStrain(Sneuber,Stress,Strain,
     &           iexit)
      SAVE
C     Given Load_Amp (could be FEA Stress Ampl.), interpolate Strain &Stress

      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
     &     SigmaxStrainAmp(250)
      integer ndata
      real FractureStress,FractureStrain

      Common/Material/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
     &                 ndata,FractureStress,FractureStrain
      iexit=0

C       We are assuming that the data has been sorted to largest life,
C       and thus smallest stress, strain etc, first.
C       If the "load", actually = Elastic_local_Stress, is below 
C       the first SelasticAmp point, then it is below
C       the fatigue limit (the first life value in the table).

      if(Sneuber .lt. SelasticAmp(1))then
C       Yes, its below. Assume straight line from (0,0)
        Strain=(Sneuber/SelasticAmp(1))*StrainAmp(1)
        Stress=(Sneuber/SelasticAmp(1))*StressAmp(1)
        return
      endif
      
C     Check if data is above Sigfp
C     Assumes: ElasMod= First Stress/ 1st strain
C     If these first points are not on the "elastic" line we could have
C     a problem here, but it shouldnt be too bad (?)
      E=StressAmp(1)/StrainAmp(1)
      if(Sneuber .ge. 
     &        SQRT(Fracturestress*Fracturestrain*E))then
C       The specimen fails in first 1/2 cycle. Damage is >= 1.
        Strain=Fracturestrain
        Stress=Fracturestress
        iexit=1
        return
      endif

C     Ok, the input appears to lie between the two extremes. Lets
C     hunt & peck.  Basically we are going to hunt our way up the stress
C     ladder.  When our point is bigger than the last rung, and smaller than
C     the next rung, we can interpolate.

      do 1000 idat=1,ndata
        if(Sneuber - SelasticAmp(idat) ) 100,200,300
C         In-elegant but fast:          -ve  0  +ve
C         When the result is -ve, then the curve point is just above
C         the point we want, and we can interpolate.

  100   continue 
C       Our stress is less than the point on the curve, we have arrived
C       Interpolate.  Try log-log interpl.  If we ever get -ve stresses
C       or strains we will definitely have to go to linear interpl.
        xslope=
     &     (ALOG10(SelasticAmp(idat)) -ALOG10(SelasticAmp(idat-1)))
     &    /(ALOG10(StressAmp(idat))-ALOG10(StressAmp(idat-1)) )
        Clog=  (ALOG10(Sneuber)-ALOG10(SelasticAmp(idat-1)) )/xslope
     &          +ALOG10(StressAmp(idat-1))
        Stress= 10.0**Clog
        xslope=
     &     (ALOG10(SelasticAmp(idat)) -ALOG10(SelasticAmp(idat-1)))
     &    /(ALOG10(StrainAmp(idat))-ALOG10(StrainAmp(idat-1)) )
        Clog=  (ALOG10(Sneuber)-ALOG10(SelasticAmp(idat-1)) )/xslope
     &          +ALOG10(StrainAmp(idat-1))
        Strain= 10.0**Clog
        return

Cx  100   continue
CxC       The stress is less than the point on the curve, we have arrived.
CxC       Interpolate.  Try log-log interpl.  If we ever get -ve stresses
CxC       or strains we will definitely have to go to linear interpl.
Cx
Cx        fraction=(ALOG10(Sneuber)-ALOG10(SelasticAmp(idat-1))  )/
Cx     &    (ALOG10(SelasticAmp(idat))-ALOG10(SelasticAmp(idat-1)))
Cx        Clog=ALOG10(StressAmp(idat-1)) + fraction *
Cx     &       (ALOG10(StressAmp(idat))-ALOG10(StressAmp(idat-1)))
Cx        Stress= 10.0**Clog

CxC       Since the Sneuber,Strain, Stress occur as triple points, we
CxC       can use the same fraction for Strain interpolation:
Cx        Clog=ALOG10(StrainAmp(idat-1)) + fraction *
Cx     &       (ALOG10(StrainAmp(idat))-ALOG10(StrainAmp(idat-1)))
Cx        Strain= 10.0**Clog
Cx        return

  200   continue
C       Sneuber is exactly equal to the curve point(idat)
C       There could be a "flat" spot in the stress-strain curve, and
C       thus the next point would have the same stress value. 
C       Check :
Cbugfix   May2005 : im didnt increment.
        im=idat
  201   im=im+1
        if(Sneuber .eq. SelasticAmp(im))then
C         yes, another is equal. is there more?
          go to 201
        endif
C       No, now the (im) point has bigger number, use the last
C       one that was equal as a match
Cbugfix   May2005 : Should be im-1
        Strain=StrainAmp(im-1)
        Stress=StressAmp(im-1)
        return

  300   continue
C       Ok, Point(idat) on the curve is still smaller, keep going
        go to 1000

 1000   continue

        write(0,*)" ERROR:saefcalc2:getload2stressStrain:",
     &            " endofDOloop"
        iexit=1
C       We should not be here.  Assume fracture strain
        Strain=Fracturestrain
        Stress=Fracturestress
        return
        end



C---------------------------------------------------------
      SUBROUTINE getSwat2Life(Smea,Cycles,iexit)
      SAVE
C     Given SmithWatsonTopper product, interpolate Cycles

      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
     &     SigmaxStrainAmp(250)
      integer ndata
      real FractureStress,FractureStrain

      Common/Material/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
     &                 ndata,FractureStress,FractureStrain
      iexit=0

C       We are assuming that the data has been sorted to largest life,
C       and thus smallest stress, strain etc, first.
C       If the Smea is below the first SigmaxStrainAmp point, then it is below
C       the fatigue limit (the first life value in the table).

C      check if data is below the fatigue limit
      if(Smea .lt. SigmaxStrainAmp(1))then
C       Yes, its below. Set it equal to - fat_limit
        Cycles= -Lifecycles(1)
        return
      endif
      
Cbug fixed here in May 19 2004:
C     Check if data is above FractureStress
      if(Smea .ge. SigmaxStrainAmp(ndata))then
C       The specimen fails in first 1/2 cycle. Damage is >= 1.
        Cycles=0.5
        iexit=1
        return
      endif

C     Ok, the input appears to lie between the two extremes. Lets
C     hunt & peck.  Basically we are going to hunt our way up the Smea value
C     ladder.  When our point is bigger than the last rung, and smaller than
C     the next rung, we can interpolate.

      do 1000 idat=1,ndata
        if(Smea - SigmaxStrainAmp(idat) ) 100,200,300
C         In-elegant but fast:          -ve  0  +ve
C         When the result is -ve, then the curve point is just above
C         the point we want, and we can interpolate.


  100   continue
C       The Smea value is less than the point on the curve, we have arrived
C       Interpolate.  Since this is Smea vs Life, use log-log interpl.
        xslope=
     &     (ALOG10(SigmaxStrainAmp(idat)) 
     &         -ALOG10(SigmaxStrainAmp(idat-1))  )
     &    /(ALOG10(Lifecycles(idat))
     &         -ALOG10(Lifecycles(idat-1)) )
        Clog=  (ALOG10(Smea)-ALOG10(SigmaxStrainAmp(idat-1)) )
     &            /xslope   +ALOG10(Lifecycles(idat-1))
        Cycles= 10.0**Clog
        return

  200   continue
C       The stress is exactly equal to the curve point(idat)
C       There could be a "flat" spot in the Smea-life curve, and
C       thus the next point would have the same Smea value. Check
Cbugfix   May2005 : im didnt increment.
        im=idat
  201   im=im+1
        if(Smea .eq. SigmaxStrainAmp(im))then
C         yes, another is equal. is there more?
          go to 201
        endif
C       No, now the (im) point has bigger Smea, use the last
C       one that was equal as a match
Cbugfix   May2005 : Should be im-1
        Cycles=Lifecycles(im-1)
        return

  300   continue
C       Ok, Point(idat) on the curve is still smaller, keep going
        go to 1000

 1000   continue
        write(0,*)" ERROR:saefcalc2:getSwat2Life.endofdoloop"
        iexit=1
C       We should not be here.  Assume fracture strain
        Cycles=0.5
        return
        end





C---------------------------------------------------------
      SUBROUTINE getloop(Stress,Straintemp,Stresstemp,npts,iexit)
      SAVE

C     Given Stress Ampl., get all stress-strain points from 0 to Stress
      real Stresstemp(250),Straintemp(250)

      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
     &     SigmaxStrainAmp(250)
      integer ndata
      real FractureStress,FractureStrain

      Common/Material/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
     &                 ndata,FractureStress,FractureStrain
      iexit=0

      Stresstemp(1)=0.
      Straintemp(1)=0
      npts=1

C       We are assuming that the data has been sorted to largest life,
C       and thus smallest stress, strain etc, first.
C       If the stress is below the first stress point, then it is below
C       the fatigue limit (the first life value in the table).

      if(Stress .lt. StressAmp(1))then
C       Yes, its below. Assume straight line from (0,0)
        npts=npts+1
        Straintemp(npts)=(Stress/StressAmp(1))*StrainAmp(1)
        Stresstemp(npts)=Stress
        go to 9000
      endif
      
C     Check if data is above Sigfp
      if(Stress .ge. Fracturestress)then
C       The specimen fails in first 1/2 cycle. Damage is >= 1.
        do 50 i=1,ndata
          npts=npts+1
          Straintemp(npts)=StrainAmp(i)
          Stresstemp(npts)=StressAmp(i)
   50   continue
        npts=npts+1
        Straintemp(npts)=Fracturestrain
        Stresstemp(npts)=Fracturestress
        iexit=1
        go to 9000
      endif

C     Ok, the input appears to lie between the two extremes. Lets
C     hunt & peck.  Basically we are going to hunt our way up the stress
C     ladder.  When our point is bigger than the last rung, and smaller than
C     the next rung, we can interpolate.

      do 1000 idat=1,ndata
        if(Stress - StressAmp(idat) ) 100,200,300
C         In-elegant but fast:        -ve  0  +ve
C         When the result is -ve, then the curve point is just above
C         the point we want, and we can interpolate.


  100   continue
C       The stress is less than the point on the curve, we have arrived
C       Interpolate.  Try log-log interpl.  If we ever get -ve stresses
C       or strains we will definitely have to go to linear interpl.
        xslope=
     &     (ALOG10(StressAmp(idat)) -ALOG10(StressAmp(idat-1)))
     &    /(ALOG10(StrainAmp(idat))-ALOG10(StrainAmp(idat-1)) )
        Clog=  (ALOG10(Stress)-ALOG10(StressAmp(idat-1)) )/xslope
     &          +ALOG10(StrainAmp(idat-1))
        npts=npts+1
        Straintemp(npts)= 10.0**Clog

        deltS=StressAmp(idat)-StressAmp(idat-1) 
        fraction=(Stress-StressAmp(idat-1))/deltS
        Stresstemp(npts)=StressAmp(idat-1) + fraction*deltS
        go to 9000

  200   continue
C       The stress is exactly equal to the curve point(idat)
C       There could be a "flat" spot in the stress-strain curve, and
C       thus the next point would have the same stress value. Check
Cbugfix   May2005 : im didnt increment
        im=idat
  201   im=im+1
        if(Stress .eq. StressAmp(im))then
C         yes, another is equal. is there more?
          go to 201
        endif
C       No, now the (im) point has bigger stress, use the last
C       one that was equal as a match
        npts=npts+1
        Straintemp(npts)=StrainAmp(im-1)
        Stresstemp(npts)=Stress
        go to 9000

  300   continue
C       Ok, Point(idat) on the curve is still smaller, keep going
        npts=npts+1
        Straintemp(npts)=StrainAmp(idat)
        Stresstemp(npts)=StressAmp(idat)

C      end of hunting loop
 1000  continue

        write(0,*)" ERROR:saefcalc2:getloop.endofdoloop"
        iexit=1
C       We should not be here.  Assume fracture strain
        npts=npts+1
        Straintemp(npts)=Fracturestrain
        Stresstemp(npts)=Fracturestress
        go to 9000

 9000   continue
        write(6,*)"#debug Plot loop for Stress Amp: ",Stress
        write(6,*)"   Amplitudes        Ampl.*2"
        do 9005 i=1,npts
        write(6,9003)Straintemp(i),Stresstemp(i),
     &   Straintemp(i)*2.0,Stresstemp(i)*2.0
 9003   format(2(1x,f7.5,1x,f6.1))
 9005   continue
        end

C==============================================================
      SUBROUTINE WR(UN,INP1)
      SAVE
C      S/R TO WRITE A LINE BUT CUT THE TRAILING BLANKS
C      Used when writing out comment lines mostly.
      CHARACTER*1 INP1(300)
      INTEGER UN
      LONG=300
C     FIND LAST NON BLANK CHAR
      DO 10 I=LONG,1,-1
        IF(INP1(I).EQ.' ')GO TO 10
C       NO? THIS IS IT
        N=I
        GO TO 20
   10 CONTINUE
C     COMPLETELY BLANK LINE
      WRITE(UN,11)
   11 FORMAT(' ')
      RETURN
   20 WRITE(UN,24)(INP1(J),J=1,N)
   24 FORMAT(300A1)
      RETURN
      END


