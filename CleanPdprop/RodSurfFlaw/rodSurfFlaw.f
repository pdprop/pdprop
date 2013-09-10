C  rodSurfFlaw.f   vers. 3.06   RoundBar Surf Crack Prop.  FAC sep 5 2013
      SAVE
C  Push-Down List crack propagation program.
C  Compile:  gfortran  -g -w -fbounds-check rodSurfFlaw.f  -o rodSurfFlaw
C  Usage:   rodSurfFlaw  scaleFactor <loadHistory >outputFile
C           also writes #crk= data to file    fadInput.rand 
C           (program will not overwrite previous copies of  fadInput.rand )
C           In BS7910:  M=fw=Mkm=Mkb=1,  Mkm,Mkb not used.

C  Program based on program "RCROK" Conle 1979 PhD thesis pg.128
C  The program is made availble to help students develop advances in crack
C  propagation software.

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
C Note that some subroutines included in this work are from GNU GPL licenced
C program:  http://fde.uwaterloo.ca/Fde/Calcs/saefcalc1.html
C
C Fork  rodSurfFlaw.f   vers. 3.06 from plateEdgeFlaw.f Sep.4 2013
C Fork  plateEdgeFlaw.f  from plateLongSurfFlaw.f  Aug. 31 2013
C       Edge is very similar except that we have a/W  rather than a/B. 
C       Also no weld feature code for Mkm and Mkb.  Also fw=M=1
C vers. 3.06   Other catch-up version fixes:
C       3.03 change "a" accumulator to  REAL*8
C            Add output of crack info to a binary file.      
C            Binary file read for  FAD  post processing.
C       2.30 Add Pm  and Pb  stress to the #crk= line printout
C       2.24 Activate SAVELEVEL code to reduce output.  
C vers. 2.23 In s/r getPeakLoads() print out "Filtered" when no changes.
C            to allow makereport1's grep  to function properly  Feb28 2013
C vers. 2.21 Correct: material file name read error.  Jan 6 2013
C Fork to  plateEdgeFlaw.f  from plateWeldflaw.f   Dec.31 2012 vers. 2.2
C       Remove S/R: readMmMb00, readMmMb90, fwSurfFlaw, readMkmMkb00,
C                   readMkmMkb90, getMmMb00a90, getMkmMkb00a90
C       Add routine getLongMmMb()

C vers. 2.1  Correction: multipy Y(sigma) by sqrt(pi*a)  Dec.29 2012
C vers. 2.0  Discretize lobj00, ldo00 lobj90,ldo90, dld00,dld90
C vers. 1.1  Divide damage or dadn  by 2.0 to make it per 1/2 cycle.
C            Add output for loads to be rainflow counted and used for
C            crack initiation analysis.
C vers. 1.00 runs ok.  
C vers. 0.91 Remove the old unused functions
C vers. 0.9 Fork: plateWeldflaw.f and plateWeldflaw+ss.f
C           In this version we have eliminated the local stress-strain stuff.
C           Thus it does the memory thing using deltaK00 and deltaK90 with no
C           tracking of the local ss.  Nov.17 2012
C           Local ss stuff is commented out with "Css"
C vers. 0.8 Create two parallel push-down lists. One for 00 and one for 90 deg.
C           crack tip.  Each has its own deltaK, stress, strain etc.
C           Each is P.D. list counted using its deltaK.
C vers. 0.6 replace Mkm and Mkb read s/r's. Add getMm and getMk s/r's.
C              (extract from testgetMmMkm.f test prog.)
C       Add fw  read and interpolate s/r. subr_getfw.f -> fwSurfFlaw()
C       Add peak pick  s/r. getPeakLoads()
C       Add getStress2Strain(), getLoad2StressStrain()
C       De-activate old functions from thesis version: XKP(),SMITH(),FLD(),DET()
C vers. 0.5 places read and interpolation for Mkm and Mkb into S/Rs. 
C           also adds read and interp. for Mm  and Mb


C In general:
C  1. Counts using Nominal "Load" = deltaK = fn{stsm(),stsb(),a,c,L...}
C  2. Computes Epsilon from Load.   ( not in this version) 
C  3. Computes Sigma from Epsilon   ( not in this version)
C  4. Computed damage as crack length

C  Material behaviour is stabilized.  No cyclic mean stress relaxation
C  or cyclic hardening or softening is performed in this version.  
C  Also no "disappearing" memory of prior deformation due to crack prop. is done.

C  Explanation of primary variables used=
C    tlimL, climL = 1 dimensional arrays of vectors containing
C       the push-ddown list loads.  TLIML is the tensile, CLIML is compr.
C  tlstr, clstr = the p.d. list associated Tens and Compr. strains
C  tlsts, clsts =               "          "               stresses
C  tldam, cldam =               "          "               1/2 cycle fat. damage
C  tlCracklength(), clCracklength() = the crack length when these occurred
C                   will be use to "disappear" memory events due to crack prop.
C  nptt,nptc =  push-down list pointers.
C  lorg, eorg, sorg  The load, strain, stress origin of a given 1/2 cycle
C  lobj, eobj, sobj         "           "      end point        "    "
C  dld, de,ds  =     the change in "    "  during a 1/2 cycle
C  totdam  =  total damage. In this program is crack length
C  ef = monotonic fracture strain
C  maxpd = maximum no of points possible in PD list
C  maxpoints = max no. of pts in Load history
C  maxnio = max no. reversals at which output is requested



C     These are used to store the equal spaced Snom. fitted data file.
      real StressAmp(1500), StrainAmp(1500), SnominalAmp(1500)
      integer*4 nmatdata,maxmatdata
C      real FractureStress,FractureStrain,Emod
      logical debugMat
      common/MATERIAL/ StressAmp,StrainAmp,SnominalAmp,
     &      snominalInterval,nmatdata,maxmatdata,debugMat
C     & FractureStress,FractureStrain,Emod

      real*4  deltaK(50),dadn(50),logdeltaK(50),logdadn(50),
     &        logdiffdK(50),logdiffdadn(50)
      logical debugDadn
      common/DADN/ deltaK,dadn,logdeltaK,logdadn,logdiffdK,logdiffdadn,
     &             ndadn,maxdadn,debugDadn


C     Matrices for storing Mm and Mb
      real*4 storeRBarMm(1000),storeRBarMb(1000)
      real*4 deltaRBarao2r,ao2rmin,ao2rmax
      integer maxstore
      common /LONGMMDATA/storeRBarMm,storeRBarMb,deltaRBarao2r,ao2rmin,
     &        ao2rmax,maxstore
      logical debugMm,lactivateMmMb

CC     Matrices for storing Mkm and Mkb
C      real*4 storeRBarMkm(1000),storeRBarMkb(1000)
C      real*4 deltaRBarzob,zobmin,zobmax,wLoB
C      integer maxstore2
C      real*4 Lweld
C      common /LONGMKDATA/storeRBarMkm,storeRBarMkb,deltaRBarzob,zobmin,
C     &        zobmax,maxstore2,wLoB
C      logical debugMk,lactivateMkmMkb

C     Matrix for storing fw
      logical debugfw,lactivatefw

      logical logiExist  !used in opening  fadInput.rand   file

C     These are the push-down list members:
      real*4 tlimL90(4000),climL90(4000)
Css     &  ,tlstr00(4000), clstr00(4000),tlsts00(4000),clsts00(4000)
Css     &  ,tlstr90(4000), clstr90(4000),tlsts90(4000),clsts90(4000)
     & ,tldam90(4000),cldam90(4000)
      integer itlnrev90(4000),iclnrev90(4000)
C     The damage (crack length) counters need to be real*8 because we
C     may be trying to add a very small crack increment to a large crack length.
      real*8 tltotCrk90(4000),cltotCrk90(4000)
      real*8 dtotdam90, daminc90, damold90
C      The reason for pushing the tltotCrk and tlnrev numbers are for future
C      disappearing memory due to crack extension. 
C
C      We could also reduce the no. of list entries by "binning"  the deltaK
C      values into, say 2000 finite values or bins. Note that its not a good
C      idea to bin (or discretize) other things such as stress-strain curves etc.
C      That would be like Wetzel's line segment memory, which causes all sorts of
C      problems near the bin "edges", mainly because of the Neuber solutions
C      interacting with the stress-strain solutions.
C      We will try the deltaK binning in a later version.

      character*300  inp300,jnp300,Cwebpage ! used to read in lines as chars.
      character*1    inpone(300),jnpone(300)
      character*10   inpten(30)
      equivalence (inpone(1), inpten(1), inp300)
      equivalence (jnpone(1),jnp300)

      character*10 name1, name2, stressunits
      character*10 strainunits, lifeunits, sorttype
      character*30 names30(10), firstfield, ctail
      character*30 cfiletype, Cdatatype
      character*80 argv, matfile, fwfile, kprimefile, dadnfile
      character*80 matname
      character*80 probType,histfile,dadnType,dadnTableFile

      integer*4 iargc,nargc

c     Load history storage.  If changing dimensions, also change same 
C                             in S/R getPeakLoads()
C      integer*4  sae(10000)
C     Stress membrane, bending, and total storage:
      real*4 stsm(5000),stsb(5000),ststot(5000)
      real*4 ststotmax,ststotmin,ststotwindow
      integer*4 nloads
      logical debugLoads
      common /LOADS/ stsm,stsb,ststot, ststotmax,ststotmin,ststotwindow,
     &               nloads,debugLoads


      real*4 ldo90,lobj90

C     Set up the discretized  deltaK table: basically the range of deltaK
C     is divided into  ndiscMax points, equally spaced by discDKinterval.
C     After each deltaK  values is computed, it is then discretized to
C     reduce the number of P.D. list entries, which normally gets very
C     large as the Mkm, Mkb factors become smaller as the cracks propagate
C     away from the weld or notch area.  Since we have discrete values of
C     deltaK we may as well discretize dadn also and, later, the 
C     corresponding values of local stress and strain.
      real*4 discDK(2000),discDadn(2000)   ! we don't really need discDK()
      real*4 discDKinterval
      integer ndiscMax        ! max storage for disc. values

c     Reversals at which PDList output occurs
      integer*4 nio(23)/2,3,4,5,8,70,100,200,500,1000,2000,5000,
     &    10000,20000,50000,100000,200000,300000,500000,1000000,
     &    5000000,10000000,100000000/
c     After last number in nio() put out data every nio(maxnio)
c     Note that as the crack propagates away from the Notch stress 
c     raiser that the PDLists may be VERY large: = maxpd
c     Thus we will need to place a limit on the no. PDlist printed out.

      logical prin
      logical stabl
      common /STAB/ stabl

      pi=3.1415927

      debugMat=.false.     !use to debug the Neuber and stress-strain portion
      debugDadn=.false.    ! use for the dadn  section
C      debugfw=.false.      ! use for the Finite width correct.
      debugMm=.false.      ! use for Mm00,Mb00,Mm90,Mb90 read and interpolate sections
C      debugMk=.false.      ! use for Mkm,Mkb,Mkm,Mkb, 00 and 90 read and interp. sections
      debugLoads=.false.   ! use to debug load read and peakpick section.

C     Storage area limits. Change this if you change above real* and int*
      maxmatdata=1500 ! Stress,strain,life store max
      maxdadn=50     ! da/dn data store max
      ndiscMax=2000   ! store for discretized dadn()

      maxloads=5000   ! max load storage
      maxpd=4000     ! Push down list storage max
      m=maxpd
      maxnio=23       ! Max store for Rev data output
      nioLast=nio(maxnio)

C     Test if -Wuninitilized  warning works in gfortran:
C      i=ifix(x)
C     Works!
    
C      Initilize some variables.
C      Origin and zero
      eorg=0.0
      sorg=0.0
      nact=0

      nptc90=0
      nptt90=0
      do 20 i=1,maxpd
        tlimL90(i)=0.    !peak in tensile nominal load or deltaK
        climL90(i)=0.    !peak in compr.  nominal load or deltaK
Css        tlstr90(i)=0. !peak in tens. strain
Css        clstr90(i)=0. !  "     comp.   "
Css        tlsts90(i)=0. !peak in tens. stress
Css        clsts90(i)=0. !  "     comp.   "
        tldam90(i)=0.    !damage or delta a for that ramp
        cldam90(i)=0.
        tltotCrk90(i)=0. ! total damage at that point
        cltotCrk90(i)=0.
        itlnrev90(i)=0.   ! rev. at which the above occured.
        iclnrev90(i)=0.

   20 continue



C---------------------------  Run time input data------------------
  184 continue
      write(6,185)
      write(0,185)
  185 format("# rodSurfFlaw.f vers. 3.06"/
     & "#Usage: rodSurfFlaw  scale <histfile  >outfile"/)

      nargc = iargc()
C     Note that in HP fortran the arg numbers must be changed by -1
C     and that iargc probably includes the "rodSurfFlaw" as an arg.
      if( nargc .ne. 1)then
        write(0,*)" rodSurfFlaw:  usage ERROR"
        write(6,*)" rodSurfFlaw:  usage ERROR"
        write(0,186)
        write(6,186)
  186  format(
     &   "#Usage:      rodSurfFlaw Scale <histfile ",
     &   " >outputFile"//
     &   "# Where Scale is the multiplier applied to all history",
     &   " points.(both the Membrane and Bending load/stresses) "/
     &   "# NOTE!: The above factor will be applied to the load"/
     &   "# history effectively AFTER the factors "/
     &   "# in the ""pdprop.env"" file are applied "
     &   "# Stopping now."
     &   )
        stop
      endif

C       The first arg is the history multiplication factor
        jvect=1
        call getarg(jvect,argv)
        read(argv,*,err= 178)scaleValue
        write(6,*)"#scaleValue= ",scaleValue
        if(scaleValue .eq. 0.0)go to 178

C        jvect=jvect+1
C        call getarg(jvect,argv)
C        read(argv,*,err=179)xmeanAdd
C        write(6,*)"#xmeanAdd= ",xmeanAdd
C       Both scaleValue and xmeanAdd have been read in successfully
        go to 190

C     Bad arguments in command line
  178 write(0,*)"# ERROR: bad scaleValue argument=",argv
      stop
C  179 write(0,*)"# ERROR: bad meanAdd argument= ",argv
C      stop

  190 continue
      write(6,191)
      write(0,191)
  191 format("#Opening pdprop.env file...")
C     Initilize the things to be read to checkable items
C     to make sure they have been entered.
      probType= " "
      iactivateMmMb= -1
      iactivateMkmMkb= -1
      iactivatefw= -1
C     The calling script will read the material file name and place the
C         results into the file "matfile"
      matfile= "matfile"

C     The history file will be read from standard input (device 5)
      histfile= "stdin"
      xmagfactorm= 0.    ! Set no data read flag
      xmagfactorb= 0.
      xmeanAddm=-1.0e20  ! Hopefully no one will ever use this :(
      xmeanAddb=-1.0e20  !   (unlikely since its basically MPa)
      
C     The calling script will read dadnType and create or copy to "dadnTable" file
      dadnType= " "
      dadnTablefile= "dadnTable"

C      Bthick= -1.
C      Width=  -1.
      radius= -1.

      azero= -1.
C      czero= -1.
C      Lweld= -1.
      maxHistReps= -1
      blockSkip= -9999.
      isavelevel= 0    !this is the default level if none is specified


C     Remove the following later if not needed
C     Older rcrock version items:
C      maxrevs=0
C      crackStart= -1.0
C      matfile= " "
C      fwfile= " "
C      xfwmult= 0.0
C      kprimefile= " "
C      xkpmult= 0.0
C      history will be read from standard input in this version.
C      histfile= " "
C      histScale= 0.
C      histMeanAdd= -9999.

C -----------   Open and read in the pdprop.env  file
C     In this file all lines should begin with a #  or are blank lines
      open(unit=10,file="pdprop.env")
      ninput=0
  800 continue
c     Loop back to here for next input line.

      read(10,"(a300)",end=380)inp300
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
Cdebug        write(0,*)"Shifted a comment: ",(inpone(j),j=1,inew)
C            Ok, now its a nice comment. Go play with it.
             go to 820
          else
C         The first non blank is not a #
          write(0,*)"# Skipped garbage in *.env file: line no.= ",ninput
          write(0,*)"# Text is= ",inp300
          write(6,*)"# Skipped garbage in *.env file: line no.= ",ninput
          write(6,*)"# Text is= ",inp300
          go to 800
          endif
  805    continue


      endif
  820 continue
C     First char was a #
C     pdprop.env file input line has a # in 1st col.  
C     See if it has a keyword.

      if(inpone(1).ne."#")then
C       something is bad in program
        write(0,*)" ERROR 820, sorry prog. messed up call ? admin"
        stop
      endif

C     Ok, its a nice comment.  Figure out if its a special tag.
      read(inp300,*)firstfield


      if(firstfield .eq."#TYPE=" .or.
     &   firstfield .eq."#Type=" .or.
     &   firstfield .eq."#type=" )then
         read(inp300,*) firstfield, probType
         write(6,833)probType
  833    format("#TYPE= ",a80)
         go to 800
      endif

      if(firstfield .eq."#ACTIVATE_MmMb=" .or.
     &   firstfield .eq."#Activate_MmMb=" .or.
     &   firstfield .eq."#activate_MmMb=" )then
         read(inp300,*) firstfield, iactivateMmMb
         write(6,835)iactivateMmMb
  835    format("#iactivateMmMb= ",i3)
         lactivateMmMb=.true.
         if(iactivateMmMb.eq.0)lactivateMmMb=.false. ! turn off
         go to 800
      endif

C      if(firstfield .eq."#ACTIVATE_MkmMkb=" .or.
C     &   firstfield .eq."#Activate_MkmMkb=" .or.
C     &   firstfield .eq."#activate_MkmMkb=" )then
C         read(inp300,*) firstfield, iactivateMkmMkb
C         write(6,837)iactivateMkmMkb
C  837    format("#iactivateMkmMkb= ",i3)
C         lactivateMkmMkb=.true.
C         if(iactivateMkmMkb.eq.0)lactivateMkmMkb=.false. ! turn off
C         go to 800
C      endif

C      if(firstfield .eq."#ACTIVATE_fw=" .or.
C     &   firstfield .eq."#ACTIVATE_FW=" .or.
C     &   firstfield .eq."#ACTIVATE_Fw=" .or.
C     &   firstfield .eq."#Activate_fw=" .or.
C     &   firstfield .eq."#activate_fw=" )then
C         read(inp300,*) firstfield, iactivatefw
C         write(6,838)iactivatefw
C  838    format("#iactivatefw= ",i3)
C         lactivatefw=.true.
C         if(iactivatefw.eq.0)lactivatefw=.false. ! turn off
C         go to 800
C      endif


C      if(firstfield .eq."#B=" .or.
C     &   firstfield .eq."#b=" )then
C         read(inp300,*) firstfield, Bthick
C         write(6,844)Bthick
C  844    format("#B= ",E14.7," # Thickness, mm.")
C         go to 800
C      endif

C      if(firstfield .eq."#W=" .or.
C     &   firstfield .eq."#w=" )then
C         read(inp300,*) firstfield, Width
C         write(6,845)Width
C  845    format("#W= ",E14.7," # width, mm.")
C         go to 800
C      endif

      if(firstfield .eq."#ri=" .or.
     &   firstfield .eq."#Ri" .or.
     &   firstfield .eq."#RI" )then
         read(inp300,*) firstfield, radius
         write(6,846)radius
  846    format("#ri= ",E14.7," # pipe radius, mm.")
         go to 800
      endif

      if(firstfield .eq."#azero=" .or.
     &   firstfield .eq."#AZERO=" )then
         read(inp300,*) firstfield, azero
         write(6,847)azero
  847    format("#azero= ",E14.7," # initial crack depth, mm.")
         go to 800
      endif

C      if(firstfield .eq."#czero=" .or.
C     &   firstfield .eq."#Czero=" .or.
C     &   firstfield .eq."#CZERO=" )then
C         read(inp300,*) firstfield, czero
C         write(6,848)czero
C  848    format("#czero= ",E14.7," # 1/2 initial crack width, mm.")
C         go to 800
C      endif

C      if(firstfield .eq."#L=" .or.
C     &   firstfield .eq."#l=" )then
C         read(inp300,*) firstfield, Lweld
C         write(6,849)Lweld
C  849    format("#L= ",E14.7," # weld feature width")
C         go to 800
C      endif

      if(firstfield .eq."#MATERIAL=" .or.
     &   firstfield .eq."#Material=" .or.
     &   firstfield .eq."#material=" )then
         read(inp300,*) firstfield, matname
         write(6,850)matname
  850    format("#MATERIAL= ",a80)
         go to 800
      endif
C
      if(firstfield .eq."#MAGFACTOR_m=" .or.
     &   firstfield .eq."#Magfactor_m=" .or.
     &   firstfield .eq."#magfactor_m=" )then
         read(inp300,*) firstfield, xmagfactorm
         write(6,852)xmagfactorm
  852    format("#MAGFACTOR_m= ",e14.7)
         go to 800
      endif

      if(firstfield .eq."#MAGFACTOR_b=" .or.
     &   firstfield .eq."#Magfactor_b=" .or.
     &   firstfield .eq."#magfactor_b=" )then
         read(inp300,*) firstfield, xmagfactorb
         write(6,854)xmagfactorb
  854    format("#MAGFACTOR_b= ",e14.7)
         go to 800
      endif

      if(firstfield .eq."#MEANADD_m=" .or.
     &   firstfield .eq."#Meanadd_m=" .or.
     &   firstfield .eq."#meanadd_m=" )then
         read(inp300,*) firstfield, xmeanAddm
         write(6,855)xmeanAddm
  855    format("#MEANADD_m= ",e14.7)
         go to 800
      endif

      if(firstfield .eq."#MEANADD_b=" .or.
     &   firstfield .eq."#Meanadd_b=" .or.
     &   firstfield .eq."#meanadd_b=" )then
         read(inp300,*) firstfield, xmeanAddb
         write(6,856)xmeanAddb
  856    format("#MEANADD_b= ",e14.7)
         go to 800
      endif
C
      if(firstfield .eq."#MAXREPS=" .or.
     &   firstfield .eq."#Maxreps=" .or.
     &   firstfield .eq."#maxreps=" )then
         read(inp300,*) firstfield, maxHistReps
         write(6,891)maxHistReps
  891    format("#MAXREPS= ",i20)
         go to 800
      endif


C     This is just a label check. The data itself will be lifted from the
C     standard file given in char variable: dadnTableFile
      if(firstfield .eq."#DADN=" .or.
     &   firstfield .eq."#Dadn=" .or.
     &   firstfield .eq."#dadn=" )then
         read(inp300,*) firstfield, dadnType
         write(6,892)dadnType
  892    format("#DADN= ",a80)
         go to 800
      endif

      if(firstfield .eq."#BLOCKSKIP=" .or.
     &   firstfield .eq."#Blockskip=" .or.
     &   firstfield .eq."#blockskip=" )then
         read(inp300,*) firstfield, blockskip
         write(6,893)blockskip
  893    format("#BLOCKSKIP= ",E14.7)
         go to 800
      endif

      if(firstfield .eq."#SAVELEVEL=" .or.
     &   firstfield .eq."#Savelevel=" .or.
     &   firstfield .eq."#savelevel=" )then
         read(inp300,*) firstfield, isavelevel
         write(6,894)isavelevel
  894    format("#SAVELEVEL= ",i2)
         go to 800
      endif


  350 continue
C     None of the above?  Then it must be a plain old
C     comment line.   Write it out too.
C       Count backwards and see where the last char is
        do 360 i=1,300
          j=300-(i-1)
          if(inpone(j).ne." ")then
C           found last char
            lastloc=j
            go to 362
          endif
  360   continue

  362   continue
      write(6,"(300a1)")(inpone(i),i=1,lastloc)
C     Go read another line
      go to 800

C     End of pdprop.env file reached
  380 continue
C     Check if critical items have been read in.
      close(unit=10)
      istop=0
      if(probType .eq. " ")then
        write(0,*)"#ERROR: #Type= not found."
        write(6,*)"#ERROR: #Type= not found."
        istop=1
      endif
      if(iactivteMmMb .eq. -1)then
        write(0,*)"ERROR: #ACTIVATE_MmMb= not found"
        write(6,*)"ERROR: #ACTIVATE_MmMb= not found"
        istop=1
      endif
C      if(iactivteMkmMkb .eq. -1)then
C        write(0,*)"ERROR: #ACTIVATE_MkmMkb= not found"
C        write(6,*)"ERROR: #ACTIVATE_MkmMkb= not found"
C        istop=1
C      endif
C      if(iactivtefw .eq. -1)then
C        write(0,*)"ERROR: #ACTIVATE_fw= not found"
C        write(6,*)"ERROR: #ACTIVATE_fw= not found"
C        istop=1
C      endif
      if(iactivateMmMb .eq. 0 )then
        write(0,383)
        write(6,383)
  383   format("#WARNING!!!  WARNING!!!  Neither Mm, Mb ",
     &         " are activated. Result: Mm=Mb=1 ")
C       this will not stop the program.
      endif
C      if(iactivatefw .eq. 0)then
C        write(0,384)
C        write(6,384)
C  384   format("#WARNING!!! WARNING!!! fw deactivated. Result fw=1.0")
C      endif

C      if(Bthick .eq. -1.0)then
C        write(0,*)"ERROR:  #B=  not found"
C        write(6,*)"ERROR:  #B=  not found"
C        istop=1
C      endif
C      if(Width .eq. -1.0)then
C        write(0,*)"ERROR:  #W=  not found"
C        write(6,*)"ERROR:  #W=  not found"
C        istop=1
C      endif
      if(radius .eq. -1.0)then
        write(0,*)"ERROR:  #ri=  not found"
        write(6,*)"ERROR:  #ri=  not found"
        istop=1
      endif
      if(azero .eq. -1.0)then
        write(0,*)"ERROR:  #azero=  not found"
        write(6,*)"ERROR:  #azero=  not found"
        istop=1
      endif
C      if(czero .eq. -1.0)then
C        write(0,*)"ERROR:  #czero=  not found"
C        write(6,*)"ERROR:  #czero=  not found"
C        istop=1
C      endif
C      if(Lweld .eq. -1.0)then
C        write(0,*)"ERROR:  #L=  not found"
C        write(6,*)"ERROR:  #L=  not found"
C        istop=1
C      endif
C 
      if(xmagfactorm .eq. 0.0)then
        write(0,*)"ERROR:  #MAGFACTOR_m=  not found"
        write(6,*)"ERROR:  #MAGFACTOR_m=  not found"
        istop=1
      endif
      if(xmagfactorb .eq. 0.0)then
        write(0,*)"ERROR:  #MAGFACTOR_b=  not found"
        write(6,*)"ERROR:  #MAGFACTOR_b=  not found"
        istop=1
      endif
      if(xmeanAddm .eq. -1.0e20)then
        write(0,*)"ERROR:  #MEANDADD_m=  not found"
        write(6,*)"ERROR:  #MEANDADD_m=  not found"
        istop=1
      endif
      if(xmeanAddb .eq. -1.0e20)then
        write(0,*)"ERROR:  #MEANDADD_b=  not found"
        write(6,*)"ERROR:  #MEANDADD_b=  not found"
        istop=1
      endif

      if(maxHistReps .eq. -1)then
        write(0,*)"#ERROR: #MAXREPS=  not found"
        write(6,*)"#ERROR: #MAXREPS=  not found"
        istop=1
      endif
      if(blockSkip .eq. -9999. )then
        write(0,*)"#ERROR: #BLOCKSKIP=  not found"
        write(6,*)"#ERROR: #BLOCKSKIP=  not found"
        istop=1
      endif

      if(istop.eq. 1)then
         write(0,*)"# Stopping..."
         write(6,*)"# Stopping..."
         stop
      endif




C------------- Read in the Material Fitted Fatigue file Table-------------------

      snominalInterval= 0.   ! The delta between the points
      write(0,701)
      write(6,701)
  701 format("# #Opening matfile: Snominal, Stain, Stress Table ..."/
     &       "# #filename= matfile")
      open(unit=10,file=matfile)

      ninput=0
  700 continue
c     Loop back to here for next input line.

      read(10,"(a300)",end=750)inp300
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
  715       format("# ERROR matfile: input line no. ",I5,
     &      " not a # and not a number. Fix program rdMatfile.f ")
            stop
         endif
C        Ok, its a number
         nmatdata=nmatdata+1
         if(nmatdata .gt. maxmatdata)then
           write(0,716)nmatdata
           write(6,716)nmatdata
  716      format("# Error too many data points:",i5,
     &        " recompile rodSurfFlaw.f or reduce matfile data")
           stop
         endif
         read(inp300,*)SnominalAmp(nmatdata),StrainAmp(nmatdata), 
     &                 StressAmp(nmatdata)
         write(6,720)SnominalAmp(nmatdata),StrainAmp(nmatdata),
     &                 StressAmp(nmatdata),nmatdata
  720    format("#matfile ",f10.3,1x,f8.5,1x,f10.3,1x,i5)
         go to 700

      endif
  730    continue
C        Look for the interval comment line:
         read(inp300,*)firstfield
         if(firstfield.eq."#SnominalINTERVAL=")then
           read(inp300,*)firstfield,snominalInterval
           write(6,732)snominalInterval
  732      format("#matfile #Found SnominalINTERVAL= ",e14.7)
         endif
C          The interval is read in e14.7  but the Snominal data is f8.3  which
C          may cause a round-off problem if we ever go from strain or stress to Snom.
         go to 700

  750    continue   ! end of file comes here
C        All data is in. Close file
         close(unit=10)

C        Check if the nominal interval was found
         if(snominalInterval.eq.0.)then   !was not read in. Compute it instead.:
C          Assuming the first point was 0 0  0
C          then we should have about nmatdata= 1001 points in these memories.
C          Thus the interval is:
           snominalInterval= SnominalAmp(nmatdata)/(nmatdata-1)
C          This interval will be used in the interpolation process when finding
C          stress and strain given Snominal.
           write(6,752)snominalInterval
  752      format("#matfile #Computed #SnominalInterval= ",E14.7)
         endif


C--------------------  Read in the dadn data table -------------------
C       Table e.g.:
C        #NAME= A36merged  data from Newman, Haddad, Klingerman
C        #Data digitized from merged graph in file a36+1015dadn-2.png
C        #
C        #MPa*Sqrt(mm)  dadn mm/cycle
C        150.216  9.62054e-08
C        176.983  4.5623e-07
C        220.235  1.16017e-06
C        287.484  3.22409e-06
C        433.167  1.06976e-05
C        763.741  7.55681e-05
C        1240.59  0.000852041
C        1471.68  0.0033073
C        1675.69  0.0107468

C        Read in the table and convert to log-log  to store and use.
C        Computations will be done using log log co-ord interpolation.
      write(0,601)
      write(6,601)
  601 format("# Opening DeltaK   dadn   Table, file=dadnTable")
      open(unit=10,file=dadnTableFile)

      ninput=0
      ndadn=0
  600 continue     !Loop back to here for next input line.
      read(10,"(a300)",end=650)inp300
      ninput=ninput+1

      if(inp300.eq." ")then    ! Check for blank line
C        write(6,"(a1)")" "
        go to 600
      endif

      if(inpone(1).ne."#")then
C        We may have a data line, or someone screwed up and put the # later in line.
C        See if 1st char is a # later in line
         do 605 i=1,300
           if(inpone(i).eq." ") go to 605
C           No? 1st non-blank found, if # its not a data line
            loc=i
           if(inpone(i).eq."#") then
C            Shift the whole mess over to begining of field
             inew=0
             do 603 j=loc,300
              inew=inew+1
              inpone(inew)=inpone(j)
  603        continue
Cdebug        write(0,*)"Shifted a comment: ",(inpone(j),j=1,inew)
C            Ok, now its a nice comment. Go play with it.
             go to 630
           else
C            The first non blank is not a #
             go to 610
           endif
  605      continue
           go to 600  !assume garbage in line.

  610      continue
C        1st non blank is not a #  Check if its a number.
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
            write(0,615)ninput
            write(6,615)ninput
  615       format(" ERROR dadnfile: input line no. ",I5,
     &      " not a # and not a number. Fix data in file dadnTable ")
            stop
         endif
C        Ok, its a number
         ndadn=ndadn+1
         if(ndadn .gt. maxdadn)then
           write(0,616)ndadn
           write(6,616)ndadn
  616      format(" Error too many data points:",i5,
     &       " recompile rodSurfFlaw.f or reduce dadnTable data")
           stop
         endif
         read(inp300,*)deltaK(ndadn),dadn(ndadn)
         logdeltaK(ndadn)= ALOG10(deltaK(ndadn))
         logdadn(ndadn)=ALOG10(dadn(ndadn))
         if(ndadn .gt.1)then  ! compute the diffs between points.
C          This will help speed up the interpolation calculation
           logdiffdk(ndadn)=logdeltaK(ndadn)-logdeltaK(ndadn-1)
           logdiffdadn(ndadn)=logdadn(ndadn)-logdadn(ndadn-1)
         endif
         write(6,620)deltaK(ndadn),dadn(ndadn),
     &              logdeltaK(ndadn),logdadn(ndadn),
     &              logdiffdk(ndadn),logdiffdadn(ndadn), ndadn
  620    format("#dadnTable ",6(e14.7,1x),i5)
         go to 600

      endif
  630    continue
C        All comment lines are just written out, or deleted
         go to 600

  650    continue   ! end of file comes here
C        All data is in. Close file
         logdiffdk(1)=0.  ! these will not be used, but set to 0. anyway
         logdiffdadn(1)=0.
         close(unit=10)
         write(6,652)ndadn
         write(0,652)ndadn
  652    format("#dadn data input completed. ndadn= ",i5)



C-------------- Initilize the Mm,Mb, Mkm,Mkb storage -----------------------------

      ivec=0  ! =0 initilize,  =1  interpolate
      ao2r=0.2   !  set to something
      if(lactivateMmMb)then
         call getRBarMmMb(ivec,iret,ao2r,xMm,xMb)
         if(iret.ne.0)then
           write(0,*)"# Error from getRBarMmMb() Init., Stopping..."
           stop
         endif
      endif

C      wLoB=Lweld/Bthick  !wLoB is in common
C      ao2r=0.2   !  set to something
C      if(lactivateMkmMkb)then
C         call getRBarMkmMkb(ivec,iret,ao2r,xMkm,xMkb)
C         if(iret.ne.0)then
C           write(0,*)"#Error: getRBarMkmMkb() Initilize. Stopping."
C           stop
C         endif
C      endif

      izout=1

C     initilize the functions
      nrev=0
Cold      ds=det(0.,nrev)
Cold      de=fld(0.,nrev)
Cold      x=SMITH(0., 0.,0., 0.,0., 0.,0)
Cold      x=xkp(0., 0)
Cold      de=0.0
c
C------------------Get Load history file------------------------------
C     Expected format:
C     #Comment 1
C     #Comment line 2
C     #time     MembraneLoad    BendingLoad
C       0.0       xxx.xxx         yyy.yyy
C       z.z       xxx.xxx         yyy.yyy

C     The time column info is not really used.
C     Membrane and Bending loads will be altered by the #MAGFACTORm  etc lines
C     in the pdprop.env file and the calling args.

      write(0,901)
      write(6,901)
  901 format("# Getting Load history file from std.input ")
C      open(unit=10,file=histfile)

      ninput=0
      nloads=0
  900 continue     !Loop back to here for next input line.
      read(5,"(a300)",end=950)inp300
      ninput=ninput+1

      if(inp300.eq." ")then    ! Check for blank line
C        write(6,"(a1)")" "
        go to 900
      endif

      if(inpone(1).ne."#")then
C        We may have a data line, or someone screwed up and put the # later in line.
C        See if 1st char is a # later in line
         do 905 i=1,300
           if(inpone(i).eq." ") go to 905
C           No? 1st non-blank found, if # its not a data line
            loc=i
           if(inpone(i).eq."#") then
C            Shift the whole mess over to begining of field
             inew=0
             do 903 j=loc,300
              inew=inew+1
              inpone(inew)=inpone(j)
  903        continue
Cdebug        write(0,*)"Shifted a comment: ",(inpone(j),j=1,inew)
C            Ok, now its a nice comment. Go play with it.
             go to 930
           else
C            The first non blank is not a #
             go to 910
           endif
  905      continue
           go to 900  !assume garbage in line.

  910      continue
C        1st non blank is not a #  Check if its a number.
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
            write(0,915)ninput
            write(6,915)ninput
  915       format(" ERROR stdin: input line no. ",I5,
     &      " not a # and not a number. Fix data in load history file")
            stop
         endif
C        Ok, its a number
         nloads=nloads+1
         if(nloads .gt. maxloads)then
           write(0,916)nloads
           write(6,916)nloads
  916      format(" Error too many data points:",i5,
     &      "recompile rodSurfFlaw.f or reduce load history data")
           stop
         endif
         read(inp300,*,err=925)xtime,xmembrane,xbending
C        Adjust to mag. and mean shift in *.env file
         xmembrane=xmembrane*xmagfactorm + xmeanAddm
         xbending= xbending*xmagfactorb  + xmeanAddb
C        Now scale by arg1  in the command that calls this program
         stsm(nloads)= xmembrane*scaleValue
         stsb(nloads)= xbending*scaleValue
C        Find the total stress for Rainflow counting
         ststot1= stsm(nloads) + stsb(nloads)
         if(nloads .eq. 1)then
            ststotmax=ststot1
            ststotmin=ststot1
         endif
         if(ststot1 .gt. ststotmax) ststotmax=ststot1
         if(ststot1 .lt. ststotmin) ststotmin=ststot1
         ststot(nloads)=ststot1
         write(6,920)stsm(nloads),stsb(nloads),ststot(nloads),nloads
  920    format("#history.org ",3(f7.1,1x),i6)
         go to 900

  925    write(0,926)nloads,inp300
         write(6,926)nloads,inp300
  926    format("#READ ERROR: history file: data line no.= ",i9/
     &   "#LINE= ",a300/ "# Stopping now.")
         stop

      endif
  930    continue
C        Yes, first char was a #
C        All comment lines are just written out, or deleted
         go to 900

  950    continue    ! end of file comes here
         if(nloads.eq.0)then
            write(6,951)
            write(0,951)
  951       format("#Error: histfile: contains no data !")
            stop
         endif
C        All data is in. Close file
C         close(unit=10)
         write(6,952)nloads,ststotmax,ststotmin
         write(0,952)nloads,ststotmax,ststotmin
  952    format("#history #Data input completed. nloads= ",i5/
     &   "#history #StotMax= ",f7.1/"#history #StotMin= ",f7.1/
     &   "#history #Where Stot = Smembrane + Sbending"//)
         

C     stsm(5000),stsb(5000), ststotal(5000) are the history stores
C-----------------We must now eliminate non-reversal points------------------
C     created by the addition of Smembrane and Sbending
C     i.e.: If Sm and Sb are not proportional (viewed over time) then when they
C     are superimposed to get the "total" stresss near the crack, some non-reversal
C     points may be in the history. These must be eliminated before applying the
C     push-down list material memory rules.  Also there may be created some small
C     1/2 cycles that are not worth running. (OK, its a value judgement) This
C     elimination window is computed as a %tage of max and min stress:

      ststotwindow= 0.02*(ststotmax-ststotmin)
      write(0,954)ststotwindow
      write(6,954)ststotwindow
  954 format("# Eliminating non-reversal points from history..."/
     &"#history #Elimination window = 2% of",
     &" (StressMax-StressMin)=",f7.2/ 
     &"#history #Any 1/2 cylce smaller than this will be eliminated...")

      ivec=1
      iret=0
      call getPeakLoads(ivec,iret)
      if(iret.ne.0)then
        write(0,*)"#ERROR: getPeakLoads() iret=",iret," Stopping."
        stop
      endif
      

C------------------ Create the discretized calculation tables------------------
      xdeltaK=0.
      discDK(1)=0.
      discDadn(1)=0.
      discDKinterval=deltaK(ndadn)/(ndiscMax-1) ! deltaKmax =deltaK(ndadn)
      i=1
      write(6,440)xdeltaK,discDadn(i),i

      do 450 i=2,ndiscMax
        xdeltaK=xdeltaK+discDKinterval
        discDK(i)=xdeltaK
        call getCracks(iret,xdeltaK,ydadn)
C           an error here in getCracks will stop program in getCracks
        discDadn(i)=ydadn
        write(6,440)xdeltaK,ydadn,i
  440   format("#discrete ",e14.7,1x,e14.7,i6)
  450 continue
      deltaKmax=xdeltaK
      write(0,455)ndiscMax
      write(6,455)ndiscMax
  455 format("#discrete #Discretized deltaK vs dadn done. no. pts=",i6)

C----------------------binary output file-------------------------------------
C     Open the binary (direct access) record file for data output
C     Get a filename for the standard input stream (doesnt work)
C      inquire(unit=5, NAME=stdinFileName)
      write(0,460)
      write(6,460)
  460 format(/"# Opening random access output file:  fadInput.rand ...")
      inquire(file="fadInput.rand", EXIST= logiExist)
      if(logiExist)then
        write(0,*)"#ERROR: an old copy of file  fadInput.rand  exists."
        write(0,*)"#   You need to rename it or remove it before we"
        write(0,*)"#   can run the simulation.  Stopping now..."

        write(6,*)"#ERROR: an old copy of file  fadInput.rand  exists."
        write(6,*)"#   You need to rename it or remove it before we"
        write(6,*)"#   can run the simulation.  Stopping now..."
        stop
      endif

C     Ok,  previous copies do not exist. open it.
      open(unit=60, file="fadInput.rand", access="direct",
     &     form= "unformatted", status= "new", recl= 36 )
      write(0,462)
      write(6,462)
  462 format("#Random Access output file: fadInput.rand   opened.")
      nrecord=1  !this will be incremented by 1 upon 1st use.
C                 The 1st record is written at the end of test.



C---------------------- Start Cycling--------------------------------------------
C     Load or nominal stress is used to model material memory. 
C     Kmax and Kmin are used to run the memory.
C     Program flow and logic are pretty much the same as concepts of
C     original program  rcrock.f from thesis.
      nblk=1
      nrev=0


      lobj90=0.
      ldo90=0.
      eobj90=0.
      sobj90=0.

      iupdown90=+1


      dtotdam90=azero  !crack length "a" for 90deg depth crack point
      damold90=0.
      daminc90=0.
C     note that there are two damage accumulators:  totdam90  and dtotdam90
C     They both contain the same numbers but dtotdam90   is REAL*8  in size
C     in order to better accumulate very small damage numbers  when the 
C     total damage is large.

      if(.not.lactivateMmMb)then ! See if Mm and Mb are shutdown
        xMm90=1.0
        xMb90=1.0
      endif

C      if(.not.lactivateMkmMkb)then ! See if Mkm and Mkb are shutdown
C        xMkm90=1.0
C        xMkb90=1.0
C      endif
       xMkm90=1.0
       xMkb90=1.0

C      if(.not.lactivatefw)xfw=1.0  ! See if FiniteWidth Corr. is shutdown
       xfw=1.0    ! See BS7910 2005  M.3.5

      write(6,119)   ! write out the header for crack data
  119 format("#crk= #   NREV    a 90      NBLK    NACT    ",
     &       "lobj90    xMm90,xMb90, nptt90")
C      write(6,120) nrev,totdam90,nblk,nact,lobj90,xMm90,xMb90,
C     &             nptt90

 3000 continue  ! =================Top of Cycling loop====================
      ldo90=lobj90
Css      eo90=eobj90
Css      so90=sobj90
      
      nrev=nrev+1
      nact=nact+1
      if(nptt90.eq.maxpd .or. nptc90.eq.maxpd )then
         write(0,118)maxpd,nptt90,nptc90
         write(6,118)maxpd,nptt90,nptc90
 118     format("#Error: one or more of the PushDown list counters >",
     &   "maxpd= ",i6," : nptt90=",i6," nptc90=",i6
     &   )
         stop
      endif

cccccccccccccccccccccccccccccccccccccccccccccc  Block Repeat
C     The user must make certain that end of block is  compatible with 
C     the begin of block.  The program will start at 0 stress(load) but
C     the history does not need to include this start point. The last
C     point in the history and the first and 2nd points MUST form reversals.
C    E.g:     2nd_last_pt    last_pt   1st_pt   2nd_pt
C                 +100        -90       90       -90       is allowed
C                 +100        +50       0         +50      is NOT allowed
      if(nrev .eq. 1 .and. stsm(1) .eq.0. .and. 
     &   stsb(1) .eq. 0.0) nact=2  !zero range, skip this half cycle
      if(nact .le. nloads) go to 3001
C       end of block, increment the blk count and reset nact=1
        nblk=nblk+1
C       daminc and damold s  are used to skip a block if damage is same
        daminc90=dtotdam90 - damold90
        damold90=dtotdam90
        nact=1
 3001 continue

Cold      lnomStress= ststot(nact)
Cold      lobj= lobj*xkp(totdam,nrev)
      ivec=1   !get interpolate data if req.
      totdam90=sngl(dtotdam90)

      ao2r=totdam90/(2.0*radius)
      if(ao2r .ge. 0.6)then
         write(0,3002)totdam90,radius,ao2r
         write(6,3002)totdam90,radius,ao2r
 3002    format("#END: a/2r >= 0.6 : a= ",e14.7," radius= ",e14.7,
     &          " a/2r= ",e14.7)
         goto 9000   !stop
      endif

      if(lactivateMmMb)then ! if not, all are equal to 1.0 (set above)
        call getRBarMmMb(ivec,iret,ao2r,xMm90,xMb90)
        if(iret .ne. 0)then
          goto 9000  !stop
        endif
      endif

C      if(lactivateMkmMkb)then ! if not, all are equal to 1.0 (set above)
C        call getRBarMkmMkb(ivec,iret,ao2r,xMkm90,xMkb90)
C        if(xMkm90 .lt. 1.0)xMkm90=1.0
C        if(xMkb90 .lt. 1.0)xMkb90=1.0
C        if(iret .ne. 0)then
C          goto 9000  !stop
C        endif
C      endif

C        call fwSurfFlaw(ivec,iret,cow,ao2r,xfw)

C     Apply factors to stresses  as per BSstandard:
C     Assume: Bulge M=1,   km=1  ktm=1  ktb=1
      stsMembrane=stsm(nact)
      stsBending=stsb(nact)
      rootPiA=sqrt(pi*totdam90)
C      lobj90=xfw*(xMkm90*xMm90*stsMembrane + xMkb90*xMb90*stsBending )
      lobj90=(xMm90*stsMembrane + xMb90*stsBending )
      lobj90=lobj90*rootPiA
C     discretize  lobj90 and lobj00
      lobj90= float(ifix(lobj90/discDKinterval)) *discDKinterval
    
Cdebug      write(6,*)"#db: lobj90,lobj00,xfw,Mkm90,Mm90 = ",
Cdebug     &           lobj90,lobj00,xfw,xMkm90,xMm90
      if(lobj90.eq.ldo90 )then
C       Skip, its a nothing reversal.  It should probably be a delta
C       check.  (but what if one is zero, but the other is not?)
        write(0,3007)nrev,nblk,nact
        write(6,3007)nrev,nblk,nact
 3007     format("#Warning: Same target as previous lobj90=ldo90 "
     &    ," at nrev=",i9," nblk=",i9," nact=",i9/
     &     "# Skipping this stress...")
        go to 3000
      endif

C     Crack length variables may have to be real*8  because adding
C     a very small da increment may not show up if crack is big.
      totdam90=sngl(dtotdam90)
      if(isavelevel .gt. 0)then
        write(6,120) nrev,totdam90,nblk,nact,lobj90,xMm90,xMb90,
     &             nptt90
  120   format("#crk=",i9,1x,e14.7,i8,1x,i6,1x,f7.2,
     &             4(1x,f8.5),i4)
      endif

C     We always write to the binary out file. RecSize is 12*4 =48
        nrecord=nrecord+1
        write(60,rec=nrecord)nrev,totdam90,nblk,nact,
     &             lobj90,xMm90,xMb90, 
     &             stsMembrane,stsBending

      prin=.false.
      if(nrev .eq. nio(izout)) prin =.true.
      if(nrev .gt. maxHistReps) prin=.true.
      if(.not. prin) go to 3050
      if(nptt90 .eq. 0  .or. nptc90 .eq. 0) go to 3050
      izout=izout+1
      if(izout.gt.maxnio)then
        izout=maxnio
        nio(izout)=nio(izout)*2
      endif
      if(nptt90 .gt. 10 .or. nptc90 .gt. 10 .or. nblk.gt.10) go to 3050

C     Often too much output.  Turn this off if isavelevel <3
      if(isavelevel .lt. 3) go to 3050


      write(6,121) nrev,totdam90,nblk,nact
  121 format("#PD Before REV= ",i9,"  Cracks: a= ",E14.7
     &  ,3X,"NBLK=",i6," NACT= ",I6
     & /"#PD: ",4x,"CLDAM90",3x,"CLIML90" )
      do 3039 ipr=1,nptc90
 3039 write(6,122) ipr,cldam90(ipr),climL90(ipr)
  122 format("#PD:",1x,i2,1x,E12.5,2x,F9.2)
   
      write(6,123)
  123 format("#PD:",5x,"TLDAM90",4x,"TLIML90")
      do 3045 ipr=1,nptt90
      write(6,122) ipr,tldam90(ipr),tlimL90(ipr)
 3045 continue

 3050 continue
C ----------------------------------------------------------------
C
C     Check for MAX. history repeats
      if(nblk .gt. maxHistReps)then
        write(0,130) nblk,totdam90,nrev
        write(6,130) nblk,totdam90,nrev
  130   format("# Max no. of History Reps. Reached: nblk= ",i10/
     &      "#TOTDAM90= ",E14.7," nrev= ",i9)
        goto 9000  !stop
      endif

      dld90=abs(lobj90-ldo90)
Cdebug      write(6,*)"#stsm(),stsb(),nact=",stsm(nact),stsb(nact),nact
Cdebug      write(6,*)"#dld90,lobj90,ldo90=",dld90,lobj90,ldo90,nrev
C      if(dld90.le. 0.05)then   !check for no actual half cycle
C     There is a slight possibility that one could have a 0 ramp in either
C     the  00 and/or the 90 deg direction, but we are not going to program 
C     for this now.

C     Going up or down ?   -----------------Crack Direction 90 deg------------
      if(lobj90 .lt. ldo90) go to 2100

C************ Going UP,  Tensile Direction **************************
C     Are we on the monotonic curve?
 1100 if(nptt90 .eq. 0) go to 1500
c     No. Has a exceedence  occurred?
 1110 if(lobj90 .gt. tlimL90(nptt90)) go to 1120

c     Is this Same as previous amplitude?
 1130 if(lobj90 .eq. tlimL90(nptt90) .and. nptt90 .ne. 1) go to 1135
      go to 1350

c     Is a loop being closed without the monotonic?
 1120 if(nptt90 .ne. 1) go to 1400

c     Is closure in connection  with the monotonic?
 1150 if(nptc90 .eq. 2) go to 1170

c     Check for impossible
 1160 if(nptc90 .eq. 1) go to 1165
      idump=1160
      go to 9999

C     Same load level as previous level.
 1135 continue
Css      dam90=SMITH( (lobj-ldo),tlsts(nptt90),clsts(nptc90),tlstr(nptt90),
Css     &            clstr(nptc90),totdam,nrev)
C     lobj90 and ldo90 exist. They are the deltaK values. Use them
      dld90=lobj90-ldo90
      if(dld90.le.0.)then
        write(6,*)"#Error: 1135:dld90,lobj90,ldo90 = ",
     &            dld90,lobj90,ldo90,nrev
      endif
Cvers1      call getCracks(iret,dld90,dam90)
Cvers1      if(iret .ne. 0)goto 9000  !   Error returned. goto stop
      jpoint=dld90/discDKinterval+1
      if(jpoint.gt.ndiscMax)then
       write(0,1136)nrev,nblk,nact,dld90,deltaKmax
       write(6,1136)nrev,nblk,nact,dld90,deltaKmax
 1136  format("#Fracture: dld90.gt.deltaKmax :",i9,1x,i6,1x,i6,
     &         2(1x,e14.7))
       go to 9000   !stop
      endif
      dam90=discDadn(jpoint)/2.0
      nptc90=nptc90-1
      dtotdam90=dtotdam90+dble(dam90)
c     The closed loop is erased and the point of rev. is already there
Css      eobj90=tlstr90(nptt90)
Css      sobj90=tlsts90(nptt90)
      tldam90(nptt90)   =dam90
      tltotCrk90(nptt90)=dtotdam90
      itlnrev90(nptt90) =nrev
      go to 5000  ! Ramp is done go to 00 code

C     An unmatched 1/2 cycle is being forgotten and a return to the 
c     monotonic curve is occurring. Count the unmatched 1/2 cycle.
 1165 continue
Css      dam90=SMITH ( (tlimL90(1)-climL90(1)), tlsts90(1), clsts90(1),
Css     &         tlstr90(1),clstr90(1),totdam90,nrev)
      dld90=tlimL90(1)-climL90(1)
      if(dld90.eq.0.)then
        write(6,*)"#Error: 1165:dld90=0 tlimL90(1),climL90(1) = ",
     &            tlimL90(1),climL90(1),nrev
      endif
Cvers1            call getCracks(iret,dld90,dam90)
Cvers1            if(iret .ne. 0)goto 9000  !   Error returned. goto stop
      jpoint=dld90/discDKinterval+1
      if(jpoint.gt.ndiscMax)then
       write(0,1136)nrev,nblk,nact,dld90,deltaKmax
       write(6,1136)nrev,nblk,nact,dld90,deltaKmax
       go to 9000  ! stop ------------------- ! Fracture
      endif
      dam90=discDadn(jpoint)/2.0
      dtotdam90=dtotdam90+dble(dam90)
      nptc90=0
      nptt90=0
      go to 1500

C     A loop is being closed and a return to the monotonic 
c     curve is occurring.  Damage is the same as the other half of
c     the cycle being closed.
 1170 continue
Css      dam90=SMITH( (tlimL90(1)-climL90(2)),
Css     &    tlsts90(1),clsts90(2),
Css     &    tlstr90(1),clstr90(2), totdam,nrev)
      dld90=tlimL90(1)-climL90(2)
      if(dld90.eq.0.)then
        write(6,*)"#Error: 1170:dld90=0 tlimL90(1),climL90(2)= ",
     &            tlimL90(1),climL90(2),nrev
      endif
Cvers1            call getCracks(iret,dld90,dam90)
Cvers1            if(iret .ne. 0)goto 9000  !   Error returned. goto stop
      jpoint=dld90/discDKinterval+1
      if(jpoint.gt.ndiscMax)then
       write(0,1136)nrev,nblk,nact,dld90,deltaKmax
       write(6,1136)nrev,nblk,nact,dld90,deltaKmax
       go to 9000  !stop  ------------------ Fracture
      endif
      dam90=discDadn(jpoint)/2.0
      dtotdam90=dtotdam90+dble(dam90)
c     Eliminate the closed loop.
      nptc90=0
      nptt90=0
      go to 1500

C     A new entry in the P.D. list is being made.
 1350 continue
      dld90=abs(lobj90-climL90(nptc90))
Css      de90=FLD(dld,nrev)
Css      ds90=DET(de,nrev)
Css      eobj90=eo90+de90
Css      sobj90=so90+ds90
Css      dam90=SMITH( dld90,sobj90,so90,eobj90,eo90,totdam90,nrev)
      if(dld90.eq.0.)then
        write(6,*)"#Error: 1350:dld90=0 lobj90,climL90(nptc90)= ",
     &            lobj90,climL90(nptc90),nrev
      endif
Cvers1            call getCracks(iret,dld90,dam90)
Cvers1            if(iret .ne. 0)goto 9000  !   Error returned. goto stop
      jpoint=dld90/discDKinterval+1
      if(jpoint.gt.ndiscMax)then
       write(0,1136)nrev,nblk,nact,dld90,deltaKmax
       write(6,1136)nrev,nblk,nact,dld90,deltaKmax
       go to 9000   !stop  ---------------------- Fracture!
      endif
      dam90=discDadn(jpoint)/2.0
      dtotdam90=dtotdam90+dble(dam90)
c     Enter the rev in the P.D. list
      nptt90=nptt90+1
      tlimL90(nptt90)=lobj90
Css      tlstr90(nptt90)=eobj90
Css      tlsts90(nptt90)=sobj90
      tldam90(nptt90)   =dam90
      tltotCrk90(nptt90)=dtotdam90
      itlnrev90(nptt90) =nrev
      go to 5000  ! Ramp is done go to 00 code

C     Deformation is occurring on the monotonic curve again.
 1500 continue
      dld90=abs(lobj90)
Css      de90=FLD(dld90*2.0,nrev)
Css      ds90=DET(de90,nrev)
Css      eobj90=de90/2.0
Css      sobj90=ds90/2.0
c     Subtract damage of previous use of monotonic curve (=0 in cyc )
      dtotdam90=dtotdam90-dble(tldam90(1) )
      go to 4000

C     A loop is being closed, count the remaining half of the
c     loop and then eliminate the loop.
 1400 continue
      dld90=tlimL90(nptt90)-climL90(nptc90)
Css      dam90=SMITH( (tlimL90(nptt90)-climL90(nptc90)), 
Css     &      tlsts90(nptt90),clsts90(nptc90),
Css     &      tlstr90(nptt90),clstr90(nptc90),totdam90,nrev)
      if(dld90.eq.0.)then
       write(6,*)"#Error:1400:dld90=0 tlimL90(nptt90),climL90(nptc90)=",
     &            tlimL90(nptt90),climL90(nptc90),nrev
      endif
Cvers1            call getCracks(iret,dld90,dam90)
Cvers1            if(iret .ne. 0)goto 9000  !   Error returned. goto stop
      jpoint=dld90/discDKinterval+1
      if(jpoint.gt.ndiscMax)then
       write(0,1136)nrev,nblk,nact,dld90,deltaKmax
       write(6,1136)nrev,nblk,nact,dld90,deltaKmax
       go to 9000  !stop --------------------- Fracture !
      endif
      dam90=discDadn(jpoint)/2.0
      dtotdam90=dtotdam90+dble(dam90)
      nptc90=nptc90-1
C     Subtract the old half cycle damage of the Re-incurred
c     stress-strain path.
      dtotdam90=dtotdam90-dble(tldam90(nptt90) )
      nptt90=nptt90-1
c     Reset the local stress-strain origin.
Css      eo90=clstr90(nptc90)
Css      so90=clsts90(nptc90)
      ldo90=climL90(nptc90)
      go to 1100

C*****  Going Down, Compressive Direction ****************************

C     Are we on the monotonic curve?
 2100 if(nptc90 .eq. 0) go to 2500

c     No. Has a exceedence  occurred?
 2110 if(lobj90 .lt. climL90(nptc90)) go to 2120

c     Is this Same as previous amplitude?
 2130 if(lobj90 .eq. climL90(nptc90) .and. nptc90 .ne. 1) go to 2135
      go to 2350

c     Is a loop being closed without the monotonic?
 2120 if(nptc90 .ne. 1) go to 2400

c     Is closure in connection  with the monotonic?
 2150 if(nptt90 .eq. 2) go to 2170

c     Check for impossible
 2160 if(nptt90 .eq. 1) go to 2165
      idump=2160
      go to 9999

C     Same load level as previous level.
 2135 continue
Css      dam90=SMITH( (ldo90-lobj90),tlsts90(nptt90),clsts90(nptc90),
Css     &      tlstr90(nptt90),clstr90(nptc90),totdam90,nrev)
      dld90=ldo90-lobj90
      if(dld90.eq.0.)then
       write(6,*)"#Error:2135:dld90=0 ldo90,lobj90=",
     &            ldo90,lobj90,nrev
      endif
Cvers1            call getCracks(iret,dld90,dam90)
Cvers1            if(iret .ne. 0)goto 9000  !   Error returned. goto stop
      jpoint=dld90/discDKinterval+1
      if(jpoint.gt.ndiscMax)then
       write(0,1136)nrev,nblk,nact,dld90,deltaKmax
       write(6,1136)nrev,nblk,nact,dld90,deltaKmax
       go to 9000  !stop --------------------- Fracture !
      endif
      dam90=discDadn(jpoint)/2.0
      dtotdam90=dtotdam90+dble(dam90)
c     The closed loop is erased and the point of rev. is already there
      nptt90=nptt90-1
Css      eobj90=clstr90(nptc90)
Css      sobj90=clsts90(nptc90)
      cldam90(nptc90)   =dam90
      cltotCrk90(nptc90)=dtotdam90
      iclnrev90(nptc90) =nrev
      go to 5000  ! Ramp is done go to 00 code

C     An unmatched 1/2 cycle is being forgotten and a return to the
c     monotonic curve is occurring. Count the unmatched 1/2 cycle.
 2165 continue
Css      dam90=SMITH( (tlimL90(1)-climL90(1)),tlsts90(1),clsts90(1),
Css     &        tlstr90(1),clstr90(1),totdam90,nrev)
      dld90=tlimL90(1)-climL90(1)
      if(dld90.eq.0.)then
       write(6,*)"#Error:2165:dld90=0 tlimL90(1),climL90(1)=",
     &            tlimL90(1),climL90(1),nrev
      endif
Cvers1            call getCracks(iret,dld90,dam90)
Cvers1            if(iret .ne. 0)goto 9000  !   Error returned. goto stop
      jpoint=dld90/discDKinterval+1
      if(jpoint.gt.ndiscMax)then
       write(0,1136)nrev,nblk,nact,dld90,deltaKmax
       write(6,1136)nrev,nblk,nact,dld90,deltaKmax
       go to 9000  !stop --------------------- Fracture !
      endif
      dam90=discDadn(jpoint)/2.0
      dtotdam90=dtotdam90+dble(dam90)
      nptc90=0
      nptt90=0
      go to 2500

C     A loop is being closed and a return to the monotonic
c     curve is occurring.  Damage is the same as the other half of
c     the cycle being closed.
 2170 continue
Css      dam90=SMITH( (tlimL90(2)-climL90(1)),tlsts90(2),
Css     &       clsts90(1),tlstr90(2),clstr90(1),totdam90,nrev)
      dld90=tlimL90(2)-climL90(1)
      if(dld90.eq.0.)then
       write(6,*)"#Error:2170:dld90=0 tlimL90(2),climL90(1)=",
     &            tlimL90(2),climL90(1),nrev
      endif
Cvers1            call getCracks(iret,dld90,dam90)
Cvers1            if(iret .ne. 0)goto 9000  !   Error returned. goto stop
      jpoint=dld90/discDKinterval+1
      if(jpoint.gt.ndiscMax)then
       write(0,1136)nrev,nblk,nact,dld90,deltaKmax
       write(6,1136)nrev,nblk,nact,dld90,deltaKmax
       go to 9000  !stop --------------------- Fracture !
      endif
      dam90=discDadn(jpoint)/2.0
      dtotdam90=dtotdam90+dble(dam90)
      nptc90=0
      nptt90=0
      go to 2500

C     A new entry in the P.D. list is being made.
 2350 continue
      dld90=abs(lobj90-tlimL90(nptt90))
Css      de90=FLD(dld90,nrev)
Css      ds90=DET(de90,nrev)
Css      eobj90=eo90-de90
Css      sobj90=so90-ds90
Css      dam90=SMITH(dld90,so90,sobj90,eo90,eobj90,totdam90,nrev)
      if(dld90.eq.0.)then
       write(6,*)"#Error:2350:dld90=0 lobj90,tlimL90(nptt90)=",
     &            lobj90,tlimL90(nptt90),nrev
      endif
Cvers1            call getCracks(iret,dld90,dam90)
Cvers1            if(iret .ne. 0)goto 9000  !   Error returned. goto stop
      jpoint=dld90/discDKinterval+1
      if(jpoint.gt.ndiscMax)then
       write(0,1136)nrev,nblk,nact,dld90,deltaKmax
       write(6,1136)nrev,nblk,nact,dld90,deltaKmax
       go to  9000  !stop --------------------- Fracture !
      endif
      dam90=discDadn(jpoint)/2.0
      dtotdam90=dtotdam90+dble(dam90)
      nptc90=nptc90+1
      if(jpoint.gt.ndiscMax)then
       write(0,1136)nrev,nblk,nact,dld90,deltaKmax
       write(6,1136)nrev,nblk,nact,dld90,deltaKmax
       go to  9000  !stop --------------------- Fracture !
      endif
      climL90(nptc90)=lobj90
Css      clstr90(nptc90)=eobj90
Css      clsts90(nptc90)=sobj90
      cldam90(nptc90)   =dam90
      cltotCrk90(nptc90)=dtotdam90
      iclnrev90(nptc90) =nrev
      go to 5000  ! Ramp is done go to 00 code

C     Deformation is occurring on the monotonic curve again. --------
 2500 continue
      dld90=abs(lobj90)
Css      de90=FLD(dld90*2.0,nrev)
Css      ds90=DET(de90,nrev)
Css      eobj90=-de90/2.0
Css      sobj90=-ds90/2.0
c     Subtract damage of previous use of monotonic curve (=0 in cyc )
      dtotdam90=dtotdam90 - dble(tldam90(1) )

C     Add the present damage.  Both Tens. and Comp. comes here.
 4000 continue
Css      dsx90=abs(sobj90)
Css      dex90=abs(eobj90)
C     dsx and dex are amplitudes
C?    dsamp=dsx/2.0
C?    deamp=dex/2.0
C?      The above was probably done to reduce Monotonic damage? or
C?      to solve the question of what is a monot. half cycle?  is it
C?      a range  or an amplitude?
C     Counting damage for monotonic curve usage is probably not a good idea.
C     For example, what if monot. usage is only on compression side?
Css      dsamp90=dsx90
Css      deamp90=dex90
Css      dam90= SMITH(dld90,dsamp90,-dsamp90,deamp90,-deamp90,
Css     &       totdam90,nrev)
      if(dld90.eq.0.)then
        write(6,*)"#Error: 4000:dld90=0 lobj90= ",
     &            lobj90,nrev
      endif
Cvers1            call getCracks(iret,dld90,dam90)
Cvers1            if(iret .ne. 0)goto 9000  !   Error returned. goto stop
      jpoint=dld90/discDKinterval+1
      if(jpoint.gt.ndiscMax)then
       write(0,1136)nrev,nblk,nact,dld90,deltaKmax
       write(6,1136)nrev,nblk,nact,dld90,deltaKmax
       go to 9000  !stop ------------------  Fracture!
      endif
      dam90=discDadn(jpoint)/2.0
      dtotdam90=dtotdam90+dble(dam90)
c     Set the appropriate push-down list arrays for monotonic.
      climL90(1)=-abs(lobj90)
      tlimL90(1)=abs(lobj90)
      nptt90=1
      nptc90=1
Css      tlstr90(1) =abs(eobj90)
Css      clstr90(1) =-abs(eobj90)
Css      tlsts90(1) =abs(sobj90)
Css      clsts90(1) =-abs(sobj90)
      tldam90(1)   =dam90
      tltotCrk90(1)=dtotdam90
      itlnrev90(1) =nrev
      cldam90(1)   =dam90
      cltotCrk90(1)=dtotdam90
      iclnrev90(1) =nrev
      
      go to 5000    ! half cycle is done, go to  00 code

C     A loop is being closed, count the remaining half of the
c     loop and then eliminate the loop.
 2400 continue
Css      dam90=SMITH( (tlimL90(nptt90)-climL90(nptc90)), 
Css     &       tlsts90(nptt90),clsts90(nptc90),
Css     &       tlstr90(nptt90),clstr90(nptc90),totdam90,nrev)
      dld90=tlimL90(nptt90)-climL90(nptc90)
      if(dld90.eq.0.)then
       write(6,*)"#Error:2400:dld90=0 tlimL90(nptt90),climL90(nptc90)=",
     &            tlimL90(nptt90)-climL90(nptc90),nrev
      endif
Cvers1            call getCracks(iret,dld90,dam90)
Cvers1            if(iret .ne. 0)goto 9000  !   Error returned. goto stop
      jpoint=dld90/discDKinterval+1
      if(jpoint.gt.ndiscMax)then
       write(0,1136)nrev,nblk,nact,dld90,deltaKmax
       write(6,1136)nrev,nblk,nact,dld90,deltaKmax
       go to 9000  !stop  ----------------------- Fracture!
      endif
      dam90=discDadn(jpoint)/2.0
      dtotdam90=dtotdam90+dble(dam90)
      nptt90=nptt90-1
C     Subtract the old half cycle damage of the Re-incurred
c     stress-strain path.
      dtotdam90=dtotdam90-dble(cldam90(nptc90) )
      nptc90=nptc90-1
c     Reset the local stress-strain origin.
Css      so90=tlsts90(nptt90)
Css      eo90=tlstr90(nptt90)
      ldo90=tlimL90(nptt90)
      go to 2100     !    continue the ramp.



 5000 continue
      go to 3000    !go back to begin of overall program cycle loop

 9000 continue
C      write out the -ve number of the last rec.  that was written.
C      All the other variables in this rec are endof test values.
       ndummy=-nrecord  !make -ve to make it unique.
       xdummy=0.
       write(60,rec=1 )ndummy,xdummy,nblk,nact,
     &            lobj90,xMm90,xMb90,
     &            stsMembrane,stsBending
      close(unit=60)

       write(0,9050)nrev,totdam90,nblk,nact,nrecord
       write(6,9050)nrev,totdam90,nblk,nact,nrecord
 9050  format("#Last: nrev= ",i10," a= ",e14.7,
     &        "   nblk= ",i10," nact= ",i10," nrecord= ",i10/)
      stop

 9999 write(6,9998) idump
 9998 format("# *** Error: near statement label no. ",i5/)
      stop
      end

C===============================================================
      SUBROUTINE getCracks(iret,dld,dam)
C     s/r to take dld, the deltaK range, and compute dam= crack increment
C     from the tables of  dadn vs deltak.
C     The table used is a log  version of dadn vs deltaK,  and desired 
C     points are interpolated using the log-log  points of the table.

      real*4  deltaK(50),dadn(50),logdeltaK(50),logdadn(50),
     &        logdiffdK(50),logdiffdadn(50)
      logical debugDadn
      common/DADN/ deltaK,dadn,logdeltaK,logdadn,logdiffdK,logdiffdadn,
     &             ndadn,maxdadn,debugDadn

C     Note that this routine is not called each reversal. The simulation
C     uses the damage tables to compute damage, normally.

      iret=0  ! return and error no. if non-zero
      if(dld.le.0.)then
        write(0,100)dld   !  This should not happen, try to compensate.
        write(6,100)dld
  100   format("#Error: call to getCracks() deltaK .le. 0,  deltaK=",
     &         e14.7)
        dam=0.
        iret=-1
        return   !a return will allow output file to be ended properly.
      endif
      dldlog=ALOG10(dld)

      n=1   ! point to the first table point.
      if(dldlog .lt. logdeltaK(n) )then  ! smaller than the first point?
C       Assume we are below the threshold
        dam=0
        return
      endif
C     Ok, we are equal to or above the first point
  300 continue
      n=n+1
      if(n.gt.ndadn)then
        write(0,350)dld,deltaK(n-1)
        write(6,350)dld,deltaK(n-1)
  350   format("#Fracture: getCracks():  deltaK = ",e14.7/
     &        "# is bigger than largest entry in dadn table: ",e14.7)
        iret=999
        return  !return to a stop. Must write last binary record.
      endif

      if(dldlog .gt. logdeltaK(n))go to 300  ! we are in the desired interval

C     No?  Ok, we are below the n data point, time to interpolate
      frac=(dldlog-logdeltaK(n-1))/logdiffdk(n)
      dam=10**(frac*logdiffdadn(n) + logdadn(n-1) )
C     Divide by 2 to make it per 1/2 cycle
      dam=dam/2.0
      return
      end
      

C================================================================
      REAL FUNCTION XKP(L,nrev)
C***** Compute Stress concentration due to :
C      (A) Notch effect
C      (B) Finite width correction

      SAVE
      real  L,fwks(100),fwkl(100)
      logical stabl
      common /STAB/ stabl
      stabl=.false.
      if(nrev .eq. 0)go to 100
      if(L .ge. xLong) go to 20
      if(L .ge. xLmin) go to 10

c     Short crack, interpolate with 'fwks'
      ip=ifix(L/xinc)+1
      xlast=(ip-1)*xinc
      ipnext=ip+1
      xkp=fwks(ip)+((L-xlast)/xinc) * (fwks(ipnext)-fwks(ip))
      if(xkp .lt. xkinit)stabl=.true.
      return

C     long crack, use "fwkl"
   10 ip=ifix((L-xLmin)/x2inc)+1
      xlast=xLmin+(ip-1)*x2inc
      ipnext=ip+1
      xkp= fwkl(ip)+((L-xlast)/x2inc)*(fwkl(ipnext)-fwkl(ip))
      if(xkp .lt. xkinit) stabl=.true.
      return

   20 write(6,117)L,nrev
      write(0,117)L,nrev
  117 format("# Crack too Long, L= ",E14.7," nrev=",I10)
      xkp=999.0
      return

c     Initilize prior to program run
  100 write(6,102)
  102 format("#Enter six constants A0, A1, ...A5 for K' polyn")
      read(5,*) a0,a1,a2,a3,a4,a5
  103 write(6,104) a0,a1,a2,a3,a4,a5
  104 format("#" 3(1x,e14.7)/3(1x,e14.7)/"#   OK? (1=yes,0=no):")
      read(5,*,err=103)iyes
      if(iyes .ne. 1) go to 100

  105 write(6,106)
  106 format(" Enter B(1/2 plate width), and D(1/2 notch width)")
      read(5,*,err=105)B,D
      write(6,*)B,D
      
  110 write(6,111)
  111 format("#Enter Max. L/D value of K' polyn. :")
      read(5,*,err=110) xLdm
      write(6,*)xLdm

C     Commence  fw*K' computation  '
      xLo1=xLdm*D
      xinc=xLo1/100.0
      xL1=0.0
      write(6,116) B,D,xLo1,xLdm
  116 format("# B= ",f6.4," D= ",f6.4,/"# Short crack up to L= ",f7.5,
     &  " at L/D max= ",f7.4
     &  //"#FWK #     L     L/D      K'      FW      FW*K'short")
      do 200 i=1,100
        xLd=xL1/D
        if(xLd .eq. 0)then
          sks=a0
        else 
          sks=a0-a1*xLd+a2*xLd**2 -a3*xLd**3 +a4*xLd**4 -A5*xLd**5
        endif
        fw=sqrt(1.0/(cos(3.141593*(xL1+D)/(2.0*B) )))
        fwks(i)=sks*fw
        write(6,114) xL1, xLd,sks, fw, fwks(i)
  114   format("#FWK ",f8.5,1x,f6.3,1x,f9.3,1x,f8.5,1x,f9.3)
C       Increment cracks
        xL1=xL1+xinc
  200 continue

      x2inc=(B-xLo1-D+xinc)/100.0
c     Long crack calcs. overlap short by length of xLo1
      xLmin=xLo1-xinc
      xL2=xLmin
      write(6,217)
  217 format("#FWK #  L      L/D    L/B     K'    FW    FW*K'long")
      do 210 i=1,100
C       Compute larger intervals (long cracks)
        xLd2=xL2/D
        xLb=xL2/B
        sks2=sqrt( (xL2+D)/xL2)
        fw2=sqrt(1.0/(cos(3.141593*(xL2+D)/(2.0*B) )))
        fwkl(i)=sks2*fw2
        write(6,219)  xL2,xLd2,xLb,sks2,fw2,fwkl(i)
  219   format("#FWK ",F8.5,5(1x,f6.3))
C       Increment cracks
        xL2=xL2+x2inc
  210 continue

      xlong=xL2-x2inc*2.0
      write(6,115)
  115 format("# Center crack & Notch assumed !"/"    Enter xkinit :")
      read(5,*)xkinit
      xkp=0.0
      return
      end


C================================================================
      REAL FUNCTION SMITH( skfw,smax,smin,emax,emin,crackL,nrev)
C      Changed mean stress correction aug 21/78
C      The function SMITH= CRACK Increment per Reversal
C     skfw = delta S * K'' * FW  
C     smax,smin = max and min local stress of hysteresis loop
c     emax,emin = max and min local strain
c     crackL = present crack length
c     nrev = no. of reversals. =0 causes initilization of function

      SAVE
      integer iovers
      logical stabl
      common/STAB/ stabl
      smith=0.
      if(nrev .eq.0) go to 200
      if(smax .le. 0) return

      de=emax-emin
      ds=smax-smin
      samp=ds/2.0

c     Mean stress code
      if( .not. stabl) go to 10
c     Consider the strain range of the portion of the loop that is 
c     above zero load:
      if(smin .lt. 0.0) de=de+smin/Emod
      xsmin=smin
      if(smin .lt. 0.0) xsmin=0.0
      dk=emod *de*2 * sqrt(smax/(smax-xsmin))
      go to 11

C     Non-stable, compressive plasticity accounted for
   10 dk= emod * de * sqrt(smax/samp)
   11 dk = dk * sqrt(3.141593 * (crackL+crack0))

C     Check for plastic zone size correction
      if(ds .le. Z) dk = dk * sqrt(xk1/(xk1-3.141593 * ds**2))
      if(iovers .eq. 1) go to 50
    5 if( dk .lt. dkth)return
      dko2=dk/2.0
      if(dko2 .ge. xkc)go to 20
      SMITH= xkc* A * (dk-dkth)**xm/(xkc-dk/2.0)
c     Divide by 2 to make it per half cycle or reversal
      SMITH=SMITH/2.0
      return

   20 write(6,103)dko2,xkc,crackl,nrev
  103 format("# Fracture, DK/2= ",f8.3," .GE. Ko= ",f8.3," L, NREV= ",
     &       F8.4,i10)
      SMITH=9999.9
      return

C     Kth Eliminated
   50 continue
      dko2=dk/2.0
      if(dko2 .GE. xkc) go to 20
      SMITH= (A*(dk**xm) )
      SMITH=SMITH/2.0
      return

C     Initilize constants
  200 write(6,101)
  101 format("# Enter E, Kth, Kc, L0,  A, M, K1, Z, iovers  :")
      read(5,*,err=200) emod,dkth,xkc,crack0, A,xm,xk1, z,iovers
      write(6,102) emod,dkth,xkc,crack0, A,xm,xk1, z,iovers
  102 format(1x,f7.0,f6.2, 1x,f5.1,f7.4,E14.7,F5.2,f8.0,f9.0,i3/
     &     "# OK ? (enter 1 or 0):")
      read(5,*,err=200)iyes
      if(iyes.ne.1) go to 200
      return
      end


C================================================================
      REAL FUNCTION FLD(dld,nrev)
C      dld is = FW * K'' * deltaS
c      This function computes the local strain range from the 
c      nominal stress range (or load range)
C         Material: G40.21-50A crack model data. Load=FW*K''*ds June 5/79

      SAVE
      real str(17)/ 0.,0.,0.,  0.0020,0.00285,0.0040,  0.00523,
     &  0.00650, 0.00796, 0.00945, 0.01118, 0.01318, 0.01540,
     &  0.01765, 0.01990, 0.02220, 0.02450/, sinc/20.0/,emod/30000./
      integer non/17/
      if(nrev .eq. 0) go to 50
      if(dld .ge. 320.0) go to 30
      if(dld .gt. 60.0)  go to 20
C     If not then elastic conditions exist
      FLD=dld/emod
      return

C     The nominal Stress-local Strain data exists, comput the local
C     strain by interpolation.
   20 continue
      ipoint=ifix(dld/sinc)+1
      dlast=(ipoint-1)*sinc
   25 ipnext=ipoint+1
      FLD= str(ipoint)+( (dld-dlast)/sinc)*(str(ipnext)-str(ipoint) )
      return

   30 write(6,106)dld,nrev
C     The stress is too large for interpolation.
  106 format("# ***Nominal Stress of ",f7.2," in REV. ",i10," is too"/
     &    T5,"large for interp. using existing data. Equation used."/
     & )
      FLD=3.1475266E-06 * (dld)**1.536107
      return

   50 so=0.0
      write(6,104)
  104 format(//T3,"Nominal Stress - Local Strain Curve"/)
      do 52 i=1,non
         write(6,103) i,so,str(i)
  103    format(T2,i4,f8.2,2x,f8.5)
         so=so+sinc
   52 continue
      FLD=0.0
      return
      end
      

C================================================================
      REAL FUNCTION DET(de,nrev)
C       Compute the local stress range from the local strain range

      SAVE
      integer non/26/
      real sts(26)/0.0, 30.0, 60.0, 76.5, 86.0, 92.2, 97.4, 101.5,
     &    105.7, 109.4, 112.8, 115.5, 117.6, 120.4, 122.8, 124.9,
     &    127.0, 129.0, 131.0, 132.3, 133.6, 134.8, 136.1, 137.8, 
     &    139.0, 140.0/
     &    ,einc/0.0010/, emod/30000.0/

c     Initilize?
      if(nrev .eq. 0) go to 50
c     The stress data exists, compute the local stress by 
c     interpolating between the points.
      if(de .gt. 0.0020) go to 20
c     If not then elastic conditions exist.
      det=emod*de
      return

   20 continue
      if(de .gt. 0.0250) go to 30
c     No?  then interpolate:
      ipoint=ifix(de/einc)+1
      dlast=(ipoint-1)*einc
   25 ipnext=ipoint+1
      det=sts(ipoint) +
     &    ((de-dlast)/einc) * (sts(ipnext)-sts(ipoint))
      return

   30 continue
      write(6,110) de,nrev
  110 format("#*** Strain Rg. of ",f8.5," in REV. ",i10," is too large"/
     &      T5,"for interpolation, Used  DS=A*DE**B equation.")
      det=361.6735*de**0.2572866
      return

   50 continue
      write(6,104)
  104 format(///T5,"# Smooth Specimen Stress-Strain curve"//
     & T8,"STRAIN", T17,"STRESS"/
     & )
      ex=0.
      do 60 i=1,non
         write(6,55)i,ex,sts(i)
   55    format("#SS    "i3,1x,f8.5,1x,f8.1)
         ex=ex+einc
   60 continue

      DET=0.0
      return
      end


C================================================================================
      SUBROUTINE getStress2Strain(Stress,Strain,iexit)
      SAVE
C     Given Stress Ampl., interpolate Strain
C---------------------------------------------------------------------------
C Subroutine from GNU GPL software: fde.uwaterloo.ca/Fde/Calcs/saefcalc1.html
C  Alterations have been made to data storage.

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

C  From saefcalc1.f:
C      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
C     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
C     &     SigmaxStrainAmp(250)
C      integer ndata
C      real FractureStress,FractureStrain
C      Common/Material/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
C     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
C     &                 ndata,FractureStress,FractureStrain

C     These are used to store the equal spaced Snom. fitted data file.
      real StressAmp(1500), StrainAmp(1500), SnominalAmp(1500)
      integer*4 nmatdata,maxmatdata
C      real FractureStress,FractureStrain,Emod
      logical debugMat
      common/MATERIAL/ StressAmp,StrainAmp,SnominalAmp,
     &      snominalInterval,nmatdata,maxmatdata,debugMat

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
      if(Stress .ge. StressAmp(nmatdata))then
C       The specimen fails in first 1/2 cycle. Damage is >= 1.
        Strain=StrainAmp(nmatdata)
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
C     Change log-log interpolation to linear.  We have 1000 points.
        xslope=
     &     (StressAmp(idat) -StressAmp(idat-1))
     &    /(StrainAmp(idat)-StrainAmp(idat-1) )
        C=  (Stress-StressAmp(idat-1) )/xslope
     &          +StrainAmp(idat-1)
C         Strain= 10.0**Clog
        Strain=C
        return

  200   continue
C       The stress is exactly equal to the curve point(idat)
C       There could be a "flat" spot in the stress-strain curve, and
C       thus the next point would have the same stress value. Check
Cbugfix   May2005 : im was not incrementing
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

        write(0,*)"ERROR:rodSurfFlaw:getStress2Strain.endofdoloop"
        iexit=1
C       We should not be here.  Assume fracture strain or largest
        Strain=StrainAmp(nmatdata)
        return
        end


C-================================================================================
      SUBROUTINE getLoad2StressStrain(Sneuber,Stress,Strain,
     &           iexit)
      SAVE
C     Given Load_Amp (could be FEA Stress Ampl.), interpolate Strain &Stress
C---------------------------------------------------------------------------
C Subroutine from GNU GPL software: fde.uwaterloo.ca/Fde/Calcs/saefcalc1.html

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

C  From saefcalc1.f:
C      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
C     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
C     &     SigmaxStrainAmp(250)
C      integer ndata
C      real FractureStress,FractureStrain
C      Common/Material/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
C     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
C     &                 ndata,FractureStress,FractureStrain

C     Storage for the equal spaced Snom.interval fitted data table in "matfile".
      real StressAmp(1500), StrainAmp(1500), SnominalAmp(1500)
      integer*4 nmatdata,maxmatdata
C      real FractureStress,FractureStrain,Emod
      logical debugMat
      common/MATERIAL/ StressAmp,StrainAmp,SnominalAmp,
     &      snominalInterval,nmatdata,maxmatdata,debugMat

      iexit=0

C       We are assuming that the data has been sorted to largest life,
C       and thus smallest stress, strain etc, first.
C       If the "load", actually = Elastic_local_Stress, is below 
C       the first SelasticAmp point, then it is below
C       the fatigue limit (the first life value in the table).

      if(Sneuber .lt. SnominalAmp(1))then  ! SnominalAmp(1) is usually=0.
C       Yes, its below. Assume straight line from (0,0) t
        Strain=(Sneuber/SnominalAmp(1))*StrainAmp(1)
        Stress=(Sneuber/SnominalAmp(1))*StressAmp(1)
        return
      endif
      
C     Check if data is above Sigfp  or largest point.
C     Assumes: ElasMod= First Stress/ 1st strain
C     If these first points are not on the "elastic" line we could have
C     a problem here, but it shouldnt be too bad (?)

      if(Sneuber .ge. SnominalAmp(nmatdata))then
C       The specimen fails in first 1/2 cycle. Damage is >= 1.
        Strain=StrainAmp(nmatdata)
        Stress=StressAmp(nmatdata)
        write(0,50)Sneuber,SnominalAmp(nmatdata)
        write(6,50)Sneuber,SnominalAmp(nmatdata)
   50   format("#Fracture:getLoad2StressStrain(): Requested Snominal=",
     &  f7.1," Max. value in matfile= ",f7.1," = Fracture")
        iexit=999                   !  this is the only return signal
        return
      endif

C     Ok, the input appears to lie between the two extremes. Lets
C     hunt & peck.  Basically we are going to hunt our way up the stress
C     ladder.  When our point is bigger than the last rung, and smaller than
C     the next rung, we can interpolate.

C     Since the SnominalAmp()  is in equal sized intervals we can jump
C     directly to the two straddle points.
       idat= ifix(Sneuber /snominalInterval)+1

C       Our stress is less than the point on the curve, we have arrived
C       Interpolate.  Try log-log interpl.  If we ever get -ve stresses
C       or strains we will definitely have to go to linear interpl.
C    Since we have about 1000 points, switch to linear interpolation:
        xslope=
     &     ( SnominalAmp(idat) -SnominalAmp(idat-1) )
     &    /( StressAmp(idat)-StressAmp(idat-1) )
        C=  ( Sneuber-SnominalAmp(idat-1) )/xslope
     &          +StressAmp(idat-1)
C        Stress= 10.0**Clog
        Stress=C
        xslope=
     &     (SnominalAmp(idat) -SnominalAmp(idat-1))
     &    /(StrainAmp(idat)-StrainAmp(idat-1) )
        C=  (Sneuber-SnominalAmp(idat-1) )/xslope
     &          +StrainAmp(idat-1)
C        Strain= 10.0**Clog
         Strain=C
         if(debugMat)then
           write(0,95)SnominalAmp(idat-1),Sneuber,SnominalAmp(idat)
           write(6,95)SnominalAmp(idat-1),Sneuber,SnominalAmp(idat)
   95      format("#matfile #Snominal: ",3(f8.2,1x))

           write(0,96)StressAmp(idat-1),Stress,StressAmp(idat)
           write(6,96)StressAmp(idat-1),Stress,StressAmp(idat)
   96      format("#matfile #Stress  : ",3(f8.2,1x))

           write(0,97)StrainAmp(idat-1),Strain,StrainAmp(idat)
           write(6,97)StrainAmp(idat-1),Strain,StrainAmp(idat)
   97      format("#matfile #Strain  : ",3(f8.5,1x))
         endif

        return
        end



C============================================================================
      SUBROUTINE getloop(Stress,Straintemp,Stresstemp,npts,iexit)
      SAVE
C     Given Stress Ampl., get all stress-strain points from 0 to Stress
C---------------------------------------------------------------------------
C Subroutine from GNU GPL software: fde.uwaterloo.ca/Fde/Calcs/saefcalc1.html

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
      real Stresstemp(250),Straintemp(250)

C  From saefcalc1.f:
C      real StressAmp(250), StrainAmp(250), SelasticAmp(250),
C     &     Lifecycles(250), PlstrainAmp(250), ElstrainAmp(250),
C     &     SigmaxStrainAmp(250)
C      integer ndata
C      real FractureStress,FractureStrain
C      Common/Material/ StressAmp,StrainAmp,SelasticAmp,Lifecycles,
C     &                 PlstrainAmp,ElstrainAmp,SigmaxStrainAmp,
C     &                 ndata,FractureStress,FractureStrain

C     Storage for the equal spaced Snom.interval fitted data table in "matfile".
      real StressAmp(1500), StrainAmp(1500), SnominalAmp(1500)
      integer*4 nmatdata,maxmatdata
C      real FractureStress,FractureStrain,Emod
      logical debugMat
      common/MATERIAL/ StressAmp,StrainAmp,SnominalAmp,
     &      snominalInterval,nmatdata,maxmatdata,debugMat

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
      
C     Check if data is above max available
      if(Stress .ge. StressAmp(nmatdata))then
C       The specimen probably fails in first 1/2 cycle. Damage is >= 1.
        do 50 i=1,ndata
          npts=npts+1
          Straintemp(npts)=StrainAmp(i)
          Stresstemp(npts)=StressAmp(i)
   50   continue
        npts=npts+1
        Straintemp(npts)=StrainAmp(nmatdata)
        Stresstemp(npts)=StressAmp(nmatdata)
        iexit=1
        go to 9000
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
Cbugfix   May2005 : im was not incrementing. fixed:
        im=idat
  201   im=im+1
        if(Stress .eq. StressAmp(im))then
C         yes, another is equal. is there more?
          go to 201
        endif
C       No, now the (im) point has bigger stress, use the last
C       one that was equal as a match
Cbugfix   May2005 : Should be im-1
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

        write(0,*)" ERROR:rodSurfFlaw:getloop.endofdoloop"
        iexit=1
C       We should not be here.  Assume fracture strain
        npts=npts+1
        Straintemp(npts)=StrainAmp(nmatdata)
        Stresstemp(npts)=StressAmp(nmatdata)
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





C======================================================================
      SUBROUTINE getPeakLoads(ivec,iret)
      SAVE
C     ivec= function instruction, usually=1
C     iret= return value 
C     s/r to look at the total stress and remove any non-reversal points
C     from the load history.   A list of removals is formed and then the
C     designated loads are removed from the memory.  In this program we are
C     focused on the total stress  ststot() to decide removal or keep.

c     Load history storage.  If changing dimensions, also change same 
C                             in mainline
C     Stress membrane, bending, and total storage:
      real*4 stsm(5000),stsb(5000),ststot(5000)
      real*4 ststotmax,ststotmin,ststotwindow
      integer*4 nloads
      logical debugLoads
      common /LOADS/ stsm,stsb,ststot, ststotmax,ststotmin,ststotwindow,
     &               nloads,debugLoads

      integer*4 iremove(5000)   ! list of points to be removed in 2nd pass of data

      ngood=0                 ! records the no. of good points
      n=0                     ! points to the data being examined
      nbad=0                  ! counts the no. of points being eliminated

C     Get the first point.
      n=n+1
      stsold=ststot(1)
      nstsold=1
      ngood=ngood+1

C     check if the first and last point in the history are same
      if(ststot(1).eq.ststot(nloads) )then
        write(0,90)
        write(6,90)
   90   format("#history #First and last pts in history are same."/
     &         "#history #Eliminating the 1st history point now")
        nbad=nbad+1
        iremove(nbad)=1
C       Save a copy for rainflow file created near end of program
        ststot1=ststot(1)
      endif


  100 continue
      n=n+1
      stsnext=ststot(n)    !get the 2nd point
      if(stsnext .eq. stsold)then ! remove the new point
        nbad=nbad+1
        iremove(nbad)=n
        go to 100
      endif
C     1st and 2nd pts have been read in
C     2nd pt is not equal to first
      iupdown=+1
      if(stsnext .lt. stsold)iupdown=-1
      stsold=stsnext
      nstsold=n

C     ok, first point and direction is set. Keep going
 2000 continue
      n=n+1
      if(n .gt. nloads)go to 9000
      stsnext=ststot(n)
C     If its equal to the old one we can eliminate the new
      if(stsnext .eq. stsold)then
        nbad=nbad+1
        iremove(nbad)=n
        go to 2000
      endif

      if(stsnext .lt. stsold)go to 3000  ! Potentially it is a reversal point, downwards
C     If not, then things are going in the same direction  if iupdown =1
C     New is going upwards--------------------------------------------------
C     If the same direction we can discard stsold  
      if(iupdown .eq. +1)then
C        Same direction but further, discard old pt
         nbad=nbad+1
         iremove(nbad)=nstsold
         stsold=stsnext
         nstsold=n
         go to 2000 ! get next point
       endif

C      iupdown must have been = -1,  thus previous direction was in compression.
C      If we get to this point then it is a potential reversal.
C      See if its size is bigger than the allowable window
       if((stsnext-stsold) .lt. ststotwindow)then  ! too small, skip new pt
         nbad=nbad+1
         iremove(nbad)=n
         go to 2000
       endif

C      Ok, it is bigger than the window. Must be a reversal
       ngood=ngood+1
       iupdown=+1
       stsold=stsnext
       nstsold=n
       go to 2000

C     New is going downwards -------------------------------------------------------
 3000 continue
C     the new point is going downwards.  Check the old direction
      if(iupdown .eq. -1)then  ! same direction
C       Discard the old point, and replace it with the new pt.
         nbad=nbad+1
         iremove(nbad)=nstsold
         stsold=stsnext
         nstsold=n
         go to 2000
      endif

C     iupdown was +1, previous direction was into tension
C     Potential change in direction.  Is it bigger than the window
      if((stsold-stsnext) .lt. ststotwindow)then
C       too small, eliminate the new point
         nbad=nbad+1
         iremove(nbad)=n
         go to 2000
      endif

C     Ok, reversal is bigger than the window. Its a reversal
      ngood=ngood+1
      iupdown=-1
      stsold=stsnext
      nstsold=n
      go to 2000

 9000 continue
      if(nbad.eq.0)then    ! no points eliminated
        write(6,9002)
 9002   format("#history #No points needed to be eliminated."/)
        ngood=nloads
        go to 9799  ! print out the filtered history and rainflow
      endif

      write(6,*)"#history #PPICK found ",ngood,"good pts, and ",nbad,
     &          " to be eliminated. List: "
      do 9010 i=1,nbad
         write(6,9009)iremove(i)
 9009    format("#history #eliminate pt.",i7)
 9010 continue

C     ------------------selction is done.  Must now eliminate and update list

C     We first need to sort the list of elimination peaks. Use ugly sort.
C      write(6,9480)
C 9480 format("#history #Sorted Eliminated point nos.:")
C      do 9500 i=1,nbad-1
CC       Find the smallest and put into position iremove(i)
C        ismallest=iremove(i)
C        ismallindex=i
C        do 9300 j=i+1,nbad
C           if(iremove(j).lt.ismallest)then
CC            Found a smaller one
C             ismallest=iremove(j)
C             jsmallest=j
C           endif
C 9300   continue
CC     Now interchange the two
C      itemp=iremove(i)
C      iremove(i)=iremove(j)
C      iremove(j)=itemp
C      write(6,9009)iremove(i)
C 9500 continue

C     The removal list is now sorted to the smallest peak no. first
C     Now go throught the stress storage and take out those to be eliminated
C     Start at the first point to be eliminated
      ielim=1       !Points to the active point in iremove()
      iput=1        !Points to next good save location, initially=1
      do 9600 i=1,nloads
        jremove=iremove(ielim)
        if(i .eq. jremove)then
C         Yes this point will be eliminated
C         increment the pointer to next elimination point
          ielim=ielim+1
          jremove=iremove(ielim)
          go to 9600
        endif
C       No?  then save the point
C       This point does not need elimination, save it
        stsm(iput)=stsm(i)
        stsb(iput)=stsb(i)
        ststot(iput)=ststot(i)
        iput=iput+1
        write(6,9590)stsm(i),stsb(i),ststot(i),i
 9590   format("#history #Filtered ",3(f7.1,1x),i6)
 9600 continue

C     Check results with a full print
 9799 continue
      do 9800 i=1,ngood
        write(6,9790)stsm(i),stsb(i),ststot(i),i
 9790   format("#history #Filteredck ",3(f7.1,1x),i6)
 9800 continue

C     Also write out a file for rainflow count + StrainStrainLife
      open(unit=10,file="loads4rain.out")
C     Compute the boundaries of the rainflow matrix
      extra=(ststotmax-ststotmin)/64.
      xmax=ststotmax+extra
      xmin=ststotmin-extra
      write(10,9838)xmax,xmin
 9838 format("#MAX= ",e14.7/"#MIN= ",e14.7/"#BEGIN")
      write(10,9839)ststot1
 9839 format(" 0 ",f7.1)
      do 9850 i=1,ngood
        write(10,9840)i,ststot(i)
 9840   format(i9,1x,f7.1,1x)
 9850 continue
      close(unit=10)
      write(0,9854)
      write(6,9854)
 9854 format("#Wrote ststot loads for rainflow: loads4rain.out")

      nloads=ngood
      return
      end

C==============================================================

      SUBROUTINE getRBarMmMb(ivec,iret,ao2r,xMm,xMb)
      SAVE
C     ivec=0  initilize storages, ivec=1  compute xMm and xMb
C     Round Bar (rod)  Surface crack
C     We have 1000 points to define the curve, so don't bother interpolating,
C     just take nearest value as a solution to   ao2r  input.
C     Use equations  M.37 from  BS7910 2005 pg 229

      real*4 storeRBarMm(1000),storeRBarMb(1000)
      real*4 deltaRBarao2r,ao2rmin,ao2rmax
      integer maxstore
      common /LONGMMDATA/storeRBarMm,storeRBarMb,deltaRBarao2r,ao2rmin,
     &        ao2rmax,maxstore

      iret=0
      if(ivec.eq.1)goto 500

      if(ivec.lt.0 .or. ivec .gt.2)then
        write(0,110)ivec
        write(6,110)ivec
  110   format("#getRBarMmMb Error: ivec=",i6," unknown."/"#Stopping")
        stop
      endif

C     ivec == 0   initilize storage
      pi=3.1415927
      term1= 1.84/pi
      ao2rmax= 0.6
      ao2rmin= 0.0  !  0.0 is not really allowed
      maxstore=1000
C     if there are maxstore  storage points, then delta must be 
      deltaRBarao2r= (ao2rmax-ao2rmin)/(float(maxstore-1) )
      n=1
      ao2r=ao2rmin
      write(6,290)
  290 format("#RbarMmMb= #  a/B      Mm       Mb        n")

  300 continue
C     According to BS7910  eq. M.37 page 229
      piao4r= pi*(ao2r/2.0)
      sinterm=sin(piao4r)
      costerm=cos(piao4r)
C      tanterm=tan(piao4r)
      tanterm=sinterm/costerm
      g= (term1 *sqrt(tanterm/piao4r) )/costerm
      
      storeRBarMm(n)=g * (0.752+ 2.02*ao2r+ 0.37*( 1-sinterm)**3 )
      storeRBarMb(n)=g * (0.923+0.199*(1-sinterm)**4 )
      write(6,320)ao2r,storeRBarMm(n),storeRBarMb(n),n
  320 format("#RbarMmMb= ",f8.6,1x,f8.6,1x,f8.6,1x,i5)
      n=n+1
      if(n.gt.maxstore)return
C     keep going
      ao2r=ao2r+deltaRBarao2r
      goto 300

  500 continue   !------------------------Run time table lookup--------
C     Get a solution
      if(ao2r.ge.ao2rmax)then
        write(0,520)ao2r,ao2rmax
        write(6,520)ao2r,ao2rmax
  520   format("#RbarMmMb= #END:  a/2r= ",f8.6," .ge.  a/2r max= ",f8.6)
        iret=999
        return   !stop
      endif
      if(ao2r.le.ao2rmin)then    !zero is not good either.
        write(0,530)ao2r,ao2rmin
        write(6,530)ao2r,ao2rmin
  530   format("#RbarMmMb= #ERROR: a/2r= ",f8.6," .le. a/2r min= ",f8.6)
        itret=999
        return  !stop
      endif

      n= ifix(ao2r/deltaRBarao2r)+1
C     When there are a lot of data points, interpolation doesnt add much,
C     so don't bother. Just take nearest.
      xMm= storeRBarMm(n)
      xMb= storeRBarMb(n)
      if(ivec .eq.2)then   !debug output
        write(0,550)ao2r,xMm,xMb
        write(6,550)ao2r,xMm,xMb
  550   format("#RbarMmMb= ",3(1x,f8.6),"#debug a/2r xMm xMb")
      endif
      return
      end

C===========================================================================
C      SUBROUTINE getRBarMkmMkb(ivec,iret,zob,xMkm,xMkb)
C      BS7910 does not presently account for welds in an Rbar 
