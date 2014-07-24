! -*- F90 -*-


!   LHAGLUE Interface to LHAPDF library of modern parton
!    density functions (PDF) with uncertainties
! 
! Authors for v4: Dimitri Bourilkov, Craig Group, Mike Whalley
! 
! Authors for v3: Dimitri Bourilkov, Craig Group, Mike Whalley
! 
! Author for v1 and v2: Dimitri Bourilkov  bourilkov@mailaps.org
!                       University of Florida
! 
! HERWIG interface by Dimitri Bourilkov and Craig Group
! 
! New numbering scheme and upgrade for LHAPDF v2.1
! by Dimitri Bourilkov and Mike Whalley
! 
! For more information, or when you cite this interface, currently
! the official reference is:
! D.Bourilkov, "Study of Parton Density Function Uncertainties with
! LHAPDF and PYTHIA at LHC", hep-ph/0305126.
! 
! The official LHAPDF page is:
! 
!    http://durpdg.dur.ac.uk/lhapdf/index.html 
! 
! The interface contains four subroutines (similar to PDFLIB).
! It can be used seamlessly by Monte Carlo generators 
! interfaced to PDFLIB or in stand-alone mode.
! 
!     For initialization (called once)
! 
!     PDFSET(PARM,VALUE)
! 
!     For the proton/pion structure functions
! 
!     STRUCTM(X,Q,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
! 
!     For the photon structure functions
! 
!     STRUCTP(X,Q2,P2,IP2,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
! 
!     For statistics ON structure functions (under/over-flows)
! 
!     PDFSTA
! 
! This interface can be invoked in 3 ways depending
! on the value of PARM(1) provided by the user when
! calling PDFSET(PARM,VALUE):
! 
!     For PYTHIA:         PARM(1).EQ.'NPTYPE'
!       (this is set automatically by PYTHIA)
! 
!     For HERWIG:         PARM(1).EQ.'HWLHAPDF'
!       (set by the USER e.g. in the main program like this:
!           AUTPDF(1) = 'HWLHAPDF'
!           AUTPDF(2) = 'HWLHAPDF'                         )
! 
!     For Stand-alone:    PARM(1).EQ.'DEFAULT'
!       (can be used for PDF studies or when interfacing
!        new generators)
! 
! The LHAPDF set/member is selected depending on the value of:
! 
!         PYTHIA:   ABS(MSTP(51)) - proton
!                   ABS(MSTP(53)) - pion
!                   ABS(MSTP(55)) - photon
! 
!         HERWIG:   ABS(INT(VALUE(1)))
! 
!    STAND-ALONE:   ABS(INT(VALUE(1)))
! 
! 
!         CONTROL switches:
!        ==================
! 
!      THE LOCATION OF THE LHAPDF LIBRARY HAS TO BE SPECIFIED
!      AS DESCRIBED BELOW (the rest is optional)
! 
!      if the user does nothing, sensible defaults
!      are active; to change the behaviour, the corresponding
!      values of LHAPARM() should be set to the values given below
! 
!   Location of the LHAPDF library of PDFs (pathname):
!      uses common block /LHAPDFC/
! 
!    If the user does nothing => default = subdir PDFsets of the 
!                               current directory (can be real subdir
!                               OR a soft link to the real location)
!    If the user sets LHAPATH => supplied by the USER who defines the
!                         path in common block COMMON/LHAPDFC/LHAPATH
!                         BEFORE calling PDFSET
! 
!    Other controls:
!    ===============
!      use common block /LHACONTROL/
! 
!   Collect statistics on under/over-flow requests for PDFs
!   outside their validity ranges in X and Q**2
!   (call PDFSTA at end of run to print it out)
! 
!       LHAPARM(16).EQ.'NOSTAT' => No statistics (faster)
!       LHAPARM(16).NE.'NOSTAT' => Default: collect statistics
! 
!   Option to use the values for the strong coupling alpha_s
!   as computed in LHAPDF in the MC generator
!   (to ensure uniformity between the MC generator and the PDF set)
!   WARNING: implemented ONLY for PYTHIA in LHAPDFv4
! 
!       LHAPARM(17).EQ.'LHAPDF' => Use alpha_s from LHAPDF
!       LHAPARM(17).NE.'LHAPDF' => Default (same as LHAPDF v1/v3)
! 
!   Extrapolation of PDFs outside LHAPDF validity range given by
!   [Xmin,Xmax] and [Q2min,Q2max]; DEFAULT => PDFs "frozen" at the
!   boundaries
! 
!       LHAPARM(18).EQ.'EXTRAPOLATE' => Extrapolate PDFs on OWN RISK
!                           WARNING: Crazy values can be returned
! 
!   Printout of initialization information in PDFSET (by default)
! 
!       LHAPARM(19).EQ.'SILENT' => No printout (silent mode)
!       LHAPARM(19).EQ.'LOWKEY' => Print 5 times (almost silent mode)
! 
!
!*********************************************************************
!
! $Id: lhaglue.f 1448 2013-09-24 15:03:20Z whalley $
!
! $Log$
! Revision 1.7  2005/12/02 14:50:54  whalley
! Changes for new CTEQ code/AB sets
!
! Revision 1.6  2005/10/18 15:35:48  whalley
! fix to allow LHAPATH to be user defined as well as lhapdf-config
!
! Revision 1.5  2005/10/18 11:47:48  whalley
! Change to only set LHAPATH once per run
!
! Revision 1.1.1.2  1996/10/30 08:29:06  cernlib
! Version 7.04
!
! Revision 1.1.1.1  1996/04/12 15:29:26  plothow
! Version 7.01
!
! v5.0  06-Oct-2005  Major change to allow multiset-initializations 
! v4.0  28-Apr-2005  PDFSTA routine; option to use Alfa_s from LHAPDF
! v4.0  21-Mar-2005  Photon/pion/new p PDFs, updated for LHAPDF v4
! v3.1  26-Apr-2004  New numbering scheme, updated for LHAPDF v2/v3
! v3.0  23-Jan-2004  HERWIG interface added
! v2.0  20-Sep-2003  PDFLIB style adopted
! v1.0  05-Mar-2003  First working version from PYTHIA to LHAPDF v1
! 
! interface to LHAPDF library
!*********************************************************************


! PDFSET
! Initialization for use of parton distributions
!  according to the LHAPDF interface.
! 
! v4.0  28-Apr-2005  Option to use Alfa_s from LHAPDF
! v4.0  21-Mar-2005  Photon/pion/new p PDFs, updated for LHAPDF v4
! v3.1  26-Apr-2004  New numbering scheme
! v3.0  23-Jan-2004  HERWIG interface added
! 
! Interface to LHAPDF library
subroutine pdfset(parm,value)
  ! Double precision and integer declarations.
  implicit double precision(a-h, o-z)
  implicit integer(i-n)
  ! Common blocks
  include 'commonlhapdf.inc'
  include 'commonlhasets.inc'
  include 'commonlhapdfc.inc'
  include 'commonlhacontrol.inc'
  include 'commonlhaglsta.inc'
  ! additions for multiset use
  double precision xxmin(nmxset),xxmax(nmxset),qq2min(nmxset),qq2max(nmxset)
  save xxmin,xxmax,qq2min,qq2max
  ! Commonblocks.
  common/pydat1/mstu(200),paru(200),mstj(200),parj(200)
  save /pydat1/
  common/pypars/mstp(200),parp(200),msti(200),pari(200)
  save /pypars/
  ! following 2 for earlier Pythia versions
  common/ludat1/mstu5(200),paru5(200),mstj5(200),parj5(200)
  save /ludat1/
  ! Interface to LHAPDFLIB.
  double precision qcdlha4, qcdlha5
  integer nfllha
  common/lhapdfr/qcdlha4, qcdlha5, nfllha
  save /lhapdfr/
  integer lhaextrp
  common/lhapdfe/lhaextrp
  save /lhapdfe/
  integer lhasilent
  common/lhasilent/lhasilent
  save /lhasilent/
  ! Interface to PDFLIB.
  common/w50511/ nptypepdfl,ngrouppdfl,nsetpdfl,modepdfl,nflpdfl,lopdfl,tmaspdfl
  save /w50511/
  double precision tmaspdfl
  ! Interface to PDFLIB.
  common/w50512/qcdl4,qcdl5
  save /w50512/
  double precision qcdl4,qcdl5
  ! Interface to PDFLIB.
  common/w50513/xmin,xmax,q2min,q2max
  save /w50513/
  double precision xmin,xmax,q2min,q2max
  ! Local arrays and character variables (NOT USED here DB)
  character*20 parm(20)
  double precision value(20)
  integer lhapathlen
  integer :: lhainput = 1
  !integer lhaselect
  integer lhaprint
  integer lhaonce
  integer lhafive
  save lhaonce
  save lhafive
  data lhaonce/0/
  data lhafive/0/
  logical first
  data first/.true./
  character*512 dirpath
  save first

  ! Initialise common blocks  
  call commoninit()

  if (first) then
     call getdirpath(dirpath)
     first = .FALSE.
  endif

  ! Init
  lhaextrp = 0
  if(lhaparm(18).EQ.'EXTRAPOLATE') then  ! Extrapolate PDFs on own risk
     lhaextrp = 1
  endif
  lhasilent = 0
  if (lhaparm(19).EQ.'SILENT') then    !  No printout (silent MODE)
     lhasilent = 1
  elseif (lhaparm(19).EQ.'LOWKEY') then ! print 5 times (lowkey mode)
     if (lhafive .lt. 6) then
        lhafive = lhafive + 1
     else
        lhasilent = 1
     endif
  endif
  if (parm(1).EQ.'NPTYPE') then        !  pythia
     if(mstp(181).ge.6) then
        lhaprint = mstu(11)
     else
        lhaprint = mstu5(11)
     endif
     if(value(1) .eq. 1) then         !   nucleon
        lhainput = abs(mstp(51))
     elseif(value(1) .eq. 2) then     !   pion
        lhainput = abs(mstp(53))
     elseif(value(1) .eq. 3) then     !   photon
        lhainput = abs(mstp(55))
     endif
     if (lhasilent .ne. 1) print *,'==== PYTHIA WILL USE LHAPDF ===='
  elseif(parm(1).EQ.'HWLHAPDF') then  !  herwig
     lhainput = abs(int(value(1)))
     if(lhaonce.eq.lhainput) return
     if(lhasilent .ne. 1) print *,'==== HERWIG WILL USE LHAPDF ===='
     lhaprint = 6
     lhaonce = lhainput
  elseif(parm(1).EQ.'DEFAULT') then  !  stand-alone
     lhainput = abs(int(value(1)))
     if(lhaonce.eq.lhainput) return
     if(lhasilent .ne. 1) print *,'==== STAND-ALONE LHAGLUE MODE TO USE LHAPDF ===='
     lhaprint = 6
     lhaonce = lhainput
  else
     print *,'== UNKNOWN LHAPDF INTERFACE CALL! STOP EXECUTION! =='
     stop
  endif
  ! Initialize parton distributions: LHAPDFLIB.
  lhapathlen=index(lhapath,' ') - 1
  lhaset = lhainput
  xmin = 1.0D-6      ! X_min for current PDF set
  xmax = 1.0D0       ! X_max for current PDF set
  q2min = 1.0D0**2   ! Q**2_min scale for current PDF set [GeV]
  q2max = 1.0D5**2   ! Q**2_max scale for current PDF set [GeV]
  ! 
  ! Protons
  ! 
  ! CTEQ Family
  if (((lhainput.ge.10000).and.(lhainput .le. 10999)).or.((lhainput.ge.19000).and.(lhainput .le. 19999))) then
     q2max = 1.0d08
     if ((lhainput .ge. 10000) .and. (lhainput .le. 10040)) then
        lhaset = 10000
        lhaname=lhapath(1:lhapathlen)//'/cteq6.LHpdf'
        q2min = 1.69d0
     elseif((lhainput .ge. 10041) .and. (lhainput .le. 10041)) then
        lhaset = 10041
        lhaname=lhapath(1:lhapathlen)//'/cteq6l.LHpdf'
        q2min = 1.69d0
     elseif((lhainput .ge. 10042) .and. (lhainput .le. 10042)) then
        lhaset = 10042
        lhaname=lhapath(1:lhapathlen)//'/cteq6ll.LHpdf'
        q2min = 1.69d0
     elseif((lhainput .ge. 10050) .and. (lhainput .le. 10090)) then
        lhaset = 10050
        lhaname=lhapath(1:lhapathlen)//'/cteq6mE.LHgrid'
        q2min = 1.69d0
     elseif((lhainput .ge. 10100) .and. (lhainput .le. 10140)) then
        lhaset = 10100
        lhaname=lhapath(1:lhapathlen)//'/cteq61.LHpdf'
        q2min = 1.69d0
     elseif((lhainput .ge. 10150) .and. (lhainput .le. 10190)) then
        lhaset = 10150
        lhaname=lhapath(1:lhapathlen)//'/cteq61.LHgrid'
        q2min = 1.69d0
     elseif((lhainput .ge. 10250) .and. (lhainput .le. 10269)) then
        lhaset = 10250
        lhaname=lhapath(1:lhapathlen)//'/cteq6AB.LHgrid'
        q2min = 1.69d0
     elseif((lhainput .ge. 10350) .and. (lhainput .le. 10390)) then
        lhaset = 10350
        lhaname=lhapath(1:lhapathlen)//'/cteq65.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-7
     elseif((lhainput .ge. 10450) .and. (lhainput .le. 10456)) then
        lhaset = 10450
        lhaname=lhapath(1:lhapathlen)//'/cteq65c.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-7
     elseif((lhainput .ge. 10460) .and. (lhainput .le. 10467)) then
        lhaset = 10460
        lhaname=lhapath(1:lhapathlen)//'/cteq65s.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-7
     elseif((lhainput .ge. 10550) .and. (lhainput .le. 10594)) then
        lhaset = 10550
        lhaname=lhapath(1:lhapathlen)//'/cteq66.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10595) .and. (lhainput .le. 10599)) then
        lhaset = 10595
        lhaname=lhapath(1:lhapathlen)//'/cteq66alphas.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10650) .and. (lhainput .le. 10653)) then
        lhaset = 10650
        lhaname=lhapath(1:lhapathlen)//'/cteq66c.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10660) .and. (lhainput .le. 10660)) then
        lhaset = 10660
        lhaname=lhapath(1:lhapathlen)//'/cteq66a0.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10661) .and. (lhainput .le. 10661)) then
        lhaset = 10661
        lhaname=lhapath(1:lhapathlen)//'/cteq66a1.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10662) .and. (lhainput .le. 10662)) then
        lhaset = 10662
        lhaname=lhapath(1:lhapathlen)//'/cteq66a2.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10663) .and. (lhainput .le. 10663)) then
        lhaset = 10663
        lhaname=lhapath(1:lhapathlen)//'/cteq66a3.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10670) .and. (lhainput .le. 10677)) then
        lhaset = 10670
        lhaname=lhapath(1:lhapathlen)//'/cteq6lg.LHgrid'
        q2min = 1.69d0
     elseif((lhainput .ge. 10770) .and. (lhainput .le. 10770)) then
        lhaset = 10770
        lhaname=lhapath(1:lhapathlen)//'/CT09MCS.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10771) .and. (lhainput .le. 10771)) then
        lhaset = 10771
        lhaname=lhapath(1:lhapathlen)//'/CT09MC1.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10772) .and. (lhainput .le. 10772)) then
        lhaset = 10772
        lhaname=lhapath(1:lhapathlen)//'/CT09MC2.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10800) .and. (lhainput .le. 10852)) then
        lhaset = 10800
        lhaname=lhapath(1:lhapathlen)//'/CT10.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10860) .and. (lhainput .le. 10870)) then
        lhaset = 10860
        lhaname=lhapath(1:lhapathlen)//'/CT10as.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10900) .and. (lhainput .le. 10952)) then
        lhaset = 10900
        lhaname=lhapath(1:lhapathlen)//'/CT10w.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10960) .and. (lhainput .le. 10970)) then
        lhaset = 10960
        lhaname=lhapath(1:lhapathlen)//'/CT10was.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10980) .and. (lhainput .le. 10980)) then
        lhaset = 10980
        lhaname=lhapath(1:lhapathlen)//'/CT10f3.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10981) .and. (lhainput .le. 10981)) then
        lhaset = 10981
        lhaname=lhapath(1:lhapathlen)//'/CT10f4.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10982) .and. (lhainput .le. 10982)) then
        lhaset = 10982
        lhaname=lhapath(1:lhapathlen)//'/CT10wf3.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 10983) .and. (lhainput .le. 10983)) then
        lhaset = 10983
        lhaname=lhapath(1:lhapathlen)//'/CT10wf4.LHgrid'
        q2min = 1.69d0
        q2max = 1.0d10
        xmin =  1.0d-8
     elseif((lhainput .ge. 19050) .and. (lhainput .le. 19050)) then
        lhaset = 19050
        lhaname=lhapath(1:lhapathlen)//'/cteq5m.LHgrid'
        xmin=1.0d-5
     elseif((lhainput .ge. 19051) .and. (lhainput .le. 19051)) then
        lhaset = 19051
        lhaname=lhapath(1:lhapathlen)//'/cteq5m1.LHgrid'
        xmin=1.0d-5
     elseif((lhainput .ge. 19053) .and. (lhainput .le. 19053)) then
        lhaset = 19053
        lhaname=lhapath(1:lhapathlen)//'/cteq5f3.LHgrid'
        xmin=1.0d-5
     elseif((lhainput .ge. 19054) .and. (lhainput .le. 19054)) then
        lhaset = 19054
        lhaname=lhapath(1:lhapathlen)//'/cteq5f4.LHgrid'
        xmin=1.0d-5
     elseif((lhainput .ge. 19060) .and. (lhainput .le. 19060)) then
        lhaset = 19060
        lhaname=lhapath(1:lhapathlen)//'/cteq5d.LHgrid'
        xmin=1.0d-5
     elseif((lhainput .ge. 19070) .and. (lhainput .le. 19070)) then
        lhaset = 19070
        lhaname=lhapath(1:lhapathlen)//'/cteq5l.LHgrid'
        xmin=1.0d-5
     elseif((lhainput .ge. 19150) .and. (lhainput .le. 19150)) then
        lhaset = 19150
        lhaname=lhapath(1:lhapathlen)//'/cteq4m.LHgrid'
        q2min = 2.56d0
        xmin=1.0d-5
     elseif((lhainput .ge. 19160) .and. (lhainput .le. 19160)) then
        lhaset = 19160
        lhaname=lhapath(1:lhapathlen)//'/cteq4d.LHgrid'
        q2min = 2.56d0
        xmin=1.0d-5
     elseif((lhainput .ge. 19170) .and. (lhainput .le. 19170)) then
        lhaset = 19170
        lhaname=lhapath(1:lhapathlen)//'/cteq4l.LHgrid'
        q2min = 2.56d0
        xmin=1.0d-5
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     
   ! New CT10 (ct12 format) family)
  elseif((lhainput .ge. 11000) .and. (lhainput .le. 11280)) then 
     q2min = 1.69d0
     q2max = 1.0d10
     xmin =  1.0d-8
 !ct10nlo
     if ((lhainput .ge. 11000) .and. (lhainput .le. 11052)) then
        lhaset = 11000
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo.LHgrid'
 !ct10nlo_as_xxxx
     elseif ((lhainput .ge. 11062) .and. (lhainput .le. 11062)) then
        lhaset = 11062
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_as_0112.LHgrid'
     elseif ((lhainput .ge. 11063) .and. (lhainput .le. 11063)) then
        lhaset = 11063
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_as_0113.LHgrid'
     elseif ((lhainput .ge. 11064) .and. (lhainput .le. 11064)) then
        lhaset = 11064
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_as_0114.LHgrid'
     elseif ((lhainput .ge. 11065) .and. (lhainput .le. 11065)) then
        lhaset = 11065
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_as_0115.LHgrid'
     elseif ((lhainput .ge. 11066) .and. (lhainput .le. 11066)) then
        lhaset = 11066
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_as_0116.LHgrid'
     elseif ((lhainput .ge. 11067) .and. (lhainput .le. 11067)) then
        lhaset = 11067
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_as_0117.LHgrid'
     elseif ((lhainput .ge. 11068) .and. (lhainput .le. 11068)) then
        lhaset = 11068
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_as_0118.LHgrid'
     elseif ((lhainput .ge. 11069) .and. (lhainput .le. 11069)) then
        lhaset = 11069
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_as_0119.LHgrid'
     elseif ((lhainput .ge. 11070) .and. (lhainput .le. 11070)) then
        lhaset = 11070
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_as_0120.LHgrid'
     elseif ((lhainput .ge. 11071) .and. (lhainput .le. 11071)) then
        lhaset = 11071
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_as_0121.LHgrid'
     elseif ((lhainput .ge. 11072) .and. (lhainput .le. 11072)) then
        lhaset = 11072
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_as_0122.LHgrid'
     elseif ((lhainput .ge. 11073) .and. (lhainput .le. 11073)) then
        lhaset = 11073
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_as_0123.LHgrid'
     elseif ((lhainput .ge. 11074) .and. (lhainput .le. 11074)) then
        lhaset = 11074
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_as_0124.LHgrid'
     elseif ((lhainput .ge. 11075) .and. (lhainput .le. 11075)) then
        lhaset = 11075
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_as_0125.LHgrid'
     elseif ((lhainput .ge. 11076) .and. (lhainput .le. 11076)) then
        lhaset = 11076
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_as_0126.LHgrid'
     elseif ((lhainput .ge. 11077) .and. (lhainput .le. 11077)) then
        lhaset = 11077
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_as_0127.LHgrid'
 !ct10nlo_nf3/4
     elseif ((lhainput .ge. 11080) .and. (lhainput .le. 11081)) then
        lhaset = 11080
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_nf3.LHgrid'
     elseif ((lhainput .ge. 11082) .and. (lhainput .le. 11083)) then
        lhaset = 11082
        lhaname=lhapath(1:lhapathlen)//'/CT10nlo_nf4.LHgrid'
 !ct10wnlo
     elseif ((lhainput .ge. 11100) .and. (lhainput .le. 11152)) then
        lhaset = 11100
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo.LHgrid'
 !ct10wnlo_as_xxxx
     elseif ((lhainput .ge. 11162) .and. (lhainput .le. 11162)) then
        lhaset = 11162
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_as_0112.LHgrid'
     elseif ((lhainput .ge. 11163) .and. (lhainput .le. 11163)) then
        lhaset = 11163
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_as_0113.LHgrid'
     elseif ((lhainput .ge. 11164) .and. (lhainput .le. 11164)) then
        lhaset = 11164
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_as_0114.LHgrid'
     elseif ((lhainput .ge. 11165) .and. (lhainput .le. 11165)) then
        lhaset = 11165
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_as_0115.LHgrid'
     elseif ((lhainput .ge. 11166) .and. (lhainput .le. 11166)) then
        lhaset = 11166
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_as_0116.LHgrid'
     elseif ((lhainput .ge. 11167) .and. (lhainput .le. 11167)) then
        lhaset = 11167
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_as_0117.LHgrid'
     elseif ((lhainput .ge. 11168) .and. (lhainput .le. 11168)) then
        lhaset = 11168
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_as_0118.LHgrid'
     elseif ((lhainput .ge. 11169) .and. (lhainput .le. 11169)) then
        lhaset = 11169
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_as_0119.LHgrid'
     elseif ((lhainput .ge. 11170) .and. (lhainput .le. 11170)) then
        lhaset = 11170
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_as_0120.LHgrid'
     elseif ((lhainput .ge. 11171) .and. (lhainput .le. 11171)) then
        lhaset = 11171
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_as_0121.LHgrid'
     elseif ((lhainput .ge. 11172) .and. (lhainput .le. 11172)) then
        lhaset = 11172
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_as_0122.LHgrid'
     elseif ((lhainput .ge. 11173) .and. (lhainput .le. 11173)) then
        lhaset = 11173
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_as_0123.LHgrid'
     elseif ((lhainput .ge. 11174) .and. (lhainput .le. 11174)) then
        lhaset = 11174
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_as_0124.LHgrid'
     elseif ((lhainput .ge. 11175) .and. (lhainput .le. 11175)) then
        lhaset = 11175
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_as_0125.LHgrid'
     elseif ((lhainput .ge. 11176) .and. (lhainput .le. 11176)) then
        lhaset = 11176
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_as_0126.LHgrid'
     elseif ((lhainput .ge. 11177) .and. (lhainput .le. 11177)) then
        lhaset = 11177
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_as_0127.LHgrid'
 !ct10wnlo_nf3/4
     elseif ((lhainput .ge. 11180) .and. (lhainput .le. 11181)) then
        lhaset = 11180
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_nf3.LHgrid'
     elseif ((lhainput .ge. 11182) .and. (lhainput .le. 11183)) then
        lhaset = 11182
        lhaname=lhapath(1:lhapathlen)//'/CT10wnlo_nf4.LHgrid'
 !ct10nnlo
     elseif ((lhainput .ge. 11200) .and. (lhainput .le. 11250)) then
        lhaset = 11200
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo.LHgrid'
 !ct10nnlo_as_xxxx
     elseif ((lhainput .ge. 11260) .and. (lhainput .le. 11260)) then
        lhaset = 11260
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0110.LHgrid'
     elseif ((lhainput .ge. 11261) .and. (lhainput .le. 11261)) then
        lhaset = 11261
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0111.LHgrid'
     elseif ((lhainput .ge. 11262) .and. (lhainput .le. 11262)) then
        lhaset = 11262
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0112.LHgrid'
     elseif ((lhainput .ge. 11263) .and. (lhainput .le. 11263)) then
        lhaset = 11263
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0113.LHgrid'
     elseif ((lhainput .ge. 11264) .and. (lhainput .le. 11264)) then
        lhaset = 11264
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0114.LHgrid'
     elseif ((lhainput .ge. 11265) .and. (lhainput .le. 11265)) then
        lhaset = 11265
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0115.LHgrid'
     elseif ((lhainput .ge. 11266) .and. (lhainput .le. 11266)) then
        lhaset = 11266
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0116.LHgrid'
     elseif ((lhainput .ge. 11267) .and. (lhainput .le. 11267)) then
        lhaset = 11267
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0117.LHgrid'
     elseif ((lhainput .ge. 11268) .and. (lhainput .le. 11268)) then
        lhaset = 11268
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0118.LHgrid'
     elseif ((lhainput .ge. 11269) .and. (lhainput .le. 11269)) then
        lhaset = 11269
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0119.LHgrid'
     elseif ((lhainput .ge. 11270) .and. (lhainput .le. 11270)) then
        lhaset = 11270
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0120.LHgrid'
     elseif ((lhainput .ge. 11271) .and. (lhainput .le. 11271)) then
        lhaset = 11271
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0121.LHgrid'
     elseif ((lhainput .ge. 11272) .and. (lhainput .le. 11272)) then
        lhaset = 11272
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0122.LHgrid'
     elseif ((lhainput .ge. 11273) .and. (lhainput .le. 11273)) then
        lhaset = 11273
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0123.LHgrid'
     elseif ((lhainput .ge. 11274) .and. (lhainput .le. 11274)) then
        lhaset = 11274
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0124.LHgrid'
     elseif ((lhainput .ge. 11275) .and. (lhainput .le. 11275)) then
        lhaset = 11275
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0125.LHgrid'
     elseif ((lhainput .ge. 11276) .and. (lhainput .le. 11276)) then
        lhaset = 11276
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0126.LHgrid'
     elseif ((lhainput .ge. 11277) .and. (lhainput .le. 11277)) then
        lhaset = 11277
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0127.LHgrid'
     elseif ((lhainput .ge. 11278) .and. (lhainput .le. 11278)) then
        lhaset = 11278
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0128.LHgrid'
     elseif ((lhainput .ge. 11279) .and. (lhainput .le. 11279)) then
        lhaset = 11279
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0129.LHgrid'
     elseif ((lhainput .ge. 11280) .and. (lhainput .le. 11280)) then
        lhaset = 11280
        lhaname=lhapath(1:lhapathlen)//'/CT10nnlo_as_0130.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif

  ! CJ12  (family)
  elseif((lhainput .ge. 12000) .and. (lhainput .le. 12238)) then 
     q2min = 1.69d0
     q2max = 1.0d10
     xmin =  1.0d-6
 !cj12min
     if ((lhainput .ge. 12000) .and. (lhainput .le. 12038)) then
        lhaset = 12000
        lhaname=lhapath(1:lhapathlen)//'/CJ12min.LHgrid'
 !cj12mid
     elseif ((lhainput .ge. 12100) .and. (lhainput .le. 12138)) then
        lhaset = 12100
        lhaname=lhapath(1:lhapathlen)//'/CJ12mid.LHgrid'
 !cj12max
     elseif ((lhainput .ge. 12200) .and. (lhainput .le. 12238)) then
        lhaset = 12200
        lhaname=lhapath(1:lhapathlen)//'/CJ12max.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif

     ! MRST Family
  elseif((lhainput .ge. 20000) .and. (lhainput .le. 20999)) then
     q2min = 1.25d0
     q2max = 1.0d07
     xmin = 1.0d-5
     if((lhainput .ge. 20000) .and. (lhainput .le. 20004)) then
        lhaset = 20000
        lhaname=lhapath(1:lhapathlen)//'/MRST2001nlo.LHpdf'
     elseif((lhainput .ge. 20050) .and. (lhainput .le. 20054)) then
        lhaset = 20050
        lhaname=lhapath(1:lhapathlen)//'/MRST2001nlo.LHgrid'
     elseif((lhainput .ge. 20060) .and. (lhainput .le. 20061)) then
        lhaset = 20060
        lhaname=lhapath(1:lhapathlen)//'/MRST2001lo.LHgrid'
     elseif((lhainput .ge. 20070) .and. (lhainput .le. 20074)) then
        lhaset = 20070
        lhaname=lhapath(1:lhapathlen)//'/MRST2001nnlo.LHgrid'
     elseif((lhainput .ge. 20100) .and. (lhainput .le. 20130)) then
        lhaset = 20100
        lhaname=lhapath(1:lhapathlen)//'/MRST2001E.LHpdf'
     elseif((lhainput .ge. 20150) .and. (lhainput .le. 20180)) then
        lhaset = 20150
        lhaname=lhapath(1:lhapathlen)//'/MRST2001E.LHgrid'
     elseif((lhainput .ge. 20200) .and. (lhainput .le. 20201)) then
        lhaset = 20200
        lhaname=lhapath(1:lhapathlen)//'/MRST2002nlo.LHpdf'
     elseif((lhainput .ge. 20250) .and. (lhainput .le. 20251)) then
        lhaset = 20250
        lhaname=lhapath(1:lhapathlen)//'/MRST2002nlo.LHgrid'
     elseif((lhainput .ge. 20270) .and. (lhainput .le. 20271)) then
        lhaset = 20270
        lhaname=lhapath(1:lhapathlen)//'/MRST2002nnlo.LHgrid'
     elseif((lhainput .ge. 20300) .and. (lhainput .le. 20301)) then
        lhaset = 20300
        lhaname=lhapath(1:lhapathlen)//'/MRST2003cnlo.LHpdf'
        q2min = 10.0d0
        xmin = 1.0d-3
     elseif((lhainput .ge. 20350) .and. (lhainput .le. 20351)) then
        lhaset = 20350
        lhaname=lhapath(1:lhapathlen)//'/MRST2003cnlo.LHgrid'
        q2min = 10.0d0
        xmin = 1.0d-3
     elseif((lhainput .ge. 20370) .and. (lhainput .le. 20371)) then
        lhaset = 20370
        lhaname=lhapath(1:lhapathlen)//'/MRST2003cnnlo.LHgrid'
        q2min = 7.0d0
        xmin = 1.0d-3
     elseif((lhainput .ge. 20400) .and. (lhainput .le. 20401)) then
        lhaset = 20400
        lhaname=lhapath(1:lhapathlen)//'/MRST2004nlo.LHpdf'
     elseif((lhainput .ge. 20406) .and. (lhainput .le. 20407)) then
        lhaset = 20406
        lhaname=lhapath(1:lhapathlen)//'/MRST2004FF3nlo.LHpdf'
     elseif((lhainput .ge. 20408) .and. (lhainput .le. 20409)) then
        lhaset = 20408
        lhaname=lhapath(1:lhapathlen)//'/MRST2004FF4nlo.LHpdf'
     elseif((lhainput .ge. 20450) .and. (lhainput .le. 20451)) then
        lhaset = 20450
        lhaname=lhapath(1:lhapathlen)//'/MRST2004nlo.LHgrid'
     elseif((lhainput .ge. 20452) .and. (lhainput .le. 20453)) then
        lhaset = 20452
        lhaname=lhapath(1:lhapathlen)//'/MRST2004FF3lo.LHgrid'
     elseif((lhainput .ge. 20454) .and. (lhainput .le. 20455)) then
        lhaset = 20454
        lhaname=lhapath(1:lhapathlen)//'/MRST2004FF4lo.LHgrid'
     elseif((lhainput .ge. 20456) .and. (lhainput .le. 20457)) then
        lhaset = 20456
        lhaname=lhapath(1:lhapathlen)//'/MRST2004FF3nlo.LHgrid'
     elseif((lhainput .ge. 20458) .and. (lhainput .le. 20459)) then
        lhaset = 20458
        lhaname=lhapath(1:lhapathlen)//'/MRST2004FF4nlo.LHgrid'
     elseif((lhainput .ge. 20460) .and. (lhainput .le. 20462)) then
        lhaset = 20460
        lhaname=lhapath(1:lhapathlen)//'/MRST2004qed.LHgrid'
     elseif((lhainput .ge. 20470) .and. (lhainput .le. 20471)) then
        lhaset = 20470
        lhaname=lhapath(1:lhapathlen)//'/MRST2004nnlo.LHgrid'
     elseif((lhainput .ge. 20550) .and. (lhainput .le. 20580)) then
        lhaset = 20550
        lhaname=lhapath(1:lhapathlen)//'/MRST2006nnlo.LHgrid'
        q2min = 1.0d0
        q2max = 1.0d09
        xmin = 1.0d-6
     elseif((lhainput .ge. 20650) .and. (lhainput .le. 20650)) then
        lhaset = 20650
        lhaname=lhapath(1:lhapathlen)//'/MRST2007lomod.LHgrid'
     elseif((lhainput .ge. 20651) .and. (lhainput .le. 20651)) then
        lhaset = 20651
        lhaname=lhapath(1:lhapathlen)//'/MRSTMCal.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! next are MSTW
  elseif((lhainput .ge. 21000) .and. (lhainput .le. 23896)) then
     q2min = 1.0d0
     q2max = 1.0d09
     xmin = 1.0d-6
     if((lhainput .ge. 21000) .and. (lhainput .le. 21040)) then
        lhaset = 21000
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008lo68cl.LHgrid'
     elseif((lhainput .ge. 21041) .and. (lhainput .le. 21080)) then
        lhaset = 21040
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008lo90cl.LHgrid'
     elseif((lhainput .ge. 21100) .and. (lhainput .le. 21140)) then
        lhaset = 21100
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo68cl.LHgrid'
     elseif((lhainput .ge. 21141) .and. (lhainput .le. 21180)) then
         lhaset = 21140
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo90cl.LHgrid'
     elseif((lhainput .ge. 21200) .and. (lhainput .le. 21240)) then
        lhaset = 21200
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo68cl.LHgrid'
     elseif((lhainput .ge. 21241) .and. (lhainput .le. 21280)) then
        lhaset = 21240
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo90cl.LHgrid'
!
     elseif((lhainput .ge. 22000) .and. (lhainput .le. 22021)) then
        lhaset = 22000
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo_asmzrange.LHgrid'
     elseif((lhainput .ge. 22100) .and. (lhainput .le. 22140)) then
        lhaset = 22100
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo68cl_asmz+68cl.LHgrid'
     elseif((lhainput .ge. 22150) .and. (lhainput .le. 22190)) then
        lhaset = 22150
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo68cl_asmz-68cl.LHgrid'
     elseif((lhainput .ge. 22200) .and. (lhainput .le. 22240)) then
        lhaset = 22200
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo68cl_asmz+68clhalf.LHgrid'
     elseif((lhainput .ge. 22250) .and. (lhainput .le. 22290)) then
        lhaset = 22250
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo68cl_asmz-68clhalf.LHgrid'
     elseif((lhainput .ge. 22300) .and. (lhainput .le. 22340)) then
        lhaset = 22300
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo90cl_asmz+90cl.LHgrid'
     elseif((lhainput .ge. 22350) .and. (lhainput .le. 22390)) then
        lhaset = 22350
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo90cl_asmz-90cl.LHgrid'
     elseif((lhainput .ge. 22400) .and. (lhainput .le. 22440)) then
        lhaset = 22400
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo90cl_asmz+90clhalf.LHgrid'
     elseif((lhainput .ge. 22450) .and. (lhainput .le. 22490)) then
        lhaset = 22450
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo90cl_asmz-90clhalf.LHgrid'
!
     elseif((lhainput .ge. 22500) .and. (lhainput .le. 22521)) then
        lhaset = 22500
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo_asmzrange.LHgrid'
     elseif((lhainput .ge. 22600) .and. (lhainput .le. 22640)) then
        lhaset = 22600
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo68cl_asmz+68cl.LHgrid'
     elseif((lhainput .ge. 22650) .and. (lhainput .le. 22690)) then
        lhaset = 22650
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo68cl_asmz-68cl.LHgrid'
     elseif((lhainput .ge. 22700) .and. (lhainput .le. 22740)) then
        lhaset = 22700
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo68cl_asmz+68clhalf.LHgrid'
     elseif((lhainput .ge. 22750) .and. (lhainput .le. 22790)) then
        lhaset = 22750
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo68cl_asmz-68clhalf.LHgrid'
     elseif((lhainput .ge. 22800) .and. (lhainput .le. 22840)) then
        lhaset = 22800
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo90cl_asmz+90cl.LHgrid'
     elseif((lhainput .ge. 22850) .and. (lhainput .le. 22890)) then
        lhaset = 22850
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo90cl_asmz-90cl.LHgrid'
     elseif((lhainput .ge. 22900) .and. (lhainput .le. 22940)) then
        lhaset = 22900
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo90cl_asmz+90clhalf.LHgrid'
     elseif((lhainput .ge. 22950) .and. (lhainput .le. 22990)) then
        lhaset = 22950
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo90cl_asmz-90clhalf.LHgrid'
!
     elseif((lhainput .ge. 23000) .and. (lhainput .le. 23040)) then
        lhaset = 23000
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008lo68cl_nf3.LHgrid'
     elseif((lhainput .ge. 23041) .and. (lhainput .le. 23080)) then
        lhaset = 23041
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008lo90cl_nf3.LHgrid'
     elseif((lhainput .ge. 23100) .and. (lhainput .le. 23140)) then
        lhaset = 23100
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008lo68cl_nf4.LHgrid'
     elseif((lhainput .ge. 23141) .and. (lhainput .le. 23180)) then
        lhaset = 23141
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008lo90cl_nf4.LHgrid'
!
     elseif((lhainput .ge. 23200) .and. (lhainput .le. 23240)) then
        lhaset = 23200
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo68cl_nf3.LHgrid'
     elseif((lhainput .ge. 23241) .and. (lhainput .le. 23280)) then
        lhaset = 23241
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo90cl_nf3.LHgrid'
     elseif((lhainput .ge. 23300) .and. (lhainput .le. 23340)) then
        lhaset = 23300
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo68cl_nf4.LHgrid'
     elseif((lhainput .ge. 23341) .and. (lhainput .le. 23380)) then
        lhaset = 23341
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo90cl_nf4.LHgrid'
!
     elseif((lhainput .ge. 23400) .and. (lhainput .le. 23414)) then
        lhaset = 23400
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo_mcrange.LHgrid'
     elseif((lhainput .ge. 23420) .and. (lhainput .le. 23434)) then
        lhaset = 23420
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo_mcrange_nf3.LHgrid'
     elseif((lhainput .ge. 23440) .and. (lhainput .le. 23454)) then
        lhaset = 23440
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo_mcrange_fixasmz.LHgrid'
     elseif((lhainput .ge. 23460) .and. (lhainput .le. 23474)) then
        lhaset = 23460
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo_mcrange_fixasmz_nf3.LHgrid'
     elseif((lhainput .ge. 23480) .and. (lhainput .le. 23486)) then
        lhaset = 23480
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo_mbrange.LHgrid'
     elseif((lhainput .ge. 23490) .and. (lhainput .le. 23496)) then
        lhaset = 23490
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nlo_mbrange_nf4.LHgrid'
!
     elseif((lhainput .ge. 23500) .and. (lhainput .le. 23540)) then
        lhaset = 23500
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo68cl_nf3.LHgrid'
     elseif((lhainput .ge. 23541) .and. (lhainput .le. 23580)) then
        lhaset = 23541
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo90cl_nf3.LHgrid'
     elseif((lhainput .ge. 23600) .and. (lhainput .le. 23640)) then
        lhaset = 23600
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo68cl_nf4.LHgrid'
     elseif((lhainput .ge. 23641) .and. (lhainput .le. 23680)) then
        lhaset = 23641
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo90cl_nf4.LHgrid'
!
     elseif((lhainput .ge. 23700) .and. (lhainput .le. 23714)) then
        lhaset = 23700
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo_mcrange.LHgrid'
     elseif((lhainput .ge. 23720) .and. (lhainput .le. 23734)) then
        lhaset = 23720
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo_mcrange_nf3.LHgrid'
     elseif((lhainput .ge. 23740) .and. (lhainput .le. 23754)) then
        lhaset = 23740
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo_mcrange_fixasmz.LHgrid'
     elseif((lhainput .ge. 23760) .and. (lhainput .le. 23774)) then
        lhaset = 23760
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo_mcrange_fixasmz_nf3.LHgrid'
     elseif((lhainput .ge. 23780) .and. (lhainput .le. 23786)) then
        lhaset = 23780
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo_mbrange.LHgrid'
     elseif((lhainput .ge. 23790) .and. (lhainput .le. 23796)) then
        lhaset = 23790
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008nnlo_mbrange_nf4.LHgrid'
     elseif((lhainput .ge. 23800) .and. (lhainput .le. 23846)) then
        lhaset = 23800
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008CPdeutnlo68cl.LHgrid'
     elseif((lhainput .ge. 23850) .and. (lhainput .le. 23896)) then
        lhaset = 23850
        lhaname=lhapath(1:lhapathlen)//'/MSTW2008CPdeutnnlo68cl.LHgrid'
!
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! old MRST98 sets
  elseif((lhainput .ge. 29000) .and. (lhainput .le. 29999)) then
     q2min = 1.25d0
     q2max = 1.0d07
     xmin = 1.0d-5
     if((lhainput .ge. 29000) .and. (lhainput .le. 29003)) then
        lhaset = 29000
        lhaname=lhapath(1:lhapathlen)//'/MRST98.LHpdf'
     elseif((lhainput .ge. 29040) .and. (lhainput .le. 29045)) then
        lhaset = 29040
        lhaname=lhapath(1:lhapathlen)//'/MRST98lo.LHgrid'
     elseif((lhainput .ge. 29050) .and. (lhainput .le. 29055)) then
        lhaset = 29050
        lhaname=lhapath(1:lhapathlen)//'/MRST98nlo.LHgrid'
     elseif((lhainput .ge. 29060) .and. (lhainput .le. 29065)) then
        lhaset = 29060
        lhaname=lhapath(1:lhapathlen)//'/MRST98dis.LHgrid'
     elseif((lhainput .ge. 29070) .and. (lhainput .le. 29071)) then
        lhaset = 29070
        lhaname=lhapath(1:lhapathlen)//'/MRST98ht.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! Fermi Family
  elseif((lhainput .ge. 30000) .and. (lhainput .le. 39999)) then
     if((lhainput .ge. 30100) .and. (lhainput .le. 30200)) then
        lhaset = 30100
        lhaname=lhapath(1:lhapathlen)//'/Fermi2002_100.LHpdf'
     elseif((lhainput .ge. 31000) .and. (lhainput .le. 32000)) then
        lhaset = 31000
        lhaname=lhapath(1:lhapathlen)//'/Fermi2002_1000.LHpdf'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
   ! a02m Family 
  elseif((lhainput .ge. 40350) .and. (lhainput .le. 40567)) then
     xmin = 1.0d-7
     q2min = 0.8d0
     q2max = 2.0d08
     if((lhainput .ge. 40350) .and. (lhainput .le. 40367)) then
        lhaset = 40350
        lhaname=lhapath(1:lhapathlen)//'/a02m_lo.LHgrid'
     elseif((lhainput .ge. 40450) .and. (lhainput .le. 40467)) then
        lhaset = 40450
        lhaname=lhapath(1:lhapathlen)//'/a02m_nlo.LHgrid'
     elseif((lhainput .ge. 40550) .and. (lhainput .le. 40567)) then
        lhaset = 40550 
        lhaname=lhapath(1:lhapathlen)//'/a02m_nnlo.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
   ! abkm09 Family 
  elseif((lhainput .ge. 40650) .and. (lhainput .le. 40975)) then
     xmin = 1.0d-7
     q2min = 0.8d0
     q2max = 2.0d08
     if((lhainput .ge. 40650) .and. (lhainput .le. 40675)) then
        lhaset = 40650 
        lhaname=lhapath(1:lhapathlen)//'/abkm09_3_nlo.LHgrid'
     elseif((lhainput .ge. 40750) .and. (lhainput .le. 40775)) then
        lhaset = 40750 
        lhaname=lhapath(1:lhapathlen)//'/abkm09_3_nnlo.LHgrid'
     elseif((lhainput .ge. 40780) .and. (lhainput .le. 40805)) then
        lhaset = 40780 
        lhaname=lhapath(1:lhapathlen)//'/abkm09_4_nlo.LHgrid'
     elseif((lhainput .ge. 40810) .and. (lhainput .le. 40835)) then
        lhaset = 40810 
        lhaname=lhapath(1:lhapathlen)//'/abkm09_4_nnlo.LHgrid'
     elseif((lhainput .ge. 40850) .and. (lhainput .le. 40875)) then
        lhaset = 40850 
        lhaname=lhapath(1:lhapathlen)//'/abkm09_5_nlo.LHgrid'
     elseif((lhainput .ge. 40950) .and. (lhainput .le. 40975)) then
        lhaset = 40950 
        lhaname=lhapath(1:lhapathlen)//'/abkm09_5_nnlo.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
   ! abm11 Family 
  elseif((lhainput .ge. 42000) .and. (lhainput .le. 42246)) then
     xmin = 1.0d-7
     q2min = 0.8d0
     q2max = 2.0d08
     if((lhainput .ge. 42000) .and. (lhainput .le. 42028)) then
        lhaset = 42000 
        lhaname=lhapath(1:lhapathlen)//'/abm11_3n_nlo.LHgrid'
     elseif((lhainput .ge. 42030) .and. (lhainput .le. 42058)) then
        lhaset = 42030 
        lhaname=lhapath(1:lhapathlen)//'/abm11_4n_nlo.LHgrid'
     elseif((lhainput .ge. 42060) .and. (lhainput .le. 42088)) then
        lhaset = 42060 
        lhaname=lhapath(1:lhapathlen)//'/abm11_5n_nlo.LHgrid'
     elseif((lhainput .ge. 42100) .and. (lhainput .le. 42128)) then
        lhaset = 42100 
        lhaname=lhapath(1:lhapathlen)//'/abm11_3n_nnlo.LHgrid'
     elseif((lhainput .ge. 42130) .and. (lhainput .le. 42158)) then
        lhaset = 42130 
        lhaname=lhapath(1:lhapathlen)//'/abm11_4n_nnlo.LHgrid'
     elseif((lhainput .ge. 42160) .and. (lhainput .le. 42188)) then
        lhaset = 42160 
        lhaname=lhapath(1:lhapathlen)//'/abm11_5n_nnlo.LHgrid'
     elseif((lhainput .ge. 42200) .and. (lhainput .le. 42220)) then
        lhaset = 42200 
        lhaname=lhapath(1:lhapathlen)//'/abm11_5n_as_nlo.LHgrid'
     elseif((lhainput .ge. 42230) .and. (lhainput .le. 42246)) then
        lhaset = 42230 
        lhaname=lhapath(1:lhapathlen)//'/abm11_5n_as_nnlo.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! Alekhin Family
  elseif((lhainput .ge. 40000) .and. (lhainput .le. 41999)) then
     if((lhainput .ge. 40100) .and. (lhainput .le. 40200)) then
        lhaset = 40100
        lhaname=lhapath(1:lhapathlen)//'/Alekhin_100.LHpdf'
     elseif((lhainput .ge. 41000) .and. (lhainput .le. 41999)) then
        lhaset = 41000
        lhaname=lhapath(1:lhapathlen)//'/Alekhin_1000.LHpdf'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! Botje Family
  elseif((lhainput .ge. 50000) .and. (lhainput .le. 59999)) then
     if((lhainput .ge. 50100) .and. (lhainput .le. 50200)) then
        lhaset = 50100
        lhaname=lhapath(1:lhapathlen)//'/Botje_100.LHpdf'
     elseif((lhainput .ge. 51000) .and. (lhainput .le. 51999)) then
        lhaset = 51000
        lhaname=lhapath(1:lhapathlen)//'/Botje_1000.LHpdf'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! ZEUS (+ ATLAS) Family
  elseif((lhainput .ge. 60000) .and. (lhainput .le. 69999)) then
     q2min = 0.3d0
     q2max = 2.0d05
     if((lhainput .ge. 60000) .and. (lhainput .le. 60022)) then
        lhaset = 60000
        lhaname=lhapath(1:lhapathlen)//'/ZEUS2002_TR.LHpdf'
     elseif((lhainput .ge. 60100) .and. (lhainput .le. 60122)) then
        lhaset = 60100
        lhaname=lhapath(1:lhapathlen)//'/ZEUS2002_ZM.LHpdf'
     elseif((lhainput .ge. 60200) .and. (lhainput .le. 60222)) then
        lhaset = 60200
        lhaname=lhapath(1:lhapathlen)//'/ZEUS2002_FF.LHpdf'
     elseif((lhainput .ge. 60300) .and. (lhainput .le. 60322)) then
        lhaset = 60300
        lhaname=lhapath(1:lhapathlen)//'/ZEUS2005_ZJ.LHpdf'
     elseif((lhainput .ge. 60400) .and. (lhainput .le. 60422)) then
        xmin=1.0d-6
        xmax=1.0d0
        q2min=1.0d0
        q2max=1.0d09
        lhaset = 60400
        lhaname=lhapath(1:lhapathlen)//'/HERAPDF01.LHpdf'
     elseif((lhainput .ge. 60430) .and. (lhainput .le. 60444)) then
        xmin=1.0d-6
        xmax=1.0d0
        q2min=1.0d0
        q2max=2.0d08
        lhaset = 60430
        lhaname=lhapath(1:lhapathlen)//'/HERAPDF01.LHgrid'
     elseif((lhainput .ge. 60500) .and. (lhainput .le. 60520)) then
        xmin=1.0d-6
        xmax=1.0d0
        q2min=1.0d0
        q2max=2.0d08
        lhaset = 60500
        lhaname=lhapath(1:lhapathlen)//'/HERAPDF10_EIG.LHgrid'
     elseif((lhainput .ge. 60530) .and. (lhainput .le. 60543)) then
        xmin=1.0d-6
        xmax=1.0d0
        q2min=1.0d0
        q2max=2.0d08
        lhaset = 60530
        lhaname=lhapath(1:lhapathlen)//'/HERAPDF10_VAR.LHgrid'
     elseif((lhainput .ge. 60550) .and. (lhainput .le. 60561)) then
        xmin=1.0d-6
        xmax=1.0d0
        q2min=1.0d0
        q2max=2.0d08
        lhaset = 60550
        lhaname=lhapath(1:lhapathlen)//'/HERAPDF10_ALPHAS.LHgrid'
     elseif((lhainput .ge. 60600) .and. (lhainput .le. 60628)) then
        xmin=1.0d-8
        xmax=1.0d0
        q2min=1.0d0
        q2max=1.0d09
        lhaset = 60600
        lhaname=lhapath(1:lhapathlen)//'/HERAPDF15NNLO_EIG.LHgrid'
     elseif((lhainput .ge. 60630) .and. (lhainput .le. 60640)) then
        xmin=1.0d-8
        xmax=1.0d0
        q2min=1.0d0
        q2max=1.0d09
        lhaset = 60630
        lhaname=lhapath(1:lhapathlen)//'/HERAPDF15NNLO_VAR.LHgrid'
     elseif((lhainput .ge. 60650) .and. (lhainput .le. 60661)) then
        xmin=1.0d-8
        xmax=1.0d0
        q2min=1.0d0
        q2max=1.0d09
        lhaset = 60650
        lhaname=lhapath(1:lhapathlen)//'/HERAPDF15NNLO_ALPHAS.LHgrid'
     elseif((lhainput .ge. 60700) .and. (lhainput .le. 60720)) then
        xmin=1.0d-8
        xmax=1.0d0
        q2min=1.0d0
        q2max=1.0d09
        lhaset = 60700
        lhaname=lhapath(1:lhapathlen)//'/HERAPDF15NLO_EIG.LHgrid'
     elseif((lhainput .ge. 60730) .and. (lhainput .le. 60742)) then
        xmin=1.0d-8
        xmax=1.0d0
        q2min=1.0d0
        q2max=1.0d09
        lhaset = 60730
        lhaname=lhapath(1:lhapathlen)//'/HERAPDF15NLO_VAR.LHgrid'
     elseif((lhainput .ge. 60750) .and. (lhainput .le. 60761)) then
        xmin=1.0d-8
        xmax=1.0d0
        q2min=1.0d0
        q2max=1.0d09
        lhaset = 60750
        lhaname=lhapath(1:lhapathlen)//'/HERAPDF15NLO_ALPHAS.LHgrid'
     elseif((lhainput .ge. 60800) .and. (lhainput .le. 60824)) then
        xmin=1.0d-8
        xmax=1.0d0
        q2min=1.0d0
        q2max=1.0d09
        lhaset = 60800
        lhaname=lhapath(1:lhapathlen)//'/LHECNLO_EIG.LHgrid'
     elseif((lhainput .ge. 65000) .and. (lhainput .le. 65030)) then
        xmin=1.0d-7
        xmax=1.0d0
        q2min=1.0d0
        q2max=1.0d09
        lhaset = 65000
        lhaname=lhapath(1:lhapathlen)//'/ATLAS-epWZ12-EIG.LHgrid'
     elseif((lhainput .ge. 65040) .and. (lhainput .le. 65051)) then
        xmin=1.0d-7
        xmax=1.0d0
        q2min=1.0d0
        q2max=1.0d09
        lhaset = 65040
        lhaname=lhapath(1:lhapathlen)//'/ATLAS-epWZ12-VAR.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! H1 Family
  elseif((lhainput .ge. 70000) .and. (lhainput .le. 79999)) then
     q2min = 1.5d0
     q2max = 3.5d04
     xmin = 5.7d-5
     if((lhainput .ge. 70050) .and. (lhainput .le. 70050)) then
        lhaset = 70050
        lhaname=lhapath(1:lhapathlen)//'/H12000ms.LHgrid'
     elseif((lhainput .ge. 70051) .and. (lhainput .le. 70070)) then
        lhaset = 70050
        lhaname=lhapath(1:lhapathlen)//'/H12000msE.LHgrid'
     elseif((lhainput .ge. 70150) .and. (lhainput .le. 70150)) then
        lhaset = 70150
        lhaname=lhapath(1:lhapathlen)//'/H12000dis.LHgrid'
     elseif((lhainput .ge. 70151) .and. (lhainput .le. 70170)) then
        lhaset = 70150
        lhaname=lhapath(1:lhapathlen)//'/H12000disE.LHgrid'
     elseif((lhainput .ge. 70250) .and. (lhainput .le. 70250)) then
        lhaset = 70250
        lhaname=lhapath(1:lhapathlen)//'/H12000lo.LHgrid'
     elseif((lhainput .ge. 70251) .and. (lhainput .le. 70270)) then
        lhaset = 70250
        lhaname=lhapath(1:lhapathlen)//'/H12000loE.LHgrid'
        ! Temporarily removed on returning to original H12000 files
        ! elseif((lhainput .ge. 70350) .and. (lhainput .le. 70350)) then
        ! lhaset = 70350
        ! lhaname=lhapath(1:lhapathlen)//'/H12000lo2.LHgrid'
        ! elseif((lhainput .ge. 70351) .and. (lhainput .le. 70370)) then
        ! lhaset = 70350
        ! lhaname=lhapath(1:lhapathlen)//'/H12000lo2E.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! GRV/GJR Family
  elseif((lhainput .ge. 80000) .and. (lhainput .le. 89999)) then
     q2min = 0.8d0
     q2max = 2.0d06
     xmin = 1.0d-9
     if((lhainput .ge. 80050) .and. (lhainput .le. 80051)) then
        lhaset = 80050
        lhaname=lhapath(1:lhapathlen)//'/GRV98nlo.LHgrid'
     elseif((lhainput .ge. 80060) .and. (lhainput .le. 80060)) then
        lhaset = 80060
        lhaname=lhapath(1:lhapathlen)//'/GRV98lo.LHgrid'
     elseif((lhainput .ge. 80150) .and. (lhainput .le. 80151)) then
        q2min = 0.3d0
        q2max = 1.0d08
        lhaset = 80150
        lhaname=lhapath(1:lhapathlen)//'/GJR08lo.LHgrid'
     elseif((lhainput .ge. 80152) .and. (lhainput .le. 80152)) then
        q2min = 0.5d0
        q2max = 1.0d08
        lhaset = 80152
        lhaname=lhapath(1:lhapathlen)//'/GJR08FFdis.LHgrid'
     elseif((lhainput .ge. 80160) .and. (lhainput .le. 80186)) then
        q2min = 0.5d0
        q2max = 1.0d08
        lhaset = 80160
        lhaname=lhapath(1:lhapathlen)//'/GJR08FFnloE.LHgrid'
     elseif((lhainput .ge. 80260) .and. (lhainput .le. 80286)) then
        q2min = 0.5d0
        q2max = 1.0d08
        lhaset = 80260
        lhaname=lhapath(1:lhapathlen)//'/GJR08VFnloE.LHgrid'
     elseif((lhainput .ge. 80360) .and. (lhainput .le. 80386)) then
        q2min = 0.55d0
        q2max = 1.0d08
        lhaset = 80360
        lhaname=lhapath(1:lhapathlen)//'/JR09FFnnloE.LHgrid'
     elseif((lhainput .ge. 80460) .and. (lhainput .le. 80486)) then
        q2min = 0.55d0
        q2max = 1.0d08
        lhaset = 80460
        lhaname=lhapath(1:lhapathlen)//'/JR09VFnnloE.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
!...NNPDF  Family
  ELSEIF((LHAINPUT .GE. 90000) .AND. (LHAINPUT .LE. 98000)) THEN
     XMIN = 1.0D-9 
     Q2MIN = 2.0D0 
     Q2MAX = 1.0D08  
     IF((LHAINPUT .GE. 90000) .AND. (LHAINPUT .LE. 90100)) THEN
     	LHASET = 90000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF10_100.LHpdf'
     ELSEIF((LHAINPUT .GE. 90200) .AND. (LHAINPUT .LE. 90300))THEN
     	LHASET = 90200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF10_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 90400) .AND. (LHAINPUT .LE. 90500))THEN
     	LHASET = 90400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF11_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 90600) .AND. (LHAINPUT .LE. 90700))THEN
     	LHASET = 90600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF12_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 90800) .AND. (LHAINPUT .LE. 90900))THEN
     	LHASET = 90800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 90901) .AND. (LHAINPUT .LE. 90901))THEN
     	LHASET = 90901
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0114_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 90902) .AND. (LHAINPUT .LE. 90902))THEN
     	LHASET = 90902
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0115_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 90903) .AND. (LHAINPUT .LE. 90903))THEN
     	LHASET = 90903
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0116_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 90904) .AND. (LHAINPUT .LE. 90904))THEN
     	LHASET = 90904
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0117_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 90905) .AND. (LHAINPUT .LE. 90905))THEN
     	LHASET = 90905
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0118_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 90906) .AND. (LHAINPUT .LE. 90906))THEN
     	LHASET = 90906
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0120_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 90907) .AND. (LHAINPUT .LE. 90907))THEN
     	LHASET = 90907
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0121_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 90908) .AND. (LHAINPUT .LE. 90908))THEN
     	LHASET = 90908
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0122_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 90909) .AND. (LHAINPUT .LE. 90909))THEN
     	LHASET = 90909
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0123_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 90910) .AND. (LHAINPUT .LE. 90910))THEN
     	LHASET = 90910
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0124_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 91000) .AND. (LHAINPUT .LE. 92000)) THEN
     	LHASET = 91000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF10_1000.LHpdf'
     ELSEIF((LHAINPUT .GE. 93000) .AND. (LHAINPUT .LE. 94000))THEN
     	LHASET = 93000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF10_1000.LHgrid'		  
     ELSEIF((LHAINPUT .GE. 95000) .AND. (LHAINPUT .LE. 96000))THEN
     	LHASET = 95000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF12_1000.LHgrid'		  
     ELSEIF((LHAINPUT .GE. 97000) .AND. (LHAINPUT .LE. 98000))THEN
     	LHASET = 97000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_1000.LHgrid'		  
     ELSE
     	WRITE(LHAPRINT,5150)  LHASET
     	STOP
        ENDIF          
!...NNPDF  Family second tranche
  ELSEIF((LHAINPUT .GE. 190000) .AND. (LHAINPUT .LE. 210100)) THEN
     XMIN = 1.0D-9 
     Q2MIN = 2.0D0 
     Q2MAX = 1.0D08  
     IF((LHAINPUT .GE. 190000) .AND. (LHAINPUT .LE. 190100))THEN
     	LHASET = 190000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0114_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 190200) .AND. (LHAINPUT .LE. 190300))THEN
     	LHASET = 190200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0115_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 190400) .AND. (LHAINPUT .LE. 190500))THEN
     	LHASET = 190400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0116_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 190600) .AND. (LHAINPUT .LE. 190700))THEN
     	LHASET = 190600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0117_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 190800) .AND. (LHAINPUT .LE. 190900))THEN
     	LHASET = 190800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0118_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 191000) .AND. (LHAINPUT .LE. 191100))THEN
     	LHASET = 191000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0120_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 191200) .AND. (LHAINPUT .LE. 191300))THEN
     	LHASET = 191200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0121_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 191400) .AND. (LHAINPUT .LE. 191500))THEN
     	LHASET = 191400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0122_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 191600) .AND. (LHAINPUT .LE. 191700))THEN
     	LHASET = 191600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0123_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 191800) .AND. (LHAINPUT .LE. 191900))THEN
     	LHASET = 191800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_as_0124_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 192000) .AND. (LHAINPUT .LE. 192100))THEN
     	LHASET = 192000 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_heraold_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 192200) .AND. (LHAINPUT .LE. 192300))THEN
     	LHASET = 192200 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_dis_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 192400) .AND. (LHAINPUT .LE. 192500))THEN
     	LHASET = 192400 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_dis+dy_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 192600) .AND. (LHAINPUT .LE. 192700))THEN
     	LHASET = 192600 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF20_dis+jet_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 192800) .AND. (LHAINPUT .LE. 192900))THEN
     	LHASET = 192800 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 192901) .AND. (LHAINPUT .LE. 193901))THEN
     	LHASET = 192901 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_1000.LHgrid'
     ELSEIF((LHAINPUT .GE. 194000) .AND. (LHAINPUT .LE. 194100))THEN
     	LHASET = 194000 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_as_0114_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 194200) .AND. (LHAINPUT .LE. 194300))THEN
     	LHASET = 194200 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_as_0115_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 194400) .AND. (LHAINPUT .LE. 194500))THEN
     	LHASET = 194400 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_as_0116_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 194600) .AND. (LHAINPUT .LE. 194700))THEN
     	LHASET = 194600 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_as_0117_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 194800) .AND. (LHAINPUT .LE. 194900))THEN
     	LHASET = 194800 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_as_0118_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 195000) .AND. (LHAINPUT .LE. 195100))THEN
     	LHASET = 195000 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_as_0120_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 195200) .AND. (LHAINPUT .LE. 195300))THEN
     	LHASET = 195200 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_as_0121_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 195400) .AND. (LHAINPUT .LE. 195500))THEN
     	LHASET = 195400 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_as_0122_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 195600) .AND. (LHAINPUT .LE. 195700))THEN
     	LHASET = 195600 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_as_0123_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 195800) .AND. (LHAINPUT .LE. 195900))THEN
     	LHASET = 195800 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_as_0124_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 196000) .AND. (LHAINPUT .LE. 196100))THEN
     	LHASET = 196000 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_mc_15_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 196200) .AND. (LHAINPUT .LE. 196300))THEN
     	LHASET = 196200 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_mc_16_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 196400) .AND. (LHAINPUT .LE. 196500))THEN
     	LHASET = 196400 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_mc_17_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 196600) .AND. (LHAINPUT .LE. 196700))THEN
     	LHASET = 196600 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_mb_425_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 196800) .AND. (LHAINPUT .LE. 196900))THEN
     	LHASET = 196800 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_mb_45_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 197000) .AND. (LHAINPUT .LE. 197100))THEN
     	LHASET = 197000 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_mb_50_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 197200) .AND. (LHAINPUT .LE. 197300))THEN
     	LHASET = 197200 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_mb_525_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 197400) .AND. (LHAINPUT .LE. 197500))THEN
     	LHASET = 197400 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_FFN_NF3_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 197600) .AND. (LHAINPUT .LE. 197700))THEN
     	LHASET = 197600 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_FFN_NF4_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 197800) .AND. (LHAINPUT .LE. 197900))THEN
     	LHASET = 197800 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_FFN_NF5_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 198000) .AND. (LHAINPUT .LE. 198100))THEN
     	LHASET = 198000 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_dis_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 198101) .AND. (LHAINPUT .LE. 199101))THEN
     	LHASET = 198101 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_dis_1000.LHgrid'
     ELSEIF((LHAINPUT .GE. 199200) .AND. (LHAINPUT .LE. 199300))THEN
     	LHASET = 199200 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_dis+dy_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 199400) .AND. (LHAINPUT .LE. 199500))THEN
     	LHASET = 199400 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_dis+jet_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 200000) .AND. (LHAINPUT .LE.200100))THEN
     	LHASET = 200000 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 200200) .AND. (LHAINPUT .LE.200300))THEN
     	LHASET = 200200 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_lo_as_0119_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 200400) .AND. (LHAINPUT .LE.200500))THEN
     	LHASET = 200400 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_lo_as_0130_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 200600) .AND. (LHAINPUT .LE.200700))THEN
     	LHASET = 200600 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_lostar_as_0119_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 200800) .AND. (LHAINPUT .LE.200900))THEN
     	LHASET = 200800 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_lostar_as_0130_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 201000) .AND. (LHAINPUT .LE.202000))THEN
     	LHASET = 201000 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_1000.LHgrid'
     ELSEIF((LHAINPUT .GE. 203000) .AND. (LHAINPUT .LE.203100))THEN
     	LHASET = 203000 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_as_0114_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 203200) .AND. (LHAINPUT .LE.203300))THEN
     	LHASET = 203200 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_as_0115_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 203400) .AND. (LHAINPUT .LE.203500))THEN
     	LHASET = 203400 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_as_0116_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 203600) .AND. (LHAINPUT .LE.203700))THEN
     	LHASET = 203600 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_as_0117_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 203800) .AND. (LHAINPUT .LE.203900))THEN
     	LHASET = 203800 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_as_0118_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 204000) .AND. (LHAINPUT .LE.204100))THEN
     	LHASET = 204000 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_as_0120_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 204200) .AND. (LHAINPUT .LE.204300))THEN
     	LHASET = 204200 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_as_0121_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 204400) .AND. (LHAINPUT .LE.204500))THEN
     	LHASET = 204400 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_as_0122_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 204600) .AND. (LHAINPUT .LE.204700))THEN
     	LHASET = 204600 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_as_0123_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 204800) .AND. (LHAINPUT .LE.204900))THEN
     	LHASET = 204800 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_as_0124_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 205000) .AND. (LHAINPUT .LE.205100))THEN
     	LHASET = 205000 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_mc_15_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 205200) .AND. (LHAINPUT .LE.205300))THEN
     	LHASET = 205200 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_mc_16_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 205400) .AND. (LHAINPUT .LE.205500))THEN
     	LHASET = 205400 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_mc_17_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 206000) .AND. (LHAINPUT .LE.206100))THEN
     	LHASET = 206000 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_mb_425_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 206200) .AND. (LHAINPUT .LE.206300))THEN
     	LHASET = 206200 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_mb_45_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 206400) .AND. (LHAINPUT .LE.206500))THEN
     	LHASET = 206400 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_mb_50_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 206600) .AND. (LHAINPUT .LE.206700))THEN
     	LHASET = 206600 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_mb_525_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 207000) .AND. (LHAINPUT .LE.207100))THEN
     	LHASET = 207000 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_dis_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 207200) .AND. (LHAINPUT .LE.207300))THEN
     	LHASET = 207200 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_dis+dy_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 207400) .AND. (LHAINPUT .LE.207500))THEN
     	LHASET = 207400 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_heraonly_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 207600) .AND. (LHAINPUT .LE.207700))THEN
     	LHASET = 207600 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_collider_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 207800) .AND. (LHAINPUT .LE.207900))THEN
     	LHASET = 207800 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF21_nnlo_nf5_100.LHgrid'
     ELSEIF((LHAINPUT .GE. 210000) .AND. (LHAINPUT .LE.210100))THEN
     	LHASET = 210000 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF22_nlo_100.LHgrid'
     ELSE
     	WRITE(LHAPRINT,5150)  LHASET
     	STOP
        ENDIF          
!...NNPDF 2.3  Family 
  ELSEIF((LHAINPUT .GE. 229000) .AND. (LHAINPUT .LE. 246700)) THEN
     XMIN = 1.0D-9 
     Q2MIN = 2.0D0 
     Q2MAX = 1.0D08  

     IF((LHAINPUT .GE. 229000) .AND. (LHAINPUT .LE. 229100))THEN
     	LHASET = 229000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0114.LHgrid'
     ELSEIF((LHAINPUT .GE. 229200) .AND. (LHAINPUT .LE. 229300))THEN
     	LHASET = 229200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0115.LHgrid'
     ELSEIF((LHAINPUT .GE. 229400) .AND. (LHAINPUT .LE. 229500))THEN
     	LHASET = 229400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0116.LHgrid'
     ELSEIF((LHAINPUT .GE. 229600) .AND. (LHAINPUT .LE. 229700))THEN
     	LHASET = 229600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0117.LHgrid'
     ELSEIF((LHAINPUT .GE. 229800) .AND. (LHAINPUT .LE. 229900))THEN
     	LHASET = 229800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0118.LHgrid'
     ELSEIF((LHAINPUT .GE. 230000) .AND. (LHAINPUT .LE. 230100))THEN
     	LHASET = 230000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0119.LHgrid'
     ELSEIF((LHAINPUT .GE. 230200) .AND. (LHAINPUT .LE. 230300))THEN
     	LHASET = 230200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0120.LHgrid'
     ELSEIF((LHAINPUT .GE. 230400) .AND. (LHAINPUT .LE. 230500))THEN
     	LHASET = 230400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0121.LHgrid'
     ELSEIF((LHAINPUT .GE. 230600) .AND. (LHAINPUT .LE. 230700))THEN
     	LHASET = 230600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0122.LHgrid'
     ELSEIF((LHAINPUT .GE. 230800) .AND. (LHAINPUT .LE. 230900))THEN
     	LHASET = 230800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0123.LHgrid'
     ELSEIF((LHAINPUT .GE. 231000) .AND. (LHAINPUT .LE. 231100))THEN
     	LHASET = 231000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0124.LHgrid'

     ELSEIF((LHAINPUT .GE. 231200) .AND. (LHAINPUT .LE. 231300))THEN
     	LHASET = 231200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_as_0114.LHgrid'
     ELSEIF((LHAINPUT .GE. 231400) .AND. (LHAINPUT .LE. 231500))THEN
     	LHASET = 231400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_as_0115.LHgrid'
     ELSEIF((LHAINPUT .GE. 231600) .AND. (LHAINPUT .LE. 231700))THEN
     	LHASET = 231600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_as_0116.LHgrid'
     ELSEIF((LHAINPUT .GE. 231800) .AND. (LHAINPUT .LE. 231900))THEN
     	LHASET = 231800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_as_0117.LHgrid'
     ELSEIF((LHAINPUT .GE. 232000) .AND. (LHAINPUT .LE. 232100))THEN
     	LHASET = 232000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_as_0118.LHgrid'
     ELSEIF((LHAINPUT .GE. 232200) .AND. (LHAINPUT .LE. 232300))THEN
     	LHASET = 232200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_as_0119.LHgrid'
     ELSEIF((LHAINPUT .GE. 232400) .AND. (LHAINPUT .LE. 232500))THEN
     	LHASET = 232400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_as_0120.LHgrid'
     ELSEIF((LHAINPUT .GE. 232600) .AND. (LHAINPUT .LE. 232700))THEN
     	LHASET = 232600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_as_0121.LHgrid'
     ELSEIF((LHAINPUT .GE. 232800) .AND. (LHAINPUT .LE. 232900))THEN
     	LHASET = 232800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_as_0122.LHgrid'
     ELSEIF((LHAINPUT .GE. 233000) .AND. (LHAINPUT .LE. 233100))THEN
     	LHASET = 233000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_as_0123.LHgrid'
     ELSEIF((LHAINPUT .GE. 233200) .AND. (LHAINPUT .LE. 233300))THEN
     	LHASET = 233200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_as_0124.LHgrid'

     ELSEIF((LHAINPUT .GE. 233400) .AND. (LHAINPUT .LE. 233500))THEN
     	LHASET = 233400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_noLHC_as_0116.LHgrid'
     ELSEIF((LHAINPUT .GE. 233600) .AND. (LHAINPUT .LE. 233700))THEN
     	LHASET = 233600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_noLHC_as_0117.LHgrid'
     ELSEIF((LHAINPUT .GE. 233800) .AND. (LHAINPUT .LE. 233900))THEN
     	LHASET = 233800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_noLHC_as_0118.LHgrid'
     ELSEIF((LHAINPUT .GE. 234000) .AND. (LHAINPUT .LE. 234100))THEN
     	LHASET = 234000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_noLHC_as_0119.LHgrid'
     ELSEIF((LHAINPUT .GE. 234200) .AND. (LHAINPUT .LE. 234300))THEN
     	LHASET = 234200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_noLHC_as_0120.LHgrid'

     ELSEIF((LHAINPUT .GE. 234400) .AND. (LHAINPUT .LE. 234500))THEN
     	LHASET = 234400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_noLHC_as_0116.LHgrid'
     ELSEIF((LHAINPUT .GE. 234600) .AND. (LHAINPUT .LE. 234700))THEN
     	LHASET = 234600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_noLHC_as_0117.LHgrid'
     ELSEIF((LHAINPUT .GE. 234800) .AND. (LHAINPUT .LE. 234900))THEN
     	LHASET = 234800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_noLHC_as_0118.LHgrid'
     ELSEIF((LHAINPUT .GE. 235000) .AND. (LHAINPUT .LE. 235100))THEN
     	LHASET = 235000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_noLHC_as_0119.LHgrid'
     ELSEIF((LHAINPUT .GE. 235200) .AND. (LHAINPUT .LE. 235300))THEN
     	LHASET = 235200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_noLHC_as_0120.LHgrid'

     ELSEIF((LHAINPUT .GE. 235400) .AND. (LHAINPUT .LE. 235500))THEN
     	LHASET = 235400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_collider_as_0116.LHgrid'
     ELSEIF((LHAINPUT .GE. 235600) .AND. (LHAINPUT .LE. 235700))THEN
     	LHASET = 235600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_collider_as_0117.LHgrid'
     ELSEIF((LHAINPUT .GE. 235800) .AND. (LHAINPUT .LE. 235900))THEN
     	LHASET = 235800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_collider_as_0118.LHgrid'
     ELSEIF((LHAINPUT .GE. 236000) .AND. (LHAINPUT .LE. 236100))THEN
     	LHASET = 236000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_collider_as_0119.LHgrid'
     ELSEIF((LHAINPUT .GE. 236200) .AND. (LHAINPUT .LE. 236300))THEN
     	LHASET = 236200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_collider_as_0120.LHgrid'

     ELSEIF((LHAINPUT .GE. 236400) .AND. (LHAINPUT .LE. 236500))THEN
     	LHASET = 236400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_collider_as_0116.LHgrid'
     ELSEIF((LHAINPUT .GE. 236600) .AND. (LHAINPUT .LE. 236700))THEN
     	LHASET = 236600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_collider_as_0117.LHgrid'
     ELSEIF((LHAINPUT .GE. 236800) .AND. (LHAINPUT .LE. 236900))THEN
     	LHASET = 236800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_collider_as_0118.LHgrid'
     ELSEIF((LHAINPUT .GE. 237000) .AND. (LHAINPUT .LE. 237100))THEN
     	LHASET = 237000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_collider_as_0119.LHgrid'
     ELSEIF((LHAINPUT .GE. 237200) .AND. (LHAINPUT .LE. 237300))THEN
     	LHASET = 237200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_collider_as_0120.LHgrid'

     ELSEIF((LHAINPUT .GE. 237400) .AND. (LHAINPUT .LE. 237500))THEN
     	LHASET = 237400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF4_as_0116.LHgrid'
     ELSEIF((LHAINPUT .GE. 237600) .AND. (LHAINPUT .LE. 237700))THEN
     	LHASET = 237600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF4_as_0117.LHgrid'
     ELSEIF((LHAINPUT .GE. 237800) .AND. (LHAINPUT .LE. 237900))THEN
     	LHASET = 237800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF4_as_0118.LHgrid'
     ELSEIF((LHAINPUT .GE. 238000) .AND. (LHAINPUT .LE. 238100))THEN
     	LHASET = 238000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF4_as_0119.LHgrid'
     ELSEIF((LHAINPUT .GE. 238200) .AND. (LHAINPUT .LE. 238300))THEN
     	LHASET = 238200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF4_as_0120.LHgrid'

     ELSEIF((LHAINPUT .GE. 238400) .AND. (LHAINPUT .LE. 238500))THEN
     	LHASET = 238400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_FFN_NF4_as_0116.LHgrid'
     ELSEIF((LHAINPUT .GE. 238600) .AND. (LHAINPUT .LE. 238700))THEN
     	LHASET = 238600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_FFN_NF4_as_0117.LHgrid'
     ELSEIF((LHAINPUT .GE. 238800) .AND. (LHAINPUT .LE. 238900))THEN
     	LHASET = 238800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_FFN_NF4_as_0118.LHgrid'
     ELSEIF((LHAINPUT .GE. 239000) .AND. (LHAINPUT .LE. 239100))THEN
     	LHASET = 239000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_FFN_NF4_as_0119.LHgrid'
     ELSEIF((LHAINPUT .GE. 239200) .AND. (LHAINPUT .LE. 239300))THEN
     	LHASET = 239200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_FFN_NF4_as_0120.LHgrid'

     ELSEIF((LHAINPUT .GE. 239400) .AND. (LHAINPUT .LE. 239500))THEN
     	LHASET = 239400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF5_as_0116.LHgrid'
     ELSEIF((LHAINPUT .GE. 239600) .AND. (LHAINPUT .LE. 239700))THEN
     	LHASET = 239600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF5_as_0117.LHgrid'
     ELSEIF((LHAINPUT .GE. 239800) .AND. (LHAINPUT .LE. 239900))THEN
     	LHASET = 239800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF5_as_0118.LHgrid'
     ELSEIF((LHAINPUT .GE. 240000) .AND. (LHAINPUT .LE. 240100))THEN
     	LHASET = 240000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF5_as_0119.LHgrid'
     ELSEIF((LHAINPUT .GE. 240200) .AND. (LHAINPUT .LE. 240300))THEN
     	LHASET = 240200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF5_as_0120.LHgrid'

     ELSEIF((LHAINPUT .GE. 240400) .AND. (LHAINPUT .LE. 240500))THEN
     	LHASET = 240400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_FFN_NF5_as_0116.LHgrid'
     ELSEIF((LHAINPUT .GE. 240600) .AND. (LHAINPUT .LE. 240700))THEN
     	LHASET = 240600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_FFN_NF5_as_0117.LHgrid'
     ELSEIF((LHAINPUT .GE. 240800) .AND. (LHAINPUT .LE. 240900))THEN
     	LHASET = 240800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_FFN_NF5_as_0118.LHgrid'
     ELSEIF((LHAINPUT .GE. 241000) .AND. (LHAINPUT .LE. 241100))THEN
     	LHASET = 241000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_FFN_NF5_as_0119.LHgrid'
     ELSEIF((LHAINPUT .GE. 241200) .AND. (LHAINPUT .LE. 241300))THEN
     	LHASET = 241200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_FFN_NF5_as_0120.LHgrid'

     ELSEIF((LHAINPUT .GE. 241400) .AND. (LHAINPUT .LE. 241500))THEN
     	LHASET = 241400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0116_mc.LHgrid'
     ELSEIF((LHAINPUT .GE. 241600) .AND. (LHAINPUT .LE. 241700))THEN
     	LHASET = 241600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0117_mc.LHgrid'
     ELSEIF((LHAINPUT .GE. 241800) .AND. (LHAINPUT .LE. 241900))THEN
     	LHASET = 241800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0118_mc.LHgrid'
     ELSEIF((LHAINPUT .GE. 242000) .AND. (LHAINPUT .LE. 242100))THEN
     	LHASET = 242000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0119_mc.LHgrid'
     ELSEIF((LHAINPUT .GE. 242200) .AND. (LHAINPUT .LE. 242300))THEN
     	LHASET = 242200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0120_mc.LHgrid'

     ELSEIF((LHAINPUT .GE. 242400) .AND. (LHAINPUT .LE. 242500))THEN
     	LHASET = 242400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF4_as_0116_mc.LHgrid'
     ELSEIF((LHAINPUT .GE. 242600) .AND. (LHAINPUT .LE. 242700))THEN
     	LHASET = 242600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF4_as_0117_mc.LHgrid'
     ELSEIF((LHAINPUT .GE. 242800) .AND. (LHAINPUT .LE. 242900))THEN
     	LHASET = 242800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF4_as_0118_mc.LHgrid'
     ELSEIF((LHAINPUT .GE. 243000) .AND. (LHAINPUT .LE. 243100))THEN
     	LHASET = 243000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF4_as_0119_mc.LHgrid'
     ELSEIF((LHAINPUT .GE. 243200) .AND. (LHAINPUT .LE. 243300))THEN
     	LHASET = 243200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF4_as_0120_mc.LHgrid'

     ELSEIF((LHAINPUT .GE. 243400) .AND. (LHAINPUT .LE. 243500))THEN
     	LHASET = 243400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF5_as_0116_mc.LHgrid'
     ELSEIF((LHAINPUT .GE. 243600) .AND. (LHAINPUT .LE. 243700))THEN
     	LHASET = 243600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF5_as_0117_mc.LHgrid'
     ELSEIF((LHAINPUT .GE. 243800) .AND. (LHAINPUT .LE. 243900))THEN
     	LHASET = 243800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF5_as_0118_mc.LHgrid'
     ELSEIF((LHAINPUT .GE. 244000) .AND. (LHAINPUT .LE. 244100))THEN
     	LHASET = 244000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF5_as_0119_mc.LHgrid'
     ELSEIF((LHAINPUT .GE. 244200) .AND. (LHAINPUT .LE. 244300))THEN
     	LHASET = 244200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_FFN_NF5_as_0120_mc.LHgrid'

     ELSEIF((LHAINPUT .GE. 244400) .AND. (LHAINPUT .LE. 244500))THEN
     	LHASET = 244400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0117_qed.LHgrid'
     ELSEIF((LHAINPUT .GE. 244600) .AND. (LHAINPUT .LE. 244700))THEN
     	LHASET = 244600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0118_qed.LHgrid'
     ELSEIF((LHAINPUT .GE. 244800) .AND. (LHAINPUT .LE. 244900))THEN
     	LHASET = 244800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0119_qed.LHgrid'

     ELSEIF((LHAINPUT .GE. 245000) .AND. (LHAINPUT .LE. 245100))THEN
     	LHASET = 245000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0117_qed_neutron.LHgrid'
     ELSEIF((LHAINPUT .GE. 245200) .AND. (LHAINPUT .LE. 245300))THEN
     	LHASET = 245200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0118_qed_neutron.LHgrid'
     ELSEIF((LHAINPUT .GE. 245400) .AND. (LHAINPUT .LE. 245500))THEN
     	LHASET = 245400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nlo_as_0119_qed_neutron.LHgrid'

     ELSEIF((LHAINPUT .GE. 245600) .AND. (LHAINPUT .LE. 245700))THEN
     	LHASET = 245600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_as_0117_qed.LHgrid'
     ELSEIF((LHAINPUT .GE. 245800) .AND. (LHAINPUT .LE. 245900))THEN
     	LHASET = 245800
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_as_0118_qed.LHgrid'
     ELSEIF((LHAINPUT .GE. 246000) .AND. (LHAINPUT .LE. 246100))THEN
     	LHASET = 246000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_as_0119_qed.LHgrid'
     ELSEIF((LHAINPUT .GE. 246200) .AND. (LHAINPUT .LE. 246300))THEN
     	LHASET = 246200
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_as_0117_qed_neutron.LHgrid'
     ELSEIF((LHAINPUT .GE. 246400) .AND. (LHAINPUT .LE. 246500))THEN
     	LHASET = 246400
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_as_0118_qed_neutron.LHgrid'
     ELSEIF((LHAINPUT .GE. 246600) .AND. (LHAINPUT .LE. 246700))THEN
     	LHASET = 246600
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDF23_nnlo_as_0119_qed_neutron.LHgrid'

     ELSE
     	WRITE(LHAPRINT,5150)  LHASET
     	STOP
        ENDIF          
!...NNPDF Pol  Family 
  ELSEIF((LHAINPUT .GE. 250000) .AND. (LHAINPUT .LE. 250100)) THEN
     XMIN = 1.0D-9 
     Q2MIN = 2.0D0 
     Q2MAX = 1.0D08  
     IF((LHAINPUT .GE. 250000) .AND. (LHAINPUT .LE. 250100))THEN
     	LHASET = 250000
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/NNPDFpol10_100.LHgrid'
     ELSE
     	WRITE(LHAPRINT,5150)  LHASET
     	STOP
        ENDIF          
!...User defined  sets 
  ELSEIF((LHAINPUT .GE. 99002) .AND. (LHAINPUT .LE. 99004)) THEN
     XMIN = 1.0D-9 
     Q2MIN = 1.0D0 
     Q2MAX = 1.0D09  
     IF((LHAINPUT .GE.99002) .AND. (LHAINPUT .LE. 99002)) THEN
     	LHASET = 99002 
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/USERGRIDQ2.LHgrid'
     ELSEIF((LHAINPUT .GE. 99003) .AND. (LHAINPUT .LE. 99003))THEN
     	LHASET = 99003
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/USERGRIDQ3.LHgrid'
     ELSEIF((LHAINPUT .GE. 99004) .AND. (LHAINPUT .LE. 99004))THEN
     	LHASET = 99004
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/USERGRIDQ4.LHgrid'
     ELSE
     	WRITE(LHAPRINT,5150)  LHASET
     	STOP
     ENDIF          
!...Nuclear PDFs HKN sets 
  ELSEIF((LHAINPUT .GE. 100050) .AND. (LHAINPUT .LE. 100169)) THEN
     XMIN = 1.0D-9 
     Q2MIN = 1.0D0 
     Q2MAX = 1.0D08  
     IF((LHAINPUT .GE.100050) .AND. (LHAINPUT .LE. 100069)) THEN
     	LHASET = 100050
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/HKNlo.LHgrid'
     ELSEIF((LHAINPUT .GE. 100150) .AND. (LHAINPUT .LE. 100169))THEN
     	LHASET =100150
     	LHANAME=LHAPATH(1:LHAPATHLEN)//'/HKNnlo.LHgrid'
     ELSE
     	WRITE(LHAPRINT,5150)  LHASET
     	STOP
        ENDIF          
     ! 
     ! Pions
     ! 
     ! OW-PI Family
  elseif((lhainput .ge. 210) .and. (lhainput .le. 212)) then
     q2min = 4.0d0
     q2max = 2.0d03
     xmin = 5.0d-03
     xmax = 0.9998d0
     if((lhainput .ge. 210) .and. (lhainput .le. 212)) then
        lhaset = 210
        lhaname=lhapath(1:lhapathlen)//'/OWPI.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! SMRS-PI Family
  elseif((lhainput .ge. 230) .and. (lhainput .le. 233)) then
     q2min = 5.0d0
     q2max = 1.31d06
     xmin = 1.0d-05
     xmax = 0.9998d0
     if((lhainput .ge. 230) .and. (lhainput .le. 233)) then
        lhaset = 230
        lhaname=lhapath(1:lhapathlen)//'/SMRSPI.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! GRV-PI Family
  elseif((lhainput .ge. 250) .and. (lhainput .le. 252)) then
     q2max = 1.00d06
     xmin = 1.0d-05
     xmax = 0.9998d0
     if((lhainput .ge. 250) .and. (lhainput .le. 251)) then
        q2min = 3.0d-1
        lhaset = 250
        lhaname=lhapath(1:lhapathlen)//'/GRVPI1.LHgrid'
     elseif((lhainput .ge. 252) .and. (lhainput .le. 252)) then
        q2min = 2.5d-1
        lhaset = 252
        lhaname=lhapath(1:lhapathlen)//'/GRVPI0.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! ABFKW-PI Family
  elseif((lhainput .ge. 260) .and. (lhainput .le. 263)) then
     q2min = 2.0d0
     q2max = 1.00d08
     xmin = 1.0d-03
     xmax = 0.9998d0
     if((lhainput .ge. 260) .and. (lhainput .le. 263)) then
        lhaset = 260
        lhaname=lhapath(1:lhapathlen)//'/ABFKWPI.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! 
     ! photons
     ! 
     ! DO-G Family
  elseif((lhainput .ge. 310) .and. (lhainput .le. 312)) then
     q2min = 1.0d01
     q2max = 1.00d04
     xmin = 1.0d-05
     xmax = 0.9d0
     if((lhainput .ge. 310) .and. (lhainput .le. 311)) then
        lhaset = 310
        lhaname=lhapath(1:lhapathlen)//'/DOG0.LHgrid'
     elseif((lhainput .ge. 312) .and. (lhainput .le. 312)) then
        lhaset = 312
        lhaname=lhapath(1:lhapathlen)//'/DOG1.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! DG-G Family
  elseif((lhainput .ge. 320) .and. (lhainput .le. 324)) then
     xmin = 1.0d-05
     xmax = 0.9998d0
     lhaset = 320
     if((lhainput .ge. 320) .and. (lhainput .le. 321)) then
        q2min = 1.0d0
        q2max = 1.0d04
        ! lhaset = 320
        lhaname=lhapath(1:lhapathlen)//'/DGG.LHgrid'
     elseif((lhainput .ge. 322) .and. (lhainput .le. 322)) then
        q2min = 1.0d0
        q2max = 5.0d01
        ! lhaset = 322
        lhaname=lhapath(1:lhapathlen)//'/DGG.LHgrid'
     elseif((lhainput .ge. 323) .and. (lhainput .le. 323)) then
        q2min = 2.0d1
        q2max = 5.0d02
        ! lhaset = 323
        lhaname=lhapath(1:lhapathlen)//'/DGG.LHgrid'
     elseif((lhainput .ge. 324) .and. (lhainput .le. 324)) then
        q2min = 2.0d2
        q2max = 1.0d04
        ! lhaset = 324
        lhaname=lhapath(1:lhapathlen)//'/DGG.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! LAC/GAL-G Family
  elseif((lhainput .ge. 330) .and. (lhainput .le. 334)) then
     q2min = 4.0d00
     q2max = 1.0d05
     xmin = 1.0d-04
     xmax = 0.9998d0
     lhaset = 330
     if((lhainput .ge. 330) .and. (lhainput .le. 332)) then
        ! lhaset = 330
        lhaname=lhapath(1:lhapathlen)//'/LACG.LHgrid'
     elseif((lhainput .ge. 333) .and. (lhainput .le. 333)) then
        q2min = 1.0d00
        ! lhaset = 333
        lhaname=lhapath(1:lhapathlen)//'/LACG.LHgrid'
     elseif((lhainput .ge. 334) .and. (lhainput .le. 334)) then
        q2min = 4.0d00
        ! lhaset = 334
        lhaname=lhapath(1:lhapathlen)//'/LACG.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! GSG/GSG96-G Family
  elseif((lhainput .ge. 340) .and. (lhainput .le. 345)) then
     q2min = 5.3d00
     q2max = 1.0d08
     xmin = 5.0d-04
     xmax = 0.9998d0
     if((lhainput .ge. 340) .and. (lhainput .le. 341)) then
        lhaset = 340
        lhaname=lhapath(1:lhapathlen)//'/GSG1.LHgrid'
     elseif((lhainput .ge. 342) .and. (lhainput .le. 343)) then
        lhaset = 341
        lhaname=lhapath(1:lhapathlen)//'/GSG0.LHgrid'
     elseif((lhainput .ge. 344) .and. (lhainput .le. 344)) then
        lhaset = 344
        lhaname=lhapath(1:lhapathlen)//'/GSG961.LHgrid'
     elseif((lhainput .ge. 345) .and. (lhainput .le. 345)) then
        lhaset = 345
        lhaname=lhapath(1:lhapathlen)//'/GSG960.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! GRV-G Family
  elseif((lhainput .ge. 350) .and. (lhainput .le. 354)) then
     q2min = 3.0d-1
     q2max = 1.0d06
     xmin = 1.0d-05
     xmax = 0.9998d0
     if((lhainput .ge. 350) .and. (lhainput .le. 352)) then
        lhaset = 350
        lhaname=lhapath(1:lhapathlen)//'/GRVG1.LHgrid'
     elseif((lhainput .ge. 353) .and. (lhainput .le. 353)) then
        q2min = 2.5d-1
        lhaset = 352
        lhaname=lhapath(1:lhapathlen)//'/GRVG0.LHgrid'
     elseif((lhainput .ge. 354) .and. (lhainput .le. 354)) then
        q2min = 6.0d-1
        q2max = 5.0d04
        lhaset = 352
        lhaname=lhapath(1:lhapathlen)//'/GRVG0.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! ACFGP-G Family
  elseif((lhainput .ge. 360) .and. (lhainput .le. 363)) then
     q2min = 2.0d00
     q2max = 5.5d05
     xmin = 1.37d-03
     xmax = 0.9998d0
     if((lhainput .ge. 360) .and. (lhainput .le. 363)) then
        lhaset = 360
        lhaname=lhapath(1:lhapathlen)//'/ACFGPG.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! WHIT-G Family
  elseif((lhainput .ge. 380) .and. (lhainput .le. 386)) then
     q2min = 4.0d00
     q2max = 2.5d03
     xmin = 1.0d-03
     xmax = 0.9998d0
     if((lhainput .ge. 380) .and. (lhainput .le. 386)) then
        lhaset = 380
        lhaname=lhapath(1:lhapathlen)//'/WHITG.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! SAS-G Family
  elseif ((lhainput .ge. 390) .and. (lhainput .le. 398)) then
     q2max = 5.0d04
     xmin = 1.0d-05
     xmax = 0.9998d0
     lhaset = 390
     if ((lhainput .ge. 390) .and. (lhainput .le. 392)) then
        q2min = 3.6d-1
        ! lhaset = 390
        lhaname=lhapath(1:lhapathlen)//'/SASG.LHgrid'
     elseif((lhainput .ge. 393) .and. (lhainput .le. 394)) then
        q2min = 4.0d00
        !           lhaset = 393
        lhaname=lhapath(1:lhapathlen)//'/SASG.LHgrid'
     elseif((lhainput .ge. 395) .and. (lhainput .le. 396)) then
        q2min = 3.6d-1
        !           lhaset = 395
        lhaname=lhapath(1:lhapathlen)//'/SASG.LHgrid'
     elseif((lhainput .ge. 397) .and. (lhainput .le. 398)) then
        q2min = 4.0d00
        !           lhaset = 397
        lhaname=lhapath(1:lhapathlen)//'/SASG.LHgrid'
     else
        write(lhaprint,5150)  lhaset
        stop
     endif
     ! Unknown Family ?! Giving up
  else
     write(lhaprint,5150)  lhaset
     stop
  endif

  lhamemb=lhainput-lhaset
  ! Now work out if we have already called this set/member
  iset = 0
  do j=1,nsets
     if (lhaname.eq.lhanames(j).and. &
          lhamemb.eq.lhamembers(j)) then
        iset = j
     endif
  enddo
  if (iset.eq.0) then
     nsets=nsets+1
     if (nsets.gt.nmxset) then
        if (LHASILENT.ne.1) then
           print *, "WARNING: too many sets initialised"
           print *,"overwriting from set 1 again"
        endif
        nsets = 1
        ! stop
     endif
     iset=nsets
     lhanames(iset)=lhaname
     lhanumbers(iset)=lhainput
     lhamembers(iset)=lhamemb
     xxmin(iset)=xmin
     xxmax(iset)=xmax
     qq2min(iset)=q2min
     qq2max(iset)=q2max
     call initpdfsetm(iset,lhaname)
     call numberpdfm(iset,lhaallmem)
     if(lhasilent .ne. 1) then
        write(lhaprint,5151)
        write(lhaprint,5152) lhaname
        write(lhaprint,5153) lhaallmem
        write(lhaprint,5154)
     endif
     if ((lhamemb.lt.0) .or. (lhamemb.gt.lhaallmem)) then
        write(lhaprint,5155)  lhamemb
        write(lhaprint,5156)  lhaallmem
        stop
     endif

     ! print *,'calling initpdf',lhamemb 
     ! print *,'calling initpdfm ',iset,lhaname,lhamemb
     ! print *,'LHAGLUE .... initializing set,member ',iset,lhamemb
     call initpdfm(iset,lhamemb)
  endif
  !  the rest is done every time pdfset is called
  !print *,'setting nset to:',iset
  call setnset(iset)
  call setnmem(iset,lhamemb)
  xmin = xxmin(iset)
  xmax = xxmax(iset)
  q2min=qq2min(iset)
  q2max=qq2max(iset)
  call GetLam4M(iset,LHAMEMB,qcdl4)
  call GetLam5M(iset,LHAMEMB,qcdl5)

  QMZ = 91.1876D0
  alphasLHA = alphasPDFM(iset,QMZ)
  if(lhasilent .ne. 1) write(lhaprint,5158) alphasLHA

  if(lhaparm(17).EQ.'LHAPDF') then
     nptypepdfl = 1      ! Proton PDFs
     nflpdfl = 4
     qcdlha4 = qcdl4
     qcdlha5 = qcdl5
     if (LHASILENT .NE. 1) write(lhaprint,5159) qcdl4, qcdl5
  else
     nptypepdfl = 1      ! Proton PDFs
     nflpdfl = 4
     alambda = 0.192d0
     qcdlha4 = alambda
     qcdlha5 = alambda
     if (parm(1).EQ.'NPTYPE') then        !  PYTHIA
        qcdl4 = alambda
        qcdl5 = alambda
     endif
  endif

  ! Formats for initialization information.
5150 format(1X,'WRONG LHAPDF set number =',I12,' given! STOP EXE!')
5151 format(1X,'==============================================')
5152 format(1X,'PDFset name ',A80)
5153 format(1X,'with ',I10,' members')
5154 format(1X,'====  initialized. ===========================')
5155 format(1X,'LHAPDF problem => YOU asked for member = ',I10)
5156 format(1X,'Valid range is: 0 - ',I10,' Execution stopped.')
  !5157 format(1X,'Number of flavors for PDF is:',I4)
5158 format(1X,'Strong coupling at Mz for PDF is:',F9.5)
5159 format(1X,'Will use for PYTHIA QCDL4, QCDL5:',2F9.5)

  return
end subroutine pdfset


!********************************************************************
! -- STRUCTA
! -- copy of PDFLIB to use the eks98 nuclear correction factors

subroutine structa(x,q,a,upv,dnv,usea,dsea,str,chm,bot,top,glu)
  implicit double precision (a-h,o-z)
  character*20 lparm
  call getlhaparm(15,lparm)
  if(lparm.eq.'EPS08') then
     call eps08(x,q,a,ruv,rdv,ru,rd,rs,rc,rb,rt,rg)
  else if (lparm(1:5).eq.'EPS09') then
    if (lparm.eq.'EPS09LO') then
      iorder=1
      ipset=1
    else if (lparm.eq.'EPS09NLO') then
      iorder=2
      ipset=1
    else if (lparm(1:8).eq.'EPS09LO,') then
      iorder=1
      read(lparm(9:),*)ipset
    else if (lparm(1:9).eq.'EPS09NLO,') then
      iorder=2
      read(lparm(10:),*)ipset
    else
      iorder=2
      ipset=1
    endif
    ia=a
    call eps09(iorder,ipset,ia,x,q,ruv,rdv,ru,rd,rs,rc,rb,rg)
    rt=1.0d0
  else
     call eks98(x,q,a,ruv,rdv,ru,rd,rs,rc,rb,rt,rg)
  endif
  call structm(x,q,upv,dnv,usea,dsea,str,chm,bot,top,glu)
  upv = ruv*upv
  dnv = rdv*dnv
  usea = ru*usea
  dsea = rd*dsea
  str = rs*str
  chm = rc*chm
  bot = rb*bot
  top = rt*top
  glu = rg*glu
  return
end subroutine structa


!*********************************************************************
! STRUCTM
! Gives parton distributions according to the LHAPDF interface.
! Two evolution codes used:
!   EVLCTEQ for CTEQ PDF sets
!   QCDNUM  for Other PDF sets
! 
! Author: Dimitri Bourilkov  bourilkov@mailaps.org
! 
! v4.0  21-Mar-2005  Photon/pion/new p PDFs, updated for LHAPDF v4
! v3.0  23-Jan-2004
! 
! interface to LHAPDF library
subroutine structm(dx,dq,upv,dnv,usea,dsea,str,chm,bot,top,glu)
  ! double precision and integer declarations.
  implicit double precision(a-h, o-z)
  implicit integer(i-n)
  ! Common blocks
  include 'commonlhapdf.inc'
  include 'commonlhasets.inc'
  include 'commonlhacontrol.inc'
  include 'commonlhaglsta.inc'
  ! commonblocks.
  common/pydat1/mstu(200),paru(200),mstj(200),parj(200)
  save /pydat1/
  ! interface to lhapdflib.
  double precision qcdlha4, qcdlha5
  integer nfllha
  common/lhapdfr/qcdlha4, qcdlha5, nfllha
  save /lhapdfr/
  integer lhaextrp
  common/lhapdfe/lhaextrp
  save /lhapdfe/
  ! interface to pdflib.
  common/w50513/xmin,xmax,q2min,q2max
  save /w50513/
  double precision xmin,xmax,q2min,q2max
  ! local variables
  double precision upv,dnv,usea,dsea,str,chm,bot,top,glu
  double precision dx,dq,x,q,f(-6:6),photon,gluino

  x = dx
  q = dq
  q2 = q**2
  ! statistics
  if(lhaparm(16).ne.'NOSTAT') then
     totnum = totnum+1.d0
     if(x .lt. xmin) xminnum = xminnum+1.d0
     if(x .gt. xmax) xmaxnum = xmaxnum+1.d0
     if(q2 .lt. q2min) q2minnum = q2minnum+1.d0
     if(q2 .gt. q2max) q2maxnum = q2maxnum+1.d0
  endif

  ! range of validity e.g. 10^-6 < x < 1, q2min < q^2 extended by
  ! freezing x*f(x,q2) at borders.
  if(lhaextrp .ne. 1) then    ! safe mode == "freeze"
     xin=max(xmin,min(xmax,x))
     q=sqrt(max(0d0,q2min,min(q2max,q2)))
  else                        ! adventurous mode == own risk !
     xin=x
  endif

  call getnset(iset)
  !print *,'calling evolvepdfm:',iset

  ! fix to allow STRUCTM to work for photon PDFs (Herwig does this)
  ! set P2 = 0.0d0 and IP2 = 0
  if(lhanumbers(iset).ge.300.and.lhanumbers(iset).le.399) then  
     p2 = 0.0d0
     ip2 = 0
     call evolvepdfpm(iset,xin,q,p2,ip2,f)
  else if (lhanumbers(iset).ge.20460.and.lhanumbers(iset).le.20462) then
     call evolvepdfphotonm(iset,xin,q,f,photon)
  else if (lhanumbers(iset).ge.10670.and.lhanumbers(iset).le.10677) then
     call evolvepdfgluinom(iset,xin,q,f,gluino)
  else
     call evolvepdfm(iset,xin,q,f)
  endif
  glu = f(0)
  dsea = f(-1)
  dnv = f(1) - dsea
  usea = f(-2)
  upv = f(2) - usea
  str = f(3)
  chm = f(4)
  bot = f(5)
  top = f(6)

  return
end subroutine structm


!*********************************************************************
! STRUCTP
! Gives parton distributions according to the LHAPDF interface.
! Used for photons.
! 
! v4.0  21-Mar-2005  Photon/pion/new p PDFs, updated for LHAPDF v4
! 
! Interface to LHAPDF library
subroutine structp(dx,dq2,p2,ip2,upv,dnv,usea,dsea,str,chm,bot,top,glu)
  ! Double precision and integer declarations.
  implicit double precision(a-h, o-z)
  implicit integer(i-n)
  ! Common blocks
  include 'parmsetup.inc'
  include 'commonlhapdf.inc'
  include 'commonlhacontrol.inc'
  include 'commonlhaglsta.inc'
  ! Commonblocks.
  common/pydat1/mstu(200),paru(200),mstj(200),parj(200)
  save /pydat1/
  ! Interface to LHAPDFLIB.
  double precision qcdlha4, qcdlha5
  integer nfllha
  common/lhapdfr/qcdlha4, qcdlha5, nfllha
  save /lhapdfr/
  integer lhaextrp
  common/lhapdfe/lhaextrp
  save /lhapdfe/
  ! Interface to PDFLIB.
  common/w50513/xmin,xmax,q2min,q2max
  save /w50513/
  double precision xmin,xmax,q2min,q2max
  ! Local variables
  double precision upv,dnv,usea,dsea,str,chm,bot,top,glu
  double precision dx,dq2,q2,x,q,f(-6:6)

  x = dx
  q2 = dq2
  ! Statistics
  if(lhaparm(16).ne.'NOSTAT') then
     totnup = totnup+1.d0
     if(x .lt. xmin) xminnup = xminnup+1.d0
     if(x .gt. xmax) xmaxnup = xmaxnup+1.d0
     if(q2 .lt. q2min) q2minnup = q2minnup+1.d0
     if(q2 .gt. q2max) q2maxnup = q2maxnup+1.d0
  endif

  ! Range of validity e.g. 10^-6 < x < 1, Q2MIN < Q^2 extended by
  ! freezing x*f(x,Q2) at borders.
  q = dsqrt(q2)
  if(lhaextrp .ne. 1) then    ! safe mode == "freeze"
     xin=max(xmin,min(xmax,x))
     q=sqrt(max(0d0,q2min,min(q2max,q2)))
  else                        ! adventurous mode == OWN RISK !
     xin=x
  endif
  call getnset(iset)
  call evolvepdfpm(iset,xin,q,p2,ip2,f)
  glu = f(0)
  dsea = f(-1)
  dnv = f(1) - dsea
  usea = f(-2)
  upv = f(2) - usea
  str = f(3)
  chm = f(4)
  bot = f(5)
  top = f(6)
  return
end subroutine structp


!*********************************************************************
! PDFSTA
! For statistics ON structure functions (under/over-flows)
! 
! Author: Dimitri Bourilkov  bourilkov@mailaps.org
! 
! 
! first introduced in v4.0  28-Apr-2005 
! 
subroutine pdfsta
  ! Double precision and integer declarations.
  implicit double precision(a-h, o-z)
  implicit integer(i-n)
  ! Common blocks
  include 'commonlhaglsta.inc'
  ! Interface to LHAPDFLIB.

  print *
  print *,'===== PDFSTA statistics for PDF under/over-flows ===='
  print *
  print *,'====== STRUCTM statistics for nucleon/pion PDFs ====='
  print *
  print *,'  total # of calls ',TOTNUM
  if(totnum .gt. 0.d0) then
     percbelow = 100.d0*xminnum/totnum
     percabove = 100.d0*xmaxnum/totnum
     print *,'  X  below PDF min ',xminnum,' or ',percbelow, ' %'
     print *,'  X  above PDF max ',xmaxnum,' or ',percabove, ' %'
     percbelow = 100.d0*q2minnum/totnum
     percabove = 100.d0*q2maxnum/totnum
     print *,'  Q2 below PDF min ',q2minnum,' or ',percbelow, ' %'
     print *,'  Q2 above PDF max ',q2maxnum,' or ',percabove, ' %'
  endif
  print *
  print *,'========= STRUCTP statistics for photon PDFs ========'
  print *
  print *,'  total # of calls ',totnup
  if(totnup .gt. 0.d0) then
     percbelow = 100.d0*xminnup/totnup
     percabove = 100.d0*xmaxnup/totnup
     print *,'  X  below PDF min ',xminnup,' or ',percbelow, ' %'
     print *,'  X  above PDF max ',xmaxnup,' or ',percabove, ' %'
     percbelow = 100.d0*q2minnup/totnup
     percabove = 100.d0*q2maxnup/totnup
     print *,'  Q2 below PDF min ',q2minnup,' or ',percbelow, ' %'
     print *,'  Q2 above PDF max ',q2maxnup,' or ',percabove, ' %'
  endif
  print *
  return
end subroutine pdfsta


subroutine pftopdg(dx,dscale,dxpdf)
  !include "pdf/expdp.inc"
  double precision dx,dscale,dupv,ddnv,dusea,ddsea,dstr,dchm,dbot,dtop,dgl,dxpdf(-6:6)
  ! Call STRUCTM in PDFLIB to get flavour content
  call structm(dx,dscale,dupv,ddnv,dusea,ddsea,dstr,dchm,dbot,dtop,dgl)
  ! Convert flavour convention of PDFLIB to PDG convention
  dxpdf(0) = dgl
  dxpdf(1) = ddnv + ddsea
  dxpdf(2) = dupv + dusea
  dxpdf(3) = dstr
  dxpdf(4) = dchm
  dxpdf(5) = dbot
  dxpdf(6) = dtop
  dxpdf(-1) = ddsea
  dxpdf(-2) = dusea
  dxpdf(-3) = dstr
  dxpdf(-4) = dchm
  dxpdf(-5) = dbot
  dxpdf(-6) = dtop
  return
end subroutine pftopdg
