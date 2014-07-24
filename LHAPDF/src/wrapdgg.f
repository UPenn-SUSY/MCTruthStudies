! -*- F90 -*-


      subroutine DGGevolvep(xin,qin,p2in,ip2in,pdf) 
      include 'parmsetup.inc' 
      real*8 xin,qin,q2in,p2in,pdf(-6:6),xval(45),qcdl4,qcdl5 
      real*8 upv,dnv,usea,dsea,str,chm,bot,top,glu 
      character*16 name(nmxset) 
      integer nmem(nmxset),ndef(nmxset),mmem 
      common/NAME/name,nmem,ndef,mmem 
      integer nset,iset 
                                                                        
      save 
      call getnset(iset) 
      call getnmem(iset,imem) 
                                                                        
      if(imem.eq.1.or.imem.eq.0) then 
        call DGPHO1(xin,qin,upv,dnv,usea,dsea,str,chm,bot,glu) 
                                                                        
      elseif(imem.eq.2) then 
        call DGPHO2(xin,qin,upv,dnv,usea,dsea,str,chm,bot,glu) 
                                                                        
      elseif(imem.eq.3) then 
        call DGPHO3(xin,qin,upv,dnv,usea,dsea,str,chm,bot,glu) 
                                                                        
      elseif(imem.eq.4) then 
        call DGPHO4(xin,qin,upv,dnv,usea,dsea,str,chm,bot,glu) 
                                                                        
      else 
        CONTINUE 
      endif 
                                                                        
      pdf(-6)= 0.0d0 
      pdf(6)= 0.0d0 
      pdf(-5)= bot 
      pdf(5 )= bot 
      pdf(-4)= chm 
      pdf(4 )= chm 
      pdf(-3)= str 
      pdf(3 )= str 
      pdf(-2)= usea 
      pdf(2 )= upv 
      pdf(-1)= dsea 
      pdf(1 )= dnv 
      pdf(0 )= glu 
                                                                        
      return 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      entry DGGread(nset) 
      read(1,*)nmem(nset),ndef(nset) 
      return 
!                                                                       
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      entry DGGalfa(alfas,qalfa) 
        call getnset(iset) 
        call getnmem(iset,imem) 
        call GetOrderAsM(iset,iord) 
        call Getlam4M(iset,imem,qcdl4) 
        call Getlam5M(iset,imem,qcdl5) 
        call aspdflib(alfas,Qalfa,iord,qcdl5) 
      return 
!                                                                       
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      entry DGGinit(Eorder,Q2fit) 
      return 
!                                                                       
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      entry DGGpdf(mem) 
      call getnset(iset) 
      call setnmem(iset,mem) 
!      imem = mem                                                       
      return 
!                                                                       
 1000 format(5e13.5) 
      END                                           
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
       SUBROUTINE DGPHO1(DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DGL) 
!********************************************************************   
!*                                                                  *   
!*    Parametrization of parton distribution functions              *   
!*    in the photon (LO analysis) - full  solution of AP eq.!       *   
!*                                                                  *   
!* authors:  M.Drees and K.Grassie (DG)                             *   
!*          /Z. Phys. C28 (1985) 451/                               *   
!*                                                                  *   
!* Prepared by:                                                     *   
!*             Krzysztof Charchula, DESY                            *   
!*             bitnet: F1PCHA@DHHDESY3                              *   
!*             decnet: 13313::CHARCHULA                             *   
!*                                                                  *   
!* Modified by:                                                     *   
!*             H. Plothow-Besch/CERN-PPE                            *   
!*                                                                  *   
!********************************************************************   
!                                                                       
      implicit real*8 (a-h,o-z) 
      double precision                                                  &
     &        A(3,4,3),AT(3),                                           &
     &        B(5,4,2,3),BT(5,2),XQPOM(2),E(2),                         &
     &        DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DGL                     
      PARAMETER (ALPEM=7.29927D-3, PI=3.141592D0) 
      PARAMETER (ALAM=0.4D0) 
!...comments                                                            
!...--------------------------------------------------                  
!...         nf=3 for   1< Q2 <32  GeV2                                 
!...         nf=4 for  32< Q2 <200 GeV2                                 
!...         nf=5 for 200< Q2 <1D4 GeV2                                 
!...--------------------------------------------------                  
!                                                                       
!...initialization of gluon parameters array for DG                     
        DATA (((A(I,J,K),I=1,3),J=1,4),K=1,3)/                          &
     &    -0.20700, -0.19870,  5.1190,                                  &
     &     0.61580,  0.62570, -0.2752,                                  &
     &     1.07400,  8.35200, -6.9930,                                  &
     &     0.00000,  5.02400,  2.2980,                                  &
     &     0.8926D-2,0.0509,  -0.2313,                                  &
     &     0.65940,  0.27740,  0.1382,                                  &
     &     0.47660, -0.39060,  6.5420,                                  &
     &     0.01975, -0.32120,  0.5162,                                  &
     &     0.03197, -0.618D-2,-0.1216,                                  &
     &     1.01800,  0.94760,  0.9047,                                  &
     &     0.24610, -0.60940,  2.6530,                                  &
     &     0.02707, -0.01067,  0.2003D-2/                               
!                                                                       
!...initialization of quark parameters array for DG                     
        DATA (((B(I,J,K,1),I=1,5),J=1,4),K=1,2)/                        &
     &     2.2850,   6.0730,  -0.4202,   -0.0808,  0.0553,              &
     &    -0.0153,  -0.8132,   0.0178,    0.6346,  1.1360,              &
     &     1.33D3, -41.310,    0.9216,    1.2080,  0.9512,              &
     &     4.2190,   3.1650,   0.1800,    0.2030,  0.0116,              &
     &    16.690,    0.1760,  -0.0208,   -0.0168, -0.1986,              &
     &    -0.7916,   0.0479,   0.3386D-2, 1.3530,  1.1000,              &
     &     1.0990D3, 1.0470,   4.8530,    1.4260,  1.1360,              &
     &     4.4280,   0.0250,   0.8404,    1.2390, -0.2779/              
        DATA (((B(I,J,K,2),I=1,5),J=1,4),K=1,2)/                        &
     &    -0.3711,  -0.1717,   0.08766,  -0.8915, -0.1816,              &
     &     1.0610,   0.7815,   0.02197,   0.2857,  0.5866,              &
     &     4.7580,   1.5350,   0.10960,   2.9730,  2.4210,              &
     &    -0.0150,   0.7067D-2,0.20400,   0.1185,  0.4059,              &
     &    -0.1207,  25.000,   -0.01230,  -0.0919,  0.02015,             &
     &     1.0710,  -1.6480,   1.16200,   0.7912,  0.9869,              &
     &     1.9770,  -0.01563,  0.48240,   0.6397, -0.07036,             &
     &    -0.8625D-2,6.4380,  -0.01100,   2.3270,  0.01694/             
        DATA (((B(I,J,K,3),I=1,5),J=1,4),K=1,2)/                        &
     &    15.8,      2.742,    0.02917,  -0.0342, -0.02302,             &
     &    -0.9464,  -0.7332,   0.04657,   0.7196,  0.9229,              &
     &    -0.5,      0.7148,   0.1785,    0.7338,  0.5873,              &
     &    -0.2118,   3.287,    0.04811,   0.08139,-0.79D-4,             &
     &     6.734,   59.88,    -0.3226D-2,-0.03321, 0.1059,              &
     &    -1.008,   -2.983,    0.8432,    0.9475,  0.6954,              &
     &    -0.08594,  4.48,     0.3616,   -0.3198, -0.6663,              &
     &     0.07625,  0.9686,   0.1383D-2, 0.02132, 0.3683/              
!                                                                       
!...specification of sets                                               
       Q2 = DQ*DQ 
         IF (Q2.LT.32.0D0) NFL=3 
         IF((Q2.GE.32.0D0).AND.(Q2.LT.200.0D0)) NFL=4 
         IF (Q2.GE.200.0D0) NFL=5 
!                                                                       
!...calculations                                                        
       ALAM2=ALAM**2 
       T=LOG(Q2/ALAM2) 
       LF=NFL-2 
!                                                                       
!...gluons                                                              
        DO 11 I=1,3 
          AT(I)=A(I,1,LF)*T**A(I,2,LF)+A(I,3,LF)*T**(-A(I,4,LF)) 
   11   CONTINUE 
        POMG=AT(1)*DX**AT(2)*(1.D0-DX)**AT(3) 
        DGL=POMG*ALPEM 
!                                                                       
!...quarks                                                              
        E(1)=1.0D0 
        IF(NFL.EQ.3) THEN 
          E(2)=9.0D0 
        ELSEIF(NFL.EQ.4) THEN 
          E(2)=10.0D0 
        ELSEIF(NFL.EQ.5) THEN 
          E(2)=55.0D0/6.0D0 
        ENDIF 
        DO 13 J=1,2 
          DO 15 I=1,5 
            BTP=B(I,1,J,LF)*T**B(I,2,J,LF) 
            BT(I,J)=BTP+B(I,3,J,LF)*T**(-B(I,4,J,LF)) 
   15     CONTINUE 
   13   CONTINUE 
!                                                                       
!...singlet & non-singlet combinations                                  
        DO 17 J=1,2 
          POM1=DX*(DX*DX+(1.D0-DX)**2)/(BT(1,J)-BT(2,J)*LOG(1.D0-DX)) 
          POM2=BT(3,J)*DX**BT(4,J)*(1.D0-DX)**BT(5,J) 
          XQPOM(J)=E(J)*POM1+POM2 
   17   CONTINUE 
!                                                                       
!...quarks flavours                                                     
        IF (NFL.EQ.3) THEN 
            DUB=ALPEM*1.D0/6.D0*(XQPOM(2)+9.D0*XQPOM(1)) 
            DDB=ALPEM*1.D0/6.D0*(XQPOM(2)-9.D0/2.D0*XQPOM(1)) 
            DSB=DDB 
            DCB=0.D0 
            DBB=0.D0 
        ELSEIF (NFL.EQ.4) THEN 
            DUB=ALPEM*1.D0/8.D0*(XQPOM(2)+6.D0*XQPOM(1)) 
            DCB=DUB 
            DDB=ALPEM*1.D0/8.D0*(XQPOM(2)-6.D0*XQPOM(1)) 
            DSB=DDB 
            DBB=0.D0 
        ELSEIF (NFL.EQ.5) THEN 
            DUB=ALPEM*1.D0/10.D0*(XQPOM(2)+15.D0/2.D0*XQPOM(1)) 
            DCB=DUB 
            DDB=ALPEM*1.D0/10.D0*(XQPOM(2)-5.D0*XQPOM(1)) 
            DSB=DDB 
            DBB=DDB 
        ENDIF 
      DUV=DUB 
      DDV=DDB 
!                                                                       
      RETURN 
      END                                           
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
       SUBROUTINE DGPHO2(DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DGL) 
!********************************************************************   
!*                                                                  *   
!*    Parametrization of parton distribution functions              *   
!*    in the photon (LO analysis) - full  solution of AP eq.!       *   
!*                                                                  *   
!* authors:  M.Drees and K.Grassie (DG)                             *   
!*          /Z. Phys. C28 (1985) 451/                               *   
!*                                                                  *   
!* Prepared by:                                                     *   
!*             Krzysztof Charchula, DESY                            *   
!*             bitnet: F1PCHA@DHHDESY3                              *   
!*             decnet: 13313::CHARCHULA                             *   
!*                                                                  *   
!* Modified by:                                                     *   
!*             H. Plothow-Besch/CERN-PPE                            *   
!*                                                                  *   
!********************************************************************   
!                                                                       
      implicit real*8 (a-h,o-z) 
      double precision                                                  &
     &        A(3,4,3),AT(3),                                           &
     &        B(5,4,2,3),BT(5,2),XQPOM(2),E(2),                         &
     &        DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DGL                     
      PARAMETER (ALPEM=7.29927D-3, PI=3.141592D0) 
      PARAMETER (ALAM=0.4D0) 
!...comments                                                            
!...--------------------------------------------------                  
!...        with nf=3 (valid for   1< Q2 <50  GeV2)                     
!...--------------------------------------------------                  
!                                                                       
!...initialization of gluon parameters array for DG                     
        DATA (((A(I,J,K),I=1,3),J=1,4),K=1,3)/                          &
     &    -0.20700, -0.19870,  5.1190,                                  &
     &     0.61580,  0.62570, -0.2752,                                  &
     &     1.07400,  8.35200, -6.9930,                                  &
     &     0.00000,  5.02400,  2.2980,                                  &
     &     0.8926D-2,0.0509,  -0.2313,                                  &
     &     0.65940,  0.27740,  0.1382,                                  &
     &     0.47660, -0.39060,  6.5420,                                  &
     &     0.01975, -0.32120,  0.5162,                                  &
     &     0.03197, -0.618D-2,-0.1216,                                  &
     &     1.01800,  0.94760,  0.9047,                                  &
     &     0.24610, -0.60940,  2.6530,                                  &
     &     0.02707, -0.01067,  0.2003D-2/                               
!                                                                       
!...initialization of quark parameters array for DG                     
        DATA (((B(I,J,K,1),I=1,5),J=1,4),K=1,2)/                        &
     &     2.2850,   6.0730,  -0.4202,   -0.0808,  0.0553,              &
     &    -0.0153,  -0.8132,   0.0178,    0.6346,  1.1360,              &
     &     1.33D3, -41.310,    0.9216,    1.2080,  0.9512,              &
     &     4.2190,   3.1650,   0.1800,    0.2030,  0.0116,              &
     &    16.690,    0.1760,  -0.0208,   -0.0168, -0.1986,              &
     &    -0.7916,   0.0479,   0.3386D-2, 1.3530,  1.1000,              &
     &     1.0990D3, 1.0470,   4.8530,    1.4260,  1.1360,              &
     &     4.4280,   0.0250,   0.8404,    1.2390, -0.2779/              
        DATA (((B(I,J,K,2),I=1,5),J=1,4),K=1,2)/                        &
     &    -0.3711,  -0.1717,   0.08766,  -0.8915, -0.1816,              &
     &     1.0610,   0.7815,   0.02197,   0.2857,  0.5866,              &
     &     4.7580,   1.5350,   0.10960,   2.9730,  2.4210,              &
     &    -0.0150,   0.7067D-2,0.20400,   0.1185,  0.4059,              &
     &    -0.1207,  25.000,   -0.01230,  -0.0919,  0.02015,             &
     &     1.0710,  -1.6480,   1.16200,   0.7912,  0.9869,              &
     &     1.9770,  -0.01563,  0.48240,   0.6397, -0.07036,             &
     &    -0.8625D-2,6.4380,  -0.01100,   2.3270,  0.01694/             
        DATA (((B(I,J,K,3),I=1,5),J=1,4),K=1,2)/                        &
     &    15.8,      2.742,    0.02917,  -0.0342, -0.02302,             &
     &    -0.9464,  -0.7332,   0.04657,   0.7196,  0.9229,              &
     &    -0.5,      0.7148,   0.1785,    0.7338,  0.5873,              &
     &    -0.2118,   3.287,    0.04811,   0.08139,-0.79D-4,             &
     &     6.734,   59.88,    -0.3226D-2,-0.03321, 0.1059,              &
     &    -1.008,   -2.983,    0.8432,    0.9475,  0.6954,              &
     &    -0.08594,  4.48,     0.3616,   -0.3198, -0.6663,              &
     &     0.07625,  0.9686,   0.1383D-2, 0.02132, 0.3683/              
!                                                                       
!...specification of sets                                               
         NFL=3 
!                                                                       
!...calculations                                                        
       Q2 = DQ*DQ 
       ALAM2=ALAM**2 
       T=LOG(Q2/ALAM2) 
       LF=NFL-2 
!                                                                       
!...gluons                                                              
        DO 11 I=1,3 
          AT(I)=A(I,1,LF)*T**A(I,2,LF)+A(I,3,LF)*T**(-A(I,4,LF)) 
   11   CONTINUE 
        POMG=AT(1)*DX**AT(2)*(1.D0-DX)**AT(3) 
        DGL=POMG*ALPEM 
!                                                                       
!...quarks                                                              
        E(1)=1.D0 
        E(2)=9.D0 
        DO 13 J=1,2 
          DO 15 I=1,5 
            BTP=B(I,1,J,LF)*T**B(I,2,J,LF) 
            BT(I,J)=BTP+B(I,3,J,LF)*T**(-B(I,4,J,LF)) 
   15     CONTINUE 
   13   CONTINUE 
!                                                                       
!...singlet & non-singlet combinations                                  
        DO 17 J=1,2 
          POM1=DX*(DX*DX+(1.D0-DX)**2)/(BT(1,J)-BT(2,J)*LOG(1.D0-DX)) 
          POM2=BT(3,J)*DX**BT(4,J)*(1.D0-DX)**BT(5,J) 
          XQPOM(J)=E(J)*POM1+POM2 
   17   CONTINUE 
!                                                                       
!...quarks flavours                                                     
      DUB=ALPEM*1.D0/6.D0*(XQPOM(2)+9.D0*XQPOM(1)) 
      DUV=DUB 
      DDB=ALPEM*1.D0/6.D0*(XQPOM(2)-9.D0/2.D0*XQPOM(1)) 
      DDV=DDB 
      DSB=DDB 
      DCB=0.D0 
      DBB=0.D0 
!                                                                       
      RETURN 
      END                                           
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
       SUBROUTINE DGPHO3(DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DGL) 
!********************************************************************   
!*                                                                  *   
!*    Parametrization of parton distribution functions              *   
!*    in the photon (LO analysis) - full  solution of AP eq.!       *   
!*                                                                  *   
!* authors:  M.Drees and K.Grassie (DG)                             *   
!*          /Z. Phys. C28 (1985) 451/                               *   
!*                                                                  *   
!* Prepared by:                                                     *   
!*             Krzysztof Charchula, DESY                            *   
!*             bitnet: F1PCHA@DHHDESY3                              *   
!*             decnet: 13313::CHARCHULA                             *   
!*                                                                  *   
!* Modified by:                                                     *   
!*             H. Plothow-Besch/CERN-PPE                            *   
!*                                                                  *   
!********************************************************************   
!                                                                       
      implicit real*8 (a-h,o-z) 
      double precision                                                  &
     &        A(3,4,3),AT(3),                                           &
     &        B(5,4,2,3),BT(5,2),XQPOM(2),E(2),                         &
     &        DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DGL                     
      PARAMETER (ALPEM=7.29927D-3, PI=3.141592D0) 
      PARAMETER (ALAM=0.4D0) 
!...comments                                                            
!...--------------------------------------------------                  
!...        with nf=4 (valid for  20< Q2 <500 GeV2)                     
!...--------------------------------------------------                  
!                                                                       
!...initialization of gluon parameters array for DG                     
        DATA (((A(I,J,K),I=1,3),J=1,4),K=1,3)/                          &
     &    -0.20700, -0.19870,  5.1190,                                  &
     &     0.61580,  0.62570, -0.2752,                                  &
     &     1.07400,  8.35200, -6.9930,                                  &
     &     0.00000,  5.02400,  2.2980,                                  &
     &     0.8926D-2,0.0509,  -0.2313,                                  &
     &     0.65940,  0.27740,  0.1382,                                  &
     &     0.47660, -0.39060,  6.5420,                                  &
     &     0.01975, -0.32120,  0.5162,                                  &
     &     0.03197, -0.618D-2,-0.1216,                                  &
     &     1.01800,  0.94760,  0.9047,                                  &
     &     0.24610, -0.60940,  2.6530,                                  &
     &     0.02707, -0.01067,  0.2003D-2/                               
!                                                                       
!...initialization of quark parameters array for DG                     
        DATA (((B(I,J,K,1),I=1,5),J=1,4),K=1,2)/                        &
     &     2.2850,   6.0730,  -0.4202,   -0.0808,  0.0553,              &
     &    -0.0153,  -0.8132,   0.0178,    0.6346,  1.1360,              &
     &     1.33D3, -41.310,    0.9216,    1.2080,  0.9512,              &
     &     4.2190,   3.1650,   0.1800,    0.2030,  0.0116,              &
     &    16.690,    0.1760,  -0.0208,   -0.0168, -0.1986,              &
     &    -0.7916,   0.0479,   0.3386D-2, 1.3530,  1.1000,              &
     &     1.0990D3, 1.0470,   4.8530,    1.4260,  1.1360,              &
     &     4.4280,   0.0250,   0.8404,    1.2390, -0.2779/              
        DATA (((B(I,J,K,2),I=1,5),J=1,4),K=1,2)/                        &
     &    -0.3711,  -0.1717,   0.08766,  -0.8915, -0.1816,              &
     &     1.0610,   0.7815,   0.02197,   0.2857,  0.5866,              &
     &     4.7580,   1.5350,   0.10960,   2.9730,  2.4210,              &
     &    -0.0150,   0.7067D-2,0.20400,   0.1185,  0.4059,              &
     &    -0.1207,  25.000,   -0.01230,  -0.0919,  0.02015,             &
     &     1.0710,  -1.6480,   1.16200,   0.7912,  0.9869,              &
     &     1.9770,  -0.01563,  0.48240,   0.6397, -0.07036,             &
     &    -0.8625D-2,6.4380,  -0.01100,   2.3270,  0.01694/             
        DATA (((B(I,J,K,3),I=1,5),J=1,4),K=1,2)/                        &
     &    15.8,      2.742,    0.02917,  -0.0342, -0.02302,             &
     &    -0.9464,  -0.7332,   0.04657,   0.7196,  0.9229,              &
     &    -0.5,      0.7148,   0.1785,    0.7338,  0.5873,              &
     &    -0.2118,   3.287,    0.04811,   0.08139,-0.79D-4,             &
     &     6.734,   59.88,    -0.3226D-2,-0.03321, 0.1059,              &
     &    -1.008,   -2.983,    0.8432,    0.9475,  0.6954,              &
     &    -0.08594,  4.48,     0.3616,   -0.3198, -0.6663,              &
     &     0.07625,  0.9686,   0.1383D-2, 0.02132, 0.3683/              
!                                                                       
!...specification of sets                                               
         NFL=4 
!                                                                       
!...calculations                                                        
       Q2 = DQ*DQ 
       ALAM2=ALAM**2 
       T=LOG(Q2/ALAM2) 
       LF=NFL-2 
!                                                                       
!...gluons                                                              
        DO 11 I=1,3 
          AT(I)=A(I,1,LF)*T**A(I,2,LF)+A(I,3,LF)*T**(-A(I,4,LF)) 
   11   CONTINUE 
        POMG=AT(1)*DX**AT(2)*(1.D0-DX)**AT(3) 
        DGL=POMG*ALPEM 
!                                                                       
!...quarks                                                              
        E(1)=1.D0 
        E(2)=10.D0 
        DO 13 J=1,2 
          DO 15 I=1,5 
            BTP=B(I,1,J,LF)*T**B(I,2,J,LF) 
            BT(I,J)=BTP+B(I,3,J,LF)*T**(-B(I,4,J,LF)) 
   15     CONTINUE 
   13   CONTINUE 
!                                                                       
!...singlet & non-singlet combinations                                  
        DO 17 J=1,2 
          POM1=DX*(DX*DX+(1.D0-DX)**2)/(BT(1,J)-BT(2,J)*LOG(1.D0-DX)) 
          POM2=BT(3,J)*DX**BT(4,J)*(1.D0-DX)**BT(5,J) 
          XQPOM(J)=E(J)*POM1+POM2 
   17   CONTINUE 
!                                                                       
!...quarks flavours                                                     
      DUB=ALPEM*1.D0/8.D0*(XQPOM(2)+6.D0*XQPOM(1)) 
      DUV=DUB 
      DDB=ALPEM*1.D0/8.D0*(XQPOM(2)-6.D0*XQPOM(1)) 
      DDV=DDB 
      DSB=DDB 
      DCB=DUB 
      DBB=0.D0 
!                                                                       
      RETURN 
      END                                           
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
       SUBROUTINE DGPHO4(DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DGL) 
!********************************************************************   
!*                                                                  *   
!*    Parametrization of parton distribution functions              *   
!*    in the photon (LO analysis) - full  solution of AP eq.!       *   
!*                                                                  *   
!* authors:  M.Drees and K.Grassie (DG)                             *   
!*          /Z. Phys. C28 (1985) 451/                               *   
!*                                                                  *   
!* Prepared by:                                                     *   
!*             Krzysztof Charchula, DESY                            *   
!*             bitnet: F1PCHA@DHHDESY3                              *   
!*             decnet: 13313::CHARCHULA                             *   
!*                                                                  *   
!* Modified by:                                                     *   
!*             H. Plothow-Besch/CERN-PPE                            *   
!*                                                                  *   
!********************************************************************   
!                                                                       
      implicit real*8 (a-h,o-z) 
      double precision                                                  &
     &        A(3,4,3),AT(3),                                           &
     &        B(5,4,2,3),BT(5,2),XQPOM(2),E(2),                         &
     &        DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DGL                     
      PARAMETER (ALPEM=7.29927D-3, PI=3.141592D0) 
      PARAMETER (ALAM=0.4D0) 
!...comments                                                            
!...--------------------------------------------------                  
!...        with nf=5 (valid for 200< Q2 <1D4 GeV2)                     
!...--------------------------------------------------                  
!                                                                       
!...initialization of gluon parameters array for DG                     
        DATA (((A(I,J,K),I=1,3),J=1,4),K=1,3)/                          &
     &    -0.20700, -0.19870,  5.1190,                                  &
     &     0.61580,  0.62570, -0.2752,                                  &
     &     1.07400,  8.35200, -6.9930,                                  &
     &     0.00000,  5.02400,  2.2980,                                  &
     &     0.8926D-2,0.0509,  -0.2313,                                  &
     &     0.65940,  0.27740,  0.1382,                                  &
     &     0.47660, -0.39060,  6.5420,                                  &
     &     0.01975, -0.32120,  0.5162,                                  &
     &     0.03197, -0.618D-2,-0.1216,                                  &
     &     1.01800,  0.94760,  0.9047,                                  &
     &     0.24610, -0.60940,  2.6530,                                  &
     &     0.02707, -0.01067,  0.2003D-2/                               
!                                                                       
!...initialization of quark parameters array for DG                     
        DATA (((B(I,J,K,1),I=1,5),J=1,4),K=1,2)/                        &
     &     2.2850,   6.0730,  -0.4202,   -0.0808,  0.0553,              &
     &    -0.0153,  -0.8132,   0.0178,    0.6346,  1.1360,              &
     &     1.33D3, -41.310,    0.9216,    1.2080,  0.9512,              &
     &     4.2190,   3.1650,   0.1800,    0.2030,  0.0116,              &
     &    16.690,    0.1760,  -0.0208,   -0.0168, -0.1986,              &
     &    -0.7916,   0.0479,   0.3386D-2, 1.3530,  1.1000,              &
     &     1.0990D3, 1.0470,   4.8530,    1.4260,  1.1360,              &
     &     4.4280,   0.0250,   0.8404,    1.2390, -0.2779/              
        DATA (((B(I,J,K,2),I=1,5),J=1,4),K=1,2)/                        &
     &    -0.3711,  -0.1717,   0.08766,  -0.8915, -0.1816,              &
     &     1.0610,   0.7815,   0.02197,   0.2857,  0.5866,              &
     &     4.7580,   1.5350,   0.10960,   2.9730,  2.4210,              &
     &    -0.0150,   0.7067D-2,0.20400,   0.1185,  0.4059,              &
     &    -0.1207,  25.000,   -0.01230,  -0.0919,  0.02015,             &
     &     1.0710,  -1.6480,   1.16200,   0.7912,  0.9869,              &
     &     1.9770,  -0.01563,  0.48240,   0.6397, -0.07036,             &
     &    -0.8625D-2,6.4380,  -0.01100,   2.3270,  0.01694/             
        DATA (((B(I,J,K,3),I=1,5),J=1,4),K=1,2)/                        &
     &    15.8,      2.742,    0.02917,  -0.0342, -0.02302,             &
     &    -0.9464,  -0.7332,   0.04657,   0.7196,  0.9229,              &
     &    -0.5,      0.7148,   0.1785,    0.7338,  0.5873,              &
     &    -0.2118,   3.287,    0.04811,   0.08139,-0.79D-4,             &
     &     6.734,   59.88,    -0.3226D-2,-0.03321, 0.1059,              &
     &    -1.008,   -2.983,    0.8432,    0.9475,  0.6954,              &
     &    -0.08594,  4.48,     0.3616,   -0.3198, -0.6663,              &
     &     0.07625,  0.9686,   0.1383D-2, 0.02132, 0.3683/              
!                                                                       
!...specification of sets                                               
         NFL=5 
!                                                                       
!...calculations                                                        
       Q2 = DQ*DQ 
       ALAM2=ALAM**2 
       T=LOG(Q2/ALAM2) 
       LF=NFL-2 
!                                                                       
!...gluons                                                              
        DO 11 I=1,3 
          AT(I)=A(I,1,LF)*T**A(I,2,LF)+A(I,3,LF)*T**(-A(I,4,LF)) 
   11   CONTINUE 
        POMG=AT(1)*DX**AT(2)*(1.D0-DX)**AT(3) 
        DGL=POMG*ALPEM 
!                                                                       
!...quarks                                                              
        E(1)=1.D0 
        E(2)=55.D0/6.D0 
        DO 13 J=1,2 
          DO 15 I=1,5 
            BTP=B(I,1,J,LF)*T**B(I,2,J,LF) 
            BT(I,J)=BTP+B(I,3,J,LF)*T**(-B(I,4,J,LF)) 
   15     CONTINUE 
   13   CONTINUE 
!                                                                       
!...singlet & non-singlet combinations                                  
        DO 17 J=1,2 
          POM1=DX*(DX*DX+(1.D0-DX)**2)/(BT(1,J)-BT(2,J)*LOG(1.D0-DX)) 
          POM2=BT(3,J)*DX**BT(4,J)*(1.D0-DX)**BT(5,J) 
          XQPOM(J)=E(J)*POM1+POM2 
   17   CONTINUE 
!                                                                       
!...quarks flavours                                                     
      DUB=ALPEM*1.D0/10.D0*(XQPOM(2)+15.D0/2.D0*XQPOM(1)) 
      DUV=DUB 
      DCB=DUB 
      DDB=ALPEM*1.D0/10.D0*(XQPOM(2)-5.D0*XQPOM(1)) 
      DDV=DDB 
      DSB=DDB 
      DCB=DUB 
      DBB=DDB 
!                                                                       
      RETURN 
      END                                           
