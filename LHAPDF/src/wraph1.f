! -*- F90 -*-


      subroutine H1evolve(xin,qin,pdf) 
      implicit real*8 (a-h,o-z) 
!******************************************************                 
!  done on  13/07/04 at  10.04.39                                       
! evolution has been made starting at q2_input =      4.000             
!  available for :       1.500 <= q2 <= 1000000.000                     
!            and :  0.000057 <= x <=  0.906052                          
!*                                                                      
!   for x outside limits, the closest limit                             
!        is assumed : f2(x>xmax,q2)=f2(xmax,q2)                         
!                     f2(x<xmin,q2)=f2(xmin,q2)                         
!  for q2 outside limits, the closest limit                             
!        is assumed : f2(x,q2>q2max)=f2(x,q2max)                        
!                     f2(x,q2<q2min)=f2(x,q2min)                        
!*                                                                      
!  comments, etc... to C. Pascaud or F. Zomer                           
!*******************************************************                
      include 'parmsetup.inc' 
      PARAMETER(n_bin_q2=86) 
      PARAMETER(n_bin_x=100) 
      REAL*4 xl_bin(n_bin_x),q2l_bin(n_bin_q2) 
      PARAMETER(ngrid=20) 
      REAL*4 f(0:ngrid,8,n_bin_x,n_bin_q2),val(8) 
      real*8 pdf(-6:6) 
      real*4 q2in,x,y 
      character*16 name(nmxset) 
      integer nmem(nmxset),ndef(nmxset),mmem 
      common/NAME/name,nmem,ndef,mmem 
      double precision gridx(nmxgridx),gridq(nmxgridq)
      integer ngridx,ngridq,jx,jq
      integer nset,iset 
      save 
!                                                                       
!      enddo                                                            
      call getnset(iset) 
      call getnmem(iset,imem) 
!                                                                       
      q2in = qin*qin 
      x=log(xin) 
      y=log(q2in) 
      DO i=2,n_bin_x 
        IF(x.LT.xl_bin(i))  goto 1 
        IF(xl_bin(i).ge.0.)  goto 1 
      ENDDO 
      i=n_bin_x 
    1 i=i-1 
      DO j=2,n_bin_q2 
        IF(y.LT.q2l_bin(j))  GOTO 2 
      ENDDO 
      j=n_bin_q2 
    2 j=j-1 
      dx=xl_bin(i+1)-xl_bin(i) 
      xd=(x-xl_bin(i))/dx 
      dy=q2l_bin(j+1)-q2l_bin(j) 
      yd=(y-q2l_bin(j))/dy 
!                                                                       
      do k=1,8 
      val(k)=f(imem,k,i,j)+xd*(f(imem,k,i+1,j)-f(imem,k,i,j))           &
     &+yd*(f(imem,k,i,j+1)-f(imem,k,i,j))                               &
     &+xd*yd*(f(imem,k,i+1,j+1)+f(imem,k,i,j)                           &
     &-f(imem,k,i+1,j)-f(imem,k,i,j+1))                                 
      enddo 
      pdf(-6) = 0.0d0 
       pdf(6) = 0.0d0 
      pdf(-5) = val(7) 
       pdf(5) = val(7) 
      pdf(-4) = val(6) 
       pdf(4) = val(6) 
      pdf(-3) = val(5) 
       pdf(3) = val(5) 
      pdf(-2) = val(4) 
       pdf(2) = val(3)+val(4) 
      pdf(-1) = val(2) 
       pdf(1) = val(1)+val(2) 
       pdf(0) = val(8) 
      return 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      
      entry H1getgrid(nset,ngridx,ngridq,gridx,gridq)
     
      ngridx=n_bin_x
      do jx=1,ngridx
          gridx(jx)=exp(xl_bin(jx))
      enddo
      ngridq=n_bin_q2
      do jq=1,ngridq
          gridq(jq)=exp(q2l_bin(jq))
      enddo
       
      return
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                       
      entry H1read(nset) 
!                                                                       
      read(1,*)nmem(nset),ndef(nset) 
      read(1,1000)xl_bin 
      read(1,1000)q2l_bin 
        do i=1,n_bin_q2 
          q2l_bin(i)=log(q2l_bin(i)) 
        enddo 
      do nm = 0,nmem(nset) 
        do jval = 1,8 
      read(1,1000)((f(nm,jval,nx,nq2),nx=1,n_bin_x),nq2=1,n_bin_q2) 
        enddo 
      enddo 
      return 
!                                                                       
      entry H1alfa(alfas,qalfa) 
      call alphah1(alfas,Qalfa) 
      return 
!                                                                       
      entry H1init(Eorder,Q2fit) 
      return 
!                                                                       
      entry H1pdf(mem) 
      call getnset(iset) 
      call setnmem(iset,mem) 
!      imem = mem                                                       
      return 
!                                                                       
 1000 format(5e13.5) 
      END                                           
!                                                                       
      subroutine alphah1(alpha,Qin) 
      implicit real*8 (a-h,o-z) 
        call getnset(nset) 
        call GetOrderAsM(nset,iord) 
        if(iord.eq.1) then 
          call alphah1nlo(alpha,Qin) 
        elseif(iord.eq.0) then 
          call alphah1lo(alpha,Qin) 
        else 
          print *,'iord = ',iord 
          stop 
        endif 
      return 
      END                                           
!                                                                       
      subroutine alphah1nlo(alpha,Qin) 
      implicit real*8 (a-h,o-z) 
!***************************************************                    
!  done on  13/07/04 at  09.10.39                                       
! evolution has been made starting at q2_input =      4.000             
!  available for :       1.500 <= q2 <= 1000000.000                     
!  for q2 outside limits, the closest limit                             
!        is assumed : f2(x,q2>q2max)=f2(x,q2max)                        
!                     f2(x,q2<q2min)=f2(x,q2min)                        
!*                                                                      
!  comments, etc... to C. Pascaud or F. Zomer                           
!**************************************************                     
      PARAMETER(n_bin_q2=102) 
      dimension q2l_bin(n_bin_q2) 
      dimension f(n_bin_q2) 
      data q2l_bin/                                                     &
     &1.500000E+00,1.600000E+00,1.700000E+00,1.800000E+00,1.900000E+00, &
     &1.959902E+00,1.960000E+00,1.960098E+00,2.000000E+00,2.100000E+00, &
     &2.200000E+00,2.300000E+00,2.400000E+00,2.500000E+00,3.000000E+00, &
     &3.500000E+00,4.000000E+00,4.500000E+00,5.000000E+00,6.000000E+00, &
     &7.000000E+00,8.000000E+00,9.000000E+00,1.000000E+01,1.500000E+01, &
     &2.000000E+01,2.024899E+01,2.025000E+01,2.025101E+01,2.500000E+01, &
     &3.000000E+01,3.500000E+01,4.000000E+01,4.500000E+01,5.000000E+01, &
     &5.500000E+01,6.000000E+01,6.500000E+01,7.000000E+01,7.500000E+01, &
     &8.000000E+01,8.500000E+01,9.000000E+01,9.500000E+01,1.000000E+02, &
     &1.500000E+02,2.000000E+02,2.500000E+02,3.000000E+02,3.500000E+02, &
     &4.000000E+02,4.500000E+02,5.000000E+02,5.500000E+02,6.000000E+02, &
     &6.500000E+02,7.000000E+02,7.500000E+02,8.000000E+02,8.500000E+02, &
     &9.000000E+02,9.500000E+02,1.000000E+03,1.500000E+03,2.000000E+03, &
     &2.500000E+03,3.000000E+03,3.500000E+03,4.000000E+03,4.500000E+03, &
     &5.000000E+03,5.500000E+03,6.000000E+03,6.500000E+03,7.000000E+03, &
     &7.500000E+03,8.000000E+03,8.500000E+03,9.000000E+03,9.500000E+03, &
     &1.000000E+04,1.500000E+04,2.000000E+04,2.500000E+04,3.000000E+04, &
     &3.500000E+04,4.000000E+04,4.500000E+04,5.000000E+04,5.500000E+04, &
     &6.000000E+04,6.500000E+04,7.000000E+04,7.500000E+04,8.000000E+04, &
     &8.500000E+04,9.000000E+04,9.500000E+04,1.000000E+05,1.500000E+05, &
     &2.000000E+05,1.000000E+06/                                        
      data f/                                                           &
     &3.935326E-01,3.849873E-01,3.773198E-01,3.703884E-01,3.640814E-01, &
     &3.605647E-01,3.605591E-01,3.605540E-01,3.585220E-01,3.537029E-01, &
     &3.492356E-01,3.450786E-01,3.411967E-01,3.375601E-01,3.222805E-01, &
     &3.104663E-01,3.009529E-01,2.930611E-01,2.863649E-01,2.755124E-01, &
     &2.669931E-01,2.600508E-01,2.542359E-01,2.492618E-01,2.318890E-01, &
     &2.210261E-01,2.205828E-01,2.205810E-01,2.205794E-01,2.139838E-01, &
     &2.085977E-01,2.042587E-01,2.006485E-01,1.975722E-01,1.949020E-01, &
     &1.925502E-01,1.904538E-01,1.885667E-01,1.868536E-01,1.852875E-01, &
     &1.838468E-01,1.825145E-01,1.812765E-01,1.801214E-01,1.790394E-01, &
     &1.709373E-01,1.656324E-01,1.617456E-01,1.587067E-01,1.562275E-01, &
     &1.541435E-01,1.523521E-01,1.507856E-01,1.493968E-01,1.481517E-01, &
     &1.470250E-01,1.459974E-01,1.450540E-01,1.441827E-01,1.433740E-01, &
     &1.426200E-01,1.419142E-01,1.412513E-01,1.362263E-01,1.328778E-01, &
     &1.303943E-01,1.284347E-01,1.268244E-01,1.254625E-01,1.242858E-01, &
     &1.232522E-01,1.223323E-01,1.215046E-01,1.207533E-01,1.200661E-01, &
     &1.194337E-01,1.188480E-01,1.183032E-01,1.177942E-01,1.173169E-01, &
     &1.168677E-01,1.134367E-01,1.111245E-01,1.093963E-01,1.080244E-01, &
     &1.068916E-01,1.059298E-01,1.050959E-01,1.043613E-01,1.037057E-01, &
     &1.031144E-01,1.025766E-01,1.020837E-01,1.016292E-01,1.012077E-01, &
     &1.008150E-01,1.004476E-01,1.001026E-01,9.977747E-02,9.728142E-02, &
     &9.558620E-02,8.711093E-02/                                        
      data init/0/ 
      if(init.eq.0) then 
        do i=1,n_bin_q2 
         q2l_bin(i) = log(q2l_bin(i)) 
        enddo 
        init=1 
      endif 
!                                                                       
      q2in = qin*qin 
!                                                                       
      y = log(q2in) 
      do j = 2, n_bin_q2 
        if (y.lt.q2l_bin(j))  goto 2 
      enddo 
      j=n_bin_q2 
    2 j=j-1 
!                                                                       
      dy = q2l_bin(j+1) - q2l_bin(j) 
      yd = (y - q2l_bin(j)) / dy 
      alpha = f(j) + yd*(f(j+1)-f(j)) 
!                                                                       
      return 
      END                                           
                                                                        
!                                                                       
      subroutine alphah1lo(alpha,Qin) 
      implicit real*8 (a-h,o-z) 
!  done on  13/07/04 at  11.32.26                                       
! evolution has been made starting at q2_input =      4.000             
!  available for :       1.500 <= q2 <= 1000000.000                     
!  for q2 outside limits, the closest limit                             
!        is assumed : f2(x,q2>q2max)=f2(x,q2max)                        
!                     f2(x,q2<q2min)=f2(x,q2min)                        
!*                                                                      
!  comments, etc... to C. Pascaud or F. Zomer                           
!**************************************************                     
      PARAMETER(n_bin_q2=102) 
      dimension q2l_bin(n_bin_q2) 
      dimension f(n_bin_q2) 
      data q2l_bin/                                                     &
     &1.500000E+00,1.600000E+00,1.700000E+00,1.800000E+00,1.900000E+00, &
     &1.959902E+00,1.960000E+00,1.960098E+00,2.000000E+00,2.100000E+00, &
     &2.200000E+00,2.300000E+00,2.400000E+00,2.500000E+00,3.000000E+00, &
     &3.500000E+00,4.000000E+00,4.500000E+00,5.000000E+00,6.000000E+00, &
     &7.000000E+00,8.000000E+00,9.000000E+00,1.000000E+01,1.500000E+01, &
     &2.000000E+01,2.024899E+01,2.025000E+01,2.025101E+01,2.500000E+01, &
     &3.000000E+01,3.500000E+01,4.000000E+01,4.500000E+01,5.000000E+01, &
     &5.500000E+01,6.000000E+01,6.500000E+01,7.000000E+01,7.500000E+01, &
     &8.000000E+01,8.500000E+01,9.000000E+01,9.500000E+01,1.000000E+02, &
     &1.500000E+02,2.000000E+02,2.500000E+02,3.000000E+02,3.500000E+02, &
     &4.000000E+02,4.500000E+02,5.000000E+02,5.500000E+02,6.000000E+02, &
     &6.500000E+02,7.000000E+02,7.500000E+02,8.000000E+02,8.500000E+02, &
     &9.000000E+02,9.500000E+02,1.000000E+03,1.500000E+03,2.000000E+03, &
     &2.500000E+03,3.000000E+03,3.500000E+03,4.000000E+03,4.500000E+03, &
     &5.000000E+03,5.500000E+03,6.000000E+03,6.500000E+03,7.000000E+03, &
     &7.500000E+03,8.000000E+03,8.500000E+03,9.000000E+03,9.500000E+03, &
     &1.000000E+04,1.500000E+04,2.000000E+04,2.500000E+04,3.000000E+04, &
     &3.500000E+04,4.000000E+04,4.500000E+04,5.000000E+04,5.500000E+04, &
     &6.000000E+04,6.500000E+04,7.000000E+04,7.500000E+04,8.000000E+04, &
     &8.500000E+04,9.000000E+04,9.500000E+04,1.000000E+05,1.500000E+05, &
     &2.000000E+05,1.000000E+06/                                        
      data f/                                                           &
     &4.395646E-01,4.308115E-01,4.229009E-01,4.157042E-01,4.091185E-01, &
     &4.054310E-01,4.054251E-01,4.054197E-01,4.032349E-01,3.980418E-01, &
     &3.932134E-01,3.887078E-01,3.844897E-01,3.805290E-01,3.637916E-01, &
     &3.507479E-01,3.401822E-01,3.313772E-01,3.238784E-01,3.116737E-01, &
     &3.020502E-01,2.941818E-01,2.875740E-01,2.819097E-01,2.620464E-01, &
     &2.495700E-01,2.490599E-01,2.490579E-01,2.490560E-01,2.413308E-01, &
     &2.350218E-01,2.299395E-01,2.257114E-01,2.221089E-01,2.189825E-01, &
     &2.162291E-01,2.137753E-01,2.115667E-01,2.095621E-01,2.077297E-01, &
     &2.060444E-01,2.044861E-01,2.030382E-01,2.016875E-01,2.004225E-01, &
     &1.909551E-01,1.847628E-01,1.802294E-01,1.766873E-01,1.737993E-01, &
     &1.713728E-01,1.692881E-01,1.674658E-01,1.658508E-01,1.644033E-01, &
     &1.630939E-01,1.619001E-01,1.608043E-01,1.597925E-01,1.588537E-01, &
     &1.579785E-01,1.571596E-01,1.563904E-01,1.505657E-01,1.466892E-01, &
     &1.438172E-01,1.415527E-01,1.396930E-01,1.381211E-01,1.367637E-01, &
     &1.355719E-01,1.345115E-01,1.335578E-01,1.326924E-01,1.319011E-01, &
     &1.311728E-01,1.304988E-01,1.298719E-01,1.292864E-01,1.287374E-01, &
     &1.282208E-01,1.242789E-01,1.216259E-01,1.196449E-01,1.180735E-01, &
     &1.167767E-01,1.156763E-01,1.147226E-01,1.138828E-01,1.131336E-01, &
     &1.124583E-01,1.118440E-01,1.112813E-01,1.107625E-01,1.102815E-01, &
     &1.098335E-01,1.094145E-01,1.090210E-01,1.086503E-01,1.058065E-01, &
     &1.038775E-01,9.426288E-02/                                        
      data init/0/ 
      if(init.eq.0) then 
        do i=1,n_bin_q2 
         q2l_bin(i) = log(q2l_bin(i)) 
        enddo 
        init=1 
      endif 
                                                                        
      q2in = qin*qin 
!                                                                       
      y = log(q2in) 
      do j = 2, n_bin_q2 
        if (y.lt.q2l_bin(j))  goto 2 
      enddo 
      j=n_bin_q2 
    2 j=j-1 
!                                                                       
      dy = q2l_bin(j+1) - q2l_bin(j) 
      yd = (y - q2l_bin(j)) / dy 
      alpha = f(j) + yd*(f(j+1)-f(j)) 
!                                                                       
      return 
      END                                           
