! -*- F90 -*-


subroutine GetOrderPDF(order)
  implicit none
  integer nset,order
  double precision Q,Q2
  nset = 1
  call GetOrderPDFM(nset,order)
  return

  entry GetRenFac(Q)
  nset = 1
  call GetRenFacM(nset,Q)
  return

  entry GetQ2fit(Q2)
  nset = 1
  call GetQ2fitM(nset,Q2)
  return
end subroutine GetOrderPDF

subroutine InitEvolve(nset)
  implicit none
  character*16 s2
  character*64 ctemp
  double precision Q,Q2,Q2fit,muR
  integer nset,order,EvlOrd,nFlav,nf,i1,i2,i3
  
  save EvlOrd,Q2fit,muR,nFlav
  
  EvlOrd=-1
  read(1,'(a)') ctemp
  i1=index(ctemp,',')
  i2=index(ctemp(i1+1:),',')
  i3=index(ctemp(i1+i2+1:),',')
  if(i3>0) then
    read(ctemp,*) s2,Q2fit,muR,nFlav
  else
    read(ctemp,*) s2,Q2fit,muR
    nFlav=-1
  endif
  if (index(s2,'lo').eq.1) EvlOrd=0
  if (index(s2,'nlo').eq.1) EvlOrd=1
  if (index(s2,'nnlo').eq.1) EvlOrd=2
  if (EvlOrd.lt.0) then
     write(*,*) 'File description error:'
     write(*,*) 'Unknown PDF evolution order ',s2
     stop
  endif
  if (muR.ne.1.0d0) then
     write(*,*) '***********************************************'
     write(*,*) '* Note than the renormalization scale is      *'
     write(*,*) '* unequal to the factorization scale for this *'
     write(*,*) '* particular PDF set.                         *'
     write(*,*) '* See manual for proper use.                  *'
     write(*,*) '***********************************************'
  endif
  ! print *,'calling readevolve', nset
  call readevolve(nset)
  return
  
  entry InitEvolveCode(nset)
  !print *,'calling initevolution', nset,Evlord,Q2fit
  call initevolution(nset,EvlOrd,Q2fit)
  return
  
  entry GetOrderPDFM(nset,order)
  order=EvlOrd
  return
  
  entry GetRenFacM(nset,Q)
  Q=muR
  return

  entry GetQ2fitM(nset,Q2)
  Q2=Q2fit
  return  

  entry GetNFlavM(nset,nf)
  nf=nFlav
  return  

end subroutine InitEvolve
