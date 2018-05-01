! Contour the force(s) over a domain of k and omegais.
! For psi and Omegac having some values.

program main
  use shiftmode 
  integer, parameter ::   nk=21
  real :: kik(nk)
  integer, parameter :: no=20
  real :: oi(no)
  real, dimension(nk,no) :: Ftnp,Fpnp,cworka,Fsum
  integer :: icl,icsw,istable
  real :: zclv(20)
  character*30 string
  
  istable=0
  psi=.1

  call initialize
  akmax=0.32             ! range of k/sqrt(psi)
!  akmax=0.7             ! range of k/sqrt(psi)
  akmin=0.002            ! Lowest k plotted.
  omegaimax=0.08         ! range of omega/sqrt(psi)
  omegacmax=.8           ! Maximum fraction of Omegac/omega_b.
  
  call pfset(3)
  call pltinit(0.,akmax,0.,omegaimax)
  call axis
  call axis2
  call axlabels('!Bk!@/!A)y!@','!Ag!@/!A)y!@')
  string='!Ay!@='
  call fwrite(psi,iwidth,2, string(lentrim(string)+1:))
  call legendline(.8,.9,258,string)
  
  nioc=24                ! Number of omega cases.
  nfac=4                 ! Number devoted to logarithmic variation
  !  facdiv=0.5
  ! This might be better for lower psi.
  !  omegacmax=1.04         ! Maximum fraction of Omegac/omega_b.
  !  facdiv=.64/1.04       ! Fraction of range that is logarithmic.
  OmObdiv=.48        ! Divide of log and linear
  dOm=.02            ! Spacing of linear contours
  ! single contour case for testing
!  OmObdiv=.66
!  nioc=1
!  nfac=1
  do ioc=1,nioc
     OmOb=OmObdiv/2**(nfac-ioc)
  !   if(ioc.gt.nfac)OmOb=OmObdiv+  &
  !        (ioc-nfac)*(omegacmax-OmObdiv)/(nioc-nfac) ! For linear
     if(ioc.gt.nfac)OmOb=OmObdiv+  &
          (ioc-nfac)*dOm  ! For linear
     
     Omegac=OmOb*sqrt(psi)/2.
     write(*,*)'OmOb=',Omob
    
  do ik=1,nk   ! k-scan
     kik(ik)=max(akmin,(ik-1)*akmax/(nk-1))
     k=kik(ik)*sqrt(psi)     
     do io=1,no  ! Iterate over omega.
        ! Non-uniform \propto f(1+f)/2 so that at f=1, spacing greater.
        oi(io)=io*omegaimax/no    *(1.+float(io)/no)/2.
        omega=(0.,0.)+sqm1*oi(io)*sqrt(psi)
        so=abs(omega**2)
        call SumHarmonics()
        Fpnp(ik,io)=real(Fpasstotal)/so
        Ftnp(ik,io)=real(Ftraptotal)/so
        !        Fsum(ik,io)=Fpnp(ik,io)+Ftnp(ik,io) ! Without F_E.
        Fsum(ik,io)=Fpnp(ik,io)+Ftnp(ik,io)-(psi*k)**2*128./315.
     enddo
  enddo
  
  call minmax2(Fsum,nk,nk,no,Fmin,Fmax)
  call color(mod(ioc-1,12)+1)
  if(Fmax.lt.0)then
     call color(15)
     istable=istable+1
     if(istable.le.1)then
        string='Stable !AW!@/!Aw!@!db!d='
        call fwrite(OmOb,iwidth,2,string(lentrim(string)+1:))
        call legendline(.1,.04,258,string)
     else
        exit  ! From the iterations.
     endif
  else
  icl=1
  zclv(1)=OmOb
  icsw=1
  ! Fool contour into plotting the zero level but labeling it with OmOb.
  ! Make the labels upright by using -Fsum
     call contourl(-Fsum+OmOb,cworka,nk,nk,no,zclv,icl,kik,oi,icsw)
  endif
  call accisflush()
  enddo
  write(*,*)'Finished'
  call pltend()
  
end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
