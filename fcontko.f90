! Contour the force(s) over a domain of k and omegais.
! For psi and Omegac having some values.
include 'fhgfunc.f'

program fcontko
  use shiftmode 
  use acpath
  integer, parameter ::   nk=21
  real :: kik(nk)
  integer, parameter :: no=20
  real :: oi(no)
  real, dimension(nk,no) :: Ftnp,Fpnp,cworka,Fsum,FE
  integer :: icl,icsw,istable
  real :: zclv(20)
  character*30 string
  
  psi=0.1
  Ty=1.
  call initialize
  akmin=0.002            ! Lowest k plotted.
  akmax=0.32/sqrt(Ty)    ! range of k/sqrt(psi)
  slopelk=growthlk(psi,beta,Ty) ! Get the low-k analytic gamma-slope.
  omegaimax=slopelk*akmax/4. ! Estimated range of omega/sqrt(psi)
  omegacmax=.8           ! Maximum fraction of Omegac/omega_b.
  nioc=20                ! Number of Omegac cases.
  nfac=4                 ! Number devoted to logarithmic variation
  OmObdiv=.48            ! Pivot from log to linear
  dOm=.02                ! Spacing of linear contours
  
  istable=0
  call pfset(3)
  call pltinit(0.,akmax,0.,omegaimax)
  call axis
  call axis2
  call axlabels('!Bk!@/!A)y!@','!Ag!@/!A)y!@')
  string='!Ay!@='
  call fwrite(psi,iwidth,2, string(lentrim(string)+1:))
  call legendline(.8,.9,258,string)
  string='!BT!dy!d!@='
  call fwrite(Ty,iwidth,2, string(lentrim(string)+1:))
  call legendline(.8,.83,258,string)
  call polyline((/0.,.1/)*akmax,(/0.,.1/)*slopelk*akmax,2)
  do ioc=1,nioc
     OmOb=OmObdiv/2**(nfac-ioc)
     if(ioc.gt.nfac)OmOb=OmObdiv+  &
          (ioc-nfac)*dOm  ! For linear
     
     Omegac=OmOb*sqrt(psi)/2.
     write(*,*)'OmOb=',Omob
    
     do ik=1,nk   ! k-scan
        kik(ik)=max(akmin,(ik-1)*akmax/(nk-1))
        k=kik(ik)*sqrt(psi)     
        do io=1,no  ! Iterate over omegai.
           ! Non-uniform \propto f(1+f)/2 so that at f=1, spacing greater.
           oi(io)=io*omegaimax/no    *(1.+float(io)/no)/2.
           omega=(0.,0.)+sqm1*oi(io)*sqrt(psi)
           so=abs(omega**2)
           call SumHarmonics()
           Fpnp(ik,io)=real(Fpasstotal)
           Ftnp(ik,io)=real(Ftraptotal)
           FE(ik,io)=(psi*k)**2*128./315.
           Fsum(ik,io)=(Fpnp(ik,io)+Ftnp(ik,io)-FE(ik,io))
!           write(*,*)oi(io),Fpnp(ik,io),Ftnp(ik,io),FE(ik,io)
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
     ! Testing the population of the contour arrays.
     write(*,*)'iarray=',iarray   !,' cv=',cvacpath
     write(*,'(i4,2f8.4)')(i,xcarray(i),ycarray(i),i=1,iarray)
     
  enddo
  write(*,*)'Finished'
  call pltend()
  
end program fcontko
