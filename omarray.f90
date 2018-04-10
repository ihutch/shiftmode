! Given psi calculate the forces for a range of k and omegas.
! Generalized to allow Omegac to be set too. 

program main
  use shiftmode 
  integer, parameter ::   nk=3
  real :: kik(nk),omik(nk),omikm(nk)
  complex :: Fcpassing(nk),Ftrapped(nk)
  real :: Fmpassing(nk),Fmtrapped(nk)
  integer, parameter :: np=20
!  integer, parameter :: np=2
  real, dimension(np,nk) :: Ftnp,Fpnp,omi,Fmtnp,Fmpnp
!  character*10 :: string
!  real :: omega0,omega1,oval
  logical :: lcompare=.true.
  external forcebalance
  psi=.1
  Omegac=.09

  call initialize
! k-scan
  akmax=0.1*sqrt(psi)
  akmin=akmax/20.
!  write(*,*)'   k      omegai   Fpassing                Ftrapped'
  do ik=1,nk
     k=max(akmin,(ik-1.)*akmax/(max(nk-1.,1.)))
     kik(ik)=k
     omegaimax=1.5*min(k,akmax/2.)
     omik(ik)=0.
     omikm(ik)=0.

     do ip=1,np  ! Iterate over omega.
        omega=(0.,0.)+sqm1*(np+0.5-ip)*omegaimax/np
        so=abs(omega**2)
        omi(ip,ik)=imag(omega)
        if(Omegac.gt.0.)then
           call SumHarmonics()
           Fmpassing(ik)=real(Fpasstotal)
           Fmtrapped(ik)=real(Ftraptotal)
           Fmtnp(ip,ik)=real(Fmtrapped(ik))/so
           Fmpnp(ip,ik)=real(Fmpassing(ik))/so
           if(lcompare)write(*,'(f6.4,3es12.3,i4)')k,imag(omega) &
                ,Fmpnp(ip,ik),Fmtnp(ip,ik),nharmonics
        endif

        if(lcompare.or.Omegac.eq.0)then
           call FpVyint()
           call FtVyint()
           Fcpassing(ik)=Fpasstotal
           Ftrapped(ik)=Ftraptotal
           Ftnp(ip,ik)=real(Ftrapped(ik))/so
           Fpnp(ip,ik)=real(Fcpassing(ik))/so
           write(*,'(f6.4,5es12.3)')k,imag(omega),Fpnp(ip,ik),Ftnp(ip,ik)
        else
           ! Always put a plausible calculation into the non-magnetic
           Ftnp(ip,ik)=Fmtnp(ip,ik)
           Fpnp(ip,ik)=Fmpnp(ip,ik)
        endif

        if(ip.gt.1)then
           Fip=Ftnp(ip,ik)+Fpnp(ip,ik)
           Fipm=Ftnp(ip-1,ik)+Fpnp(ip-1,ik)
           if(Fip*Fipm.le.0)then
              omik(ik)=(Fipm*omi(ip,ik)+Fip*omi(ip-1,ik))/(Fip+Fipm)
              write(*,'(a,f8.4,a,f8.4)')'k=',k,' Intercept  at omegai=',omik(ik)
           endif
           Fip=Fmtnp(ip,ik)+Fmpnp(ip,ik)
           Fipm=Fmtnp(ip-1,ik)+Fmpnp(ip-1,ik)
           if(Fip*Fipm.le.0)then
              omikm(ik)=(Fipm*omi(ip,ik)+Fip*omi(ip-1,ik))/(Fip+Fipm)
              write(*,'(a,f8.4,a,f8.4,i3)')'k=',k,' Interceptm at omegai=',omikm(ik),nharmonics
           endif
        endif

     enddo
  enddo

!  write(*,*)'omik',omik
!  write(*,*)'kik ',kik

  ik=nk
  call pfset(3)
  call multiframe(2,1,2)
  call pltinit(0.,omi(1,nk),-.2,.2)
!  call pltinit(0.,omi(1,nk),-1.,1.)
  call axis()
  call axlabels(' ','Normalized Force')
  call winset(.true.)
  call legendline(.8,.9,258,'F!dp!d')
  call legendline(.5,.1,258,'F!dt!d')
  do ik=nk,1,-1
     call color(ik)
        call dashset(2)
        call polyline(omi(:,ik),Ftnp(:,ik),np)
        call polyline(omi(:,ik),Fpnp(:,ik),np)
        call dashset(0)
     if(Omega.ne.0.and.lcompare)then
        call polyline(omi(:,ik),Fmtnp(:,ik),np)
        call polyline(omi(:,ik),Fmpnp(:,ik),np)
     endif
  enddo
  call color(15)

  call pltinit(0.,omi(1,nk),-3.,3.)
  call axis()
  call axlabels('!Aw!B!di!d!@','Normalized Force')
  call winset(.true.)
  do ik=nk,1,-1
     call dashset(2)
     call color(ik)
     call polyline(omi(:,ik),Ftnp(:,ik),np)
     call polyline(omi(:,ik),Fpnp(:,ik),np)
     call dashset(0)
     if(Omegac.ne.0..and.lcompare)then
        call polyline(omi(:,ik),Fmtnp(:,ik),np)
        call polyline(omi(:,ik),Fmpnp(:,ik),np)
     endif
  enddo
  call pltend()
  call multiframe(0,0,0)

  call autoplot(kik,omikm,nk)
  call axlabels('k','omegai')
  call dashset(1)
  call polyline(kik,omik,nk)
  call dashset(0)
  call pltend()

  stop
   
end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function forcebalance(omega0,omega1,frac,omega2)
  use shiftmode
  implicit none
  real, intent(in) :: omega0,omega1,frac
  real, intent(out) :: omega2
  omega2=(1-frac)*omega0+frac*omega1
  omega=complex(0.,omega2)
  call FpVyint()
  call FtVyint()
  forcebalance=real(Ftraptotal)+real(Fpasstotal)
  write(*,*)frac,omega2,forcebalance
end function forcebalance

include 'bisectroot.f' 
