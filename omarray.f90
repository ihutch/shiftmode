! Given psi and k calculate the forces for a range of omegas.

!include 'fhgfunc.f'

program main
  use shiftmode 
  integer, parameter ::   nk=11
  real :: kik(nk),omik(nk)
  complex :: Fcpassing(nk),Ftrapped(nk)
  integer, parameter :: np=20
!  integer, parameter :: np=2
  real, dimension(np,nk) :: Ftnp,Fpnp,omi
  character*10 :: string
  real :: omega0,omega1,oval
  external forcebalance
  psi=.02
        call initialize
! k-scan
  akmax=0.3*sqrt(psi)
  akmin=akmax/200.
  write(*,*)'   k      omegai   Fpassing                Ftrapped'
  do ik=1,nk
     k=akmin+(ik-1.)*(akmax-akmin)/(max(nk-1.,1.))
     kik(ik)=k
     omegaimax=1.3*min(k,akmax/2.)
     omik(ik)=0.
     omega0=akmin
     omega1=2.*k
     nbi=10
     call bisectroot(omega0,omega1,forcebalance,nbi,frac,oval)
     write(*,*)'oval=',oval
     omega=complex(0.,oval)
     so=abs(omega**2)
!     Ftnp(ip,ik)=real(Ftraptotal)/so
!     Fpnp(ip,ik)=real(Fcpassing(ik))/so

     if(.true.)then
     do ip=1,np  ! Iterate over omega.
        omega=(0.,0.)+sqm1*(np+1-ip)*omegaimax/np
        so=abs(omega**2)
!        write(*,*)'nx, ne, nvy,   xL,    pL,   omegar,  omegai,     k     psi   beta'
!        write(*,'(3i4,7f8.4)')nx,ne,nvy,xL,pL,real(omega),imag(omega),k,psi,beta
        call dentcalc2()
        call FtVyint()
        Fcpassing(ik)=Fpasstotal
        Ftrapped(ik)=Ftraptotal
        write(*,'(f6.4,5es12.3)')k,imag(omega),Fcpassing(ik),Ftraptotal
        omi(ip,ik)=imag(omega)
        Ftnp(ip,ik)=real(Ftraptotal)/so
        Fpnp(ip,ik)=real(Fcpassing(ik))/so
        if(ip.gt.1)then
           if((Ftnp(ip,ik)+Fpnp(ip,ik))  &
                *(Ftnp(ip-1,ik)+Fpnp(ip-1,ik)).le.0)then
              omik(ik)=omi(ip,ik)
              write(*,*)ip,Ftnp(ip,ik),Ftnp(ip-1,ik)
           endif
        endif
     enddo
     endif
  enddo

!  write(*,*)'omik',omik
!  write(*,*)'kik ',kik

  ik=nk
  call pfset(3)
  call multiframe(2,1,2)
  call pltinit(0.,omi(1,nk),-.05,.05)
  call axis()
  call axlabels('omega','Normalized Force')
  call winset(.true.)
  do ik=nk,1,-1
     call color(ik)
     call polyline(omi(:,ik),Ftnp(:,ik),np)
     call polyline(omi(:,ik),Fpnp(:,ik),np)
  enddo
  call color(15)

  call pltinit(0.,omi(1,nk),-1.,1.)
  call axis()
  call axlabels('omega','Normalized Force')
  call winset(.true.)
  do ik=nk,1,-1
     call color(ik)
     call polyline(omi(:,ik),Ftnp(:,ik),np)
     call polyline(omi(:,ik),Fpnp(:,ik),np)
  enddo
  call pltend()
  call multiframe(0,0,0)

  call autoplot(kik,omik,nk)
  call pltend()

  stop
  call pfset(3)
  call pltinit(0.,np*psistep,-np*psistep*1.3,np*psistep*4.)
  call axis
  call axlabels('!Ay!@','Force !BF!dt!d , F!dp!d!@ (/!Aw!@!u2!u)')
  do ik=1,nk
     call polyline(psinp,Ftnp(:,ik),np)
     call polyline(psinp,Fpnp(:,ik),np)
     if(ik.eq.1)then
        call jdrwstr(wx2nx(psinp(np/2)),wy2ny(Ftnp(np/2,ik)),'Trapped',-1.)
     elseif(ik.eq.2)then
        call jdrwstr(wx2nx(psinp(np/2)),wy2ny(Fpnp(np/2,ik)),'Passing',-1.)
     endif
     call fwrite(kik(ik),iwdth,4,string)
     call jdrwstr(wx2nx(psinp(2*np/3)),wy2ny(Ftnp(2*np/3,ik)),string,1.)
  enddo
   
end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function forcebalance(omega0,omega1,frac,omega2)
  use shiftmode
  implicit none
  real, intent(in) :: omega0,omega1,frac
  real, intent(out) :: omega2
  omega2=(1-frac)*omega0+frac*omega1
  omega=complex(0.,omega2)
  call dentcalc2()  
  call FtVyint()
  forcebalance=real(Ftraptotal)+real(Fpasstotal)
  write(*,*)frac,omega2,forcebalance
end function forcebalance

include 'bisectroot.f' 
