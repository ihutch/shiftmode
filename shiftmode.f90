!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module shiftmode
  ! We use x in place of z, because x is real.
  integer, parameter :: nx=20, ne=20, nvy=20
  ! Position arrays
  real :: dx,k=1
  real :: x(nx),phi(nx),tau(nx),v(nx)
  complex :: omega=(0.,0.2),omegad,sqm1=(0.,1.)
  complex :: Lt(nx),phiut(nx),dent(nx),denttemp(nx),priorcontrib(nx)
  complex :: ft1(nx),ft2(nx),ft3(nx),ft1int(nx),ft2int(nx),ft3int(nx)
  ! Energy arrays
  real :: Emax=4
  real :: de,E(ne),fe0(ne),fe0de(ne)
  complex :: fte(nx,ne)
  ! vy arrays
  real :: vy(nvy),fy(nvy),fywy(nvy)
  real :: Ty=1.,vymax=4.
  real :: dvy,xL=30.,psi=.5,pL=4.
  integer :: idebug=2
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initialize
    ! x-arrays
    dx=xL*2./(nx-1.)
    x=dx*(/(i-1,i=1,nx)/)-xL
    phi=psi/cosh(x/pL)
    ! vy-arrays
    dvy=vymax*2./(nvy-1.)
    vy=dvy*(/(i-1,i=1,nvy)/)-vymax
    fy=exp(-vy**2/(2.*Ty))
    fywy=-fy/Ty
  end subroutine initialize
  !--------------------------------------------
  subroutine tauxcalc(vinf,omegad)
    ! integrate dx/v to get v(x), tau(x) and Lt(x,t=0)
    real vinf
    complex omegad,Ltint
    tau(1)=0.
    v(1)=sqrt(vinf**2+2.*phi(1))
    Ltint=(v(1)-vinf)*sqm1
    phiut(1)=v(1)-vinf
!    phiut(1)=0.
!    write(*,*)v(1),vinf,phiut(1)
    ! \int_0^tau exp(-i*omegad*tau) dtau to begin with.  
    do i=2,nx
       v(i)=sqrt(vinf**2+2.*phi(i))
       tau(i)=tau(i-1)+dx*(v(i-1)+v(i))/(v(i-1)*v(i))
       Lt(i)=Ltint
       Ltint=(v(i)-vinf)*exp(-sqm1*omegad*tau(i))*sqm1
       Lt(i)=(Lt(i)+Ltint)/2. *(tau(i)-tau(i-1))
       ! Add accel term and adjust tau at end to t=0.
       phiut(i)=v(i)-vinf+ Lt(i)*exp(sqm1*tau(i)*omegad)
    enddo

  end subroutine tauxcalc
  !--------------------------------------------
  subroutine ftcalc(vinf,fe,dfe) ! Integrate over vy to get ft-parallel
    ft1int=0. ; ft2int=0. ; ft3int=0.
    ! Integrate (histogram) over vy to get the three ft terms.
    ! Total fte is omega*ft1int-k*ft2int+k*ft3int. Done in dentcalc.
    do i=1,nvy
       omegad=omega-k*vy(i)
       call tauxcalc(vinf,omegad) ! Integrate along past orbit.
       ! Maybe synthesize phiut explicitly here.
       ! Maybe store things as a function of vy in order to reuse.
       ft1=dfe*fy(i)*phiut
       ft2=dfe*fy(i)*vy(i)*phiut
       ft3=fe*fy(i)*vy(i)*phiut
       ft1int=ft1int+ft1*dvy
       ft2int=ft2int+ft2*dvy
       ft3int=ft3int+ft3*dvy
    enddo
  end subroutine ftcalc
  !--------------------------------------------
  subroutine dentcalc ! Integrate over energy to get passing n-tilde.
    ! On the way, calculate the f(x,E) array. Only non-adiabatic.
    vinfmax=sqrt(2.*Emax)
    dvinf=vinfmax/ne
    priorcontrib=0.
    do i=1,ne
       vinf=i*dvinf
       E(i)=vinf**2/2.
       fe0(i)=exp(-E(i)) ! Unit height Maxwellian with unit temperature.
       fe0de(i)=fe0(i)
       call ftcalc(vinf,fe0(i),fe0de(i))
       ! ftcalc calls tauxcalc which calculates v(nx) for this vinf.
       fte(:,i)=omega*ft1int+k*(-ft2int+ft3int)*sqm1
       if(idebug.gt.1)write(*,'(i4,10f8.4)')i,E(i),vinf,fe0(i),omega,k

       if(i.eq.1.and.idebug.gt.1)then
          write(*,*)'  i    E(i)   vinf   fe0(i)  omegar  omegai    k'
          write(*,*)'ft1int'
          write(*,'(10f8.4)')ft1int
          write(*,*)'ft2int'
          write(*,'(10f8.4)')ft2int
          write(*,*)'ft3int'
          write(*,'(10f8.4)')ft3int
          write(*,*)'fte(:,i)'
          write(*,'(10f8.4)')fte(:,i)
          write(*,*)'phiut'
          write(*,'(10f8.4)')phiut
       endif
       denttemp=fte(:,i)*vinf/v
       dent=dent+(priorcontrib+denttemp)*dvinf
       priorcontrib=denttemp       
    enddo
    ! Must multiply by m_e\Delta to get the non-adiabatic n-tilde.
  end subroutine dentcalc

end module shiftmode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
  use shiftmode 
  write(*,*)'nx, ne, nvy,   xL,    pL,   omegar,  omegai,     k'
  write(*,'(3i4,5f8.4)')nx,ne,nvy,xL,pL,real(omega),imag(omega),k
  call initialize
  call tauxcalc(.2,(.2,0.))
  if(idebug.gt.0)then
     write(*,*)'phi'
     write(*,'(10f8.4)')phi
     write(*,*)'v'
     write(*,'(10f8.4)')v
     write(*,*)'tau'
     write(*,'(10f8.3)')tau
     write(*,*)'Lt'
     write(*,'(10f8.4)')real(Lt)
     write(*,*)'fy'
     write(*,'(10f8.4)')fy
  endif
  
  call dentcalc
  write(*,*)'dent'
  write(*,'(10f8.4)')dent

  call multiframe(2,1,2)
  call autoplot(x,real(dent),nx)
  call boxtitle('dent real')
  call autoplot(x,imag(dent),nx)
  call pltend()
  
  
end program main

! A problem when tau*omegad is large arises for the slowest passing
! particles leading to a NAN. We need a systematic way to fix this.
! I think the way is to rescale the tau_0 when tau becomes too large.
