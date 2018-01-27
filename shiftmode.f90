!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module shiftmode
  ! We use x in place of z, because x is real.
  integer, parameter :: nx=100, ne=50, nvy=50
  real :: xL=29.,Emax=2.,vymax=4.
  real :: psi=.1,pL=4.,k=.0, Ty=1.
  complex :: omega=(0.,.01)
  integer :: idebug=1
  ! Position arrays
  real :: dx
  real :: x(nx),phi(nx),phiprime(nx),tau(nx),v(nx),work(nx)
  complex :: omegad,sqm1=(0.,1.)
  complex :: Lt(nx),phiut(nx),dent(nx),denttemp(nx),priorcontrib(nx)
  complex :: ft1(nx),ft2(nx),ft3(nx),ft1int(nx),ft2int(nx),ft3int(nx)
  ! Energy arrays
  real :: de,E(ne),fe0(ne),fe0de(ne)
  complex :: fte(nx,ne)
  ! vy arrays
  real :: vy(nvy),fy(nvy),fywy(nvy)
  real :: dvy
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initialize
    ! x-arrays
    dx=xL*2./(nx-1.)
    x=dx*(/(i-1,i=1,nx)/)-xL
    phi=psi/cosh(x/pL)
    phiprime=-psi*sinh(x/pL)/cosh(x/pL)**2
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
    scalefactor=1.e10
    taustep=alog(scalefactor)/max(imag(omegad),1.e-5)
    tau(1)=0.
    v(1)=sqrt(vinf**2+2.*phi(1))
    Ltint=(v(1)-vinf)*sqm1
    phiut(1)=v(1)-vinf
    ! \int_0^tau exp(-i*omegad*tau) dtau to begin with.  
    do i=2,nx
       v(i)=sqrt(vinf**2+2.*phi(i))
       tau(i)=tau(i-1)+dx*(v(i-1)+v(i))/(v(i-1)*v(i))
       dtau=(tau(i)-tau(i-1))
       if(tau(i).gt.taustep)then  !Rescale to avoid overflow.
          tau(i)=tau(i)-taustep
          Ltint=Ltint*exp(sqm1*omegad*taustep)
!          write(*,*)'scaletau',i,taustep,tau(i),Ltint
       endif
       Lt(i)=Ltint
       Ltint=(v(i)-vinf)*exp(-sqm1*omegad*tau(i))*sqm1
       Lt(i)=(Lt(i)+Ltint)/2. *dtau
       ! Add accel term and adjust tau at end to t=0.
       phiut(i)=v(i)-vinf+ sqm1*Lt(i)*exp(sqm1*tau(i)*omegad)
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
  subroutine dentcalc ! Integrate over v_parallel to get passing n-tilde.
    ! On the way, calculate the f(x,E) array. Only non-adiabatic.
    vinfmax=sqrt(2.*Emax)
    dvinf=vinfmax/ne
    do i=1,ne
       vinf=i*dvinf
       E(i)=vinf**2/2.
       fe0(i)=exp(-E(i)) ! Unit height Maxwellian with unit temperature.
       fe0de(i)=fe0(i)
       call ftcalc(vinf,fe0(i),fe0de(i))
       ! ftcalc calls tauxcalc which calculates v(nx) for this vinf.
       ! we might wish to store it.
       fte(:,i)=(omega*ft1int+k*(-ft2int+ft3int))*sqm1
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
!       if(i.eq.1)priorcontrib=0.
       if(i.eq.1)priorcontrib=fte(:,i)*vinf/v
       denttemp=fte(:,i)*vinf/v
       dent=dent+0.5*(priorcontrib+denttemp)*dvinf
       priorcontrib=denttemp       
    enddo
    ! Must multiply by m_e\Delta to get the non-adiabatic n-tilde.
    ! Then integrate times d\phi/dx to get the momentum.
  end subroutine dentcalc
  !--------------------------------------------
  subroutine plotfte
    call multiframe(2,1,2)
    call autoplot(x,real(fte(:,1)),nx)
    call axlabels('','Real(fte)')
    do i=2,ne
       call color(mod(i-2,14)+1)
       call polyline(x,real(fte(:,i)),nx)
    enddo
    call color(15)
    call autoplot(x,imag(fte(:,1)),nx)
    call axlabels('x','Imag(fte)')
    do i=2,ne
       call color(mod(i-2,14)+1)
       call polyline(x,imag(fte(:,i)),nx)
    enddo
    call pltend()
  end subroutine plotfte
  !--------------------------------------------
  subroutine plotdent  
    call multiframe(3,1,0)
    call autoplot(x,real(dent),nx)
    call axis2()
    call axlabels('','!ER!@(n-tilde)')
    call autoplot(x,imag(dent),nx)
    call axis2()
    call axlabels('','!EI!@(n-tilde)')
    call autoplot(x,phi,nx)
    call winset(.true.)
    call dashset(2)
    call polyline(x,phiprime+psi/2.,nx)
    call winset(.false.)
    call axis2()
    call axlabels('!Bx!@','!Af!@')
    call pltend()
  end subroutine plotdent
  !--------------------------------------------
end module shiftmode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
  use shiftmode 
  write(*,*)'nx, ne, nvy,   xL,    pL,   omegar,  omegai,     k'
  write(*,'(3i4,5f8.4)')nx,ne,nvy,xL,pL,real(omega),imag(omega),k
  call initialize
  if(idebug.gt.0)then
     call tauxcalc(.2,(.2,0.))
     write(*,*)'phi'
     write(*,'(10f8.4)')phi
     write(*,*)'v'
     write(*,'(10f8.4)')v
     write(*,*)'tau'
     write(*,'(10f8.3)')tau
     write(*,*)'real(Lt)'
     write(*,'(10f8.4)')real(Lt)
     write(*,*)'fy'
     write(*,'(10f8.4)')fy
  endif
  
  call dentcalc
  call plotfte
  call plotdent

end program main

