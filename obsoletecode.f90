  subroutine dentcalc ! Integrate over v_parallel to get passing n-tilde.
    ! This is the outermost integral. It calls inner integrals (dvy(dtau)).
    ! On the way, calculate the f(x,E) array. Only non-adiabatic.
    ! Also calculate the adiabatic density perturbation, denad.
    ! So far only positive velocity direction. 
    vinfmax=sqrt(2.*Emax)      ! The uppermost orbit
    dvinf=vinfmax/ne           ! Use steps of constant vinf spacing.
    denad=0.                   ! adiabatic density 
    dent=0.                    ! n-tilde
    vinf=0.
    priorcontrib=0.
    priorda=0.
    do i=1,ne                  ! Integrate over ne orbits (energies)
       vinf=i*dvinf
       E(i)=vinf**2/2.
       fe0(i)=exp(-E(i)) ! Unit height Maxwellian with unit temperature.
       fe0de(i)=-fe0(i)  ! Derivative wrt E is same.
       ! ftildecalc calls tauxcalc which also calculates v(nx) for this vinf.
       call ftildecalc(vinf,fe0(i),fe0de(i))
       fte(:,i)=(omega*ft1int+k*(-ft2int+ft3int))*sqm1  ! tilde-f_x(x,E)
       if(.false.)then
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
       endif
       denttemp=fte(:,i)*vinf/v  
       dent=dent+0.5*(priorcontrib+denttemp)*dvinf  ! Non-adiabatic
       priorcontrib=denttemp
! (-dphi/dx=\hat\phi)*(df/dE dv) for \Delta=1. ! Adiabatic not needed.
       denadtemp=-phiprime*fe0de(i)*vinf/v
       denad=denad+0.5*(priorda+denadtemp)*dvinf 
       priorda=denadtemp
! end of adiabatic section.
    enddo
  end subroutine dentcalc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine passingdiags
    ! Main passing calculation  
    call ftildecalc(vinf,1.,1.) ! Separatrix ft for plotting.
    fte(:,0)=(omega*ft1int+k*(-ft2int+ft3int))*sqm1  ! tilde-f_x(x,E)
    call dentcalc
    call plotfte
    call plotdent
    priorcontrib=dent
    call dentcalc2
    write(*,*)'Histogrammed calculation'
    call multiframe(2,1,2)
    call autoplot(x,real(dent),nx)
    call axlabels('x','Real(tilde-n)')
    call color(3)
    call polyline(x,real(priorcontrib),nx)
    call color(15)
    call autoplot(x,imag(dent),nx)
    call axlabels('x','Imag(tilde-n)')
    call color(3)
    call polyline(x,imag(priorcontrib),nx)
    call color(15)
    call pltend()
  end subroutine passingdiags
