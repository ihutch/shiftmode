!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine passingdiags
    ! Main passing calculation  
    call ftildecalc(vinf,1.,1.) ! Separatrix ft for plotting.
    fte(:,0)=(omega*ft1int+k*(-ft2int+ft3int))*sqm1  ! tilde-f_x(x,E)
    call FpVyint()
    call plotfte
    call plotdent
    priorcontrib=dent
    call FpVyint()
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

  !--------------------------------------------
  subroutine plotfte
    call multiframe(2,1,2)
    call autoplot(x,real(fte(:,0)),nx)
    call axlabels('','Real(fte)')
    do i=1,ne
!       call color(mod(i-2,14)+1)
       call color(i)
       call polyline(x,real(fte(:,i)),nx)
    enddo
    call color(15)
    call autoplot(x,imag(fte(:,0)),nx)
    call axlabels('x','Imag(fte)')
    do i=1,ne
!       call color(mod(i-2,14)+1)
       call color(i)
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
!    call autoplot(x,imag(dent),nx)
!    call axis2()
    !    call axlabels('','!EI!@(n-tilde)')
    call autoplot(x,denad,nx)
    call axis2()
    call axlabels('','Adiabatic n')
    call winset(.true.)
    call dashset(2)
!    call polyline(x,real(dent)+denad,nx)
    call winset(.false.)
    call dashset(0)
    call autoplot(x,phi,nx)
    call axis2()
    call axlabels('!Bx!@','!Af!@')
    call pltend()
  end subroutine plotdent
  !--------------------------------------------
