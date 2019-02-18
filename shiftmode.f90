!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module shiftmode
  ! We use x in place of z, because x is real.
!  integer, parameter :: nx=50, ne=100, nvy=50  !standard
!  integer, parameter :: nx=20, ne=20, nvy=20   !low resolution
!  integer, parameter :: nx=50, ne=400, nvy=30  ! fcontko highres
  integer, parameter :: nx=300, ne=800, nvy=30  ! Scratch
  real, parameter :: pi=3.1415926, sq2pi=sqrt(2.*3.1415926)
  real :: psi=.1,pL=4.,k=.01, Ty=1.         ! psi, sech4width, k, Ty
  real :: xL=20.,Emax=4.,vymnorm=4.,vymax   ! Hole length, Energy, v_y
  real :: vdrift=0.                         ! Shift of f(v) in hole frame.
  real :: beta                              ! inverse hole parallel temp.
  complex :: omega=(0.0,.01)                ! complex frequency
  integer :: idebug=0
  complex :: omegad,sqm1=(0.,1.),Ftraptotal,Fpasstotal,Fsum
  ! Position arrays
  real :: dx
  real :: x(nx),phi(nx),phiprime(nx),tau(nx),v(nx)
  real :: xt(nx),phit(nx),phitprime(nx)      ! Values in trapped region
  complex :: Lt(nx),CapPhi(nx)               ! tilde L,CapPhi
  ! Energy arrays
  integer :: iwpow=2
  real :: de,E(ne),vinfarray(ne),dvinf,fe0(ne),fe0de(ne) ! Passing orbits
  real :: vinfneg(ne),fe0neg(ne),fe0deneg(ne)! Negative velocity passing.
  real :: Eplus(ne),Eminus(ne)
  real :: Wt(ne),Wtscaled(ne),vpsiarray(ne)  ! Trapped orbit energies etc.
  real :: xlen(ne),tbe(ne)                   ! Trapped orbit length, period.
  complex :: fte(nx,0:ne)
  complex :: Ftrap(ne)
  ! vy arrays
  real :: vy(nvy),fy(nvy),fywy(nvy),dtaumax
  real :: dvy
  complex :: Ftrapvy(nvy),Fpassvy(nvy)
  ! Magnetic field information
  real :: Omegac=0.
  integer :: nharmonics
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initialize
    ! x-arrays
    dx=xL*2./(nx-1.)                            ! Step size in x
    x=dx*(/(i-1,i=1,nx)/)-xL                    ! x-position array
    phi=psi/cosh(x/pL)**4                       ! Potential
    phiprime=-psi*sinh(x/pL)/cosh(x/pL)**5*4/pL ! phi x-gradient
    ! vinf arrays
    vinfmax=sqrt(2.*Emax)      ! The uppermost orbit
    dvinf=vinfmax/ne           ! Use steps of constant vinf spacing.
    vinfarray=dvinf*(/(i-0.5,i=1,ne)/)  ! vinf array 
    vinfneg=-vinfarray                  ! negative velocity
    E=vinfarray**2/2.                   ! Energy in hole frame.
    Eplus=(vinfarray-vdrift)**2/2.      ! Energy in Maxwellian rest frame.
    Eminus=(vinfneg-vdrift)**2/2.
    fe0=exp(-Eplus)/sq2pi     ! Normalized Maxwellian with unit temperature.
    fe0de=-fe0*(1-vdrift/vinfarray)     ! Derivative wrt E.
    fe0neg=exp(-Eminus)/sq2pi ! Normalized Maxwellian with unit temperature.
    fe0deneg=-fe0neg*(1-vdrift/vinfneg) ! Derivative wrt E.
!    write(*,'(3f8.4)')(E(i),fe0de(i),fe0deneg(i),i=1,ne)
    ! vpsi arrays
    Wt=-psi*(/((i-0.5)/ne,i=1,ne)/)**iwpow
    Wtscaled=(-Wt)**(1./iwpow)
    vpsiarray=sqrt(2.*(psi+Wt))
    ! vy-arrays
    vymax=vymnorm*sqrt(Ty)
    dvy=vymax*2./(nvy-1.)                       ! vy-step
    vy=dvy*(/(i-1,i=1,nvy)/)-vymax              ! vy array.
    fy=exp(-vy**2/(2.*Ty))/sqrt(2.*pi*Ty)       ! fy Normalized Maxwellian
    fywy=-fy/Ty                                 ! gradient
    dtaumax=0.
    beta=-1.-(15./16.)*sqrt(pi/psi) 
!      beta=beta*(.8-.9/abs(beta))    ! Ad hoc correction density
!      beta=beta*(1.-.75/abs(beta))    ! Ad hoc correction force.
  end subroutine initialize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!Trapped particle routines. Can only be called after initialization.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine dFdvpsidvy(vpsi,trapforce,tb,xlent)
    ! Calculate the past integral of tau vs x and L(tau) and
    ! phut(taut) vs x for a trapped particle of velocity vpsi at peak
    ! of potential. Return trapforce and orbit period tb These
    ! quantities are the differential force density with respect to vy
    ! and vpsi. And so have to be integrated f(vpsi,vy) dvpsi dvy. The
    ! contributions to force are proportional to df/dW_parallel and
    ! df/dW_perp evaluated at these orbit energies, now done externally
    ! xlent returns the orbit length extreme, used for diagnostic purposes.
    real :: vpsi,tb,xlent
    complex :: trapforce

    complex :: exptau,Ltbb2,exptb,exptbb2,sumfactor,Ltfactor
    Wj=vpsi**2/2-psi
    if(Wj.gt.0)then
       write(*,*)'Positive energy in trapped particle code'
       return
    endif

    
    ! New Alternative Determine the starting node istart and remesh.
    call orbitend(Wj,xm)
    xlent=xm

    ! 2.99 makes istart only just inside the orbit.
    dxt=2.*xlent/(nx-2.99)  ! make array run from -xlen-r*dxt to xlen+r*dxt
    xt=dxt*((/(i-1,i=1,nx)/)-(nx-1.)/2.)         ! x-position array
    phit=psi/cosh(xt/pL)**4                        ! Potential
    phitprime=-psi*sinh(xt/pL)/cosh(xt/pL)**5*4/pL ! phi x-gradient
    istart=2  ! The first position that is inside the orbit.
!    write(*,*)phitprime(1:2),phitprime(nx-1:nx)

   
    iend=nx+1-istart
    if(istart.gt.nx/2)stop 'No steps on remeshed trapped orbit'    
    tau=0.
    v=0.
    Lt=0.
    C=0.
    CapPhi=0.
    ! Integrate to get tau and L(tau)
    v(istart)=sqrt(2.*(phit(istart)+Wj))
    tau(istart)=v(istart)/phitprime(istart) ! To get to v at accel phiprime.
    Ltfactor=exp(sqm1*omegad*tau(istart))
    Lt(istart)=-0.5*v(istart)*(1.-Ltfactor)/(omegad)
    do i=istart+1,iend
       v(i)=sqrt(2.*(phit(i)+Wj))
       if(v(i).lt.0.7*vpsi)then    ! Use phiprime to determine step
          dtau=(v(i)-v(i-1))*2./(phitprime(i)+phitprime(i-1))
       else                        ! Use v to determine step. 
          dtau=(xt(i)-xt(i-1))*0.5*(v(i-1)+v(i))/(v(i-1)*v(i))
       endif
       tau(i)=tau(i-1)+dtau
       Ltfactor=exp(sqm1*omegad*dtau) ! Current exponential
       vmean=(v(i)+v(i-1))/2.
!       vmean=2./(1./v(i)+1./v(i-1))  ! makes no difference
       Lt(i)=Ltfactor*Lt(i-1)-vmean*(1.-Ltfactor)/(omegad)   ! New integral
    enddo
    tbb2=tau(iend)+tau(istart)              ! Use symmetry on end point
    tau(iend+1)=tbb2
    tb=2.*tbb2
    Ltfactor=exp(sqm1*omegad*(tbb2-tau(iend))) ! Current exponential
    vmean=v(iend)/2.
    Ltbb2=Ltfactor*Lt(iend)-vmean*(1.-Ltfactor)/(omegad)
    Lt(iend+1)=Ltbb2
! Calculate CapPhi from L(i), and on the way integrate force round the orbit
    exptb=exp(sqm1*omegad*tb)
    exptbb2=exp(sqm1*omegad*tbb2)
    sumfactor=1./(1.-exptb) ! Not now used.
    trapforce=0.
    CapPhi=0.
! Form CapPhi as omegad*Lt, and the contribution to trapforce as
! \int sqm1*phitprime*(omega*dfe ...)*CapPhi  dtau vpsi
    do i=istart,iend
       exptau=exp(sqm1*omegad*tau(i))
       CapPhi(i)=omegad*((1.-exptb)*Lt(i)+exptau*(exptbb2-1.)*Ltbb2)
! The dtau to be applied to this subsequent tau integral. 
       dtau=tau(i)-tau(i-1)
       trapforce=trapforce+sqm1*  &
       (CapPhi(i)*phitprime(i)+CapPhi(i-1)*phitprime(i-1))/2.*dtau*vpsi
    enddo
  end subroutine dFdvpsidvy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine FtEint(Ftotal,dfperpdWperp,fperp)
    ! Integrate over fe (trapped). Wj=vpsi^2/2-psi. So vpsi=sqrt(2(psi+Wj))
    ! We must use cells that fill the Wj range 0 to -psi.
    ! New version puts evaluations at ends of integration steps.
    complex :: Ftotal,dFdvpsi,dFdvprev,exptb,exptbprev
    if(idebug.eq.-2)then
       write(*,*)'FtEint',omegad,beta,dfperpdWperp
       write(*,*)'  vpsi   fe*s2pi            Ftrap       ' &
         ,'       tb     dvpsi    dFtotal'
       call pltinit(-12.,12.,-8.,8.)
       call axis()
    endif
!    iwpow=2         ! Equal W spacing is ipow=1
    Ftotal=0.
    Wjprev=0.
    feprev=1/sqrt(2.*pi)
    dfeprev=-beta*feprev*fperp
    dfeperpprev=feprev*dfperpdWperp
    vpsiprev=sqrt(2.*psi)
    do i=1,ne-1
       Wj=psi*(-(float(i)/ne)**iwpow)
       fe=exp(-beta*Wj)/sqrt(2.*pi) ! Normalized f_\parallel
       dfe=-beta*fe*fperp           ! df_||/dW_||
       dfeperp=fe*dfperpdWperp      ! df/dW_perp
       vpsi=sqrt(2.*(psi+Wj))
       dvpsi=vpsiprev-vpsi
       ! calculate the force dFdvpsi for this vpsi and dvy element:
       call dFdvpsidvy(vpsi,dFdvpsi,tbe(i),xlen(i))
       if(.not.(abs(dFdvpsi).ge.0))then
          write(*,*)'dFdvpsi NaN?',dFdvpsi,Wj,vpsi,dvpsi
          write(*,'(a,8f8.4)')' omega,k,Omegac=',omega,k,Omegac
          stop
       endif
       exptb=exp(sqm1*omegad*tbe(i))
!       write(*,*)i,'exptb',exptb
       if(i.eq.1)then
          dFdvprev=dFdvpsi
          exptbprev=exptb
       endif
       if(idebug.eq.-2)then
          write(*,'(2f8.4,a,2es12.4,a,f8.3,f9.5,es12.4)')vpsi&
               ,fe*sqrt(2.*pi) &
               ,' (',Ftrap(i),')',tbe(i),dvpsi,real(Ftrap(i))*dvpsi
       endif
       ! Ftrap(i) becomes the contribution to F from this step.
       Ftrap(i)=0.5*(dFdvpsi*(omegad*dfe-(omegad-omega)*dfeperp) &
        /(1-exptb)&
        +dFdvprev*(omegad*dfeprev-(omegad-omega)*dfeperpprev) &
        /(1-exptbprev)  &
        )*dvpsi
       ! Multiply by 2. to account for \pm v_\psi.
       Ftotal=Ftotal+2.*Ftrap(i)       ! Add to Ftotal integral.
       vpsiprev=vpsi
       feprev=fe
       dFdvprev=dFdvpsi
       dfeprev=dfe
       dfeperpprev=dfeperp
       exptbprev=exptb
    enddo
    ! Calculate end by extrapolation.
    Ftrap(ne)=Ftrap(ne-1)+0.5*(Ftrap(ne-1)-Ftrap(ne-2))/dvpsi*vpsi
    Ftotal=Ftotal+2.*Ftrap(i)*vpsi
  end subroutine FtEint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine FtEintOld(Ftotal,dfperpdWperp,fperp)
  ! Here we need to iterate over fe (trapped). Wj=vpsi^2/2-psi. So
  ! vpsi=sqrt(2(psi+Wj)) We must use cells that fill the Wj range 0 to
  ! -psi. Old version with evaluations in centers of cells.
    complex :: Ftotal
    Ftotal=0.
    if(idebug.eq.-2)then
       write(*,*)'FtEint',omegad,beta,dfperpdWperp
       write(*,*)'  vpsi   fe*s2pi            Ftrap       ' &
         ,'       tb     dvpsi    dFtotal'
       call pltinit(-12.,12.,-8.,8.)
       call axis()
    endif
!    iwpow=2         ! Equal W spacing is ipow=1
    do i=1,ne
       Wj=psi*(-((i-0.5)/ne)**iwpow)
       fe=exp(-beta*Wj)/sqrt(2.*pi) ! Normalized f_\parallel
       dfe=-beta*fe*fperp           ! df_||/dW_||
!       dfe=-beta*fe           ! df_||/dW_|| old error
       dfeperp=fe*dfperpdWperp      ! df/dW_perp
       vpsi=sqrt(2.*(psi+Wj))
       dvpsi=-sqrt(2.*psi*(1.-(float(i)/ne)**iwpow)) &
            +sqrt(2.*psi*(1.-(float(i-1)/ne)**iwpow))
       ! calculate the force Ftrap for this dvpsi and dvy element:
       call dFdvpsidvy(vpsi,Ftrap(i),tbe(i),xlen(i))
       if(.not.(abs(Ftrap(i)).ge.0))then
          write(*,*)'Ftrap NAN?',Ftrap(i),Wj,vpsi!,dvpsi
          write(*,'(a,8f8.4)')' omega,k,Omegac=',omega,k,Omegac
          stop
       endif
       if(idebug.eq.-2)then
          write(*,'(2f8.4,a,2es12.4,a,f8.3,f9.5,es12.4)')vpsi&
               ,fe*sqrt(2.*pi) &
               ,' (',Ftrap(i),')',tbe(i),dvpsi,real(Ftrap(i))*dvpsi
       endif
       ! Ftrap(i) becomes the contribution to F from this energy.
       Ftrap(i)=Ftrap(i)*(omegad*dfe-(omegad-omega)*dfeperp)*dvpsi
       Ftrap(i)=Ftrap(i)/(1-exp(sqm1*omegad*tbe(i))) ! From subfactor
       ! Multiply by 2. to account for \pm v_\psi.
       Ftotal=Ftotal+2.*Ftrap(i)       ! Add to Ftotal integral.
    enddo
  end subroutine FtEintOld
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine FtVyint()
    ! Integrate (histogram) over vy to get Ftraptotal.
    Ftraptotal=0.
    do i=1,nvy
       omegad=omega-k*vy(i)
       call FtEint(Ftrapvy(i),fywy(i),fy(i))
       Ftraptotal=Ftraptotal+Ftrapvy(i)*dvy
       if(idebug.ne.0)write(*,'(a,i4,es10.2,a,f8.4,2e12.4)')&
            'i,R(omegad)',i,real(omegad),  &
            '  vy,R(Ftrapvy)',vy(i),real(Ftrapvy(i)),real(Ftraptotal)
    enddo
  end subroutine FtVyint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Passing particle routines (New integration order)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dFdvinfdvy(vinf,passforce)  
  ! Calculate the past integral of tau vs x and L(tau) and
  ! CapPhi(taut) vs x for passing particles of velocity vinf.
  ! Return the differential force density with respect to vy
  ! and vinf. This has to be integrated f(vinf,vy) dvinf dvy. The
  ! contributions to force are proportional to df/dW_parallel and
  ! df/dW_perp evaluated at these orbit energies but that scaling 
  ! is now done externally. vinf and v are positive in this calculation.
  ! So one should pass it vinfarray NOT vinfneg.   
    real    :: vinf
    complex :: passforce
    complex :: expdtau
    tau(1)=0.
    thisdtaumax=0.
    v(1)=sqrt(vinf**2+2.*phi(1))
    Lt(1)=0.                                          ! Lt integrated
    CapPhi(1)=v(1)-vinf
    CapPhi(1)=0.                                       ! More consistent.
    passforce=0.
! Do \int_-\infty^t exp(-i*omegad*tau) dtau, passing: Lt(i). And CapPhi(i)
! And Force integration.
    do i=2,nx
       v(i)=sqrt(vinf**2+2.*phi(i))
       dtau=dx*0.5*(v(i-1)+v(i))/(v(i-1)*v(i))        ! tau increment varies
       if(dtau.gt.thisdtaumax)thisdtaumax=dtau
       tau(i)=tau(i-1)+dtau
       expdtau=exp(sqm1*omegad*dtau)
       vmean=(v(i)+v(i-1))/2.-vinf
       Lt(i)=expdtau*Lt(i-1)-vmean*(1.-expdtau)/(omegad)   ! New integral
       CapPhi(i)=omegad*Lt(i)
       ! Histogram version.
       passforce=passforce+sqm1*(CapPhi(i)*phiprime(i))*(vinf/v(i))*dx
    enddo
    if(idebug.gt.1)then
       write(*,'(a,2es12.4,a,f8.4,a,2es12.4)')'omegad',omegad,' vin',vinf,' CapPhi(n/2)',CapPhi(nx/2)
    endif
  end subroutine dFdvinfdvy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine FpEint(Ftotal,dfperpdWperp,fperp)
    ! Integrate over vinf (v_parallel) to get passing force for specified v_y,
    ! passing particles.
    complex :: Ftotal
    real :: dfperpdWperp,fperp
    real :: dfe,dfeperp
    complex :: passforce
    Ftotal=0.
    do i=1,ne                  ! Integrate over ne orbits (energies)
       dfe=fe0de(i)*fperp
       dfeperp=fe0(i)*dfperpdWperp
       call dFdvinfdvy(vinfarray(i),passforce)  
       ! No factor q_e has been applied to CapPhi, so no q_e factor here.
       Ftotal=Ftotal+passforce*dvinf*(omegad*dfe-(omegad-omega)*dfeperp) 
    enddo
    if(vdrift.eq.0.)then
! For stationary holes negative velocities just double the effective force
       Ftotal=2.*Ftotal
    else                   ! Need to integrate over negative velocity.
       do i=1,ne
          dfe=fe0deneg(i)*fperp
          dfeperp=fe0neg(i)*dfperpdWperp
          call dFdvinfdvy(vinfarray(i),passforce)  ! Pass postive vinf.
          Ftotal=Ftotal+passforce*dvinf*(omegad*dfe-(omegad-omega)*dfeperp) 
       enddo
    endif
  end subroutine FpEint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine FpVyint()
    ! Integrate (histogram) over vy to get total force.
    Fpasstotal=0.
    do i=1,nvy
       omegad=omega-k*vy(i)
       call FpEint(Fpassvy(i),fywy(i),fy(i))
       Fpasstotal=Fpasstotal+Fpassvy(i)*dvy
       if(idebug.ne.0)write(*,'(a,i4,es10.2,a,f8.4,2e12.4)')&
            'i,R(omegad)',i,real(omegad),  &
            '  vy,R(Fpassvy)',vy(i),real(Fpassvy(i)),real(Fpasstotal)
    enddo
  end subroutine FpVyint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routines for magnetized plasma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine SumHarmonics()
    implicit none
    ! Sum harmonic contributions to obtain Forces for specified
    ! k, Omegac. Does the job of FpVyint, FtVyint. 
    real :: EIm(0:nvy),xit
    integer :: m,ncalc,ifirst
    data ifirst/0/
! How many harmonics do we need? Regard vymax as the velocity relative
! to the thermal perpendicular speed. But don't allow less than +-4.
    nharmonics=9999    ! Default too large.
    if(Omegac.gt.0)nharmonics=max(4,nint(abs(k*vymax)/Omegac))
    if(.not.nharmonics.le.nvy-1)then   ! B-Field too low.
       if(mod(ifirst,20).eq.0)then
          write(*,'(a,i4,a,f6.4,a)')'Too many magnetized harmonics', &
               nharmonics,' Ignoring Omega_c=',Omegac,' Zero-B calculation.'
       endif
       ifirst=ifirst+1
       call FpVyint
       call FtVyint
       return
    endif
! Truly magnetized case:
    xit=k*sqrt(Ty)/Omegac
! Calculate the Integer[0.] exp*I[2] Bessel functions 0 to nharmonics
    call RIBESL(xit**2,0.,nharmonics+1,2,EIm,ncalc)
    if(.not.ncalc.eq.nharmonics+1)then ! All orders not calculated correctly.  
       write(*,'(a,i3,a,i3,a)')'Bessel functions',ncalc,' of',nharmonics+1, &
       ' are precise. Others may be irrelevant.' 
    endif
! m=0 always used.! fy is Maxwellian hence the fywy,fy.
    omegad=omega
    call FpEint(Fpassvy(1),-1./Ty,1.)
    call FtEint(Ftrapvy(1),-1./Ty,1.)
    Fpasstotal=Fpassvy(1)*EIm(0)
    Ftraptotal=Ftrapvy(1)*EIm(0)
    do m=1,nharmonics
       omegad=omega+m*Omegac
       call FpEint(Fpassvy(m+1),-1./Ty,1.)
       call FtEint(Ftrapvy(m+1),-1./Ty,1.)
!       if(.false.)then
       if(real(omega).eq.0)then   !Short cut.
          Fpasstotal=Fpasstotal+2*real(Fpassvy(m+1))*EIm(m)
          Ftraptotal=Ftraptotal+2*real(Ftrapvy(m+1))*EIm(m)
       else      ! Full sum over plus and minus m.
          Fpasstotal=Fpasstotal+Fpassvy(m+1)*EIm(m)
          Ftraptotal=Ftraptotal+Ftrapvy(m+1)*EIm(m)
          omegad=omega-m*Omegac
          call FpEint(Fpassvy(m+1),-1./Ty,1.)
          call FtEint(Ftrapvy(m+1),-1./Ty,1.)
          Fpasstotal=Fpasstotal+Fpassvy(m+1)*EIm(m)
          Ftraptotal=Ftraptotal+Ftrapvy(m+1)*EIm(m)
       endif
    enddo
  end subroutine SumHarmonics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine orbitend(Wj,x0)
    ! Find by bisection the turning point of the orbit whose energy is Wj.
    ! That is where phi=-Wj. But we return the value below it: x0.
    real :: Wj,xm,x1,x0,phim,phi0,phi1
    nbi=20
    x0=0
    x1=xL
    phi0=psi/cosh(x0/pL)**4+Wj
    phi1=psi/cosh(x1/pL)**4+Wj
    if(sign(1.,phi1).eq.sign(1.,phi0))then
       write(*,*)'orbitend energies do not bracket zero',phi0,phi1
       stop
    endif
    do i=1,nbi
       xm=(x1+x0)/2.
       phim=psi/cosh(xm/pL)**4+Wj
! Which value to replace.
       if(sign(1.,phim).eq.sign(1.,phi0))then
          if(x0.eq.xm)exit ! should never happen
          phi0=phim
          x0=xm
       else
          if(x1.eq.xm)exit
          phi1=phim
          x1=xm
       endif
    enddo
  end subroutine orbitend
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module shiftmode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
