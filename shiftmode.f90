!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module shiftmode
  ! We use x in place of z, because x is real.
!  integer, parameter :: nx=50, ne=100, nvy=50  !standard
!  integer, parameter :: nx=20, ne=20, nvy=20   !low resolution
!  integer, parameter :: nx=50, ne=400, nvy=30  ! fcontko highres
  integer, parameter :: nx=100, ne=800, nvy=30  ! Low omegai
!  integer, parameter :: nx=300, ne=400, nvy=30  ! Low omegai
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
  complex :: PlotCapPhi(nx)
  integer :: istart,iend
  ! Energy arrays
  integer :: iwpow=2
  real :: de,E(ne),vinfarray(ne),dvinf,fe0(ne),fe0de(ne) ! Passing orbits
  real :: vinfneg(ne),fe0neg(ne),fe0deneg(ne)! Negative velocity passing.
  real :: Eplus(ne),Eminus(ne)
  real :: Wt(ne),Wtscaled(ne),vpsiarray(ne)  ! Trapped orbit energies etc.
  real :: xlen(ne),tbe(ne)                   ! Trapped orbit length, period.
  complex :: fte(nx,0:ne)
  complex :: Ftrap(ne)
  complex :: omegab(0:ne),Fnonres(0:ne)
  ! vy arrays
  real :: vy(nvy),fy(nvy),fywy(nvy)
!  real :: dtaumax  ! seems unused
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
    ! vinf passing arrays
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
    ! vpsi trapped arrays
    Wt=-psi*(/(float(i)/ne,i=1,ne)/)**iwpow
    Wtscaled=(-Wt)**(1./iwpow)
    vpsiarray=sqrt(2.*(psi+Wt))
    beta=-1.-(15./16.)*sqrt(pi/psi) 
    ! vy-arrays
    vymax=vymnorm*sqrt(Ty)
    dvy=vymax*2./(nvy-1.)                       ! vy-step
    vy=dvy*(/(i-1,i=1,nvy)/)-vymax              ! vy array.
    fy=exp(-vy**2/(2.*Ty))/sqrt(2.*pi*Ty)       ! fy Normalized Maxwellian
    fywy=-fy/Ty                                 ! gradient
!    dtaumax=0.
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
    if(.true.)then !Rescale the x to put more points near the end
       xt=sign(1.,xt)*xlent*(1.-(abs(abs(xt/xlent)-1.))**2)
       xt(1)=2*xt(2)-xt(3)
       xt(nx)=2*xt(nx-1)-xt(nx-2)
    endif
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
! This is the integral 0-tb: giving CapPhi:
       CapPhi(i)=omegad*((1.-exptb)*Lt(i)+exptau*(exptbb2-1.)*Ltbb2)
! The dtau to be applied to this subsequent tau integral. 
       dtau=tau(i)-tau(i-1)
! This gives the integral 0-tb/2 of CapPhi without resonant muliplier
! or the f' terms so it needs to be multiplied by (2)(res-mult)(f'terms):
       trapforce=trapforce+sqm1*  &
       (CapPhi(i)*phitprime(i)+CapPhi(i-1)*phitprime(i-1))/2.*dtau*vpsi
    enddo
  end subroutine dFdvpsidvy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine FtEint(Ftotal,dfperpdWperp,fperp)
    ! Integrate over fe (trapped). Wj=vpsi^2/2-psi. So vpsi=sqrt(2(psi+Wj))
    ! We must use cells that fill the Wj range 0 to -psi.
    ! New version puts evaluations at ends of integration steps.
    complex :: Ftotal,dFdvpsi,dFdvprev,exptb,exptbprev,cdvpsi,dob
    real :: obi
    character string*100
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
!    dfeprev=-beta*feprev*fperp  ! Old Schamel approx
    dfeprev=(15./(16.*sqrt(2.*psi))-1./sqrt(2.*pi))*fperp ! New exact analytic
    dfeperpprev=feprev*dfperpdWperp
    vpsiprev=sqrt(2.*psi)
    omegab(0)=0.
    Fnonres(0)=0.                 !Don't add zero energy point.
    exptbprev=0.                  !Silence warnings
    do i=1,ne-1
       Wj=psi*(-(float(i)/ne)**iwpow)
       ! Old Schamel approximate form.
!       fe=exp(-beta*Wj)/sqrt(2.*pi) ! Normalized f_\parallel
!       dfe=-beta*fe*fperp           ! df_||/dW_||
       ! New exact sech^4 hole form.
       sqWj=sqrt(-Wj)
       fe=((2./pi)*sqWj+(15./16.)*Wj/sqrt(psi)+experfcc(sqWj)/sqrt(pi))/sqrt(2.)
       dfe=((15./16.)/sqrt(psi)-experfcc(sqWj)/sqrt(pi))/sqrt(2.)*fperp
       ! End of choice
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
       omegab(i)=2.*pi/tbe(i)
       if(idebug.eq.-3  &
!            .and.mod(i,ne/4).eq.1 &
!            .and.abs(-Wj/real(omega**2)-1.).lt..15 &
            .and.abs(real(omegab(i))/real(omega)-1.).lt..05 &
            .and.omega.eq.omegad)then
          !Diagnostic of CapPhi
          PlotCapPhi=CapPhi  ! Save the Phi_0 for plotting.
          if(i.eq.-10)then
             call autoinit(tau(istart:iend), &
                  real(PlotCapPhi(istart:iend)),iend-istart+1)
             call axis
          endif
          call autoplot(tau(istart:iend), &
               real(PlotCapPhi(istart:iend)),iend-istart+1)
          call legendline(.1,.1,0,' real(!AF!@)')
          call dashset(1)
          call polyline(tau(istart:iend), &
               imag(PlotCapPhi(istart:iend)),iend-istart+1)
          call legendline(.1,.15,0,' imag(!AF!@)')
          call axlabels('tau','!AF!@')
          call axis2
          call dashset(2)
!          call polyline(tau(istart:iend), &
!               imag(PlotCapPhi(istart:iend)),iend-istart+1)
          call dashset(0)
          write(string,'(a,f8.5,a,e10.2,a,f8.5)') &
            'or=',real(omega),' oi=',imag(omega),' ob=',real(omegab(i))
          call boxtitle(string(1:lentrim(string)))
          call polymark(tau(istart:iend),v(istart:iend)/5.,iend-istart+1,1)
          call polymark(tau(istart:iend),xt(istart:iend)/200.,iend-istart+1,2)
          call legendline(.7,.95,1,' v/5')
          call legendline(.7,.9,1,' z/200')
          call pltend
       endif
       call pathshift(i,obi)
       if(omegad.ne.omega.and.obi.ne.0)then
!          write(*,'(a,i4,6f10.5)')'i,obi,omega,omegad',i,obi,omega,omegad
!          obi=0.
       endif
       omegab(i)=omegab(i)+sqm1*obi
       exptb=exp(sqm1*omegad*2*pi/omegab(i))
       if(.not.abs(exptb).lt.1.e10)exptb=1.e10  ! Enable divide by infinity
       if(i.eq.1)then
          dFdvprev=dFdvpsi
          exptbprev=exptb
       endif
       dob=omegab(i)-omegab(i-1)
       cdvpsi=dvpsi*(1.+sqm1*imag(dob)/real(dob))
       ! Strictly to get dFdvpsi we need to multiply by the omega f' terms
       Fnonres(i)=dFdvpsi*(omegad*dfe-(omegad-omega)*dfeperp)
       ! and correct for the imaginary shift of omegab:
       Fnonres(i)=Fnonres(i)+sqm1*real(Fnonres(i)-Fnonres(i-1))  &
            /real(omegab(i)-omegab(i-1))*obi
       if(tbe(i)*imag(omegad).lt.-4.)then ! Hack to fix giant dFdvpsi problem.
          write(*,*)'drop',tbe(i),imag(omegad),dFdvpsi
          Fnonres(i)=0.
       endif
       ! Then multiply by the resonance factor and the complex dvpsi and sum.
       Ftrap(i)=0.5*(Fnonres(i)/(1.-exptb) &
            +Fnonres(i-1)/(1.-exptbprev))*cdvpsi
       if(.not.(abs(Ftrap(i)).ge.0))then
          write(*,*)'Ftrap NAN?',i
          write(*,*)Fnonres(i),Ftrap(i)
          write(*,*)omegab(i-1),omegab(i)
          write(*,*)omegad,omegad/omegab(i)
          write(*,*)exptb,exptbprev
          stop
       endif
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
    Ftrap(ne)=Ftrap(ne-1)+0.5*(Ftrap(ne-1)-Ftrap(ne-2))
    Ftotal=Ftotal+2.*Ftrap(i)
    Fnonres(ne)=Fnonres(ne-1) !+0.5*(Fnonres(ne-1)-Fnonres(ne-2))
    ! Now Ftrap(i) is a quantity when simply summed over all ne positions
    ! and multiplied by 2 gives the total force. 
    if(.false.)then
    ! Diagnostic plots
    call autoplot(real(omegab(1:ne)),real(Fnonres(1:ne)),ne)
    call dashset(2)
    call polyline(real(omegab(1:ne)),imag(Fnonres(1:ne)),ne)
!    call dashset(0)
!    call polyline(real(omegab(1:ne)),real(1.e3*Ftrap(1:ne)),ne)
    call dashset(0)
    call polyline(real(omegab),imag(.001/(1.-exp(sqm1*omegad*2.*pi/omegab))),ne)
    call pltend()
    endif
  end subroutine FtEint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pathshift(i,obi)
! Calculate the required shift of the omegab path below the axis for this
! mesh point.     
    real :: obi
    integer :: el,ielsign
    real :: doel,doel2,dob
    real,parameter :: Sc=4,Rc=2./Sc

    ! Select the closest odd resonance such that el*omegab=omegar
    ielsign=int(sign(1.,real(omegad)))
    el=(2*int( (abs(real(omegad))/real(omegab(i))-1)/2. )+1)*ielsign
    doel=real(omegab(i))-real(omegad)/el
    doel2=real(omegab(i))-real(omegad)/(el+2*ielsign)
    if(abs(doel2).lt.abs(doel))then
       el=el+2*ielsign
       doel=doel2
    endif
    ! Calculate the required omegab imaginary part.
    dob=real(omegab(i)-omegab(i-1))
    obi=-ielsign*max(0.,dob/Rc-imag(omegad)/abs(el)) &
         *max(0.,1.-(doel/(Sc*dob))**2)
!    if(obi.ne.0.and.abs(el).lt.10)write(*,'(i4,f8.4,i5,f8.4,a,4f10.5)') &
!         i,sqrt(-Wt(i)),el,real(omegab(i))/(real(omegad)/el) &
!         ,' obi',obi,real(omegad-omega)/omegac,dob,doel
  end subroutine pathshift
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
!****************************************************************
! This is exp(X^2)*erfc(X)
      FUNCTION expERFCC(X)
      Z=ABS(X)      
      T=1./(1.+0.5*Z)
      expERFCC=T*EXP(-1.26551223+T*(1.00002368+T*(.37409196+             &
     &    T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+        &
     &    T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
      IF (X.LT.0.) expERFCC=2.*exp(z**2)-expERFCC
      END
!********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module shiftmode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
