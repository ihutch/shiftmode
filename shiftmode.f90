!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module shiftmode
  ! We use x in place of z, because x is real.
!  integer, parameter :: nx=50, ne=60, nvy=200
  integer, parameter :: nx=100, ne=30, nvy=50
  real :: xL=20.,Emax=4.,vymax=4.            ! Hole length, Energy, v_y 
  real :: psi=.1,pL=4.,k=.01, Ty=1.          ! psi, sech4width, k, Ty
  real :: beta                               ! inverse hole parallel temp.
  real :: pi=3.1415926   
  complex :: omega=(0.0,.01)                  ! complex frequency
  integer :: idebug=0
  complex :: omegad,sqm1=(0.,1.),Ftraptotal,Fpasstotal
  ! Position arrays
  real :: dx
  real :: x(nx),phi(nx),phiprime(nx),tau(nx),v(nx)
  real :: denad(nx),priorda(nx)
  real :: xt(nx),phit(nx),phitprime(nx)      ! Values in trapped region
  complex :: Lt(nx),phiut(nx),dent(nx)       ! tilde L,phiu,density
  complex :: priorcontrib(nx),phiutest(nx)
  complex :: ft1(nx),ft2(nx),ft3(nx),ft1int(nx),ft2int(nx),ft3int(nx)
  ! Energy arrays
  real :: de,E(ne),fe0(ne),fe0de(ne)
  real :: xlen(ne),tbe(ne)                  ! Trapped orbit length, period.
  complex :: fte(nx,0:ne)
  complex :: Ftrap(ne)
  ! vy arrays
  real :: vy(nvy),fy(nvy),fywy(nvy),dtaumax
  real :: dvy
  complex :: Ftrapvy(nvy)
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initialize
    ! x-arrays
    dx=xL*2./(nx-1.)                            ! Step size in x
    x=dx*(/(i-1,i=1,nx)/)-xL                    ! x-position array
    phi=psi/cosh(x/pL)**4                       ! Potential
    phiprime=-psi*sinh(x/pL)/cosh(x/pL)**5*4/pL ! phi x-gradient
    ! vy-arrays
    dvy=vymax*2./(nvy-1.)                       ! vy-step
    vy=dvy*(/(i-1,i=1,nvy)/)-vymax
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
  subroutine cxcalc(vpsi,omegad,trapforce,tb,dfe,dfeperp,xlent)
    ! Calculate the past integral of tau vs x and L(tau) and
    ! phut(taut) vs x for a trapped particle of velocity vpsi at peak
    ! of potential. Return trapforce and orbit period tb These
    ! quantities are the differential force density with respect to vy
    ! and vpsi. And so have to be integrated f(vpsi,vy) dvpsi dvy. The
    ! contributions to force are proportional to df/dW_parallel and
    ! df/dW_perp evaluated at these orbit energies, which we denote dfe
    ! and dfeperp. xlent returns the orbit length extreme, used for
    ! diagnostic purposes.
    real :: vpsi,dfe,dfeperp,tb,xlent
    complex :: omegad,trapforce

    complex :: exptau,Ltemp,Ltint,Ltbb2,exptb,exptbb2,sumfactor,Ltfactor
    Wj=vpsi**2/2-psi
    if(Wj.gt.0)then
       write(*,*)'Positive energy in trapped particle code'
       return
    endif
    vj=sqrt(-2.*Wj)
    ! Determine the starting and ending nodes istart, iend for this
    ! orbit given the previously calculated symmetric phi(x)-array.
    do istart=1,nx/2
       if(phi(istart).gt.-Wj)exit
    enddo
    ! Then remesh the orbit over just its actual x-extent.
    xlent=abs(x(istart-2))
    do j=1,3  ! Iterate the length determination if needed.
       dxt=xlent*2./(nx-1.)                           ! Step size in x
       xt=dxt*(/(i-1,i=1,nx)/)-xlent                  ! x-position array
       phit=psi/cosh(xt/pL)**4                        ! Potential
       phitprime=-psi*sinh(xt/pL)/cosh(xt/pL)**5*4/pL ! phi x-gradient
       do istart=1,nx/2
          if(phit(istart).gt.-Wj)exit
       enddo
       if(istart.lt.nx/4)then
          exit
       else
          xlent=abs(xt(istart-2))
       endif
    enddo
    iend=nx+1-istart
    if(istart.gt.nx/2)stop 'No steps on remeshed trapped orbit'    
    tau=0.
    v=0.
    Lt=0.
    C=0.
    phiut=0.
    dtauprior=0.
    ! Integrate to get tau and L(tau)
    v(istart)=sqrt(2.*(phit(istart)+Wj))
    tau(istart)=v(istart)/phitprime(istart) ! To get to v at accel phiprime.
    Ltfactor=exp(sqm1*omegad*tau(istart))
    Lt(istart)=-0.5*v(istart)*(1.-Ltfactor)/(omegad)
    Ltint=exp(-sqm1*omegad*tau(istart))
    do i=istart+1,iend
       v(i)=sqrt(2.*(phit(i)+Wj))
       if(v(i).lt.0.7*vpsi)then    ! Use phiprime to determine step
          dtau=(v(i)-v(i-1))*2./(phitprime(i)+phitprime(i-1))
       else                        ! Use v to determine step. 
          dtau=dxt*0.5*(v(i-1)+v(i))/(v(i-1)*v(i))
       endif
       tau(i)=tau(i-1)+dtau
       Ltfactor=exp(sqm1*omegad*(tau(i)-tau(i-1))) ! Current exponential
       Ltemp=Ltint
       Ltint=exp(-sqm1*omegad*tau(i)) ! Current exponential
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
! Calculate phiut from L(i), and on the way integrate force round the orbit
    exptb=exp(sqm1*omegad*tb)
    exptbb2=exp(sqm1*omegad*tbb2)
    sumfactor=1./(1.-exptb)
    trapforce=0.
    phiut=0.
! Form phiut as omegad*Lt, and the contribution to trapforce as
! \int sqm1*phitprime*(omega*dfe ...)*phiut  dtau vpsi
    do i=istart,iend
       exptau=exp(sqm1*omegad*tau(i))
       phiut(i)=omegad*((1.-exptb)*Lt(i)+exptau*(exptbb2-1.)*Ltbb2) &
            *sumfactor
! The dtau to be applied to this subsequent tau integral. 
       dtau=tau(i)-tau(i-1)
       trapforce=trapforce+sqm1*  &
       (phiut(i)*phitprime(i)+phiut(i-1)*phitprime(i-1))/2. &
            *(omegad*dfe-(omegad-omega)*dfeperp)*dtau*vpsi
    enddo
    if(idebug.lt.-2)then
       write(*,'(a,2i3,5f8.3)')'istart,iend,dxt,tbb2,sumfac' &
            ,istart,iend,dxt,tbb2,sumfactor
       write(*,*)'v'
       write(*,'(10f8.4)')v
       write(*,*)'xt'
       write(*,'(10f8.4)')xt
       write(*,*)'tau'
       write(*,'(10f8.2)')tau
       write(*,*)'imag(Lt)'
       write(*,'(10f8.5)')imag(Lt)
       write(*,*)'real(phiut)'
       write(*,'(10f8.5)')real(phiut)
    endif
  end subroutine cxcalc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine FtEint(omegad,Ftotal,dfperpdWperp,fperp)
  ! Here we need to iterate over fe (trapped). Wj=vpsi^2/2-psi. So
  ! vpsi=sqrt(2(psi+Wj)) We must use cells that fill the Wj range 0 to
  ! -psi. 
    complex :: omegad,Ftotal
    Ftotal=0.
    if(idebug.eq.-2)then
       write(*,*)'FtEint',omegad,beta,dfperpdWperp
       write(*,*)'  vpsi   fe*s2pi            Ftrap       ' &
         ,'       tb     dvpsi    dFtotal'
       call pltinit(-12.,12.,-8.,8.)
       call axis()
    endif
    iwpow=2         ! Equal W spacing is ipow=1
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
       call cxcalc(vpsi,omegad,Ftrap(i),tbe(i),dfe,dfeperp,xlen(i))
       if(.not.(abs(Ftrap(i)).ge.0))then
          write(*,*)'Ftrap NAN?',Ftrap(i),Wj,dfe,dfeperp,vpsi,dvpsi
          stop
       endif
       if(idebug.eq.-2)then
          write(*,'(2f8.4,a,2es12.4,a,f8.3,f9.5,es12.4)')vpsi&
               ,fe*sqrt(2.*pi) &
               ,' (',Ftrap(i),')',tbe(i),dvpsi,real(Ftrap(i))*dvpsi
       endif
       Ftotal=Ftotal+2.*Ftrap(i)*dvpsi  ! Add to Ftotal integral. 
       ! Here we do not multiply by vpsi because that was done in cxcalc.
       ! But we multiply by 2. to account for \pm v_\psi.
    enddo
  end subroutine FtEint
  !--------------------------------------------
  subroutine FtVyint()
    ! Integrate (histogram) over vy to get the three ft terms.
    Ftraptotal=0.
    do i=1,nvy
       omegad=omega-k*vy(i)
       call FtEint(omegad,Ftrapvy(i),fywy(i),fy(i))
       Ftraptotal=Ftraptotal+Ftrapvy(i)*dvy
       if(idebug.ne.0)write(*,'(a,i4,es10.2,a,f8.4,2e12.4)')&
            'i,R(omegad)',i,real(omegad),  &
            '  vy,Im(Ftrapvy)',vy(i),imag(Ftrapvy(i)),imag(Ftraptotal)
    enddo
  end subroutine FtVyint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Passing particle routines.
  !--------------------------------------------
  subroutine tauxcalcnew(vinf,omegad,thisdtaumax)  
    ! Passing particles
    ! integrate dx/v to get v(x), tau(x), Lt(x,t=0), phiut, at vinf,omegad
    ! return maximum dtau used for accuracy checking
    real    :: vinf,thisdtaumax
    complex :: omegad      ! Input omega'
    complex :: expdtau
    logical :: ldone=.false.
    tau(1)=0.
    thisdtaumax=0.
    v(1)=sqrt(vinf**2+2.*phi(1))
    Lt(1)=0.                                          ! Lt integrated
    phiut(1)=v(1)-vinf
    phiut(1)=0.                                       ! More consistent.
! Do \int_-\infty^t exp(-i*omegad*tau) dtau, passing: Lt(i). And phiut(i)
    do i=2,nx
       v(i)=sqrt(vinf**2+2.*phi(i))
       dtau=dx*0.5*(v(i-1)+v(i))/(v(i-1)*v(i))        ! tau increment varies
       if(dtau.gt.thisdtaumax)thisdtaumax=dtau
       tau(i)=tau(i-1)+dtau
       expdtau=exp(sqm1*omegad*dtau)
       vmean=(v(i)+v(i-1))/2.-vinf
       Lt(i)=expdtau*Lt(i-1)-vmean*(1.-expdtau)/(omegad)   ! New integral
       !       phiut(i)=omegad*Lt(i)*exp(sqm1*omegad*tau(i))
       phiut(i)=omegad*Lt(i)
    enddo
    if(idebug.gt.0.and..not.ldone)then
       write(*,'(a,2f8.4,a,f8.4,a,f8.2)')'omegad',omegad,' vinf', &
         vinf,' taustep',taustep
       write(*,*)'phi'
       write(*,'(10f8.4)')phi
       write(*,*)'v'
       write(*,'(10f8.4)')v
       write(*,*)'tau'
       write(*,'(10f8.2)')tau
       write(*,*)'real(Lt)'
       write(*,'(10f8.1)')real(Lt)
       write(*,*)'real(phiut)'
       write(*,'(10f8.2)')real(phiut)
       ldone=.true.
    endif
  end subroutine tauxcalcnew
  !--------------------------------------------
  subroutine ftildecalc(vinf,fe,dfe)   ! Integrate over vy to get ft-parallel
                   ! fe is parallel distrib, dfe is its derivative wrt E.
    ft1int=0. ; ft2int=0. ; ft3int=0.
    ! Integrate (histogram) over vy to get the three ft terms.
    ! Total fte is omega*ft1int-k*ft2int+k*ft3int. Done in dentcalc.
    ! It might be more efficient to have an array of omegad and do
    ! all the tauxcalcs at once on the array.
    do i=1,nvy
       omegad=omega-k*vy(i)
       call tauxcalcnew(vinf,omegad,dtaumax) ! Integrate along past orbit
                                          ! to get phiut, hence ft1-3
       ft1=dfe*fy(i)*phiut
       ft2=dfe*fy(i)*vy(i)*phiut
       ft3=fe*fywy(i)*vy(i)*phiut
       ft1int=ft1int+ft1*dvy
       ft2int=ft2int+ft2*dvy
       ft3int=ft3int+ft3*dvy
    enddo
  end subroutine ftildecalc
  !--------------------------------------------
  subroutine dentcalc2 
    ! Integrate over v_parallel to get passing n-tilde and force.
    ! Alternate using histogram vinf distribution.
    ! This is the outermost integral. It calls inner integrals (dvy(dtau)).
    ! On the way, calculate the f(x,E) array. Only non-adiabatic.
    ! Also calculate the adiabatic density perturbation, denad.
    ! So far only positive velocity direction. 
    complex :: Force
    sq2pi=sqrt(2.*pi)
    vinfmax=sqrt(2.*Emax)      ! The uppermost orbit
    dvinf=vinfmax/ne           ! Use steps of constant vinf spacing.
    denad=0.                   ! adiabatic density 
    dent=0.                    ! n-tilde
    vinf=0.
    do i=1,ne                  ! Integrate over ne orbits (energies)
       vinf=(i-0.5)*dvinf
       E(i)=vinf**2/2.
       fe0(i)=exp(-E(i))/sq2pi ! Normalized Maxwellian with unit temperature.
       fe0de(i)=-fe0(i)   ! Derivative wrt E is minus the same.
       ! ftildecalc calls tauxcalc which also calculates v(nx) for this vinf.
       call ftildecalc(vinf,fe0(i),fe0de(i))
       ! No factor q_e has been applied to phiut, so no q_e factor here.
       fte(:,i)=sqm1*(omega*ft1int+k*(-ft2int+ft3int))  ! tilde-f_x(x,E)
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
       dent=dent+(fte(:,i)*vinf/v)*dvinf  ! Non-adiabatic
    enddo
! For stationary holes negative velocities just double the effective ne
! (Although actually it is the force they double and non-zero hole speed
! would require this to be done differently.)
    dent=dent*2
! Integrate to get force:
    Force=0.
    do i=1,nx
       Force=Force+phiprime(i)*dent(i)*dx
    enddo
    Fpasstotal=Force
  end subroutine dentcalc2
!  include 'obsoletecode.f90'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module shiftmode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
