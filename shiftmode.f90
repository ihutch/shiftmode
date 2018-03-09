!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module shiftmode
  ! We use x in place of z, because x is real.
  integer, parameter :: nx=100, ne=50, nvy=50
  real :: xL=19.,Emax=2.,vymax=4.            ! Hole length, Energy, v_y 
  real :: psi=.05,pL=4.,k=.01, Ty=1.          ! psi, sech4width, k, Ty
  real :: beta=-5.                           ! inverse hole parallel temp.
  real :: pi=3.1415926   
  complex :: omega=(0.0,.01)                  ! complex frequency
  integer :: idebug=0
  complex :: omegad,sqm1=(0.,1.),Ftraptotal
  ! Position arrays
  real :: dx
  real :: x(nx),phi(nx),phiprime(nx),tau(nx),v(nx)
  real :: denad(nx),priorda(nx),denadtemp(nx)
  real :: xt(nx),phit(nx),phitprime(nx)      ! Values in trapped region
  complex :: Lt(nx),phiut(nx),dent(nx)       ! tilde L,phiu,density
  complex :: denttemp(nx),priorcontrib(nx),C(nx)
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
  end subroutine initialize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!Trapped particle routines. Can only be called after initialization.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine cxcalc(vpsi,omegad,trapforce,tb,dfe,dfeperp,xlent)
    ! Calculate the past integral of tau vs x and C(tau) and
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

    complex :: Ltint,Ltemp,Ctbb2,exptb,exptbb2,sumfactor,exptau
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
    C=0.
    phiut=0.
    ! Integrate to get tau and C(tau)
    v(istart)=sqrt(2.*(phit(istart)+Wj))
    tau(istart)=v(istart)/phitprime(istart) ! To get to v at accel phiprime.
    Ltint=exp(-sqm1*omegad*tau(istart))
    C(istart)=-0.5*v(istart)*(Ltint-1.)/(omegad)
    do i=istart+1,iend
       v(i)=sqrt(2.*(phit(i)+Wj))
       if(v(i).lt.0.7*vpsi)then    ! Use phiprime to determine step
          dtau=(v(i)-v(i-1))*2./(phitprime(i)+phitprime(i-1))
       else                        ! Use v to determine step. 
          dtau=dxt*0.5*(v(i-1)+v(i))/(v(i-1)*v(i))
       endif
       tau(i)=tau(i-1)+dtau
       Ltemp=Ltint
       Ltint=exp(-sqm1*omegad*tau(i)) ! Current exponential
       vmean=(v(i)+v(i-1))/2.
       C(i)=C(i-1)-vmean*(Ltint-Ltemp)/(omegad)   ! New integral
    enddo
    tbb2=tau(iend)+tau(istart)
    tau(iend+1)=tbb2
    tb=2.*tbb2
    vmean=v(iend)/2.
    Ctbb2=C(iend)-vmean*(exp(-sqm1*omegad*tbb2)-Ltint)/(omegad)
    C(iend+1)=Ctbb2
! Calculate phiut from C(i), and on the way integrate force round the orbit
    exptb=exp(sqm1*omegad*tb)
    exptbb2=exp(sqm1*omegad*tbb2)
    sumfactor=1./(1.-exptb)
    trapforce=0.
    phiut=0.
!    write(*,*)'omegad,exptb,exptbb2',omegad,exptb,exptbb2
    do i=istart,iend
       exptau=exp(sqm1*omegad*tau(i))
       phiut(i)=omegad*exptau*((1.-exptb)*C(i)+(exptb-exptbb2)*Ctbb2) &
            *sumfactor
       trapforce=trapforce+sqm1*phiut(i)*phitprime(i) &
            *(omegad*dfe-(omegad-omega)*dfeperp)*dxt
    enddo
    if(idebug.lt.-2)then
       write(*,*)'istart,xlent,dxt',istart,xlent,dxt
       write(*,*)'xt'
       write(*,'(10f8.3)')xt
       write(*,*)'istart,iend,tbb2,sq2psi',istart,iend,tbb2,sqrt(2.*psi)
       write(*,*)'sumfactor',sumfactor
       write(*,*)'v'
       write(*,'(10f8.4)')v
       write(*,*)'xt'
       write(*,'(10f8.4)')xt
       write(*,*)'tau'
       write(*,'(10f8.2)')tau
       write(*,*)'real(C)'
       write(*,'(10f8.3)')real(C)
       write(*,*)'real(phiut)'
       write(*,'(10f8.2)')real(phiut)
    endif
  end subroutine cxcalc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine FtEint(omegad,Ftotal,dfperpdWperp)
  ! Here we need to iterate over fe (trapped). Wj=vpsi^2/2-psi. So
  ! vpsi=sqrt(2(psi+Wj)) We must use cells that fill the Wj range 0 to
  ! -psi.  Equal W steps are best.
    complex :: omegad,Ftotal
    character*20 :: string
    Ftotal=0.
    if(idebug.eq.-2)then
       write(*,*)'FtEint',omegad,beta,dfperpdWperp
       write(*,*)'  vpsi   fe*s2pi            Ftrap       ' &
         ,'       tb     dvpsi    dFtotal'
       call pltinit(-8.,8.,-8.,8.)
       call axis()
    endif
    iwpow=2         ! Equal W spacing ipow=1
    do i=1,ne
       Wj=psi*(-((i-0.5)/ne)**iwpow)
       fe=exp(-beta*Wj)/sqrt(2.*pi) ! Normalized f_\parallel
       dfe=-beta*fe                 ! df_||/dW_||
       dfeperp=dfperpdWperp*fe      ! df/dW_perp
       vpsi=sqrt(2.*(psi+Wj))
       dvpsi=-sqrt(2.*psi*(1.-(float(i)/ne)**iwpow)) &
            +sqrt(2.*psi*(1.-(float(i-1)/ne)**iwpow))
       ! calculate the force Ftrap for this dvpsi and dvy element:
       call cxcalc(vpsi,omegad,Ftrap(i),tbe(i),dfe,dfeperp,xlen(i))
       if(idebug.eq.-2)then
!          write(*,*)'tau'
          write(string,'(10f8.2)')tbe(i)
          if(i.eq.0)then
!             call axlabels('x','v')
!             call autoplot(xt,real(phiut),nx)
          else
             call cyccolor(i,14)
             call polyline(xt,10.*v+2.,nx)
!             call polyline(xt,real(C),nx)
             call polyline(xt,imag(phiut),nx)
             call polyline(xt,real(phiut)-4.,nx)
             if(maxval(abs(real(phiut))).gt.3.)then
!                write(*,'(10f8.4)'),real(phiut)
                idb=idebug
                idebug=-3
!                call cxcalc(vpsi,omegad,Ftrap(i),tbe(i),dfe,dfeperp,xlen(i))
                idebug=idb
             endif
          endif
!          call legendline(-0.4,1-i*.05,258,string)
       endif
       if(idebug.eq.-2)then
          write(*,'(2f8.4,a,2es12.4,a,f8.3,f9.5,es12.4)')vpsi&
               ,fe*sqrt(2.*pi) &
               ,' (',Ftrap(i),')',tbe(i),dvpsi,real(Ftrap(i))*dvpsi
       endif
       Ftotal=Ftotal+Ftrap(i)*vpsi*dvpsi  ! Add to Ftotal integral. 
    enddo
    if(idebug.eq.-2)then 
       call pltend()
!       call autoplot(real(Ftrap),imag(Ftrap),ne)
!       call pltend()
    endif
  end subroutine FtEint
  !--------------------------------------------
  subroutine FtVyint()
    ! Integrate (histogram) over vy to get the three ft terms.
    ! Total fte is omega*ft1int-k*ft2int+k*ft3int. Done in dentcalc.
    ! It might be more efficient to have an array of omegad and do
    ! all the tauxcalcs at once on the array.
    Ftraptotal=0.
    do i=1,nvy
       omegad=omega-k*vy(i)
       call FtEint(omegad,Ftrapvy(i),fywy(i))
!       write(*,*)'dvy,Ftraptotal',dvy,Ftraptotal
       Ftraptotal=Ftraptotal+Ftrapvy(i)*dvy
       if(idebug.ne.0)write(*,'(a,i4,es10.2,a,f8.4,2e12.4)')&
            'i,R(omegad)',i,real(omegad),  &
            '  vy,Im(Ftrapvy)',vy(i),imag(Ftrapvy(i)),imag(Ftraptotal)
    enddo
  end subroutine FtVyint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Passing particle routines.
  subroutine tauxcalc(vinf,omegad,thisdtaumax)  ! Passing particles
    ! integrate dx/v to get v(x), tau(x), Lt(x,t=0), phiut, at vinf,omegad
    ! return maximum dtau used for accuracy checking
    real    :: vinf,thisdtaumax
    complex :: omegad      ! Input omega'
    complex :: Ltint,Ltemp
    logical :: ldone=.false.


    thisdtaumax=0.
    scalefactor=1.e4       ! To avoid overflows we need intermediate scaling
    taustep=alog(scalefactor)/max(imag(omegad),1.e-5) ! tau-value at which
    tau(1)=0.                                         ! to rescale.
    v(1)=sqrt(vinf**2+2.*phi(1))
    Ltint=1.                                          ! prior exp(0)
    Lt(1)=0.                                          ! Lt integrated
    phiut(1)=v(1)-vinf
! Do \int_-\infty^t exp(-i*omegad*tau) dtau, passing: Lt(i). And phiut(i)
    do i=2,nx
       v(i)=sqrt(vinf**2+2.*phi(i))
       dtau=dx*0.5*(v(i-1)+v(i))/(v(i-1)*v(i))        ! tau increment varies
       if(dtau.gt.thisdtaumax)thisdtaumax=dtau
       tau(i)=tau(i-1)+dtau
       Lt(i)=Lt(i-1)
       if(tau(i).gt.taustep)then        ! Rescale to avoid overflow.
          do
             tau(i)=tau(i)-taustep         ! Subtract from tau
             Ltint=Ltint*exp(sqm1*omegad*taustep)  ! Scale down prior values
             Lt(i)=Lt(i-1)*exp(sqm1*omegad*taustep)! by scalefactor.
             if(tau(i).le.taustep)exit
          enddo
       endif
       Ltemp=Ltint
       Ltint=exp(-sqm1*omegad*tau(i)) ! Current exponential
       vmean=(v(i)+v(i-1))/2.-vinf
       Lt(i)=Lt(i)-vmean*(Ltint-Ltemp)/(omegad)   ! New integral
       ! Adjust scales so that tau at end is t=0 in
       ! phiut and the effective integral is from -infty to 0.
       !       phiut(i)=omegad*Lt(i)*exp(sqm1*omegad*tau(i))
       phiut(i)=omegad*Lt(i)*exp(sqm1*omegad*tau(i))
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
    
  end subroutine tauxcalc
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
       call tauxcalc(vinf,omegad,dtaumax) ! Integrate along past orbit
                                          ! to get phiut, hence ft1-3
       ft1=dfe*fy(i)*phiut
       ft2=dfe*fy(i)*vy(i)*phiut
       ft3=fe*fy(i)*vy(i)*phiut
       ft1int=ft1int+ft1*dvy
       ft2int=ft2int+ft2*dvy
       ft3int=ft3int+ft3*dvy
    enddo
  end subroutine ftildecalc
  !--------------------------------------------
  subroutine dentcalc2 ! Integrate over v_parallel to get passing n-tilde.
    ! Alternate using histogram vinf distribution.
    ! This is the outermost integral. It calls inner integrals (dvy(dtau)).
    ! On the way, calculate the f(x,E) array. Only non-adiabatic.
    ! Also calculate the adiabatic density perturbation, denad.
    ! So far only positive velocity direction. 
    vinfmax=sqrt(2.*Emax)      ! The uppermost orbit
    dvinf=vinfmax/ne           ! Use steps of constant vinf spacing.
    denad=0.                   ! adiabatic density 
    dent=0.                    ! n-tilde
    vinf=0.
    do i=1,ne                  ! Integrate over ne orbits (energies)
       vinf=(i-0.5)*dvinf
       E(i)=vinf**2/2.
       fe0(i)=exp(-E(i)) ! Unit height Maxwellian with unit temperature.
       fe0de(i)=fe0(i)   ! Derivative wrt E is same.
       ! ftildecalc calls tauxcalc which also calculates v(nx) for this vinf.
       call ftildecalc(vinf,fe0(i),fe0de(i))
       fte(:,i)=(omega*ft1int+k*(-ft2int+ft3int))*sqm1  ! tilde-f_x(x,E)
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
       dent=dent+denttemp*dvinf  ! Non-adiabatic
! (-dphi/dx=\hat\phi)*(df/dE dv) for \Delta=1. ! Adiabatic not needed.
       denadtemp=-phiprime*fe0de(i)*vinf/v
       denad=denad+denadtemp*dvinf 
! end of adiabatic section.
    enddo
  end subroutine dentcalc2
!  include 'obsoletecode.f90'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module shiftmode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
  use shiftmode 
  integer, parameter ::   nk=21
  real :: kik(nk)
  complex :: Fcpassing(nk),Ftrapped(nk),Ftotal
  call initialize
  write(*,*)'nx, ne, nvy,   xL,    pL,   omegar,  omegai,     k     psi   beta'
  write(*,'(3i4,7f8.4)')nx,ne,nvy,xL,pL,real(omega),imag(omega),k,psi,beta
!  call passingdiags

  idebug=-2
!  call FtVyint()
!  write(*,*)'Ftraptotal',Ftraptotal
  if(.false.)then
     omegad=complex(.05,.003)
     call FtEint(omegad,Ftotal,1.)
     write(*,*)'k*vy,Ftotal (trapped)',real(omegad),Ftotal
  endif

  idebug=0
! k-scan
  omega=(0.0,.015)
  akmax=.04
  write(*,*)'   k     Fpassing                Ftrapped'
  do ik=1,nk
     k=(ik-1.)*akmax/(nk-1.)
     kik(ik)=k
     call dentcalc2
!     write(*,*)'dtaumax=',dtaumax,' kvymax.dtaumax=',k*vymax*dtaumax
! Integrate phi'*n-tilde dx     
     Fcpassing(ik)=0.
     do i=1,nx
        Fcpassing(ik)=Fcpassing(ik)+phiprime(i)*dent(i)*dx
     enddo
     call FtVyint()
     Ftrapped(ik)=Ftraptotal
     write(*,'(f6.4,4es12.3)')k,Fcpassing(ik),Ftraptotal
  enddo
  call multiframe(2,1,3)
  call autoplot(kik,real(Fcpassing),nk)
  call axis2()
  call axlabels('k','Fpassing')
  call polyline(kik,imag(FCpassing),nk)
  call autoplot(kik,real(Ftrapped),nk)
  call axis2()
  call axlabels('k','Ftrapped')
  call polyline(kik,imag(Ftrapped),nk)
  call dashset(2)
  call polyline(kik,real(-Fcpassing),nk)
  call pltend()
  
end program main
