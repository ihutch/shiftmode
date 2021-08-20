! This is a general calculation of the force arising from an
! linearized oscillatory displacement (shift) perturbation of a
! localized structure that gives potential ENERGY phi. The potential
! is either a hill or a valley and tends to zero at large distances,
! but its derivative passes through zero only once, at z=0, which is
! the only potential extremum (discouting infinity). When the extremum
! is a minimum, there are locally trapped particles, but when it is a
! maximum, reflected orbits are not locally trapped. There are
! therefore three types of orbit: Passing, Trapped, or Reflected,
! which must be treated differently.

! The force is the integral dz of -d\phi/dz times the non-adiabatic
! perturbed f, which is an integral over the past time d\tau of the
! past orbit. Both integrals can be expressed as time integrals.
! However, equal intervals of neither time nor space are universally
! optimal choices. Equal intervals of space are suboptimal near a
! reflection. Equal intervals of time are suboptimal for orbits of
! passing energy close to zero (because they move slowly at large z
! which is a less important region). 

! The proposed z array is 
!    zi= z1+ z2/(1+2K)[2K(i/n)+(i/n)^2]
! where repelled z1=zR, z2=(zm-zR), 
! attracted is z1=0, z2=zR trapped, z1=0, z2=zm passing.
! And K=sqrt(max(0,W/psi-1)) repelled, K=-1-sqrt(max(0,-W/psi)) attacted.  
! Repelled is \psi>0, attracted \psi<0.

module shiftgen
  integer, parameter :: ngz=50,nge=100
  real, dimension(-ngz:ngz) :: zg,vg,ones=1.,phig,phigprime,taug
  complex, dimension(-ngz:ngz) :: Lg
  complex :: omegag=(1.,0.),sqm1=(0.,1.),Ftot
  real :: psig=.1,Wg,zm=10.,v0
  integer :: ivs
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine makezg(isigma)
    integer :: isigma
! Calculate the zg mesh, and on it phig, vg, phigprime for an incoming
! orbit from z=isigma*zm (v sign -isigma) or trapped orbit from its
! isigma end.
    z0=0.
! Find any reflection position zR. ! Put equal to 0 if W>max(phi).
    call orbitend(Wg,z0,isigma*zm)
!    write(*,*)Wg,psig,'Reflection?',z0
    zR=z0
    ivs=-1
    if(psig.gt.0)then ! Repelling potential
       gK=sqrt(max(0.,Wg/psig-1.))
       z1=zR; z2=isigma*zm-zR     ! Or maybe z2=zm
       if(Wg.lt.psig)ivs=1  ! Reflected orbit, all z are negative.
    elseif(psig.lt.0)then !Attracted
       gK=-1.-10.*sqrt(max(0.,-Wg/psig))
       if(zR.eq.0.)then
          z2=isigma*zm
       else
          z2=zR
       endif
       z1=0.
    else  
       write(*,*)'WARNING psi is zero'
       stop
    endif
    write(*,*)'makezg: z1,z2,gK',z1,z2,gK
    do i=0,ngz
       zi=float(i)/ngz
       zg(i)=z1+z2*((2.*gK+zi)*zi)/(1.+2.*gK)
       phig(i)=phigofz(zg(i))
       vg(i)=-ivs*sqrt(2.*max(0.,Wg-phig(i)))
       phigprime(i)=phigprimeofz(zg(i))
       if(i.gt.0)then
          zg(-i)=ivs*zg(i)
          phig(-i)=phig(i)
          vg(-i)=-ivs*vg(i)
          phigprime(-i)=ivs*phigprime(i)
       endif
    enddo
  end subroutine makezg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine orbitend(Wj,z0,zL)
    ! Find by bisection the turning point if any of the orbit whose
    ! energy is Wj lying between zL (normally negative) and z0 (=0).  
    ! That is where phi=-Wj, return z0 s.t. |z0| is just lower.
    real :: Wj,zm,z1,z0,enrgym,enrgy0,enrgy1
    nbi=20
    z0=0.
    z1=zL
    enrgy0=phigofz(z0)-Wj
    enrgy1=phigofz(z1)-Wj
    if(sign(1.,enrgy1).eq.sign(1.,enrgy0))then
!       write(*,*)'orbitend energies do not bracket zero',enrgy0,enrgy1
       return
    endif
    do i=1,nbi
!       write(*,*)i,z0,enrgy0,z1,enrgy1
       zm=(z1+z0)/2.
       enrgym=phigofz(zm)-Wj
! Which value to replace.
       if(sign(1.,enrgym).eq.sign(1.,enrgy0))then
          if(z0.eq.zm)exit ! should never happen
          enrgy0=enrgym
          z0=zm
       else
          if(z1.eq.zm)exit
          enrgy1=enrgym
          z1=zm
       endif
    enddo
  end subroutine orbitend
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine testmakezg
    use shiftmode
    complex :: dForceg
    psig=.5
    isigma=-1
    call pltinit(-zm,zm,min(0.,psig),max(1.8*psig,.8*abs(psig)))
    call axis
    call axlabels('z','W')
    nw=10
    do i=1,nw
       if(psig.gt.0)Wg=1.01*psig*i/nw
       if(psig.lt.0)Wg=psig+1.42*abs(psig)*i/nw
       write(*,*)'Wg=',Wg
!       call makezg(isigma)
       call LofW(Wg,isigma,dForceg)
       call polymark(zg(-ngz:ngz),Wg*ones(-ngz:ngz),2*ngz+1,i)
       write(*,*)'dForceg=',dForceg
       write(*,'(10f8.4)')(zg(j),j=0,ngz)
       write(*,'(10f8.4)')(vg(j),j=0,ngz)
    enddo
    call pltend
  end subroutine testmakezg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function phigofz(zval)
    phigofz=psig/cosh(zval/4.)**4
  end function phigofz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function phigprimeofz(zval)
    phigprimeofz=-psig*sinh(zval/4.)/cosh(zval/4.)**5
  end function phigprimeofz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine LofW(Wg,isigma,dForceg)
  ! Calculate the past integral of tau vs x and Lg(tau) and
  ! vs z for orbit of energy W. CapPhi is i*omegag*Lg, where
  ! Lg= \int_{-\infty}^{t} (v-v_0)*exp(-i*omegag(tau-t)dtau.
  ! For untrapped particles (W>=0), v_0=vinf=sqrt(2*W).
  ! For trapped particles (W<0), v_0=sqrt(2*(W-psig)).  
  ! The integration is done on a z-grid such that dz=v*dtau.
  ! isigma is the sign of zg at the start of orbit (v sign -isigma)
  ! The total needs to account for v-sign=+isigma as well.

  ! Integrate Lg dt (=dz/v) to get the differential
  ! force density with respect to vy and vinf when multiplied by
  ! df/dW_parallel.
    complex :: dForceg,Lgfactor,Lgb2
    
    if(Wg.ge.0.)then
       v0=sqrt(2.*Wg)
    else
       if(Wg-psig.lt.0)stop 'Wg below minimum potential ERROR'
       v0=0.
    endif
    call makezg(isigma)
    taug(-ngz)=0.
    Lg(-ngz)=0.
    dForceg=0.
    do i=-ngz+1,ngz
       phigp=0.5*(phigprime(i)+phigprime(i-1))
       vmean=0.5*(vg(i)+vg(i-1))
       if(abs(phigp).gt.abs(vg(i)))then   ! Use dtau=dv*dtau/dv
          dtau=(vg(i)-vg(i-1))/phigp
       else                               ! Use dtau=dx*dtau/dx
          dtau=(zg(i)-zg(i-1))*vmean/(vg(i-1)*vg(i))
       endif
       taug(i)=taug(i-1)+dtau
       Lgfactor=exp(sqm1*omegag*dtau) ! Current exponential
       Lg(i)=Lgfactor*Lg(i-1)-(vmean-v0)*(1.-Lgfactor)/omegag
       dForceg=dForceg+sqm1*omegag*Lg(i)*phigp*dtau    ! sign?
       if(.not.real(dForceg).lt.1.e6)write(*,*)'real(dForceg)',real(dForceg)
    enddo
    if(Wg.lt.0)then ! Trapped particle correction and reintegration.
       Lgb2=Lg(ngz)
       tbg2=taug(ngz)
       Lg=Lg+exp(sqm1*omegag*taug)*(exp(sqm1*omegag*tbg2)-1.) &
            /(1.-exp(2.*sqm1*omegag*tbg2))*Lgb2 ! Prior Trapped bounces.
       dForceg=0.
       do i=-ngz+1,ngz
          phigp=0.5*(phigprime(i)+phigprime(i-1))
          dtau=taug(i)-taug(i-1)
          dForceg=dForceg+sqm1*omegag*Lg(i)*phigp*dtau    ! sign?
       enddo
    endif
  end subroutine LofW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module shiftgen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use shiftgen
call testmakezg
end
