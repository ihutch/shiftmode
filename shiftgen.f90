! This is a general calculation of the force arising from an
! linearized oscillatory displacement (shift) perturbation of a
! localized structure that gives potential ENERGY phi. The potential
! is either a hill or a valley and tends to zero at large distances,
! but its derivative passes through zero only once, at z=0, which is
! the only potential extremum (discounting infinity). When the extremum
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

! The z array is 
!    z(i)= z1+ z2/(1+2K)[2K(i/n)+(i/n)^2]
! when repelled, z1=zR, z2=(zm-zR), where zR is reflection point. 
! When attracted z1=0, z2=zR trapped, z2=zm passing.
! And K=sqrt(max(0,W/psi-1)) repelled, K=-1-sqrt(max(0,-W/psi)) attacted.  
! Repelling hill psi>0, attracted valley psi<0.

module shiftgen
  integer, parameter :: ngz=100,nge=100
  real, dimension(-ngz:ngz) :: zg,vg,ones=1.,phig,phigprime,taug
  complex, dimension(-ngz:ngz) :: Lg,CapPhig
  complex :: omegag=(1.,0.),sqm1=(0.,1.),Ftot
  real :: psig=.1,Wg,zm=10.,v0,z0,z1,z2,zR
  integer :: ivs,iws
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
          z2=zR*(1.+1.e-4/ngz) ! Force trapped end vg to zero.
       endif
       z1=0.
    else  
       write(*,*)'ERROR psi is zero'
       stop
    endif
!    write(*,*)'makezg: z1,z2,gK',z1,z2,gK
    do i=0,ngz
       zi=float(i)/ngz
       zg(i)=ivs*(z1+z2*((2.*gK+zi)*zi)/(1.+2.*gK))
       phig(i)=phigofz(zg(i))
       vg(i)=isigma*ivs*sqrt(2.*max(0.,Wg-phig(i)))
       phigprime(i)=phigprimeofz(zg(i))
       if(i.gt.0)then
          zg(-i)=ivs*zg(i)
          phig(-i)=phig(i)
          vg(-i)=-ivs*vg(i)
          phigprime(-i)=ivs*phigprime(i)
       endif
    enddo
!    write(*,*)zg(ngz-2),zg(ngz-1),zg(ngz)
!    write(*,*)vg(ngz-2),vg(ngz-1),vg(ngz)
  end subroutine makezg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine orbitend(Wj,z0,zL)
    ! Find by bisection the turning point if any of the orbit whose
    ! energy is Wj lying between zL (normally negative) and z0 (=0).  
    ! That is where phi=-Wj, return z0 s.t. |z0| is just lower.
    real :: Wj,zm,z1,z0,enrgym,enrgy0,enrgy1
    nbi=25
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
  real function phigofz(zval)
    phigofz=psig/cosh(zval/4.)**4
  end function phigofz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function phigprimeofz(zval)
    phigprimeofz=-psig*sinh(zval/4.)/cosh(zval/4.)**5
  end function phigprimeofz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine LofW(Wg,isigma,dForceg)
  ! Calculate the past integral of tau and Lg(tau) and
  ! vs z for orbit of energy W. Where
  ! Lg= \int_{-\infty}^{t} (v-v_0)*exp(-i*omegag(tau-t)dtau.
  ! For untrapped particles (W>=0), v_0=vinf=sqrt(2*W).
  ! For trapped particles (W<0), v_0=sqrt(2*(W-psig)).  
  ! The integration is done on a z-grid such that dz=v*dtau.
  ! isigma is the sign of zg at the start of orbit (v-sign -isigma)
  ! The total needs to account for v-sign=+isigma as well.

  ! Integrate Lg dt (=dz/v) to get the differential
  ! force density with respect to vy and vinf when multiplied by
  ! df/dW_parallel.
    complex :: dForceg,Lgfactor,Lgdtau,expdtau,CapPhigdtau
    complex :: exptbb2,exptau,forcedelta

    if(Wg.ge.0.)then
       v0=-isigma*sqrt(2.*Wg)
    else
       if(Wg-psig.lt.0)stop 'Wg below minimum potential ERROR'
       v0=0.
    endif
    call makezg(isigma)
!    write(*,*)'v0,vg(-ngz)',v0,vg(-ngz),vg(-ngz+1)
    taug(-ngz)=0.
    Lg(-ngz)=0.
    dForceg=0.
    ips=int(sign(1.,psig))
    iws=0
    do i=-ngz+1,ngz
       phigp=0.5*(phigprime(i)+phigprime(i-1))
       vmean=0.5*(vg(i)+vg(i-1))
       if(zR.ne.0.and.abs(vg(i)).lt.0.5*sqrt(2.*max(Wg,Wg-psig)))then
          dtau=-((vg(i)-vg(i-1))/phigp) ! Use dtau=dv*dtau/dv
!          write(*,*)'ips',i,ips,zR,psig,Wg
          if(ips.gt.0)iws=i               ! Reflected track
       else                               ! Use dtau=dx*dtau/dx
          dtau=((zg(i)-zg(i-1))*vmean/(vg(i-1)*vg(i)))
          if(ips.le.0..or.zR.eq.0)iws=i   ! Attracted or unreflected
       endif
       if(.not.dtau.lt.1e6)write(*,*)i,'dtau=',dtau,v0,vg(i-1)
       
       if(dtau.gt.10)write(*,*)'JUMP?',phigp,vg(i),vmean
       taug(i)=taug(i-1)+dtau
       Lgfactor=exp(sqm1*omegag*dtau) ! Current exponential
       Lg(i)=Lgfactor*Lg(i-1)-(vmean-v0)*(1.-Lgfactor)/omegag
       forcedelta=0.5*(phigprime(i)+phigprime(i-1))*( &
            (Lg(i-1)-(vmean-v0)/(-sqm1*omegag))*(Lgfactor-1)/(sqm1*omegag) &
            +(vmean-v0)/(-sqm1*omegag)*dtau)*abs(v0)
       dForceg=dForceg-sqm1*omegag*forcedelta
! Simple version without step integral correction.       
!       dForceg=dForceg-sqm1*omegag* 0.5*(Lg(i)*phigprime(i)+Lg(i-1)&
!            &*phigprime(i-1))*abs(v0)*dtau
       if(.not.real(dForceg).lt.1.e6)write(*,*)'real(dForceg)',real(dForceg)
    enddo
!    write(*,*)'v0,vmean,phigprime,dtau',v0,vmean,phigprime(ngz/3),dtau,forcedelta
    if(Wg.lt.0)then ! Trapped particle correction and reintegration.
! In shiftmode the division by the resonant denominator is done
! outside the routine because it involves complicated negotiation of
! the resonance to preserve accuracy for trapped particles. 
! So trapped dForceg needs to be divided by ()
! Prior Trapped bounces, simplified expression.
!       Lg=Lg - exp(sqm1*omegag*taug)/(1.+exp(sqm1*omegag*taug(ngz)))*Lg(ngz)
!       exptb=exp(sqm1*omegag*2.*taug(ngz))
       exptbb2=exp(sqm1*omegag*taug(ngz))
       dForceg=0.
       CapPhig(-ngz)=0.
       do i=-ngz+1,ngz
          exptau=exp(sqm1*omegag*taug(i))
          dtau=taug(i)-taug(i-1)
          Lgfactor=exp(sqm1*omegag*dtau) ! Current exponential
          CapPhig(i)=omegag*(exptbb2-1.)*(-(1.+exptbb2)*Lg(i)+exptau*Lg(ngz))
          if(.false.)then  ! New Force integration
             Lgdtau=(Lg(i-1)-(vmean-v0)/(-sqm1*omegag))*(Lgfactor-1)&
                  &/(sqm1*omegag) +(vmean-v0)/(-sqm1*omegag)*dtau
             expdtau=(Lgfactor-1)*exptau/Lgfactor/(sqm1*omegag)
             CapPhigdtau=omegag*(exptbb2-1.)*(-(1.+exptbb2)*Lgdtau&
                  &+expdtau*Lg(ngz))
             dForceg=dForceg-sqm1*0.5*(phigprime(i)+phigprime(i-1)) &
                  &*CapPhigdtau*abs(vg(0))
          else ! Old:
             dForceg=dForceg -sqm1*0.5*(CapPhig(i)*phigprime(i)+CapPhig(i-1)&
                  & *phigprime(i-1))*abs(vg(0))*dtau
          endif
       enddo
    endif
  end subroutine LofW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine testLofW
    use shiftmode
    integer, parameter :: nw=10
    complex, dimension(-ngz:ngz,nw) :: Lgw
    real, dimension(-ngz:ngz,nw) :: zgw,vgw,taugw
    complex :: dForceg(nw),forcet(nw)
    integer :: iwsa(nw)
    real :: Wn(nw)
    logical :: lplotmz=.true.
    omegag=(.5,0.)
    omegad=omegag
    psig=-.5
    psi=abs(psig)
    call initialize
    isigma=-1
    if(lplotmz)call pltinit(-zm,zm,min(0.,psig),max(1.8*psig,.8*abs(psig)))
    if(lplotmz)call axis
    if(lplotmz)call axlabels('z','W')
    do i=1,nw
       if(psig.gt.0)Wg=1.1*psig*i/nw
       if(psig.lt.0)Wg=psig+1.42*abs(psig)*i/nw
       write(*,*)'Wg=',Wg,' omegag=',omegag
       Wn(i)=Wg
       call LofW(Wg,isigma,dForceg(i))
       if(lplotmz)call polymark(zg(-ngz:ngz),Wg*ones(-ngz:ngz),2*ngz+1,i)
       iwsa(i)=min(iws,ngz)
       Lgw(-ngz:ngz,i)=Lg(-ngz:ngz) ! Save for plotting.
       zgw(-ngz:ngz,i)=zg(-ngz:ngz)
       vgw(-ngz:ngz,i)=vg(-ngz:ngz)
       taugw(-ngz:ngz,i)=taug(-ngz:ngz)
!       write(*,*)'dForceg=',dForceg(i),' taug',taug(ngz)
!       write(*,'(10f8.4)')(zg(j),j=-ngz,ngz)
!       write(*,'(10f8.3)')(taug(j),j=-ngz,ngz)
!       write(*,'(10f8.4)')(vg(j),j=-ngz,ngz)
!       write(*,'(10f8.4)')(real(Lg(j)),j=-ngz,ngz)
!       write(*,'(10f8.4)')(imag(Lg(j)),j=-ngz,ngz)
!       write(*,*)'dForceg=',dForceg(i),' taug',taug(ngz)
! Testing against shiftmode.
       if(Wg.lt.0)then
          vpsi=-isigma*sqrt(2.*(Wg-psig))
          call dFdvpsidvy(vpsi,forcet(i),tb,xlent)
          tdur=tb/2
          xLend=xlent
       elseif(psig.lt.0)then
          vinf=-isigma*sqrt(2.*Wg)
          xL=zm
          call initialize          
          call dFdvinfdvy(vinf,forcet(i))
          tdur=tau(nx)
          xLend=xL
       endif
       write(*,'(a, 5f10.5)')'End position        ',zg(ngz),xLend
       write(*,'(a, 5f10.4)')'Time duration       ',taug(ngz),tdur
       write(*,'(a, 5f10.6)')'Lg     Lt           ',Lg(ngz),Lt(iend+1)
       write(*,'(a, 5f10.6)')'Force real,imag     ',dForceg(i),forcet(i)
          
    enddo
    if(lplotmz)call pltend

    if(psig.lt.0)call pltinit(-zm,zm,0.,-isigma*1.5*sqrt(2.*abs(psig)))
    if(psig.ge.0)call pltinit(-zm,zm,-1.5*sqrt(2.*abs(psig)),1.5*sqrt(2.*abs(psig)))
    call axis
    call axlabels('zg','vg')
    do i=1,nw
       call color(i)
       call polyline(zgw(:,i),vgw(:,i),2*ngz+1)
       call polyline(zgw(:,i),-isigma*taugw(:,i)/taugw(ngz,i),2*ngz+1)
!       call polymark(zgw(:,i),taugw(:,i)/taugw(ngz,i),2*ngz+1,10)
       
       call polymark(zgw(iwsa(i),i),vgw(iwsa(i),i),1,1)
    enddo
    call color(15)
    call pltend
    
    call multiframe(2,1,3)
    call minmax2(real(Lgw),2*ngz+1,2*ngz+1,nw,amin,amax)
    call pltinit(-zm,zm,amin,amax)
    call axis
    call axlabels('z','Real(Lg)')
    do i=1,nw
       call color(i)
       call polyline(zgw(-ngz:ngz,i),real(Lgw(-ngz:ngz,i)),2*ngz+1)
       call polymark(zgw(-ngz:ngz,i),real(Lgw(-ngz:ngz,i)),2*ngz+1,i)
    enddo
    call color(15)
    call minmax2(imag(Lgw),2*ngz+1,2*ngz+1,nw,amin,amax)
    call pltinit(-zm,zm,amin,amax)
    call axis
    call axlabels('z','Imag(Lg)')
    do i=1,nw
       call color(i)
       call polyline(zgw(-ngz:ngz,i),imag(Lgw(-ngz:ngz,i)),2*ngz+1)
       call polymark(zgw(-ngz:ngz,i),imag(Lgw(-ngz:ngz,i)),2*ngz+1,i)
    enddo
    call pltend
    call multiframe(0,0,0)

    call minmax(dForceg,2*nw,amin,amax)
    call pltinit(Wn(1),Wn(nw),amin,amax)
    call axis
    call axlabels('W','Force')
    call polyline(Wn,real(dForceg),nw)
    call dashset(1)
    call polyline(Wn,imag(dForceg),nw)
    call pltend
  end subroutine testLofW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module shiftgen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use shiftgen
call testLofW
end program
