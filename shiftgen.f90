! This is a general calculation of the force arising from an
! linearized oscillatory displacement (shift) perturbation of a
! localized structure that gives potential ENERGY phi. The potential,
! which is prescribed by two functions phigofz(z) and phigprimeofz(z),
! is either a hill or a valley and tends to zero at large distances,
! but its derivative passes through zero only once, at z=0, which is
! the only potential extremum (discounting infinity). When the extremum
! is a minimum, there are locally trapped particles, but when it is a
! maximum, reflected orbits are not locally trapped. There are
! therefore three types of orbit: Passing, Trapped, or Reflected,
! which must be treated differently.

! The orbit's equation of motion is dv/dt=-d\phi/dz, which means for
! species s of different mass, that the time is scaled differently: to
! omega_{ps}. Consequently, for a particular perturbation frequency
! omega, when there are multiple species the scaled value omegag must
! be set to omega/omega_{ps}, different for different mass
! species. For different charge sign, the hill peak psig must likewise
! be opposite. The length scale is normally the Debye length for some
! reference temperature, the default spatial extent is |zm|=10.  The
! units of returned force

! The force is the integral dz of -d\phi/dz times the non-adiabatic
! perturbed f, which is an integral over the past time d\tau of the
! past orbit, integrated over velocity, using total (distant) density
! of unity, and given (linearized) for a unit perturbing z-shift.

! Both space and past time integrals can be expressed as time
! integrals.  However, equal intervals of neither time nor space are
! universally optimal choices. Equal intervals of space are suboptimal
! near a reflection. Equal intervals of time are suboptimal for orbits
! of passing energy close to zero (because they move slowly at large z
! which is a less important region).  The z array is z(i)= z1+
! z2/(1+2K)[2K(i/n)+(i/n)^2] when repelled, z1=zR, z2=(zm-zR), where
! zR is reflection point.  When attracted z1=0, z2=zR trapped, z2=zm
! passing.  And K=sqrt(max(0,W/psi-1)) repelled,
! K=-1-sqrt(max(0,-W/psi)) attracted.  Repelling hill psi>0, attracted
! valley psi<0.

module shiftgen
  integer, parameter :: ngz=100,nge=50
  real, dimension(-ngz:ngz) :: zg,vg,ones=1.,phig,phigprime,taug
  complex, dimension(-ngz:ngz) :: Lg,CapPhig
  complex :: omegag=(1.,0.),sqm1=(0.,1.),Ftot
  complex :: omegabg(nge),Forcegarray(nge),Forcegp(nge),Forcegr(nge)
  real :: Wgarray(nge),Wgarrayp(nge),Wgarrayr(nge),vinfarrayp(nge)&
       &,vinfarrayr(nge),tbr(nge),tbp(nge)
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
    zR=z0
    ivs=-1
    if(psig.gt.0)then ! Repelling potential
       gK=sqrt(max(0.,Wg/psig-1.))
       z1=zR; z2=isigma*zm-zR     ! Or maybe z2=zm
       if(Wg.lt.psig)ivs=1  ! Reflected orbit, all z are same sign.
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
    if(Wg.eq.0)then
       write(*,*)'ZeroWg',z1,z2,gK,zR
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
       if(.not.phigprime(i).lt.1.e20)then
          write(*,*)'phigprime NAN',i,phigprime(i),zg(i),z1,z2,zi,zR,gK,Wg,isigma
          stop
       endif
    enddo
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
!    if(sign(1.,enrgy1).eq.sign(1.,enrgy0))then
    if(enrgy1*enrgy0.ge.0)then
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
  subroutine LofW(Wgi,isigma,dForceg)
  ! Calculate the past integral of tau and Lg(tau)
  ! vs z for orbit of energy W, starting at isigma side. Where
  ! Lg= \int_{-\infty}^{t} (v-v_0)*exp(-i*omegag(tau-t)dtau.
  ! For untrapped particles (W>=0), v_0=vinf=sqrt(2*W).
  ! For trapped particles (W<0), v_0=sqrt(2*(W-psig)).  
  ! The integration is done on a z-grid such that dz=v*dtau.
  ! isigma is the sign of zg at the start of orbit (v-sign=-isigma)
  ! The total needs to account for v-sign=+isigma as well. However,
  ! for a symmetric hole that just multiplies force by 2. 

  ! Integrate Lg dt (=dz/v) to get the differential
  ! force density with respect to vy and vinf when multiplied by
  ! df/dW_parallel.
    complex :: dForceg,Lgfactor
    complex :: exptbb2

    Wg=Wgi
    call makezg(isigma)
    if(psig.gt.0.)then                ! Repelling hill. 
       if(Wg.gt.phigofz(zm))then      ! Sufficiently positive Wg
          v0=-isigma*sqrt(2.*Wg)
          vpsig=abs(v0)
       else                           ! Negative W or orbit beyond z-domain.
          write(*,*)'Orbit beyond z. Wg=',Wg
          dForceg=0.
          Lg=0.
          return
       endif
    else                               ! Attracting valley.
       if(Wg-psig.lt.0)stop 'Wg below minimum potential ERROR'
       v0=0.
       vpsig=abs(vg(0))
    endif
!    write(*,*)'v0,vg(-ngz)',v0,vg(-ngz),vg(-ngz+1)
    taug(-ngz)=0.
    Lg(-ngz)=0.
    dForceg=0.
    ips=int(sign(1.,psig))
    iws=0
    do i=-ngz+1,ngz
       phigp=0.5*(phigprime(i)+phigprime(i-1))
       vmean=0.5*(vg(i)+vg(i-1))
       if(zR.ne.0.and.abs(vg(i)).le.0.5*sqrt(2.*max(Wg,Wg-psig)))then
          dtau=-((vg(i)-vg(i-1))/phigp) ! Use dtau=dv*dtau/dv
          if(ips.gt.0)iws=i               ! Reflected track
       else                               ! Use dtau=dx*dtau/dx
          dtau=((zg(i)-zg(i-1))*vmean/(vg(i-1)*vg(i)))
          if(ips.le.0..or.zR.eq.0)iws=i   ! Attracted or unreflected
       endif
       if(.not.dtau.lt.1e6)then
          write(*,*)i,'dtau=',dtau,v0,vg(i-1),vg(i),zR,Wg,vpsig,phigp
          stop
       endif
       if(dtau.gt.10)write(*,*)'JUMP?',phigp,vg(i),vmean
       taug(i)=taug(i-1)+dtau
       Lgfactor=exp(sqm1*omegag*dtau) ! Current exponential
       Lg(i)=Lgfactor*Lg(i-1)-(vmean-v0)*(1.-Lgfactor)/omegag
! Simple version without step integral correction.       
       dForceg=dForceg-sqm1*omegag* 0.5*(Lg(i)*phigprime(i)+Lg(i-1)&
               &*phigprime(i-1))*vpsig*dtau
       if(.not.real(dForceg).lt.1.e6)then
          write(*,*)'real(dForceg)',real(dForceg)
          write(*,*)i,Wg,phigp,vpsig,dtau,vg(i),vg(i-1),zR
          stop
       endif
    enddo
    if(Wg.lt.0)then     ! Trapped orbit. Add resonant term. 
       exptbb2=exp(sqm1*omegag*taug(ngz))
! This form is to be divided by (1-exptb) full resonant denominator.
! But it would probably be better to remove a (1.-exptbb2**2) factor and
! Divide only the resonant term by (1.+exptbb2) later.
       dForceg=dForceg*(1.-exptbb2**2) &
            + sqm1*Lg(ngz)**2*omegag**2*(1.-exptbb2)*vpsig
! In shiftmode the division by the resonant denominator is done
! outside the routine because it involves complicated negotiation of
! the resonance to preserve accuracy for trapped particles. 
    elseif(Wg.lt.psig)then  ! Reflected orbit. Add integration correction.
       dForceg=dForceg+sqm1*2.*vg(-ngz)**2*vpsig
    endif
  end subroutine LofW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine FgRepelEint(Ftotalg,isigma)
! Integrate the force over energy to obtain the full parallel distribution
! This version ignores possible dependence on v-perp (for now). 
! Just for positive (repelling) psig.
    complex Ftotalg,Ftotalpg,Ftotalrg
    Emaxg=5.
    do i=1,nge  ! Passing
       Wgarray(i)=psig+Emaxg*(i/float(nge))**2
       vinfarrayp(i)=-isigma*sqrt(2.*Wgarray(i))
       call LofW(Wgarray(i),isigma,Forcegarray(i))
       if(i.eq.1)then
          dvinf=vinfarrayp(i)-sqrt(2.*psig)
          Ftotalpg=Forcegarray(i)*dfdWpar(vinfarrayp(i))&
               &*dvinf*omegag
       else
          dvinf=vinfarrayp(i)-vinfarrayp(i-1)
          Ftotalpg=Ftotalpg+0.5*(Forcegarray(i)*dfdWpar(vinfarrayp(i)) &
               +Forcegarray(i-1)*dfdWpar(vinfarrayp(i-1)))*dvinf*omegag
       endif
       Forcegp(i)=Forcegarray(i)*omegag*dfdWpar(vinfarrayp(i))
       tbp(i)=taug(ngz)
    enddo
    Wgarrayp=Wgarray
    do i=1,nge  ! Reflected
       Wgarray(i)=psig*(1.-(i/float(nge))**2)
       vinfarrayr(i)=-isigma*sqrt(2.*Wgarray(i))
       call LofW(Wgarray(i),isigma,Forcegarray(i))
       if(i.eq.1)then
          dvinf=(vinfarrayr(i)-sqrt(2.*psig))
          Ftotalrg=Forcegarray(i)*dvinf*dfdWpar(vinfarrayr(i))&
               &**omegag
       else
          dvinf=vinfarrayr(i)-vinfarrayr(i-1)
          Ftotalrg=Ftotalrg+0.5*(Forcegarray(i)*dfdWpar(vinfarrayr(i)) &
               +Forcegarray(i-1)*dfdWpar(vinfarrayr(i-1)))*dvinf*omegag
       endif
       Forcegr(i)=Forcegarray(i)*omegag*dfdWpar(vinfarrayr(i))
       tbr(i)=taug(ngz)
    enddo
    Wgarrayr=Wgarray
    Ftotalg=Ftotalpg+Ftotalrg
  end subroutine FgRepelEint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function dfdWpar(vinf)
! The derivative of the distant distribution func wrt energy W at 
! velocity vinf, in the rest frame of the hole. Example. Maxwellian
! of temperature T. f=exp(-vinv^2/2T)/sqrt(2 pi T)
    Tinf=1.
    dfdWpar=-exp(-vinf**2/(2.*Tinf))*vinf/Tinf/sqrt(3.1415926*Tinf)
  end function dfdwpar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pathshiftg(i,obi)
! Calculate the required shift of the omegab path below the real axis
! for this ith energy mesh value, assuming the previous omegabg(i-1)
! is known.
    real :: obi
    integer :: el,ielsign
    real :: doel,doel2,dob
    real,parameter :: Sc=4,Rc=2./Sc

    ! Select the closest odd resonance such that el*omegabg=omegar
    ielsign=int(sign(1.,real(omegag)))
    el=(2*int( (abs(real(omegag))/real(omegabg(i))-1)/2. )+1)*ielsign
    doel=real(omegabg(i))-real(omegag)/el
    doel2=real(omegabg(i))-real(omegag)/(el+2*ielsign)
    if(abs(doel2).lt.abs(doel))then
       el=el+2*ielsign
       doel=doel2
    endif
! Calculate the required omegabg imaginary part, which is that el*obi 
! must be at least |el|dob/Rc below imag(omegag) if omegabg is closer to
! the real resonance than Sc(=4) times the omegabg step size.
    dob=real(omegabg(i)-omegabg(i-1))
    obi=-ielsign*max(0.,dob/Rc-imag(omegag)/abs(el)) &
         *max(0.,1.-(doel/(Sc*dob))**2)
  end subroutine pathshiftg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module shiftgen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine testLofW
    use shiftgen
    use shiftmode
    integer, parameter :: nw=60
    complex, dimension(-ngz:ngz,nw) :: Lgw
    real, dimension(-ngz:ngz,nw) :: zgw,vgw,taugw
    complex :: dForceg(nw),forcet(nw)
    integer :: iwsa(nw)
    real :: Wn(nw)
    logical :: lplotmz=.true.
    omegag=(.9,0.02)
    psig=.5
    omegad=omegag
    psi=abs(psig)
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
       if(lplotmz)call color(mod(i-1,15)+1)
       if(lplotmz)call polymark(zg(-ngz:ngz),Wg*ones(-ngz:ngz),2*ngz+1,10)
       iwsa(i)=min(iws,ngz)
       Lgw(-ngz:ngz,i)=Lg(-ngz:ngz) ! Save for plotting.
       zgw(-ngz:ngz,i)=zg(-ngz:ngz)
       vgw(-ngz:ngz,i)=vg(-ngz:ngz)
       if(Wg.gt.0.and.Wg.lt.psig)Lgw(:,i)=Lg+(vg-vg(-ngz))/omegag
       taugw(-ngz:ngz,i)=taug(-ngz:ngz)
!       write(*,*)'dForceg=',dForceg(i),' taug',taug(ngz)
!       write(*,'(10f8.4)')(zg(j),j=-ngz,ngz)
!       write(*,'(10f8.3)')(taug(j),j=-ngz,ngz)
!       write(*,'(10f8.4)')(vg(j),j=-ngz,ngz)
!       write(*,'(10f8.4)')(real(Lg(j)),j=-ngz,ngz)
!       write(*,'(10f8.4)')(imag(Lg(j)),j=-ngz,ngz)
!       write(*,*)'dForceg=',dForceg(i),' taug',taug(ngz)
! Testing against shiftmode. We must use positive velocity in the
! shiftmode calculation because it is not set up to use negative
! integration direction. The sign change is then in the print out. 
      if(Wg.lt.0)then          
          vpsig=sqrt(2.*(Wg-psig))
          call dFdvpsidvy(vpsig,forcet(i),tb,xlent)
          tdur=tb/2
          xLend=xlent
       elseif(psig.lt.0)then
          vinf=sqrt(2.*Wg)
          xL=zm
          call initialize          
          call dFdvinfdvy(vinf,forcet(i))
          tdur=tau(nx)
          xLend=xL
       endif
       write(*,'(a, 5f10.5)')'End position        ',zg(ngz),xLend
       write(*,'(a, 5f10.4)')'Time duration       ',taug(ngz),tdur
       write(*,'(a, 5f10.6)')'Lg     Lt           ',Lg(ngz),-isigma*Lt(iend+1)
       write(*,'(a, 5f10.6)')'Force real,imag     ',dForceg(i),forcet(i)
          
    enddo
    if(lplotmz)call pltend

    call pfset(3)
    call multiframe(2,1,0)
    if(psig.lt.0)call pltinit(-zm,zm,0.,-isigma*1.5*sqrt(2.*abs(psig)))
    if(psig.ge.0)call pltinit(-zm,zm,-1.5*sqrt(2.*abs(psig)),1.5*sqrt(2.*abs(psig)))
    call axis
    call axlabels('zg','vg')
    do i=1,nw
       call color(mod(i-1,15)+1)
       call polyline(zgw(:,i),vgw(:,i),2*ngz+1)
!       call polyline(zgw(:,i),-isigma*taugw(:,i)/taugw(ngz,i),2*ngz+1)
!       call polymark(zgw(:,i),taugw(:,i)/taugw(ngz,i),2*ngz+1,10)
       
       call polymark(zgw(iwsa(i),i),vgw(iwsa(i),i),1,1)
    enddo
    call color(15)

    call minmax2(taugw,2*ngz+1,2*ngz+1,nw,tmin,tmax)
    call pltinit(-zm,zm,0.,tmax)
    call axis
    call axlabels('zg','tau')
    do i=1,nw
       call color(mod(i-1,15)+1)
       call polyline(zgw(:,i),taugw(:,i),2*ngz+1)
    enddo
    call pltend

    call multiframe(2,1,3)
    call minmax2(real(Lgw),2*ngz+1,2*ngz+1,nw,amin,amax)
    call pltinit(-zm,zm,amin,amax)
    call axis
    call axlabels('z','Real(Lg)')
    do i=1,nw
       call color(mod(i-1,15)+1)
       call polyline(zgw(-ngz:ngz,i),real(Lgw(-ngz:ngz,i)),2*ngz+1)
!       call polymark(zgw(-ngz:ngz,i),real(Lgw(-ngz:ngz,i)),2*ngz+1,i)
    enddo
    call color(15)
    call minmax2(imag(Lgw),2*ngz+1,2*ngz+1,nw,amin,amax)
    call pltinit(-zm,zm,amin,amax)
    call axis
    call axlabels('z','Imag(Lg)')
    do i=1,nw
       call color(mod(i-1,15)+1)
       call polyline(zgw(-ngz:ngz,i),imag(Lgw(-ngz:ngz,i)),2*ngz+1)
!       call polymark(zgw(-ngz:ngz,i),imag(Lgw(-ngz:ngz,i)),2*ngz+1,i)
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
    call dashset(0)
    call pltend
  end subroutine testLofW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine testFrepel
    use shiftgen
    complex :: Ftotalg
    omegag=(1.,0.01)
    psig=.5
    isigma=-1    
    call FgRepelEint(Ftotalg,isigma)
    
    call pltinit(vinfarrayr(nge),vinfarrayp(nge),0.,Wgarrayp(nge))
    call axis
    call polymark(vinfarrayp,Wgarrayp,nge,1)
    call polyline(vinfarrayp,Wgarrayp,nge)
    call polyline(vinfarrayr,Wgarrayr,nge)
    call polymark(vinfarrayr,Wgarrayr,nge,2)
    call polyline(vinfarrayp,tbp/10.,nge)
    call polyline(vinfarrayr,tbr/10.,nge)
    call legendline(.2,.9,0,' tb')
    call axlabels('vinf','Wg')
    call pltend
!    write(*,*)forcegr,forcegp
    call minmax(forcegp,2*nge,pmin,pmax)
    call minmax(forcegr,2*nge,rmin,rmax)
    call pltinit(vinfarrayr(nge),vinfarrayp(nge),min(pmin,rmin),max(pmax,rmax))
    call axis
    call axlabels('vinf','force')
    call polyline(vinfarrayr,real(forcegr),nge)
    call polymark(vinfarrayr,real(forcegr),nge,1)
    call polyline(vinfarrayr,imag(forcegr),nge)
    call polymark(vinfarrayr,imag(forcegr),nge,2)
    call polyline(vinfarrayp,real(forcegp),nge)
    call polymark(vinfarrayp,real(forcegp),nge,1)
    call polyline(vinfarrayp,imag(forcegp),nge)
    call polymark(vinfarrayp,imag(forcegp),nge,2)
    call legendline(.6,.7,1,' real')
    call legendline(.6,.8,2,' imag')
    call pltend
  end subroutine testFrepel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  call testLofW
  call testFrepel
end program
