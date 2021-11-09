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

! An orbit's equation of motion is dv/dt=-d\phi/dz, which means for
! species s of different mass, that the time is scaled differently: to
! omega_{ps}. Consequently, for a particular perturbation frequency
! omega, when there are multiple species the scaled value omegag must
! be set to omega/omega_{ps}, different for different mass
! species. For different charge sign, the hill peak psig must likewise
! be opposite. The length scale is normally the Debye length for some
! reference temperature, the default spatial extent is |zm|=10.
! But for energies near zero is sometimes automatically lengthened.

! The force is the integral dz of -d\phi/dz times the non-adiabatic
! perturbed f, which is an integral over the past time d\tau of the
! orbit, integrated over velocity, using total (distant) density
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
  complex :: omegag=(1.,0.),sqm1g=(0.,1.),Ftot,dFordirect
  complex :: omegadiff,omegaonly
  real :: psig=.1,Wg,zm=10.,v0,z0,z1,z2,zR,kg=0.,Omegacg=5.
  real :: vshift=0.,Tinf=1.,Tperpg=1.
  integer :: ivs,iws
  integer, parameter :: ngz=100,nge=200,nhmax=50
  integer :: iwpowg=2,nharmonicsg
  real,parameter :: pig=3.1415926
! Spatial Arrays
  real, dimension(-ngz:ngz) :: zg,vg,phig,phigprime,taug
  complex, dimension(-ngz:ngz) :: Lg,CapPhig
! Parallel energy arrays
  complex :: omegabg(0:nge)
  complex, dimension(nge) :: Forcegarray,Forcegp,Forcegr,Forcegt
  real :: Wgarray(nge),Wgarrayp(nge),Wgarrayr(nge),vinfarrayp(nge)&
       &,vinfarrayr(nge),tbr(nge),tbp(nge),Wgarrayu(nge)
  complex :: Fnonresg(0:nge),Ftrapg(0:nge)
! Perpendicular Harmonic force arrays.
  complex, dimension(-nhmax:nhmax) :: Frg,Ftg,Fpg
! Total forces
  complex :: Ftotalmode
  complex :: Ftotalrg,Ftotalpg,Ftotalsumg
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine makezg(isigma)
    integer :: isigma
! Calculate the zg mesh, and on it phig, vg, phigprime for an incoming
! orbit from z=isigma*zm (v sign -isigma) or trapped orbit from its
! isigma end.
    z0=0.
    zmfac=1.
! Make sure we don't miss a shallow trapped orbit
1    if(Wg.lt.0..and.Wg.gt.phigofz(zmfac*isigma*zm))then
       zmfac=1.2*zmfac
       goto 1
    endif
! Find any reflection position zR. ! Put equal to 0 if W>max(phi).
    call orbitendg(Wg,z0,zmfac*isigma*zm)
    zR=z0
    ivs=-1
    if(psig.gt.0)then ! Repelling potential
       gK=sqrt(max(0.,Wg/psig-1.))
       z1=zR; z2=isigma*zm-zR
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
          write(*,*)'phigprime NAN',i,phigprime(i),zg(i),z1,z2,zi,zR,&
               gK,Wg,isigma
          stop
       endif
    enddo
  end subroutine makezg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine orbitendg(Wj,z0,zL)
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
  end subroutine orbitendg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Fdirect(Wgi,isigma,dForceg)
! Calculate the past integral of tau and force 
! vs z for orbit of energy W, starting at isigma side. Where
! for untrapped particles (W>=0), v_0=vinf=sqrt(2*W).
! For trapped particles (W<0), v_0=sqrt(2*(W-psig)).  
! The integration is done on a z-grid such that dz=v*dtau.
! isigma is the sign of zg at the start of orbit (v-sign=-isigma)
! The total needs to account for v-sign=+isigma as well. However,
! for a symmetric hole that just multiplies force by 2. 
! Force integration not using v-integration by parts
! dF/dv/xi/(omega df/dW)=(i)\int -dphi/dz|_t \int^t (-dphi/dz|_tau)
!          exp(-i omega(tau-t)) dtau |dz(t)|        where dz(t)=v(t)dt.
! But if we want the contribution per dvinf, dv/dvinf=vinf/v so that 
! dF/dvinf/xi/(omega df/dW)=(i)\int -dphi/dz|_t \int^t (-dphi/dz|_tau)
!          exp(-i omega(tau-t)) dtau |vinf dt|.
! Inner integral is CapPhi=\int^t (-dphi/dz|_tau) exp(-i omega(tau-t)) dtau 
! done as exp(-i omega(tau-t)) dtau = d[exp()]/(-i omega)
! The result needs to be multiplied by df/dWpar.
    complex :: dForceg,CPfactor,exptbb2
    
    Wg=Wgi
    call makezg(isigma) ! Sets various time array values.
    vpsig=vg(-ngz)            ! Untrapped default ...
    if(Wg.lt.0.)vpsig=vg(0)   ! Trapped particle f(v) reference.
    ips=int(sign(1.,psig))
    iws=0                     ! dtau algorithm switch index (plotting only)
    taug(-ngz)=0.
    CapPhig(-ngz)=0.
    dForceg=0.
    if(Wg.lt.phigofz(zm).and.psig.gt.0)return  ! Reflected Energy too small
    do i=-ngz+1,ngz
       phigp=0.5*(phigprime(i)+phigprime(i-1))
       vmean=0.5*(vg(i)+vg(i-1))
! Calculate dtau
       if(zR.ne.0.and.abs(vg(i)).le.0.5*sqrt(2.*max(Wg,Wg-psig)))then
          dtau=-((vg(i)-vg(i-1))/phigp) ! Use dtau=dv*dtau/dv
          if(ips.gt.0)iws=i               ! Reflected track
       else                               ! Use dtau=dx*dtau/dx
          dtau=((zg(i)-zg(i-1))*vmean/(vg(i-1)*vg(i)))
          if(ips.le.0..or.zR.eq.0)iws=i   ! Attracted or unreflected
       endif
       if(.not.dtau.lt.1e6)then           ! Test for NAN error
          write(*,*)i,'dtau=',dtau,v0,vg(i-1),vg(i),zR,Wg,vpsig,phigp
          stop
       endif
       taug(i)=taug(i-1)+dtau
       CPfactor=exp(sqm1g*omegag*dtau) ! Current exponential
       CapPhig(i)=CPfactor*CapPhig(i-1)-phigp*(1.-CPfactor)*sqm1g/omegag
       dForceg=dForceg-sqm1g*&
            0.5*(CapPhig(i)*phigprime(i)+CapPhig(i-1)*phigprime(i-1))&
            *abs(vpsig*dtau)
    enddo
    if(Wg.lt.0)then     ! Trapped orbit. Add resonant term. 
       exptbb2=exp(sqm1g*omegag*taug(ngz))
       vpsig=abs(vg(0))
! This form is to be divided by (1-exptb) full resonant denominator.
       dForceg=dForceg*(1.-exptbb2**2) &
            + sqm1g*CapPhig(ngz)**2*(1.-exptbb2)*vpsig
! But it would probably be better to remove a (1.-exptbb2) factor and
! Divide only the resonant term by (1.+exptbb2) later Half Period:
!       dForceg=dForceg*(1.+exptbb2) + sqm1g*CapPhig(ngz)**2*vpsig
! The division by the resonant denominator is done outside the routine
! because it involves complicated negotiation of the resonance to
! preserve accuracy for trapped particles.
    endif
  end subroutine Fdirect
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine FgRepelEint(Ftotalg,isigma)
! Integrate the force over energy to obtain the full parallel distribution
! This version ignores possible dependence on v-perp (for now). 
! Just for positive (repelling) psig. isigma is the entering z-sign.
! For symmetric potentials and f(v), the returned Ftotalg can simply be
! doubled to give the total force since then it is symmetric in isigma.
    complex Ftotalg
    integer, save :: idone=0
    if(idone.eq.0.and.abs(omegaonly-omegag).gt.1.e-6)then
       write(*,*)'WARNING: omegaonly and omegag are unequal.'
       write(*,*)omegaonly,omegag
       idone=idone+1
    endif
    Emaxg=4.*Tinf+vshift**2
    call FgPassingEint(Ftotalpg,isigma,Emaxg)
    call FgReflectedEint(Ftotalrg,isigma)
    Ftotalg=Ftotalpg+Ftotalrg
  end subroutine FgRepelEint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine FgPassingEint(Ftp,isigma,Emaxg)
    implicit none
    complex Ftp
    integer :: isigma,i
    real :: Emaxg
    logical :: lcompare=.false.
    integer, parameter :: ippow=3
    real :: dvinf,fvinf,dfe,dfeperp,dfepre,dfeperppre
    omegadiff=omegag-omegaonly
    do i=1,nge  ! Passing, corrected for psig sign.
       Wgarray(i)=max(psig,0.)+Emaxg*(i/float(nge))**ippow
       vinfarrayp(i)=-isigma*sqrt(2.*Wgarray(i))
       call Fdirect(Wgarray(i),isigma,Forcegarray(i))
       dfe=dfdWpar(vinfarrayp(i),fvinf) ! fparallel slope and value
       dfeperp=-fvinf/Tperpg
       if(i.eq.1)then  ! Start of integration approximated
          dvinf=abs(vinfarrayp(i))-sqrt(2.*max(psig,0.))
          Ftp=dvinf*Forcegarray(i)*(omegag*dfe-omegadiff*dfeperp)
       else
          dvinf=abs(vinfarrayp(i)-vinfarrayp(i-1))
          Ftp=Ftp+dvinf*0.5*&
               (Forcegarray(i)*(omegag*dfe-omegadiff*dfeperp) &
               +Forcegarray(i-1)*(omegag*dfepre-omegadiff*dfeperppre))
       endif
       Forcegp(i)=Forcegarray(i)*(omegag*dfe-omegadiff*dfeperp)
       tbp(i)=taug(ngz)
       if(lcompare)call Fpasscompare(i)
       dfepre=dfe
       dfeperppre=dfeperp
    enddo
    Wgarrayp=Wgarray
  end subroutine FgPassingEint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine FgReflectedEint(Ftr,isigma)
    implicit none
    complex :: Ftr
    integer :: isigma,i
    real :: dvinf,fvinf,dfe,dfeperp,dfepre,dfeperppre
    omegadiff=omegag-omegaonly
    do i=1,nge  ! Reflected
       Wgarray(i)=psig*(1.-(i/float(nge))**2)
       vinfarrayr(i)=-isigma*sqrt(2.*Wgarray(i))
       call Fdirect(Wgarray(i),isigma,Forcegarray(i))
       dfe=dfdWpar(vinfarrayp(i),fvinf) ! fparallel slope and value
       dfeperp=-fvinf/Tperpg
       if(i.eq.1)then
          dvinf=abs(vinfarrayr(i))-sqrt(2.*psig)
          Ftr=dvinf*Forcegarray(i)*(omegag*dfe-omegadiff*dfeperp)
       else
          dvinf=abs(vinfarrayr(i)-vinfarrayr(i-1))
          Ftr=Ftr+dvinf*0.5* &
               (Forcegarray(i)*(omegag*dfe-omegadiff*dfeperp) &
               +Forcegarray(i-1)*(omegag*dfepre-omegadiff*dfeperppre))
       endif
       Forcegr(i)=Forcegarray(i)*omegag*dfdWpar(vinfarrayr(i),fvinf)
       tbr(i)=taug(ngz)
       dfepre=dfe
       dfeperppre=dfeperp
    enddo
    Wgarrayr=Wgarray
  end subroutine FgReflectedEint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine FgAttractEint(Ftotalg,isigma)
    complex Ftotalg
    integer, save :: idone=0
    if(idone.eq.0.and.abs(omegaonly-omegag).gt.1.e-6)then
       write(*,*)'WARNING: omegaonly and omegag are unequal.'
       write(*,*)omegaonly,omegag
       idone=idone+1
    endif
    Emaxg=4.*Tinf+vshift**2
    call FgPassingEint(Ftotalpg,isigma,Emaxg)
    call FgTrappedEint(Ftotalrg,-1/Tperpg,1.,isigma)
 ! Hacked dfperpdWperp, fperp, because only Maxwellian perp distrib.    
    Ftotalpg=2.*Ftotalpg ! Because both passing v-directions not yet done.
    Ftotalg=Ftotalpg+Ftotalrg
!    write(*,*)'Fs',Ftotalpg,Ftotalrg
  end subroutine FgAttractEint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine FgTrappedEint(Ftotal,dfperpdWperp,fperp,isigma)
    logical :: lcompare=.false.  ! Make True Only if using shiftmode
    ! Integrate over fe (trapped). Wj=vpsi^2/2-psi. So vpsi=sqrt(2(psi+Wj))
    ! Wj range 0 to -psi.
    complex :: Ftotal,exptb,cdvpsi,dob,dFdvpsig,resdenom,resdprev
    real :: dfperpdWperp,fperp,obi
    integer :: isigma

    omegadiff=omegag-omegaonly
    Ftotal=0.
    Wjprev=0.
    Ftotalmode=0.
    dfe=dfdWptrap(Wjprev,feprev)*fperp
    dfeperpprev=feprev*dfperpdWperp
    vpsiprev=sqrt(-2.*psig)
    omegabg(0)=0.
    Fnonresg(0)=0.                !Don't add zero energy point.
    resdprev=1.
    do i=1,nge-1       ! Trapped
       Wgarray(i)=psig*((float(i)/nge)**iwpowg)
       Wj=Wgarray(i)
! Functionalized version.
       dfe=dfdWptrap(Wj,fe)*fperp
       dfeperp=fe*dfperpdWperp      ! df/dW_perp
       vpsi=sqrt(2.*(-psig+Wj))
       vinfarrayr(i)=vpsi ! reflected==trapped for attracting hill.
       dvpsi=vpsiprev-vpsi
! Calculate the force dFdvpsi for this vpsi and dvy element for one transit:
       call Fdirect(Wgarray(i),isigma,dFdvpsig)
       Forcegarray(i)=dFdvpsig
       Forcegr(i)=Forcegarray(i)*(omegag*dfe-omegadiff*dfeperp)
       omegabg(i)=2.*pig/(2.*taug(ngz))
       call pathshiftg(i,obi)
       omegabg(i)=omegabg(i)+sqm1g*obi
       exptb=exp(sqm1g*omegag*pig/omegabg(i))   ! Half period
       if(.not.abs(exptb).lt.1.e10)exptb=1.e10  ! Enable divide by infinity
       resdenom=1.-exptb**2                        !Full period version
!       resdenom=1.+exptb                        ! Half period version
       if(i.eq.1)resdprev=resdenom
       dob=omegabg(i)-omegabg(i-1)
       cdvpsi=dvpsi*(1.+sqm1g*imag(dob)/real(dob))
       ! Strictly to get dFdvpsi we need to multiply by the omega f' terms
       Fnonresg(i)=dFdvpsig*(omegag*dfe-omegadiff*dfeperp)
       ! and correct for the imaginary shift of omegabg:
       Fnonresg(i)=Fnonresg(i)+sqm1g*real(Fnonresg(i)-Fnonresg(i-1))  &
            /real(omegabg(i)-omegabg(i-1))*obi
       if(taug(ngz)*imag(omegag).lt.-2.)then!Hack fix giant dFdvpsi problem
          write(*,*)'Drop',taug(ngz),imag(omegag),dFdvpsig
          Fnonresg(i)=0. 
       endif
! Then divide by the resonance denominator multiply by complex dvpsi and sum.
       Ftrapg(i)=0.5*(Fnonresg(i)/resdenom + Fnonresg(i-1)/resdprev)*cdvpsi
    ! Now Ftrapg(i) is a quantity when simply summed over all ne positions
    ! and multiplied by 2 gives the total force. 
       if(.not.(abs(Ftrapg(i)).ge.0))then
          write(*,*)'Ftrapg NAN?',i
          write(*,*)Fnonresg(i),Ftrapg(i)
          write(*,*)omegabg(i-1),omegabg(i)
          write(*,*)omegag,omegag/omegabg(i)
          write(*,*)resdenom,resdprev
          write(*,*)dvpsi,vpsi,vpsiprev,Wj
          write(*,*)fe,dfe,dfeperp
          stop
       endif
       ! Multiply by 2. to account for \pm v_\psi.
       Ftotal=Ftotal+2.*Ftrapg(i)       ! Add to Ftotal integral.
       if(lcompare)call  Ftrapcompare(i,dfe,dfeperp,obi,resdenom&
            &,resdprev,cdvpsi,Ftotal,dFdvpsig,vpsi)
       vpsiprev=vpsi
       feprev=fe
       dfeprev=dfe
       dfeperpprev=dfeperp
       resdprev=resdenom
    enddo
    ! Calculate end by extrapolation.
    Ftrapg(nge)=Ftrapg(nge-1)+0.5*(Ftrapg(nge-1)-Ftrapg(nge-2))
    Fnonresg(nge)=Fnonresg(nge-1)
    Ftotal=Ftotal+2.*Ftrapg(nge)
    Wgarray(nge)=psig
    Wgarrayr=Wgarray
    if(lcompare)call  Ftrapcompare(i,dfe,dfeperp,obi,resdenom&
         &,resdprev,cdvpsi,Ftotal,dFdvpsig,vpsi)
  end subroutine FgTrappedEint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine FgEint(Ftotalg,isigma)
    complex :: Ftotalg
    if(psig.gt.0)then
       call FgRepelEint(Ftotalg,isigma)
    else
       call FgAttractEint(Ftotalg,isigma)
    endif
  end subroutine FgEint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine SumHarmonicsg(isigma)
    implicit none
    integer :: isigma
    ! Sum harmonic contributions to obtain Forces for specified
    ! kg, Omegacg. Does the job of v-perp integration. 
    real :: EIm(0:nhmax),xit,vymax,hnum
    integer :: m,ncalc
    real :: Oceff   ! The effective Omegacg
! The maximum needed perp velocity, such that kg*vymax/Oc=nharmonicsg
    vymax=3.5*sqrt(2.*Tperpg)
    Oceff=max(1.e-6,max(Omegacg,kg*vymax/nhmax)) ! Don't allow zero Oceff.
! How many harmonics do we actually need? Nominally:
    hnum=kg*vymax/Oceff
! If nharmonicsg is small, then use rather more for accuracy.
    nharmonicsg=min(nhmax,int(hnum*(1.+3./(hnum+1.))))
    xit=kg*sqrt(Tperpg)/Oceff
!    write(*,*)vymax,hnum,xit
! Calculate the Integer[0.] exp*I[2] Bessel functions 0 to nharmonicsg
    call RIBESL(xit**2,0.,nharmonicsg+1,2,EIm,ncalc)
    if(.not.ncalc.eq.nharmonicsg+1)then ! All orders not calculated correctly.  
       write(*,'(a,i3,a,i3,a)')'Bessel functions',ncalc,' of',nharmonicsg+1, &
       ' are precise. Others may be irrelevant.' 
    endif
! m=0 always used.! fy is Maxwellian.
! But the fywy,fy need to be fixed in inner routines.
    omegaonly=omegag
    call FgEint(Fpg(0),isigma)
!    write(*,*)'SumHarmonicsg: Oceff,nharmonicsg=',Oceff,nharmonicsg !,hnum
    write(*,'(a,e11.4,''('',2e12.4,'')'',i4)')&
         ' EIm(0),Ftt(0)  =',EIm(0),Fpg(0),nharmonicsg
    Ftotalsumg=Fpg(0)*EIm(0)
    do m=1,nharmonicsg
       omegag=omegaonly+m*Oceff
       if(abs(omegag).lt.1.e-5)omegag=omegag+sqm1g*1.e-5 ! Prevent zero.
       call FgEint(Fpg(m),isigma)
       if(real(omegaonly).eq.0)then   !Short cut.
          Ftotalsumg=Ftotalsumg+2*real(Fpg(m))*EIm(m)
       else      ! Full sum over plus and minus m.
          Ftotalsumg=Ftotalsumg+Fpg(m)*EIm(m)
          omegag=omegaonly-m*Oceff
          if(abs(omegag).lt.1.e-5)omegag=omegag+sqm1g*1.e-5
          call FgEint(Fpg(-m),isigma)
          Ftotalsumg=Ftotalsumg+Fpg(-m)*EIm(m)
       endif
       write(*,'(a,e11.4,''('',2e12.4,'')('',2e12.4,'')'')')&
            ' EIm(m),Fpg(+-m)=',EIm(m),Fpg(m),Fpg(-m)
    enddo
    omegag=omegaonly
  end subroutine SumHarmonicsg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function dfdWpar(vinf,fvinf)
! The derivative of the distant distribution func wrt energy W at 
! velocity vinf, in the rest frame of the hole. fvinf is returned the f(v). 
! Example. Maxwellian of temperature T. f=exp(-vinv^2/2T)/sqrt(2 pi T)
!    dfdWpar=-exp(-vinf**2/(2.*Tinf))/Tinf/sqrt(2.*3.1415926*Tinf)
! Two half-density Maxwellians symmetrically shifted by vshift
! f=0.5*(exp(-(vinv-vshift)^2/2T)+exp(-(vinv+vshift)^2/2T))/sqrt(2 pi T))
    e1=exp(-(vinf-vshift)**2/(2.*Tinf))
    e2=exp(-(vinf+vshift)**2/(2.*Tinf))
    fvinf=0.5*(e1+e2)/sqrt(2.*3.1415926*Tinf)
    dfdWpar=0.5*(-e1*(vinf-vshift)-e2*(vinf+vshift)) &
         &/sign(max(abs(vinf*Tinf),1.e-6),vinf)/sqrt(2.*3.1415926*Tinf)
  end function dfdWpar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! If some different form of phi is required. Replace these functions,
! e.g. with forms that call functions outside the module. 
  real function phigofz(zval)
    phigofz=psig/cosh(zval/4.)**4
  end function phigofz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function phigprimeofz(zval)
    phigprimeofz=-psig*sinh(zval/4.)/cosh(zval/4.)**5
  end function phigprimeofz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function dfdWptrap(Wj,fe)
! Return the energy derivative and the parallel distribution function
! This simple version ignores charge contribution from repelled species.
    sqWj=sqrt(-Wj)
    fe=((2./pig)*sqWj+(15./16.)*Wj/sqrt(-psig)+experfcc(sqWj)&
       /sqrt(pig))/sqrt(2.)        ! New exact sech^4 hole form.
    dfdWptrap=((15./16.)/sqrt(-psig)-experfcc(sqWj)/sqrt(pig))/sqrt(2.)
end function dfdWptrap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fvinfplot
! vshift, Tinf, defined in module shiftgen,     
    integer, parameter :: ngpl=100
    real, dimension(ngpl) :: vplinf,fplinf
    character*10 string
    vpmax=4.*sqrt(Tinf)+vshift
    do i=1,ngpl
       vplinf(i)=-vpmax+2.*i/ngpl*vpmax
       blah=dfdWpar(vplinf(i),fvinf)
       fplinf(i)=fvinf   ! This is really f, not v.
    enddo
    call autoplot(vplinf,fplinf,ngpl)
    call axlabels('v','f!di!A;!@!d')
    call fwrite(vshift,iwidth,2,string)
    call legendline(0.1,0.9,258,'v!ds!d='//string(1:iwidth))
    call pltend
  end subroutine fvinfplot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pathshiftg(i,obi)
! Calculate the required shift of the omegabg path below the real axis
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
  include 'LofW.f90'    ! Obsolete; only for testing.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module shiftgen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The standard routines for comparison are:
include 'compareroutines.f90'
! They can be replaced by dummies for compilation without shiftmode.
