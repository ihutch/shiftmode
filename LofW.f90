!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine LofW(Wgi,isigma,dForceg) !Obsolete. Retained only for testing
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
    complex :: dForceg,CPfactor
    complex :: exptbb2,ddF1,ddF2,ddF3
    
    Wg=Wgi
    v0=-isigma*sqrt(2.*Wg)
    vpsig=abs(v0)
    call makezg(isigma)
    if(psig.gt.0.)then                ! Repelling hill. 
       if(Wg.le.phigofz(zm))then      ! Negative W or orbit beyond z-domain.
!          write(*,*)'Orbit beyond z. Wg=',Wg
          dForceg=0.
          Lg=0.
          return
       endif
    else                               ! Attracting valley.
       if(Wg-psig.lt.0)stop 'Wg below minimum potential ERROR'
       if(Wg.lt.0.)then
          v0=0.
          vpsig=abs(vg(0))
       endif
    endif
!    write(*,*)'v0,vg(-ngz)',v0,vg(-ngz),vg(-ngz+1)
    taug(-ngz)=0.
    Lg(-ngz)=0.
    dForceg=0.
    dFordirect=0.
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
!       Lgfactor=exp(sqm1g*omegag*dtau) ! Current exponential
          CPfactor=exp(sqm1g*omegag*dtau) ! Current exponential
          CapPhig(i)=CPfactor*CapPhig(i-1)-phigp*(1.-CPfactor)*sqm1g/omegag
          ddF1=-sqm1g*0.5*(CapPhig(i)*phigprime(i)+CapPhig(i-1)&
               &*phigprime(i-1)) *abs(vg(-ngz)*dtau)
          ddF3=(CPfactor-1.)/(sqm1g*omegag)
          ddF2=-sqm1g*phigp*(ddF3*CapPhig(i-1)-phigp*(dtau-ddF3)&
               &*sqm1g/omegag) *abs(vg(-ngz))  ! Integral corrected.
          dFordirect=dFordirect+ddF2
          Lg(i)=CPfactor*Lg(i-1)-(vmean-v0)*(1.-CPfactor)/omegag
! Old simple version without step integral correction.
          ddF1=-sqm1g*omegag*0.5*(Lg(i)*phigprime(i)+Lg(i-1)&
               &*phigprime(i-1))*vpsig*dtau
! Integral corrected version is incorrect for reflected particles.
          ddF2=-sqm1g*phigp*(ddF3*Lg(i-1)*omegag+(vmean-v0)*(dtau-ddF3) &
               &*sqm1g/omegag) *abs(vg(-ngz))  ! Integral corrected.
          dForceg=dForceg+ddF1
!          dForceg=dForceg-sqm1g*omegag* 0.5*(Lg(i)*phigprime(i)+Lg(i-1)&
!               &*phigprime(i-1))*vpsig*dtau
!          write(*,*)dForceg,dFordirect
       if(.not.real(dForceg).lt.1.e6)then
          write(*,*)'real(dForceg)',real(dForceg)
          write(*,*)i,Wg,phigp,vpsig,dtau,vg(i),vg(i-1),zR
          stop
       endif
    enddo
    if(Wg.lt.0)then     ! Trapped orbit. Add resonant term. 
       exptbb2=exp(sqm1g*omegag*taug(ngz))
! This form is to be divided by (1-exptb) full resonant denominator.
! But it would probably be better to remove a (1.-exptbb2**2) factor and
! Divide only the resonant term by (1.+exptbb2) later.
       dForceg=dForceg*(1.-exptbb2**2) &
            + sqm1g*Lg(ngz)**2*omegag**2*(1.-exptbb2)*vpsig
! In shiftmode the division by the resonant denominator is done
! outside the routine because it involves complicated negotiation of
! the resonance to preserve accuracy for trapped particles. 
    elseif(Wg.lt.psig)then  ! Reflected orbit. Add integration correction.
       dForceg=dForceg+sqm1g*2.*vg(-ngz)**2*vpsig      !  *omegag
    endif
  end subroutine LofW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
