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

  real, dimension(0:ngz) :: zg,ones=1.
  real :: psig=.1,Wg,zm=10.
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine makezg
    z0=0.
! Find any reflection position zR. ! Put equal to 0 if W>max(phi).
    call orbitend(Wg,z0,zm)
!    write(*,*)Wg,psig,'Reflection?',z0
    zR=z0
    if(psig.gt.0)then ! Repelled
       gK=sqrt(max(0.,Wg/psig-1.))
       z1=zR; z2=zm-zR   ! Or maybe z2=zm
    elseif(psig.lt.0)then !Attracted
       gK=-1.-10.*sqrt(max(0.,-Wg/psig))
       if(zR.eq.0.)then
          z2=zm
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
    enddo
  end subroutine makezg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine orbitend(Wj,z0,zL)
    ! Find by bisection the turning point of the orbit whose energy is Wj.
    ! That is where phi=-Wj. But we return the value below it: x0.
    real :: Wj,zm,z1,z0,phim,phi0,phi1
    nbi=20
    pL=4.
    z0=0.
    z1=zL
    phi0=psig/cosh(z0/pL)**4-Wj
    phi1=psig/cosh(z1/pL)**4-Wj
    if(sign(1.,phi1).eq.sign(1.,phi0))then
!       write(*,*)'orbitend energies do not bracket zero',phi0,phi1
       return
    endif
    do i=1,nbi
!       write(*,*)i,z0,phi0,z1,phi1
       zm=(z1+z0)/2.
       phim=psig/cosh(zm/pL)**4-Wj
! Which value to replace.
       if(sign(1.,phim).eq.sign(1.,phi0))then
          if(z0.eq.zm)exit ! should never happen
          phi0=phim
          z0=zm
       else
          if(z1.eq.zm)exit
          phi1=phim
          z1=zm
       endif
    enddo
  end subroutine orbitend
!****************************************************************
  subroutine testmakezg
    psig=.5
    call pltinit(0.,zm,min(0.,psig),max(1.8*psig,.8*abs(psig)))
    call axis
    call axlabels('z','W')
    nw=10
    do i=1,nw
       if(psig.gt.0)Wg=1.5*psig*i/nw
       if(psig.lt.0)Wg=psig+1.5*abs(psig)*i/nw
       write(*,*)'Wg=',Wg
       call makezg
       call polymark(zg,Wg*ones,ngz+1,i)
       write(*,'(10f8.4)')(zg(j),j=0,ngz)
    enddo
    call pltend
  end subroutine testmakezg
!****************************************************************
  
end module shiftgen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use shiftgen
call testmakezg
end
