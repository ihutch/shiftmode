! Calculate Ft contribution for a range of vparallel and omega'.

program dFtdWs
  use shiftmode
  implicit none

  integer, parameter :: nomegad=1

  complex :: Ft(ne,nomegad),Fomega(nomegad)
  integer :: j,iwidth,ip0,lentrim
  real :: omegarmax,omegar(nomegad),omegai
  real :: pmin,pmax,pdum,xleg,yleg
  real :: dvdob(ne)
  character*30 string
  character*20 wvar
  logical :: limag=.true.
  real :: pwr
  
  psi=.16
  omegarmax=(nomegad/max(nomegad-1.,1.))*sqrt(psi)/2.    *.1
  omegai=.0002
  call initialize

  write(*,*)'case, omegar,   omegai,     Ftr       Fti'
  do j=1,nomegad
     omegar(j)=j*omegarmax/nomegad
     omegad=omegar(j)+sqm1*omegai
     omega=omegad   ! This is needed for 1-D (high Omega) case.
     ! It determines the balance of perp and parallel gradients.
     call FtEint(Fomega(j),-1./Ty,1.)
!     call FtEintOld(Fomega(j),-1./Ty,1.)
     write(*,'(i3,6f10.6)')j,omegar(j),omegai,Fomega(j)
     Ft(:,j)=Ftrap
  enddo

!  write(*,*)real(PhiInt(:,1))
!  write(*,'(10f8.4)')sWj-Wtscaled
  ip0=2
  call pfset(3)

  xleg=0.05
  yleg=0.48
  call minmax2(real(Ft),ne,ne,nomegad,pmin,pmax)
  if(limag)call minmax2(imag(Ft),ne,ne,nomegad,pdum,pmax)
  if(limag)xleg=.4
  if(limag)yleg=.8
  call pltinit(Wtscaled(ip0),Wtscaled(ne),pmin,pmax)
  call charsize(.018,.018)
  call axis
  call axis2
  call iwrite(iwpow,iwidth,string)
  wvar='(-W!d||!d)!u1/'//string(1:1)//'!u'
  call axlabels(wvar,'!BdF!dt!d/d!@'//wvar)
  call polyline((/Wtscaled(ip0),Wtscaled(ne)/),(/0.,0./),2)
  call legendline(xleg,yleg,258,'  !Aw!B!dr!d!@')
  string='!Aw!B!di!d!@='
  call fwrite(omegai,iwidth,6,string(lentrim(string)+1:))
  call legendline(xleg,.85,258,string)
  string='!Ay!@='
  call fwrite(psi,iwidth,2,string(lentrim(string)+1:))
  call legendline(xleg,.9,258,string)
  call charsize(.015,.015)
  do j=1,nomegad
     call color(j)
     call iwrite(j,iwidth,string)
     call labeline(Wtscaled(ip0),real(Ft(ip0:ne,j)),ne-ip0+1,string,iwidth)
     string(2:2)=' '
     call fwrite(omegar(j),iwidth,3,string(3:))
     call legendline(xleg,yleg-.05*j,258,string)
     if(limag)then
        call dashset(2)
        call polyline(Wtscaled(ip0),imag(Ft(ip0:ne,j)),ne-ip0+1)
        call dashset(0)
     endif
  enddo
  call color(15)
  call pltend()

  ! Plot omega_b versus sqrt(-W)=Wtscaled.
  call autoplot(Wtscaled/sqrt(psi),real(omegab)/sqrt(psi),ne)
  call axlabels('(-W/!Ay!@)!u0.5!u','!Aw!@!dbr!d/!A)y!@')
  call axis2
  string='!Ay!@='
  call fwrite(psi,iwidth,2,string(lentrim(string)+1:))
  call legendline(xleg,.9,258,string)
  call winset(.true.)
  call dashset(1)
  ! Formula
  pwr=0.45
  call polyline(Wtscaled/sqrt(psi), &
       sqrt(1/2.)/((sqrt(psi)/Wtscaled)**pwr+(sqrt(2.))**pwr-1.)**(1/pwr),ne)
  call dashset(4)
  call polyline(Wtscaled/sqrt(psi), &
       Wtscaled/sqrt(2.*psi),ne)
  call dashset(2)
  call polyline(Wtscaled/sqrt(psi), &
       (Wtscaled/(2.*sqrt(psi))),ne)
!  call dashset(3)
!  call polyline(Wtscaled/sqrt(psi), &
!       (Wtscaled/sqrt(psi))**.89/2.,ne)
  call dashset(0)
  call pltend

  if(.false.)then
  call autoplot(vpsiarray(ip0),real(omegab(ip0:)),ne-ip0)
  call axlabels('v!d!Ay!@!d','!Aw!@!dbr!d')
  call pltend()
  
  do j=ip0,ne-ip0
     dvdob(j)=(vpsiarray(j+1)**2-vpsiarray(j-1)**2)/(real(omegab(j+1))**2-real(omegab(j-1))**2)
  enddo
  call autoplot(real(omegab(ip0:)),dvdob(ip0:),ne/2)
  call axlabels('!Aw!@!dbr!d','dv!d!Ay!@!d!u2!u/d!Aw!@!dbr!d!u2!u')
  call pltend()
  endif

  call autoplot(psi-vpsiarray**2/2,4.7*real(omegab)**2.25/psi**.127,ne)
  call axis2
  call axlabels('!Ay!@-v!d!Ay!@!d!u2!u/2', &
       '!Aw!@!dbr!d!u2.25!u4.7/!Ay!@!u0.127!u')
  call legendline(xleg,.9,258,string)
  call pltend


!  call autoplot(vpsiarray(ip0:)**2,real(omegab(ip0:))**2,ne-ip0)
!  call lautoplot(vpsiarray(ip0:)**2,real(omegab(ip0:))**2,ne-ip0,.true.,.true.)
!  call axlabels('v!d!Ay!@!d!u2!u','!Aw!@!dbr!d!u2!u')  
!  call pltend
  
end program dFtdWs
