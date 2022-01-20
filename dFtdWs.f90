! Calculate Ft contribution for a range of vparallel and omega'.

program dFtdWs
  use shiftmode
  implicit none

  integer, parameter :: nomegad=9

  complex :: Ft(ne,nomegad),Fomega(nomegad),Fbroad
  integer :: j,iwidth,ip0,lentrim,ip1,i
  real :: omegarmax,omegar(nomegad),omegai
  real :: pmin,pmax,pdum,xleg,yleg
  real :: bestfit(ne)
  character*30 string
  character*20 wvar
  logical :: limag=.false.
  real :: pwr,cw,zB,fd,fn
  real :: plotx(ne),ploty(ne)
  
  psi=.16
  omegarmax=(nomegad/max(nomegad-1.,1.))*sqrt(psi)/2.    *.1
  omegarmax=.225
  omegai=.001
  call initialize

  write(*,*)'case, omegar,   omegai,     Ftr       Fti    Fbroadr  Fbroadi  '
  do j=1,nomegad
     omegar(j)=j*omegarmax/nomegad
     omegad=omegar(j)+sqm1*omegai
     omega=omegad   ! This is needed for 1-D (high Omega) case.
     ! It determines the balance of perp and parallel gradients.
     call FtEint(Fomega(j),-1./Ty,1.)
     ! Calculate the contribution from the broad part above resonance.
     Fbroad=0.
     ip1=int(2.*omegar(j)*ne/Wtscaled(ne))
     do i=ip1,ne
        Fbroad=Fbroad+Ftrap(i)
     enddo
     write(*,'(i3,9f10.6)')j,omegar(j),omegai,Fomega(j),Fbroad
 !,real(Fbroad)/imag(Fbroad)
     Ft(:,j)=Ftrap*2.*ne/Wtscaled(ne) ! Convert to true dFt/d(-W).
     ! So now we need to multiply Ft by dWtscaled and sum to integrate.
  enddo

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
  if(nomegad.gt.1)then
     call legendline(xleg,yleg,258,'  !Aw!B!dr!d!@')
  else
     call fwrite(omegar(1),iwidth,4,string)
     call legendline(xleg,yleg,258,'!Aw!B!dr!d!@='//string)
  endif
  string='!Aw!B!di!d!@='
  call fwrite(omegai,iwidth,4,string(lentrim(string)+1:))
  call legendline(xleg,.85,258,string)
  string='!Ay!@='
  call fwrite(psi,iwidth,2,string(lentrim(string)+1:))
  call legendline(xleg,.9,258,string)
  call charsize(.015,.015)
  do j=1,nomegad
     call color(j)
     if(nomegad.gt.1)then
        call iwrite(j,iwidth,string)
        call labeline(Wtscaled(ip0),real(Ft(ip0:ne,j)),ne-ip0+1,string,iwidth)
        string(2:2)=' '
        call fwrite(omegar(j),iwidth,3,string(3:))
        call legendline(xleg,yleg-.05*j,258,string)
     else
        call polyline(Wtscaled(ip0),real(Ft(ip0:ne,j)),ne-ip0+1)
        call legendline(xleg,yleg-.05,0,' real')
     endif
     if(limag.and.j.eq.1)then
        call dashset(2)
        call color(4)
        call polyline(Wtscaled(ip0),imag(Ft(ip0:ne,j)),ne-ip0+1)
  ! Scale the imaginary part to show it comes from real part*2oi/or.
        call legendline(xleg,yleg-.1,0,' imaginary')
        call color(5)
        call dashset(4)
        call polyline(Wtscaled(ip1),imag(Ft(ip1:ne,j))*.5*omegar(j)/omegai &
             ,ne-ip1)
        call legendline(xleg,yleg-.15,0,' imaginary!AXw!B!dr!d!@/2!Aw!B!di!d!@')
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
  pwr=.5
  ploty=sqrt(1/2.)/((sqrt(psi)/Wtscaled)**(pwr)+(sqrt(2.))**pwr-1.)**(1/pwr)
  call polyline(Wtscaled/sqrt(psi),ploty,ne)
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

  zB=1.5
  cw=(4*zB/3.14159)
  fn=cw/(sqrt(2.)-1.)-1.
  fd=cw*sqrt(2.)/(sqrt(2.)-1.)-1.
  bestfit=(1.+fn*Wtscaled/sqrt(2.*psi))/ &
       (1.+fd*Wtscaled/sqrt(2.*psi))
  call autoplot(Wtscaled/sqrt(psi),real(omegab(1:))/sqrt(psi)&
       /(Wtscaled/sqrt(2.*psi)),ne-1)
  call axlabels('(-W/!Ay!@)!u0.5!u','omegab/sqrt(-2W), Power Fit, Ratio')
  call dashset(2)
  call polyline(Wtscaled/sqrt(psi),real(omegab(1:))/sqrt(psi)/ploty,ne-1)
!  call polyline(Wtscaled/sqrt(psi),ploty/(Wtscaled/sqrt(2.*psi)),ne-1)
  call dashset(3)
  call polyline(Wtscaled/sqrt(psi),bestfit,ne-1)
  call dashset(0)
  call pltend

  if(.false.)then
  call autoplot(psi-vpsiarray**2/2,4.7*real(omegab(1:ne))**2.25/psi**.127,ne-1)
  call axis2
  call axlabels('!Ay!@-v!d!Ay!@!d!u2!u/2', &
       '!Aw!@!dbr!d!u2.25!u4.7/!Ay!@!u0.127!u')
  call legendline(xleg,.9,258,string)
  call pltend
  endif

  plotx=Wtscaled/psi
  ploty=4.7*real(omegab(1:ne))**2.25/psi**.127/Wtscaled**2
  call autoplot(plotx,ploty,ne-1)
  call axis2
  call axlabels('(-W/!Ay!@)!u0.5!u', &
       '!Aw!@!dbr!d!u2.25!u4.7/!Ay!@!u0.127!u/(-W/!Ay!@)')
    call legendline(xleg,.9,258,string)
  call pltend

  

!  call autoplot(vpsiarray(ip0:)**2,real(omegab(ip0:))**2,ne-ip0)
!  call lautoplot(vpsiarray(ip0:)**2,real(omegab(ip0:))**2,ne-ip0,.true.,.true.)
!  call axlabels('v!d!Ay!@!d!u2!u','!Aw!@!dbr!d!u2!u')  
!  call pltend
  
end program dFtdWs
