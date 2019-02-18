! Calculate Ft contribution for a range of vparallel and omega'.

program dFtdWs
  use shiftmode
  implicit none

  integer, parameter :: nomegad=9

  complex :: Ft(ne,nomegad),Fomega(nomegad)
  integer :: j,iwidth,ip0,lentrim
  real :: omegarmax,omegar(nomegad),omegai
  real :: pmin,pmax
  character*20 string
  character*20 wvar
  
  psi=0.16
  omegarmax=(nomegad/(nomegad-1.))*sqrt(psi)/2.
  omegai=.0001
  call initialize

  do j=1,nomegad
     omegar(j)=j*omegarmax/nomegad
     omegad=omegar(j)+sqm1*omegai
     call FtEint(Fomega(j),-1./Ty,1.)
!     call FtEintOld(Fomega(j),-1./Ty,1.)
     write(*,*)j,omegar(j),Fomega(j)
     Ft(:,j)=Ftrap
  enddo

!  write(*,*)real(PhiInt(:,1))
!  write(*,'(10f8.4)')sWj-Wtscaled
  ip0=2
  call pfset(3)

  call minmax2(real(Ft),ne,ne,nomegad,pmin,pmax)
  call pltinit(Wtscaled(ip0),Wtscaled(ne),pmin,pmax)
  call charsize(.018,.018)
  call axis
  call axis2
  call iwrite(iwpow,iwidth,string)
  wvar='(-W!d||!d)!u1/'//string(1:1)//'!u'
  write(*,*)wvar
  call axlabels(wvar,'!BdF!dt!d/d!@'//wvar)
  call polyline((/Wtscaled(ip0),Wtscaled(ne)/),(/0.,0./),2)
  call legendline(.05,.48,258,'  !Aw!B!dr!d!@')
  string='!Aw!B!di!d!@='
  call fwrite(omegai,iwidth,3,string(lentrim(string)+1:))
  call legendline(.05,.85,258,string)
  string='!Ay!@='
  call fwrite(psi,iwidth,2,string(lentrim(string)+1:))
  call legendline(.05,.9,258,string)
  call charsize(.015,.015)
  do j=1,nomegad
     call color(j)
     call iwrite(j,iwidth,string)
     call labeline(Wtscaled(ip0),real(Ft(ip0:ne,j)),ne-ip0+1,string,iwidth)
     string(2:2)=' '
     call fwrite(omegar(j),iwidth,3,string(3:))
     call legendline(.05,.48-.05*j,258,string)
  enddo
  call color(15)
  call pltend()
   
  
end program dFtdWs
