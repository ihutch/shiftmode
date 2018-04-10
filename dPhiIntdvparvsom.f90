! Calculate PhiInt for a range of vparallel and omega'.

program dphiint
  use shiftmode
  implicit none

  integer, parameter :: nomegad=5

  complex :: PhiInt(ne,nomegad),Ft(ne,nomegad),Fomega(nomegad)
  integer :: i,j,iwidth,ip0,lentrim
  real :: omegarmax,omegar(nomegad),omegai,dvpsi,vpsi(ne),Wj(ne),sWj(ne)
  real :: pmin,pmax
  character*20 string
  character*20 wvar
  logical :: lold=.false.
  
  psi=0.16
  omegarmax=(nomegad/(nomegad-1.))*sqrt(psi)/2.
  omegai=.004
  call initialize

  do j=1,nomegad
     omegar(j)=j*omegarmax/nomegad
     omegad=omegar(j)+sqm1*omegai
     if(lold)then
     iwpow=2         ! Equal W spacing is ipow=1
     do i=1,ne
        Wj(i)=-psi*((i-0.5)/ne)**iwpow
        sWj(i)=(-Wj(i))**(1./iwpow)
        vpsi(i)=sqrt(2.*(psi+Wj(i)))
        dvpsi=-sqrt(2.*psi*(1.-(float(i)/ne)**iwpow)) &
             +sqrt(2.*psi*(1.-(float(i-1)/ne)**iwpow))
        ! calculate the force PhiInt for this dvpsi and omegad element:
        call dFdvpsidvy(vpsi(i),PhiInt(i,j),tbe(i),xlen(i))
!        Ftrap(i)=Ftrap(i)*(omegad*dfe-(omegad-omega)*dfeperp)
        PhiInt(i,j)=PhiInt(i,j)*dvpsi
        ! We do not multiply by vpsi because that was done in dFdvpsidvy.
        ! But we multiply by 2. to account for \pm v_\psi. Maybe!!
     enddo
     endif
     ! Alternative call of shiftmode to do the above.
     call FtEint(Fomega(j),-1./Ty,1.)
     write(*,*)j,omegar(j),Fomega(j)
     Ft(:,j)=Ftrap
  enddo

!  write(*,*)real(PhiInt(:,1))
!  write(*,'(10f8.4)')sWj-Wtscaled
  ip0=2
  call pfset(3)
  
  if(lold)then

  call minmax2(real(PhiInt),ne,ne,nomegad,pmin,pmax)
  call pltinit(Wtscaled(ip0),Wtscaled(ne),pmin,pmax)
!  call scalewn(Wtscaled(ip0),Wtscaled(ne),pmin,pmax,.true.,.false.)
  call axis
  call iwrite(iwpow,iwidth,string)
  call axlabels('(-W!d||!d)!u1/'//string(1:1)//'!u','!AJFf!@''!Bdz!@')
  call polyline((/Wtscaled(ip0),Wtscaled(ne)/),(/0.,0./),2)
  call legendline(.05,1.-.05,258,'    !Aw!B!dr!d!@')
  do j=1,nomegad
     call color(j)
     call iwrite(j,iwidth,string)
     call labeline(Wtscaled(ip0),real(PhiInt(ip0:ne,j)),ne-ip0+1,string,iwidth)
     string(2:2)=' '
     call fwrite(omegar(j),iwidth,3,string(3:))
     call legendline(.05,.95-.05*j,258,string)
  enddo
  call color(15)
  string='!Aw!B!di!d!@='
  call fwrite(omegai,iwidth,3,string(lentrim(string)+1:))
  call legendline(.05,.05,258,string)
  string='!Ay!@='
  call fwrite(psi,iwidth,2,string(lentrim(string)+1:))
  call legendline(.05,.1,258,string)
  call pltend()

else

  call minmax2(real(Ft),ne,ne,nomegad,pmin,pmax)
  call pltinit(Wtscaled(ip0),Wtscaled(ne),pmin,pmax)
  call axis
  call iwrite(iwpow,iwidth,string)
  wvar='(-W!d||!d)!u1/'//string(1:1)//'!u'
  write(*,*)wvar
  call axlabels(wvar,'!BdF!dt!d/d!@'//wvar)
  call polyline((/Wtscaled(ip0),Wtscaled(ne)/),(/0.,0./),2)
  call legendline(.05,1.-.05,258,'    !Aw!B!dr!d!@')
  do j=1,nomegad
     call color(j)
     call iwrite(j,iwidth,string)
     call labeline(Wtscaled(ip0),real(Ft(ip0:ne,j)),ne-ip0+1,string,iwidth)
     string(2:2)=' '
     call fwrite(omegar(j),iwidth,3,string(3:))
     call legendline(.05,.95-.05*j,258,string)
  enddo
  call color(15)
  string='!Aw!B!di!d!@='
  call fwrite(omegai,iwidth,3,string(lentrim(string)+1:))
  call legendline(.05,.05,258,string)
  string='!Ay!@='
  call fwrite(psi,iwidth,2,string(lentrim(string)+1:))
  call legendline(.05,.1,258,string)
  call pltend()
   
  endif
  
end program dphiint
