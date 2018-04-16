! Calculate for an array of k and psi the trapped and passing forces.
! And F_E, and plot them for chosen omega_i, psimax, akmax, etc.

include 'fhgfunc.f'

program main
  use shiftmode 
  integer, parameter ::   nk=20
  real :: kik(nk)
  complex :: Fcpassing(nx),Ftrapped(nk)
  integer, parameter :: np=11
!  integer, parameter :: np=2
  real :: psinp(np),hnp(np),gnp(np),hmgnp(np),gjnp(np),pnp(np)
  real, dimension(np,nk) :: Ftnp,Fpnp
  character*10 :: string
  omega=(0.0,.060)
  psimax=1.
  psistep=psimax/(np-1)
  akmax=3.*imag(omega)
!  akmax=1.2*imag(omega)
  akmin=.000
  flower=-.3
  fupper=1.2
  so=-imag(omega)**2  ! equivalent.
  so=real(omega**2)   ! equivalent.
  ! so is a Negative quantity equal to omega**2, since omega_r=0.
  ! \ddot{U}=(i\omega_i)^2*\Delta= so*\Delta. Therefore normalizing to
  ! so is normalizing to a force from negative \ddot{U}, since Delta
  ! is +ve. Still puzzling.

  do ip=1,np  ! Iterate over psi.
     psi=psistep*(ip-1)
     if(ip.eq.1)psi=0.2*psistep
     psinp(ip)=psi
     call initialize
     write(*,*)'nx, ne, nvy,   xL,    pL,   omegar,  omegai,     k     psi   beta'
     write(*,'(3i4,7f8.4)')nx,ne,nvy,xL,pL,real(omega),imag(omega),k,psi,beta

! k-scan
     write(*,*)'   k     Fpassing                Ftrapped'
     do ik=1,nk
        k=akmin+(ik-1.)*(akmax-akmin)/(max(nk-1.,1.))
        kik(ik)=k
        call FpVyint()
        Fcpassing(ik)=Fpasstotal
        call FtVyint()
        Ftrapped(ik)=Ftraptotal
        write(*,'(f6.4,4es12.3)')k,Fcpassing(ik),Ftraptotal
        Ftnp(ip,ik)=real(Ftraptotal)/so
        Fpnp(ip,ik)=real(Fcpassing(ik))/so
     enddo

     call fhgfunc(psi,20.,100,pL,gave,have,gjave,pint)
  
     hnp(ip)=have
     gnp(ip)=gave
     gjnp(ip)=gjave
     pnp(ip)=-pint
     hmgnp(ip)=(have-gave)
  enddo

  call pfset(3)
  call pltinit(0.,np*psistep,-np*psistep*1.3,np*psistep*4.)
  call axis
  call axlabels('!Ay!@','Force !BF!dt!d , F!dp!d!@ (/!Aw!@!u2!u!AD!@)')
  call winset(.true.)
  do ik=1,nk
     call polyline(psinp,Ftnp(:,ik),np)
     call polyline(psinp,Fpnp(:,ik),np)
     if(ik.eq.1)then
        call jdrwstr(wx2nx(psinp(np/2)),wy2ny(Ftnp(np/2,ik)),'Trapped',-1.)
     elseif(ik.eq.2)then
        call jdrwstr(wx2nx(psinp(np/2)),wy2ny(Fpnp(np/2,ik)),'Passing',-1.)
     endif
     call fwrite(kik(ik),iwdth,4,string)
     call jdrwstr(wx2nx(psinp(2*np/3)),wy2ny(Ftnp(2*np/3,ik)),string,1.)
  enddo
  call color(3)
  call dashset(1)
  call polyline(psinp,gnp,np)
  call polyline(psinp,hmgnp,np)
  !     call polyline(psinp,gjnp,np) ! trapped force w/o boundary term
  !     call polyline(psinp,pnp,np)  ! passing force w/o boundary term
  call dashset(0)
  call color(15)
  call pltend()

  ! Force versus k.
!  call minmax(Ftnp(np,:),nk,Fmin,Fmax)
  !  call pltinit(-0.1*k,k,Fmin,Fmax)
  Fmax=Ftnp(np,1)+Fpnp(np,1)
  !  call pltinit(-0.1*k,k,-np*psistep*1.3,np*psistep*5.)
  call pltinit(-0.1*k/imag(omega),k/imag(omega),flower*Fmax,fupper*Fmax)
  call axis
  call axis2
  call axlabels('!Bkv!dt!d/!Aw!@!di!d',  &
       'Force !BF!dt!d+F!dp!d!@ (/!Aw!@!u2!u!AD!@)')
  call winset(.true.)
  call polyline((/-0.1*k/imag(omega),k/imag(omega)/),(/0.,0./),2)
       
  do ip=1,np     
     call fwrite(psinp(ip),iwidth,2,string)
     call polyline(kik/imag(omega),Ftnp(ip,:)+Fpnp(ip,:),nk)
     call jdrwstr(wx2nx(0.),wy2ny(Ftnp(ip,1)+Fpnp(ip,1)),string(1:iwidth),-1.)
     call fwrite(imag(omega),iwidth,3,string)
     call legendline(0.5,0.94,258,'!Aw!@!di!d='//string(1:iwidth))
     !     call polyline(kik,Fpnp(ip,:),nk)
     call dashset(2)
     call color(3)
     call polyline(kik/imag(omega),(psistep*max(ip-1.,0.02)*kik)**2*128/315/so,nk)
     call color(15)
     call dashset(0)
  enddo
  call jdrwstr(wx2nx(0.),wy2ny(1.1*(Ftnp(np,1)+Fpnp(np,1))),'!Ay!@=',-1.)
  call dashset(2)
  call color(3)
  call legendline(.05,.06,0,' !BF!dE!d!@ (/!Aw!@!u2!u)')
  call dashset(0)
  call pltend
  
end program main
