! Verify the coding by comparing with f,g etc. At zero k, omega_r and
! small omega_i for a range of psis.

include 'fhgfunc.f'

program main
  use shiftmode 
  integer, parameter ::   nk=6
  real :: kik(nk)
  complex :: Fcpassing(nx),Ftrapped(nk)
  integer, parameter :: np=20
!  integer, parameter :: np=2
  real :: psinp(np),hnp(np),gnp(np),hmgnp(np),gjnp(np),pnp(np)
  real, dimension(np,nk) :: Ftnp,Fpnp
  character*10 :: string
  omega=(0.0,.002)
  psistep=.01
  so=-imag(omega)**2
  do ip=1,np  ! Iterate over psi.
     psi=psistep*ip
     psinp(ip)=psi
     call initialize
     write(*,*)'nx, ne, nvy,   xL,    pL,   omegar,  omegai,     k     psi   beta'
     write(*,'(3i4,7f8.4)')nx,ne,nvy,xL,pL,real(omega),imag(omega),k,psi,beta
     !  call passingdiags

! k-scan
     akmax=.002
     akmin=.000
     write(*,*)'   k     Fpassing                Ftrapped'
     do ik=1,nk
        k=akmin+(ik-1.)*(akmax-akmin)/(max(nk-1.,1.))
        kik(ik)=k
        call dentcalc2
        ! Integrate phi'*n-tilde dx.  
        Fcpassing(ik)=0.
        do i=1,nx
           Fcpassing(ik)=Fcpassing(ik)+phiprime(i)*dent(i)*dx
        enddo
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
  call axlabels('!Ay!@','Force !BF!dt!d , F!dp!d!@ (/!Aw!@!u2!u)')
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
   
end program main
