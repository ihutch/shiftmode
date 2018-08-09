! Verify the coding by comparing with f,g etc. At zero k, omega_r and
! small omega_i for a range of psis.

include 'fhgfunc.f'

program main
  use shiftmode 
  integer, parameter ::   nk=1
  real :: kik(nk)
  complex :: Fcpassing(nx),Ftrapped(nk)
  integer, parameter :: np=20
!  integer, parameter :: np=5
  real :: psinp(np),hnp(np),gnp(np),hmgnp(np),gjnp(np),pnp(np)
  real, dimension(np,nk) :: Ftnp,Fpnp
!  character*10 :: string
  psistep=.01*20/np
  do ip=1,np  ! Iterate over psi.
     psi=psistep*ip
     psinp(ip)=psi
! k-scan disabled here because nk=1 use the minimum.
     akmax=sqrt(psi)/8.
     akmin=.000
!     write(*,*)'   k     Fpassing                Ftrapped'
     do ik=1,nk
        k=akmin+(ik-1.)*(akmax-akmin)/max(nk-1.,1.)
        kik(ik)=k
        omega=sqrt(psi)*(0.001+k/4.)*sqm1
        so=-imag(omega)**2
        so=real(omega**2)   ! the same here.
        call initialize
        write(*,*)'nx, ne, nvy,   xL,    pL,   omegar,  omegai,     k     psi   beta'
        write(*,'(3i4,7f8.4)')nx,ne,nvy,xL,pL,real(omega),imag(omega),k,psi,beta
     !  call passingdiags
        idebug=0
        call FpVyint()
        idebug=0
        Fcpassing(ik)=Fpasstotal
        call FtVyint()
        Ftrapped(ik)=Ftraptotal
        write(*,'(a,f8.4,4es12.3)')' k,Fpass,Ftrap',k,Fcpassing(ik),Ftraptotal
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
  call axis2
  call axlabels('!Ay!@','Force !BF!dt!d , F!dp!d!@ (/!Aw!@!u2!u!AD!@)')
  call winset(.true.)
  do ik=1,nk
!     write(*,*)ik,Ftnp(:,ik)
!     write(*,*)ik,Fpnp(:,ik)
     call polyline(psinp,Ftnp(:,ik),np)
     call polyline(psinp,Fpnp(:,ik),np)
     if(ik.eq.1)then
        call jdrwstr(wx2nx(psinp(np/2)),wy2ny(Ftnp(np/2,ik)), &
             '!BF!dt!d!@ Trapped',-1.)
        call jdrwstr(wx2nx(psinp(np/2)),wy2ny(Fpnp(np/2,ik)), &
             '!BF!dp!d!@ Passing',-1.)
     endif
!     call fwrite(kik(ik),iwdth,4,string)
!     call jdrwstr(wx2nx(psinp(2*np/3)),wy2ny(Ftnp(2*np/3,ik)),string,1.)
  enddo
  call legendline(.65,.58,0.,'Numerical')
  call color(3)
  call dashset(1)
  call legendline(.65,.5,0.,'Analytic')
  call polyline(psinp,gnp,np)
  call polyline(psinp,hmgnp,np)
  !     call polyline(psinp,gjnp,np) ! trapped force w/o boundary term
  !     call polyline(psinp,pnp,np)  ! passing force w/o boundary term
  call dashset(0)
  call color(15)
  call pltend()
   
end program main
