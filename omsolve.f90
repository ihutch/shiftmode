! Given psi and k calculate the forces for a range of omegas.

include 'fhgfunc.f'

program main
  use shiftmode 
  integer, parameter ::   nk=15
  real :: kik(nk),omik(nk)
!  complex :: Fcpassing(nk),Ftrapped(nk)
  !  integer, parameter :: np=7
    integer, parameter :: np=5
  real, dimension(np,nk) :: Ftnp,Fpnp
  character*20 :: string
  real :: omega0,omega1,oval,omax,komax
  external forcebalance,growthlk

  Ty=0.01
  psimax=1.
  pgmax=1. ! Silence warnings.
  ovalmax=0.2*sqrt(psi)/sqrt(Ty) ! Old Guess at peak omegai
  do ip=1,np
     decade=(ip-1)/3
     ipmod=mod(ip,3)
     if(ipmod.eq.1)then
        psi=psimax/10**decade
     elseif(ipmod.eq.2)then
        psi=0.5*psimax/10**decade
     else
        psi=0.2*psimax/10**decade
     endif
     sqpsi=sqrt(psi)
     omax=-999.
     call initialize
     if(ip.eq.1)then
        pgmax=growthlk(psi,beta,Ty)
        call pfset(3)
        pkmax=0.3/sqrt(Ty)
        ovalmax=pkmax*pgmax/4. ! Guess at peak omegai
!        call pltinit(0.,pkmax,0.,0.08/Ty**0.33)
        call pltinit(0.,pkmax,0.,pkmax*pgmax/4.)
        call axis()
        call axis2()
        call axlabels('!Bk!@/!A)y!@','!Ag!@/!A)y!@')
        call accisflush()
        xg=.75
        yg=.94
        string='!Ay!@'
        call legendline(xg,yg,257,string)
        string='!BT!dy!d='
        call fwrite(Ty,width,2,string(10:))
        call legendline(.1,yg,257,string)
        call polyline([0.,.1]*pkmax,[0.,.1]*pgmax*pkmax,2)
        call polyline([0.,.1]*pkmax,[.0625,.0625],2)
     endif
     ! k-scan
     akmax=0.3*sqrt(psi)/sqrt(Ty)
     akmin=akmax/100.
     write(*,*)'   k      omegai   Fpassing                Ftrapped'
     ikmax=nk
     do ik=1,nk
        k=akmin+(ik-1.)*(akmax-akmin)/(max(nk-1.,1.))
        !     k=max(akmin,(ik-1.)*akmax/(max(nk-1.,1.)))
        kik(ik)=k
        omega1=min(4.*k/sqrt(Ty),1.3*ovalmax)
        omega0=omega1*(k*(1-akmax)+akmin)
        if(np.gt.1)then
           aknorm=k/akmax
! NG           omega1=ovalmax/sqpsi*max(aknorm*(1.-aknorm),0.1*aknorm)
           omega0=0.1/sqpsi*ovalmax*max(aknorm*(1.-aknorm),0.03*aknorm)
!           write(*,'(a,4f8.4)')'aknorm,omega0,omega1,ovalmax' &
!                ,aknorm,omega0,omega1,ovalmax
        endif
        nbi=10
        call bisectroot(omega0,omega1,forcebalance,nbi,frac,oval)
        if(frac.gt.0)then
           write(*,'(f7.5,5es12.3)')k,oval,Fpasstotal,Ftraptotal 
           omega=complex(0.,oval)
           so=abs(omega**2)
           omik(ik)=oval
           !           write(*,*)'oval=',oval/sqpsi,' ovalmax=',ovalmax/sqpsi
           if(oval.gt.omax)then
              omax=oval
              komax=k
           endif
       else
           omik(ik)=0.
           ikmax=ik-1
           exit
        endif
     enddo
     if(oval.ne.-999.)then
        ovalmax=omax   ! Subsequently use the previous peak.
     endif

!     write(*,*)'psi,komax,omax=',psi,komax,omax
!     write(*,*)'prediction=    ',psi,sqrt(psi)/8.,sqrt(psi)/16.
!     write(*,*)'16*omax/sqrt(psi)=',16.*omax/sqrt(psi)
!     write(*,*)'omik(2)/kik(2)=',omik(2)/kik(2),' omik(3)/kik(3)=',omik(3)/kik(3)
!     write(*,*)'omax/komax=',omax/komax
     call color(ip)
     call dashset(ip)
     call polyline(kik/sqpsi,omik/sqpsi,ikmax)
     string=' '
     call fwrite(psi,width,2,string(2:))
     call legendline(xg,yg-.05*ip,0,string)
     call dashset(0)
     call polyline((/0.,.1/)*pkmax,(/0.,.1/)*growthlk(psi,beta,Ty)*pkmax,2)
     call accisflush()
  enddo
!  call drcstr('1/16')
  call pltend()
   
end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function forcebalance(omega0,omega1,frac,omega2)
  use shiftmode
  implicit none
  real, intent(in) :: omega0,omega1,frac
  real, intent(out) :: omega2
  omega2=(1-frac)*omega0+frac*omega1
  omega=complex(0.,omega2)
  call FpVyint()
  call FtVyint()
  !  forcebalance=real(Ftraptotal)+real(Fpasstotal)
  ! Include F_E :
  forcebalance=real(Ftraptotal)+real(Fpasstotal)-(psi*k)**2*128./315.
end function forcebalance

include 'bisectroot.f' 
