! Solve for omega as a function of Omegac for various k,
! and specified psi.
program osolvomk
  use iterfind
  integer, parameter :: nomax=100,nk=1
  real, dimension(nomax) :: oc,or,oi
  complex :: omegap
  integer :: isigma=-1,no=10
  real :: kp
  character*100 filename
  character*20 string

  psip=.09
  vs=1.
  kp=.1
  Ocmax=2.
  call parseos(psip,vs,kp,Ocmax,no)
  
  write(filename,'(a,2i2.2,a,i3.3,a,i3.3,a,i3.3,a,i3.3,a)')   &
       'F',no,nk, &
       'B',min(999,abs(nint(100*Ocmax))),   &
       'k',min(999,abs(nint(100*kp))),   &
       'v',min(999,abs(nint(100*vs))),  &
       'p',min(999,abs(nint(100*psip))),'.dat'

  write(*,'(a,$)')'Attempt opening  '//filename(1:lentrim(filename))
  open(12,file=filename,status='old',err=1)
  read(12,*)nof,ncol
  if(nof.ne.no)goto 1
  do i=1,nof
     read(12,*)oc(i),or(i),oi(i)
  enddo
  close(12)
  write(*,*)' Succeeded.'
  goto 2        ! Just plot the read-in data.
1  continue     ! No or bad old file. Calculate and save.
  close(12,status='delete')
  write(*,*)'  Failed. Calculating ...'

  oi=0;or=0.
  thek=kp*sqrt(psip)
  do i=1,no
     Omegacp=sqrt(psip)*Ocmax*(i-.7)/(no-.7)
     oc(i)=Omegacp
     write(*,'(a,f8.4,a,f8.4,a)')'Find for Omegac=',Omegacp,' k=',kp
     call iterfindroot(psip,vs,Omegacp,omegap,thek,isigma,nit)
     write(*,*)omegap,nit
     if(nit.gt.1.and.nit.le.niterfind)then
        or(i)=real(omegap)/sqrt(psip)
        oi(i)=imag(omegap)/sqrt(psip)
     endif
  enddo
  
  open(12,file=filename,status='new')
  write(12,*)no,2
  do i=1,no
     write(12,*)oc(i),or(i),oi(i)
  enddo
  close(12)

2 continue                  !Plot
  call charsize(.02,.02)
  call minmax(or,no,rmin,rmax)
  call minmax(oi,no,gmin,gmax)
  call pltinit(0.,oc(no),0.,1.1*max(rmax,gmax))
  call axis;call axis2
  call polyline(oc,or,no)
  call axlabels('!AW!@/!A)y!@','!Aw!@/!A)y!@')
  call legendline(.7,.15,0,'real')
  call dashset(2)
  call legendline(.7,.1,0,'imag')
  call polyline(oc,oi,no)
  call fwrite(kp,iwidth,3,string)
  call legendline(.2,1.05,258,'k/!A)y!@='//string(1:iwidth))
  call fwrite(vs,iwidth,2,string)
  call legendline(.6,1.05,258,'v!ds!d='//string(1:iwidth))
  call pltend

end program osolvomk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine parseos(psip,vsp,kp,Ocmax,no)
  real :: kp
  character*50 :: argument
  do i=1,iargc()
     call getarg(i,argument)
     if(argument(1:2).eq.'-p')read(argument(3:),*)psip
     if(argument(1:2).eq.'-k')read(argument(3:),*)kp
     if(argument(1:3).eq.'-vs')read(argument(4:),*)vsp
!     if(argument(1:3).eq.'-or')read(argument(4:),*)ormax
!     if(argument(1:3).eq.'-oi')read(argument(4:),*)oimax
     if(argument(1:3).eq.'-oc')read(argument(4:),*)Ocmax
     if(argument(1:2).eq.'-n')read(argument(3:),*)no
     if(argument(1:3).eq.'-Ti')then
        read(argument(4:),*)Ti
        write(*,'(a,f8.2)')'Setting Ti to',Ti
        call Tset(1.,Ti)
     endif
     if(argument(1:2).eq.'-c')ipfset=-3
     if(argument(1:2).eq.'-h')goto 1
  enddo
  return
1 continue
  write(*,'(a,f7.3)')'-p... set psip     [',psip
  write(*,'(a,f7.3)')'-k... set kp       [',kp
  write(*,'(a,f7.3)')'-vs.. set vs       [',vsp
  write(*,'(a,f7.3)')'-oc.. set Ocmax    [',Ocmax
  write(*,'(a,i5)')  '-n.. set Ocmax     [',no
  stop
end subroutine parseos
