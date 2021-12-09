! Solve for omega as a function of Omegac for various k,
! and specified psi.
program osolvomk
  use iterfind
  integer, parameter :: nomax=100,nkmax=20
  real, dimension(nomax,nkmax) :: or,oi
  real, dimension(nomax) :: oc
  real, dimension(nkmax) :: kparray
  complex :: omegap
  character*100 filename
  character*20 string
  integer :: isigma=-1,no=10,nk=1
  real :: kp=.1,psip=0.09,Ocmax=2.,vs=1.

  call parseos(psip,vs,kp,Ocmax,no,nk)
  
  write(filename,'(a,2i2.2,a,i3.3,a,i3.3,a,i3.3,a,i3.3,a)')   &
       'F',no,nk, &
       'B',min(999,abs(nint(100*Ocmax))),   &
       'k',min(999,abs(nint(100*kp))),   &
       'v',min(999,abs(nint(100*vs))),  &
       'p',min(999,abs(nint(100*psip))),'.omk'

  write(*,'(a,$)')'Attempt opening  '//filename(1:lentrim(filename))//' '
  open(12,file=filename,status='old',err=1)
  read(12,*)nof,ncol,nkf
  if(nof.ne.no.or.nkf.ne.nk)goto 1
  do i=1,nof
     read(12,*)oc(i),(or(i,j),oi(i,j),j=1,nk)
  enddo
  read(12,*)(kparray(j),j=1,nk)
  close(12)
  write(*,*)' Succeeded.'
  goto 2        ! Just plot the read-in data.
1  continue     ! No or bad old file. Calculate and save.
  close(12,status='delete')
  write(*,*)'  Failed. Calculating ...'

  oi=0;or=0.
  do j=1,nk
     kparray(j)=kp*j/nk
     thek=kparray(j)*sqrt(psip)
  do i=1,no
     if(j.eq.1)oc(i)=Ocmax*(i-.7)/(no-.7)
     Omegacp=sqrt(psip)*oc(i)
     write(*,'(a,f8.4,a,f8.4,a)')'Find for Omegac=',oc(i),' k=',kparray(j)
     call iterfindroot(psip,vs,Omegacp,omegap,thek,isigma,nit)
     write(*,*)omegap,nit
     if(nit.gt.1.and.nit.le.niterfind)then
        or(i,j)=real(omegap)/sqrt(psip)
        oi(i,j)=imag(omegap)/sqrt(psip)
     endif
  enddo
  enddo
  
  open(12,file=filename,status='new')
  write(12,*)no,2*nk,nk
  do i=1,no
     write(12,*)oc(i),(or(i,j),oi(i,j),j=1,nk)
  enddo
  write(12,*)(kparray(j),j=1,nk)
  close(12)

2 continue                  !Plot
  call charsize(.02,.02)
  call minmax(or(1:no,1:nk),no*nk,rmin,rmax)
  call minmax(oi(1:no,1:nk),no*nk,gmin,gmax)
  omax=1.1*max(rmax,gmax)
  call pltinit(0.,oc(no),0.,omax)
  call axis;call axis2
  call axlabels('!AW!@/!A)y!@','!Aw!@/!A)y!@')
  call fwrite(psip,iwidth,2,string)
  call legendline(.2,1.05,258,'!Ay!@='//string(1:iwidth))
  call fwrite(vs,iwidth,2,string)
  call legendline(.6,1.05,258,'v!ds!d='//string(1:iwidth))
  call color(1)
  call jdrwstr(wx2nx(oc(no)*.7),wy2ny(max(.05*omax,or(no,1)-0.05&
       &*omax)),'real',0.)
  call color(4)
  call jdrwstr(wx2nx(oc(no)*.5),wy2ny(max(.05*omax,oi(no,1)+0.05&
       &*omax)),'imag',0.)
  call winset(.true.)
  do j=1,nk
     call dashset(j-1)
     call color(1)
     call polyline(oc,or(1,j),no)
     call color(4)
     call polyline(oc,oi(1,j),no)
     call color(15)
     call fwrite(kparray(j),iwidth,2,string)
     call legendline(.73,.97-j*.05,0,' '//string(1:iwidth))
  enddo
  call legendline(.76,.95,258,'k=')
  call pltend

end program osolvomk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine parseos(psip,vsp,kp,Ocmax,no,nk)
  real :: kp
  character*50 :: argument
  integer :: ipset=0
  do i=1,iargc()
     call getarg(i,argument)
     if(argument(1:2).eq.'-p')read(argument(3:),*)psip
     if(argument(1:2).eq.'-k')read(argument(3:),*)kp
     if(argument(1:3).eq.'-vs')read(argument(4:),*)vsp
!     if(argument(1:3).eq.'-or')read(argument(4:),*)ormax
!     if(argument(1:3).eq.'-oi')read(argument(4:),*)oimax
     if(argument(1:3).eq.'-oc')read(argument(4:),*)Ocmax
     if(argument(1:3).eq.'-no')read(argument(4:),*)no
     if(argument(1:3).eq.'-nk')read(argument(4:),*)nk
     if(argument(1:3).eq.'-Ti')then
        read(argument(4:),*)Ti
        write(*,'(a,f8.2)')'Setting Ti to',Ti
        call Tset(1.,Ti)
     endif
     if(argument(1:2).eq.'-w')then
        if(lentrim(argument).gt.2)then
           read(argument(3:),*,err=2)ipset
           call pfset(ipset)
        endif
     endif
2    if(argument(1:2).eq.'-h')goto 1
  enddo
  
  return
1 continue
  write(*,'(a,f7.3)')'-p... set psip     [',psip
  write(*,'(a,f7.3)')'-k... set kp       [',kp
  write(*,'(a,f7.3)')'-vs.. set vs       [',vsp
  write(*,'(a,f7.3)')'-oc.. set Ocmax    [',Ocmax
  write(*,'(a,i5)')  '-no.. set no       [',no
  write(*,'(a,i5)')  '-nk.. set nk       [',nk
  write(*,'(a,i5)')  '-w... write file?  [',ipset
  stop
end subroutine parseos
