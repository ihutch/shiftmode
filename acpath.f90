!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module acpath
  integer, parameter :: nc=4000
  real, dimension(nc) :: xcarray,ycarray
  real :: cvacpath
  integer :: iarray
contains
  ! replacement for accis dummy
subroutine acpathdoc(imax,xc,yc,cv)
  integer imax
  real xc(imax),yc(imax),cv
!  write(*,*)'Calling local acpathdoc'
  iarray=imax
  do i=1,imax
     xcarray(i)=xc(i)
     ycarray(i)=yc(i)
  enddo
  cvacpath=cv
end subroutine acpathdoc
end module acpath
