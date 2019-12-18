! Contour the force(s) over a domain of k and omegar.
! For psi and omegai fixed and a range of Omegac values
! Do fcontkor -h for option list.

program fcontkor
  use shiftmode 
  use acpath
  implicit none
  integer, parameter ::   nk=31,no=31
  integer :: nkf,nof
  real :: kik(nk),or(no)
  real, dimension(nk,no) :: cworka,FE
  complex, dimension(nk,no) :: Forcecomplex,Fpcomplex,Ftcomplex
  real :: zclv(20)
  character*30 string,filename
  integer :: i,io,ik,ioc,iwidth,nioc,icl,icsw,nocmax,noc1
  real :: akmax,akmin,Fmax,Fmin,omegacmax,omegacmin,omegarmax,omegarmin
  real :: so,Omegacp,oi,psip,Typ
  logical :: lreadit
  integer lentrim
  external lentrim
  
  psi=.64
  Ty=1.
  akmin=0.002                     ! Lowest k/sqrt(psi) plotted.
  akmax=0.7/sqrt(Ty)              ! range of k/sqrt(psi)
  omegarmax=.55*akmax**0.75       ! Estimated range of omega/sqrt(psi)
  omegarmin=0.0
  omegacmax=0.80                  ! Maximum Omegac/sqrt(psi)=Omegacp
  omegacmin=.10                   ! Minimum ditto
  nioc=8                          ! Number of Omegac cases.
  oi=.001
! To select a different set of cases adjust these parameters:
  noc1=1                  
  nocmax=nioc           ! Default case range
!  noc1=5;nocmax=8    ! Choose a particular case range
!!!!!!!!!!!!!!! Rudimentary argument reading/setting.
  do i=1,iargc()
     call getarg(i,string)
     if(string(1:2).eq.'-p')read(string(3:),*,err=201)psi
     if(string(1:2).eq.'-i')read(string(3:),*,err=201)noc1
     if(string(1:2).eq.'-f')read(string(3:),*,err=201)nocmax
     if(string(1:2).eq.'-h')goto 203
     goto 202
201  write(*,*)'Error reading from argument:',string
203  write(*,*)'Options: -p.. set psi; -i.. set initial Oc; -f.. set final'
     stop
202  continue
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  call pfset(3)
  call pltinit(0.,akmax,0.,omegarmax*sqrt(psi))
  call axis
  call axis2
  call axlabels('!Bk!@/!A)y!@','!Aw!B!dr!d!@/!A)y!@')
  string='!Ay!@='
  call fwrite(psi,iwidth,2, string(lentrim(string)+1:))
  call legendline(.8,1.03,258,string)
  string='!Aw!B!di!d!@='
  call fwrite(oi,iwidth,3, string(lentrim(string)+1:))
  call legendline(.6,1.03,258,string)
  call legendline(0.,1.03,258,'Labels !AW!@/!A)y!@')
  !  string='!BT!dy!d!@='
!  call fwrite(Ty,iwidth,2, string(lentrim(string)+1:))
!  call legendline(.8,.83,258,string)
!  call polyline((/0.,.1/)*akmax,(/0.,.1/)*slopelk*akmax,2)
  call initialize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Loop over Omegac values
  do ioc=noc1,nocmax
     Omegacp=(Omegacmin+(Omegacmax-Omegacmin)*(ioc-1.)/(nioc-1.))
     write(*,*)'Omegacp=',Omegacp
!!!!!!!!!!!!!!!!!!!!!!!
     ! Create filename in accordance with specified parameters
     write(filename,'(a,2i2.2,a,i3.3,a,i3.3,a,i3.3,a)')   &
          'F',nk,no, &
          'Oc',min(999,abs(nint(100*Omegacp))),   &
          'i',min(999,abs(nint(1000*oi))),  &
          'p',min(999,abs(nint(100*psi))),'.arr'
!     write(*,*)'Filename=',filename
!!!!!!!!!!!!!!!!!!!!!!!
     ! Try to open the file.
     lreadit=.false.
     open(9,file=filename,status='old',form='unformatted',err=101) 
     read(9,err=101)nkf,nof
     if(nkf.ne.nk.or.nof.ne.no)then
        write(*,*)'Reading from ',filename
        write(*,*)'File array dimensions',nkf,nof,' not compatible'
        write(*,*)'Adjust allocation or delete file'
        stop
     endif
     read(9,err=101)kik,or
     read(9,err=101)psip,Omegacp,Typ,oi,akmax,omegarmax
     read(9,err=101)forcecomplex,Ftcomplex,Fpcomplex
     close(9)
     lreadit=.true.
     write(*,*)'Read forcecomplex from file: ',filename
     goto 102
101  write(*,*)'Failed to open or read from file: ', filename,'Calculating...'
102  continue

     Omegac=Omegacp*sqrt(psi)
     if(.not.lreadit)then
!!!!!!!!!!!!!!!!!! Calculate the arrays.
        do ik=1,nk   ! k-scan
           kik(ik)=max(akmin,(ik-1)*akmax/(nk-1))
           k=kik(ik)*sqrt(psi)     
           do io=1,no  ! Iterate over omegar.
              or(io)=(omegarmin+(io-1.)*(omegarmax-omegarmin)/(no-1.))*sqrt(psi)
              omega=complex(or(io),oi)
              so=abs(omega**2)
              call SumHarmonics()
              Fpcomplex(ik,io)=Fpasstotal
              Ftcomplex(ik,io)=Ftraptotal
              FE(ik,io)=(psi*k)**2*128./315.
              Forcecomplex(ik,io)=Fpasstotal+Ftraptotal-FE(ik,io)
           enddo
        enddo
!!!!!!!!!!!!!!!!!!!
     
!!!!!!!!!!!!!!!!!!!! Attempt to write but skip if file exists.
        open(9,file=filename,status='new',form='unformatted',err=103)
        write(*,*)'Opened new file: ',filename,' and writing'
        write(9)nk,no
        write(9)kik,or
        write(9)psip,Omegacp,Ty,oi,akmax,omegarmax
        write(9)forcecomplex,Ftcomplex,Fpcomplex
        goto 104
103     write(*,*)'New File: ',filename,' cannot be opened; not rewriting.'
104     close(9)
     endif

! Plot this case data.  
     call minmax2(Forcecomplex,nk,nk,no,Fmin,Fmax)
     icl=1
     zclv(1)=Omegacp
     icsw=1
     ! Fool contour into plotting the zero level but labeling it with Omegac.
     call color(2)
     if(noc1.eq.nocmax)icl=-1
     call contourl(real(Forcecomplex)+Omegacp, &
          cworka,nk,nk,no,zclv,icl,kik,or,icsw)
     if(noc1.ne.nocmax)call legendline(0.1,-.1,0,'real')
     call color(4)
     call contourl(imag(Forcecomplex)+Omegacp, &
          cworka,nk,nk,no,zclv,icl,kik,or,icsw)
     if(noc1.ne.nocmax)call legendline(.6,-.1,0,'imag')
     call accisflush()
     if(noc1.eq.nocmax)then ! Single case: plot more contours
        icsw=1
        call color(2)
        call dashset(1)
        icl=0
        zclv(1)=10
        call contourl(real(Forcecomplex), &
             cworka,nk,nk,no,zclv,icl,kik,or,icsw)
        call legendline(0.1,-.1,0,'real')
        call color(4)
        call dashset(2)
        icl=0
        zclv(1)=10
        call contourl(imag(Forcecomplex), &
             cworka,nk,nk,no,zclv,icl,kik,or,icsw)
        call legendline(.6,-.1,0,'imag')
        call dashset(0)
     endif
  enddo
  write(*,*)'Finished'
  call pltend()
  
end program fcontkor
