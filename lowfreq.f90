! Plot the imaginary part of Fp and Ft at low omegar, and near zero omegai
! as function of omegar and psi. 

program lowfreq
  use shiftmode
  integer, parameter ::   nor=41,noi=1,nunconv=3,npsi=10
  real :: or(nor),oi(noi),FE
  complex, dimension(nor,noi) ::  omegacomplex,forcecomplex
  complex, dimension(nor,noi) ::  Ftcomplex,Fpcomplex
  complex, dimension(npsi) :: Ftcpsi,Fpcpsi
  real, dimension(npsi) :: psiF
  logical :: lplot2=.true.,lreadit=.false.
  character*30 filename,argument
  real :: kp
!  logical :: lcont=.false.,lplot=.false.
  logical :: lcont=.true.,lplot=.true.
  real :: ormax,oimax
  complex :: omegap
  character*80 string

  idebug=-3
  Omegacmax=20.

  k=0.
  psip=.64
!  psip=.04
!  psip=.16
!  psip=.09
!  psip=.25
!  psip=.04
!  psip=.36
  Omegacp=Omegacmax

  if(Omegacmax.ge.5)then               ! High-B scaling.
     oimax=psip**1.5*8.e-4             ! High-B scaling
     ormax=psip**0.75*.07              ! High-B scaling
  endif
  Typ=1.

  if(noi.eq.1)then
     lplot=.false.
     lplot2=.false.
  endif
!!!!!!!!!!!!!!!!!!!!!!!
  ! Create filename in accordance with passed parameters
  write(filename,'(a,2i2.2,a,i3.3,a,i3.3,a,i3.3,a)')   &
       'F',nor,noi, &
       'Oc',min(999,abs(int(100*Omegacp))),   &
       'k',min(999,abs(int(1000*kp))),  &
       'p',min(999,abs(int(100*psip))),'.arr'
!!!!!!!!!!!!!!!!!!!!!!!
 do i=1,iargc()   ! Check cmdline for filenames.
     call getarg(i,argument)
     if(.not.argument(1:1).eq.'-')then
        filename=argument
        write(*,'(a,$)')'Specified file: '
     endif
  enddo
!!!!!!!!!!!!!!!!!!!!!!!
  ! Try to open the file.
  open(12,file=filename,status='old',form='unformatted',err=101) 
  read(12,err=101)norf,noif
  if(norf.ne.nor.or.noif.ne.noi)then
     write(*,*)'Reading from ',filename
     write(*,*)'File array dimensions',norf,nori,' not compatible'
     write(*,*)'Adjust allocation or delete file'
     stop
  endif
  read(12,err=101)or,oi
  read(12,err=101)psip,Omegacp,Typ,kp,omegap,ormax,oimax
  read(12,err=101)omegacomplex,forcecomplex,Ftcomplex,Fpcomplex
  close(12)
  lreadit=.true.
  write(*,*)'Read forcecomplex from file: ',filename
  goto 102
101 write(*,*)'Failed to open or read from file: ', filename
102 continue
!!!!!!!!!!!!!!!!!!!!!!!  
! Set shiftmode values from parameters passed or read.
  Omegac=Omegacp
  psi=psip
  sqpsi=sqrt(psi)
  tqpsi=psi**0.75
  uqpsi=psi**1.5
  qqpsi=psi**.25
  k=kp
  Ty=Typ
  FE=(psi*k)**2*128./315.
  call initialize

  if(.not.lreadit)then    ! Failed to read from file so calculate
     if(lcont)then
        write(*,*)nor,noi,Omegac,k,psi

     ! Contruct the forcecomplex matrix
        write(*,'(a)')'ior,  ioi  omegar omegai        Fpass           Ftrap       FT'
        do ior=1,nor
!           or(ior)=(ior-1)*ormax/(nor-1.)
!           if(noi.eq.1.and.ior.eq.1)or(ior)=0.5*ormax/(nor-1.) 
           ! Uniform log plot alternative.
           or(ior)=ormax*10.**(-2.5*(1.-(ior-1.)/(nor-1.))) *1.3
           dioi=0.00   ! Offset of oi(1) from zero.
           do ioi=1,noi
              oi(ioi)=(ioi-1+dioi)*oimax/(noi-1+dioi)
              if(noi.eq.1)oi(ioi)=1.*dioi*oimax
              omegacomplex(ior,ioi)=complex(or(ior),oi(ioi))
              omega=omegacomplex(ior,ioi)
              omegad=omega
              call SumHarmonics()
              forcecomplex(ior,ioi)=Fpasstotal+Ftraptotal-FE
              Ftcomplex(ior,ioi)=Ftraptotal
              Fpcomplex(ior,ioi)=Fpasstotal
              write(*,'(2i4,8f8.4)')ior,ioi,omega/uqpsi,Fpasstotal,Ftraptotal &
                   ,forcecomplex(ior,ioi)
           enddo
        enddo
        
        ! Attempt to write but skip if file exists.
        open(12,file=filename,status='new',form='unformatted',err=103)
        write(*,*)'Opened new file: ',filename,' and writing'
        write(12)nor,noi
        write(12)or,oi
        write(12)psip,Omegacp,Typ,kp,omegap,ormax,oimax
        write(12)omegacomplex,forcecomplex,Ftcomplex,Fpcomplex
        goto 104
103     write(*,*)'New File: ',filename,' cannot be opened; not rewriting.'
104     close(12)
     endif
  endif

  write(*,*)'Omegac/omegab=',2*Omegac/sqrt(psi)

  write(*,'(a,f8.4,a,f8.4,a,f8.5,a,f8.4)')'psi=',psi,  &
       ' Omegac=',Omegac,' k=',k,' Omegac/omegab=',2*Omegac/sqrt(psi)
  omega=complex(0.7*sqrt(psi)/8.,0.7*sqrt(psi)/16.*1.2)

  err=1.
  omegap=omega
  call pfset(3)
  !     call lautoplot(or/tqpsi,imag(Ftcomplex(1:,1)),nor,.true.,.true.)
  call pltinit(0.,1.,0.,1)
  ytop=imag(Ftcomplex(nor,1))
  ybot=imag(-Fpcomplex(1,1))
  xbot=or(1)/tqpsi
  xtop=or(nor)/tqpsi
  
  write(string,'(''Imaginary Forces !Ay!@='',f7.4,'' !Aw!@!di!d='',e9.2)') &
       psi,oi(noi)
  write(*,*)'noi,oi=',ioi,oi
  call fitscale(xbot,xtop,ybot,ytop,.true.,.true.)
  call polyline(or/tqpsi,imag(Ftcomplex(1:,1)),nor)
  call boxtitle(string(1:lentrim(string)))
  call axis
  call axis2
  call axlabels('!Aw!B!dr!d/!Ay!@!u3/4!u','F')
  call winset(.true.)
  call legendline(.1,.9,0,'   F!dt!d')
  call dashset(4)
  call polyline(or/tqpsi,imag(Ftcomplex(nor,1))*(or/or(nor))**4.5,nor)
  call legendline(.1,.8,0,' !Aw!@!dr!d, !Aw!@!dr!d!u3!u, !Aw!@!dr!d!u4.5!u') 
  call polyline(or/tqpsi,imag(Ftcomplex(1,1))*(or/or(1))**1,nor)
  call color(1)
  call dashset(2)
  call legendline(.1,.85,0,'  -F!dp!d')
  call polyline(or/tqpsi,imag(-Fpcomplex(1:,1)),nor)
  call dashset(4)
  call polyline(or/tqpsi,imag(-Fpcomplex(nor,1))*(or/or(nor))**3,nor)
  call dashset(0)
  call pltend

! Calculate and plot the force scaling with psi.
  write(*,*)'ipsi   psi      omega       Ftcpsi/psi^3.2   Fpcpsi/psi^3.2'
  do ipsi=1,npsi
     psi=psip*10.**(-2.*(1.-float(ipsi)/npsi))
     psiF(ipsi)=psi
     oimax=psi**1.5*8.e-4             ! High-B scaling
     ormax=psi**0.75*.07              ! High-B scaling
     sqpsi=sqrt(psi)
     tqpsi=psi**0.75
     uqpsi=psi**1.5
     qqpsi=psi**.25
     FE=(psi*k)**2*128./315.

     call initialize
     omega=complex(ormax,oimax)
     omegad=omega   ! Probably not necessary.
     call SumHarmonics()
     Ftcpsi(ipsi)=Ftraptotal
     Fpcpsi(ipsi)=Fpasstotal
     write(*,'(i3,7f8.4)')ipsi,psi,omega,Ftcpsi(ipsi)/psi**3.2,Fpcpsi(ipsi)/psi**3.2
  enddo

  call lautoplot(psiF,imag(Ftcpsi),npsi,.true.,.true.)
  call boxtitle('Imaginary Force at constant !Aw!@/!Ay!@!u0.75!u=0.07')
  call axis2
  call axlabels('!Ay!@','F')
  call legendline(.1,.9,0,'   F!dt!d')
  call dashset(2)
  call legendline(.1,.85,0,'  -F!dp!d')
  call polyline(psiF,imag(-Fpcpsi),npsi)
  call dashset(4)
  call polyline(psiF,imag(Ftcpsi(1))*(psiF/psiF(1))**3.2,npsi)
  call legendline(.1,.8,0,' !Ay!@!u3.2!u')
  call pltend


end program lowfreq
  
