! Plot the imaginary part of Fp and Ft at low omegar, and near zero omegai
! as function of omegar and psi. 

program lowfreq
  use shiftmode
  integer, parameter ::   nor=41,noi=1,nunconv=3
  real :: or(nor),oi(noi),FE
  complex, dimension(nor,noi) ::  omegacomplex,forcecomplex
  complex, dimension(nor,noi) ::  Ftcomplex,Fpcomplex
  real, dimension(nor,noi) :: cworka
  integer :: icl
  real :: zclv(20)
  logical :: lplot2=.true.,lreadit=.false.
  character*30 string,filename,argument
  real :: kp
!  logical :: lcont=.false.,lplot=.false.
  logical :: lcont=.true.,lplot=.true.
  real :: ormax,oimax
  complex :: omegap
!  character*30 string
  psip=.64
!  psip=.04
!  psip=.16
!  psip=.09
!  psip=.25
!  psip=.04
!  psip=.36
  Omegacmax=10.
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
  ! Create filename in accordance with passed parameters
  write(filename,'(a,2i2.2,a,i3.3,a,i3.3,a,i3.3,a)')   &
       'F',nor,noi, &
       'Oc',min(999,abs(int(100*Omegacp))),   &
       'k',min(999,abs(int(1000*kp))),  &
       'p',min(999,abs(int(100*psip))),'.arr'

  do i=1,iargc()   ! Check cmdline for filenames.
     call getarg(i,argument)
     if(.not.argument(1:1).eq.'-')then
        filename=argument
        write(*,'(a,$)')'Specified file: '
     endif
  enddo
!  write(*,*)filename

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
  write(*,*)'Read forcecomplex from file'
  goto 102
101 write(*,*)'Failed to open or read from file: ', filename
102 continue

  sqpsi=sqrt(psi)
  tqpsi=psi**0.75
  uqpsi=psi**1.5
  qqpsi=psi**.25
  
! Set shiftmode values from parameters passed or read.
  Omegac=Omegacp
  psi=psip
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
           or(ior)=(ior-1)*ormax/(nor-1.)
           if(noi.eq.1.and.ior.eq.1)or(ior)=0.5*ormax/(nor-1.)
           ! Uniform log plot alternative.
           or(ior)=ormax*10.**(-2.5*(1.-(ior-1.)/(nor-1.)))
           dioi=0.04   ! Offset of oi(1) from zero.
           do ioi=1,noi
              oi(ioi)=(ioi-1+dioi)*oimax/(noi-1+dioi)
              if(noi.eq.1)oi(ioi)=1.*dioi*oimax
              omegacomplex(ior,ioi)=complex(or(ior),oi(ioi))
              omega=omegacomplex(ior,ioi)
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
  if(lplot)then
     call pfset(3)
     call pltinit(0.,ormax,0.,oimax)
     call charsize(0.02,0.02)
     call axis
     call axis2
     call axlabels('!Aw!B!dr!d','!Aw!B!di!d!@')
     if(lcont)then
        icsw=1
        call color(1)
        call dashset(4)
        icl=0
        zclv(1)=20
        call contourl(real(forcecomplex),cworka,nor,nor,noi,zclv,icl,or,oi,icsw)
        call legendline(0.1,-.1,0,'real')
        icl=-1
        zclv(1)=0.
        call color(1)
        call dashset(0)
        call contourl(real(forcecomplex),cworka,nor,nor,noi,zclv,icl,or,oi,icsw)
        icl=0
        zclv(1)=20
        call color(2)
        call dashset(2)
        call contourl(imag(forcecomplex),cworka,nor,nor,noi,zclv,icl,or,oi,icsw)
        call legendline(.6,-.1,0,'imag')
        call dashset(0)
        icl=-1
        zclv(1)=0.0000
        call color(2)
        call dashset(0)
        call contourl(imag(forcecomplex),cworka,nor,nor,noi,zclv,icl,or,oi,icsw)
        call color(15)
        call fwrite(k,iwidth,3,string)
        call legendline(0.1,1.04,258,'k='//string)
        call fwrite(omegac,iwidth,3,string)
        call legendline(.45,1.04,258,'!AW!@='//string)
        call fwrite(psi,iwidth,2,string)
        call legendline(.8,1.04,258,'!Ay!@='//string)
        call dashset(0)
     endif
     call pltend
  endif

  zoif=.01  ! Iteration minimum oimag limit factor.
  nzo=0
  write(*,'(a,f8.4,a,f8.4,a,f8.5,a,f8.4)')'psi=',psi,  &
       ' Omegac=',Omegac,' k=',k,' Omegac/omegab=',2*Omegac/sqrt(psi)
  omega=complex(0.7*sqrt(psi)/8.,0.7*sqrt(psi)/16.*1.2)

  err=1.
  omegap=omega
  if(lplot2)then
     call pfset(3)
     call pltinit(0.,ormax/tqpsi,0.,oimax/uqpsi)
     call charsize(0.02,0.02)
     call axis
     call axis2
     call axlabels('!Aw!B!dr!d/!Ay!@!u3/4!u','!Aw!B!di!d!@/!Ay!@!u3/2!u')
     icsw=1
     icl=0
     zclv(1)=10
     call color(1)
     call dashset(2)
     call contourl(real(forcecomplex)/psi**2.5,cworka,nor,nor,noi,zclv,icl, &
          or/tqpsi,oi/uqpsi,icsw)
     call legendline(0.1,-.1,0,'real')
     icl=0
     zclv(1)=20
     call color(2)
     call dashset(0)
     call contourl(imag(forcecomplex)/psi**2.5,cworka,nor,nor,noi,zclv,icl, &
          or/tqpsi,oi/uqpsi,icsw)
     call legendline(.65,-.1,0,'imag')
     call dashset(0)
     call color(15)
     call fwrite(k/qqpsi,iwidth,3,string)
     call legendline(0.05,1.04,258,'k/!Ay!@!u1/4!u='//string)
     call fwrite(omegac,iwidth,2,string)
     call legendline(.45,1.04,258,'!AW!@='//string)
     call fwrite(psi,iwidth,2,string)
     call legendline(.8,1.04,258,'!Ay!@='//string)
     call color(6)
     call polymark(real(omega/tqpsi),imag(omega/uqpsi),1,3)
     call color(15)
     call pltend
  else  ! 1-d plot of Ft vs omegar
     call pfset(3)
!     call lautoplot(or/tqpsi,imag(Ftcomplex(1:,1)),nor,.true.,.true.)
     call pltinit(0.,1.,0.,1)
     ytop=imag(Ftcomplex(nor,1))
     ybot=imag(-Fpcomplex(1,1))
     xbot=or(1)/tqpsi
     xtop=or(nor)/tqpsi
     write(*,*)ybot,ytop
     call fitscale(xbot,xtop,ybot,ytop,.true.,.true.)
     call polyline(or/tqpsi,imag(Ftcomplex(1:,1)),nor)
     call axis
     call axis2
     call axlabels('!Aw!B!dr!d/!Ay!@!u3/4!u','F')
!     call axlabels('!Aw!B!dr!d/!Ay!@!u3/4!u','!Aw!B!di!d!@/!Ay!@!u3/2!u')
     call winset(.true.)
     call legendline(.1,.9,0,'   F!dt!d')
     call dashset(4)
     call polyline(or/tqpsi,imag(Ftcomplex(nor,1))*(or/or(nor))**4,nor)
     call polyline(or/tqpsi,imag(Ftcomplex(1,1))*(or/or(1))**1,nor)
     call color(1)
     call dashset(2)
     call legendline(.1,.8,0,'  -F!dp!d')
     call polyline(or/tqpsi,imag(-Fpcomplex(1:,1)),nor)
     call dashset(4)
     call polyline(or/tqpsi,imag(-Fpcomplex(nor,1))*(or/or(nor))**3,nor)
!     call polyline(or,imag(-Fpcomplex(1,1))*(or/or(1))**1,nor)
     call pltend
  endif

end program lowfreq
  
