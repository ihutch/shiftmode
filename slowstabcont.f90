! New version of omegacont.f90 to use for slow hole stability including
! the effects of ion force Fi.
! Contour the real and imaginary force(s) over a complex omega domain.
! For k, psi and Omegac having some values.
! Solve for the complex omega that makes the complex force zero.


program fomegasolve
  real :: kp,kmid
!  logical :: lcont=.false.,lplot=.false.
  logical :: lcont=.true.,lplot=.true.
  integer, parameter :: nk=1,noc=1
  complex :: omegap
  real :: Omegacarr(noc),karr(nk)
  real :: ormax,oimax
!  character*30 string

! Default parameters
  kmid=0.
  Omegacmax=10
  rangek=.5    !fractional k-range
  Typ=1.

  psip=.25
  vsin=1.25         ! Maxwellian ion component velocity shift
  ormax=.15
  oimax=.04
  call parsefoarguments(psip,vsin,ormax,oimax)
  if(vsin.le.1.05)oimax=.15
  if(vsin.le..9)oimax=.3
  
  do ik=1,nk
     if(nk.gt.1)then
        kp=kmid*(1.+rangek*(2.*(ik-1.)/max(1.,nk-1.)-1.))
     else
        kp=kmid
     endif
     karr(ik)=kp/sqrt(psip)
!     write(*,*)'kp=',kp,' karr=',karr(ik)
     do ioc=1,noc
        Omegacp=ioc*Omegacmax/noc
        Omegacarr(ioc)=Omegacp/sqrt(psip)
        call fomegacont(psip,Omegacp,Typ,kp,vsin,lcont,lplot,err&
             &,omegap,ormax,oimax)
     enddo
  enddo
end program fomegasolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine parsefoarguments(psip,vsin,ormax,oimax)
  character*20 argument
  do i=1,iargc()
     call getarg(i,argument)
     if(argument(1:2).eq.'-p')read(argument(3:),*)psip
     if(argument(1:3).eq.'-vs')read(argument(4:),*)vsin
     if(argument(1:3).eq.'-or')read(argument(4:),*)ormax
     if(argument(1:3).eq.'-h')goto 1
  enddo
  return
  1 write(*,*)'-p psi, -vs vshift, -or -oi real, imaginary omega'
  end subroutine parsefoarguments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fomegacont(psip,Omegacp,Typ,kp,vsin,lcont,lplot,err,omegap&
     &,ormax,oimax)
  use shiftmode
  use shiftgen
  complex :: omegap
  integer, parameter ::   nor=21,noi=21,nunconv=3
  real :: or(nor),oi(noi),FE,kp
  complex ::  omegacomplex(nor,noi),forcecomplex(nor,noi),Fi
  complex ::  Fpcomplex(nor,noi),Ftcomplex(nor,noi),Ficomplex(nor,noi)
  real, dimension(nor,noi) :: cworka
  integer :: icl
  real :: zclv(20)
  logical :: lcont,lplot,lplot2=.true.,lreadit=.false.,lions=.true.
  character*30 string,filename,argument
  real ormax,oimax

  ! Create filename in accordance with passed parameters
  write(filename,'(a,2i2.2,a,i3.3,a,i3.3,a,i3.3,a)')   &
       'F',nor,noi, &
       'Oc',min(999,abs(int(100*Omegacp))),   &
       'v',min(999,abs(int(100*vsin))),  &
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
  read(12,err=101)omegacomplex,forcecomplex,Fpcomplex,Ftcomplex,Ficomplex
  read(12,end=100)lions
  close(12)
100  lreadit=.true.
  write(*,*)'Read forcecomplex from file'
  goto 102
101 write(*,*)'Failed to open or read from file: ', filename
102 continue

  if(lplot)then
     vshift=vsin
     call fvinfplot
  endif
  
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
        Fi=-FE          ! Unless ionforce is called just old k effect.
        do ior=1,nor
           or(ior)=(ior-1)*ormax/(nor-1.)
           dioi=-0.04   ! Offset of oi(1) from zero.
           do ioi=1,noi
              oi(ioi)=(ioi-1+dioi)*oimax/(noi-1+dioi)
              omegacomplex(ior,ioi)=complex(or(ior),oi(ioi))
              omega=omegacomplex(ior,ioi)
              call SumHarmonics()
! Here we replace -FE with the ion force +Fi
              Fpcomplex(ior,ioi)=Fpasstotal
              Ftcomplex(ior,ioi)=Ftraptotal
              if(lions)then
                 call ionforce(Fi,omega,psip,vsin)
                 write(*,*)'-FE,Fi',-FE,Fi
              endif
              Ficomplex(ior,ioi)=Fi
              forcecomplex(ior,ioi)=Fpasstotal+Ftraptotal+Fi
              write(*,'(2i4,8f8.4)')ior,ioi,omega,Fpasstotal,Ftraptotal &
                   ,forcecomplex(ior,ioi)
           enddo
        enddo

        ! Attempt to write but skip if file exists.
        open(12,file=filename,status='new',form='unformatted',err=103)
        write(*,*)'Opened new file: ',filename,' and writing'
        write(12)nor,noi
        write(12)or,oi
        write(12)psip,Omegacp,Typ,kp,omegap,ormax,oimax
        write(12)omegacomplex,forcecomplex,Fpcomplex,Ftcomplex,Ficomplex
        write(12)lions
        goto 104
103     write(*,*)'New File: ',filename,' cannot be opened; not rewriting.'
104     close(12)
     endif
  endif

  write(*,*)'Omegac/omegab=',2*Omegac/sqrt(psi)
  if(lplot)call lplot1(or,oi,nor,noi,vsin,omegac,psi,Ftcomplex)
  call legendline(.1,.9,258,'Ftcomplex')
  call pltend
  if(lplot)call lplot1(or,oi,nor,noi,vsin,omegac,psi,Fpcomplex)
  call legendline(.1,.9,258,'Fpcomplex')
  call pltend
  if(lplot)call lplot1(or,oi,nor,noi,vsin,omegac,psi,Fpcomplex+Ftcomplex)
  call legendline(.1,.9,258,'Fp+Ft')
  call pltend
  if(lplot)call lplot1(or,oi,nor,noi,vsin,omegac,psi,Ficomplex)
  call legendline(.1,.9,258,'Ficomplex')
  call pltend
  if(lplot)call lplot1(or,oi,nor,noi,vsin,omegac,psi,forcecomplex)
  call legendline(.1,.9,258,'forcecomplex')
  call pltend

  call multiframe(2,1,3)
!  call orealplot(or,nor,vsin,omegac,psi,real(Fpcomplex(:,noi)+Ftcomplex(:,1)))
  call ocomplot(or,nor,vsin,omegac,psi,(Fpcomplex(:,noi)+Ftcomplex(:,1)))
  call legendline(.1,.9,258,'Fp+Ft at imag(!Aw!@)=0')
!  call orealplot(or,nor,vsin,omegac,psi,real(Ficomplex(:,1)))
  call ocomplot(or,nor,vsin,omegac,psi,(Ficomplex(:,1)))
  call legendline(.1,.9,258,'Fi')
  call multiframe(0,0,0)

  if(.not.lions)then     ! Find solution by iteration.
  zoif=.01  ! Iteration minimum oimag limit factor.
  nzo=0
  write(*,'(a,f8.4,a,f8.4,a,f8.5,a,f8.4)')'psi=',psi,  &
       ' Omegac=',Omegac,' k=',k,' Omegac/omegab=',2*Omegac/sqrt(psi)
  omega=complex(0.7*sqrt(psi)/8.,0.7*sqrt(psi)/16.*1.2)
  call SumHarmonics
  err=1.
  do i=1,10
     if(lplot)then
        call color(6)
! Omit intermediate step indications.
!        call polymark(real(omega),imag(omega),1,ichar('0')+i)
     endif
!     write(*,*)'omega,Fsum',omega,Fpasstotal+Ftraptotal-FE
     call complexnewton(FE,err)
     if(imag(omega).lt.zoif*oimax)then
        nzo=nzo+1
        omega=complex(real(omega),zoif*oimax)
        if(nzo.ge.nunconv)then
           write(*,'(a,i2,a)')'Uncoverged after',nzo,' negative omegai'
           err=1.
           exit
        endif
     endif
     call SumHarmonics()
     Fsum=real(Fpasstotal+Ftraptotal-FE)
     write(*,'(a,i4,5f10.6)')'i,omega,Fsum,err=',i+1,omega,Fsum,err
     if(err.lt..5e-4)exit
  enddo
  if(lplot.and.nzo.lt.2)call polymark(real(omega),imag(omega),1,'@')
  endif
  if(lplot)call pltend()
  omegap=omega
  sqpsi=sqrt(psi)
  tqpsi=psi**0.75
  uqpsi=psi**1.5
  qqpsi=psi**.25
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
  endif
  
end subroutine fomegacont
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ionforce(Fi,omega,psiin,vsin)
  use shiftgen
  complex :: Fi,omega,Ftotalg
  real :: psiin,vsin
  real, parameter :: mime=1836
  omegag=omega*sqrt(mime)
  psig=psiin
  isigma=-1
  vshift=vsin
  call FgRepelEint(Ftotalg,isigma)
  Fi=2.*Ftotalg
end subroutine ionforce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine lplot1(or,oi,nor,noi,vsin,omegac,psi,forcecomplex)
  real :: or(nor),oi(noi)
  complex ::  forcecomplex(nor,noi)
!  complex ::  Fpcomplex(nor,noi),Ftcomplex(nor,noi),Ficomplex(nor,noi)
  real, dimension(nor,noi) :: cworka
  integer :: icl
  real :: zclv(20)
  character*30 string
  call pfset(3)
  call pltinit(0.,or(nor),0.,oi(noi))
  call charsize(0.02,0.02)
  call axis
  call axis2
  call axlabels('!Aw!B!dr!d','!Aw!B!di!d!@')
  
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
!        call fwrite(k,iwidth,3,string)
!        call legendline(0.1,1.04,258,'k='//string)
  call fwrite(vsin,iwidth,3,string)
  call legendline(0.1,1.04,258,'v!ds!d='//string)
  call fwrite(omegac,iwidth,3,string)
  call legendline(.45,1.04,258,'!AW!@='//string)
  call fwrite(psi,iwidth,2,string)
  call legendline(.8,1.04,258,'!Ay!@='//string)
  call dashset(0)
end subroutine lplot1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine orealplot(or,nor,vsin,omegac,psi,orealforce)
  real ::  or(nor)
  real ::  orealforce(nor)

  call minmax(orealforce,nor,fmin,fmax)
  call pltinit(0.,or(nor),fmin,fmax)
  call axis; call axis2
  call axlabels('real(!Aw!@)','real(Force)')
  call polyline(or,orealforce,nor)
end subroutine orealplot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ocomplot(or,nor,vsin,omegac,psi,ocomforce)
  real ::  or(nor)
  complex ::  ocomforce(nor)
  call minmax(ocomforce,2*nor,fmin,fmax)
  call pltinit(0.,or(nor),fmin,fmax)
  call axis; call axis2
  call axlabels('real(!Aw!@)','Force')
  call polyline(or,real(ocomforce),nor)
  call legendline(.1,.8,0,'real')
  call dashset(1)
  call polyline(or,imag(ocomforce),nor)
  call legendline(.1,.7,0,'imag')
  call dashset(0)
end subroutine ocomplot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine complexnewton(FE,err)
  ! Take a Newton Step in finding the complex omega root of
  ! complex force by inverting the 2x2 matrix from adjacent evaluations.
  use shiftmode
  real :: FE,err
  real :: eps1=.05,J11,J12,J21,J22,domega1,domega2,det
  complex :: om1,om2,om3,f1,f2,f3,domega
  det=1.
  f1=Fpasstotal+Ftraptotal-FE
  om1=omega
!  if(err.gt.0.1)then ! Not a big win
  ! Calculate the Jacobian's coefficients.
  om2=om1+eps1*abs(real(om1))
  omega=om2
  call SumHarmonics()
  f2=Fpasstotal+Ftraptotal-FE
  !  write(*,*)om1,om2,f2
  if(.not.om2-om1.ne.0)stop 'om2-om1=0'
  J11=real(f2-f1)/real(om2-om1)
  J21=imag(f2-f1)/real(om2-om1)
  om3=om1+complex(0.,eps1*imag(om1))
  omega=om3
  call SumHarmonics()
  f3=Fpasstotal+Ftraptotal-FE
  !  write(*,*)om1,om3,f3
  if(.not.om3-om1.ne.0)stop 'om3-om1=0'
  J12=real(f3-f1)/imag(om3-om1)
  J22=imag(f3-f1)/imag(om3-om1)
  ! Solve J domega = -f to find domega
  det=J11*J22-J21*J12
  if(det.eq.0)then
     write(*,*)'Det=0'
     write(*,*)om1,om2,om3,f1,f2,f3
     write(*,'(2f10.4)')J11,J12,J21,J22
     ! Try to recover.
     det=1.
  endif
!  endif
  domega1=-(J22*real(f1)-J12*imag(f1))/det
  domega2=-(-J21*real(f1)+J11*imag(f1))/det
!  write(*,*)J11*domega1+J12*domega2,J21*domega1+J22*domega2
  domega=complex(domega1,domega2)
  omega=om1+domega
  err=abs(domega/omega)
!  write(*,*)'domega1,domega2,omega',domega1,domega2,omega
end subroutine complexnewton
