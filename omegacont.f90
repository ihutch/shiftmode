! Contour the real and imaginary force(s) over a complex omega domain.
! For k, psi and Omegac having some values.
! Solve for the complex omega that makes the complex force zero.

program fomegasolve
  real :: kp,kmid
!  logical :: lcont=.false.,lplot=.false.
  logical :: lcont=.true.,lplot=.true.
  integer, parameter :: nk=1,noc=1
  complex :: osolved(nk,noc),omegap
  real :: Omegacarr(noc),karr(nk)
  real :: ormax,oimax
  character*30 string
  psip=.64
!  psip=.04
!  psip=.16
  psip=.09
!  psip=.25
!  psip=.04
!  psip=.36
  Omegacmax=sqrt(psip)/2.      *1.125   !.875 ! .625 .685
  oimax=sqrt(psip)/16.         *3.
  ormax=2.*sqrt(psip)/16.      *.4
  ormax=Omegacmax
  kmid=sqrt(psip)/8.           *1.2
  
  Omegacmax=10                ! Values adjusted for new sech4 analytics.
  if(Omegacmax.ge.5)then               ! High-B scaling.
     oimax=psip**1.5*2.1e-3            ! High-B scaling
     ormax=psip**0.75*.1               ! High-B scaling
     kmid=psip**0.25*.16               ! High-B scaling
  endif
  rangek=.5    !fractional k-range
  Typ=1.

  do ik=1,nk
     if(nk.gt.1)then
        kp=kmid*(1.+rangek*(2.*(ik-1.)/max(1.,nk-1.)-1.))
     else
        kp=kmid
     endif
     karr(ik)=kp/sqrt(psip)
     write(*,*)'kp=',kp,' karr=',karr(ik)
     do ioc=1,noc
        Omegacp=ioc*Omegacmax/noc
        Omegacarr(ioc)=Omegacp/sqrt(psip)
        call fomegacont(psip,Omegacp,Typ,kp,lcont,lplot,err,omegap,ormax,oimax)
        if(err.lt.1.e-4)then
           osolved(ik,ioc)=omegap/sqrt(psip)
        else
           osolved(ik,ioc)=complex(1.e-10,1.e-10)
        endif
     enddo
     write(*,'(3f8.4)')(Omegacarr(ioc),osolved(ik,ioc),ioc=1,noc)
     if(ik.eq.1)then
        call pfset(3)
        call pltinit(0.,Omegacarr(noc),0.,.2)
        if(Omegacarr(noc).gt.1.) &
        call scalewn(0.,Omegacarr(noc),1.e-5,.4,.false.,.true.) ! logarithmic
        call charsize(0.02,0.02)
        call axis
        call axis2
        call axlabels('!AW!@/!A)y!@','!Aw!@/!A)y!@')
        call legendline(.2,.9,0,'!Aw!@!di!d')
        call dashset(2)
        call legendline(.2,.85,0,'!Aw!@!dr!d')
        call dashset(0)
        call accisflush
        call winset(.true.)
     endif
     call color(ik)
     call polyline(Omegacarr,imag(osolved(ik,:)),noc)
     call dashset(2)
     call polyline(Omegacarr,real(osolved(ik,:)),noc)
     call dashset(0)
!     string='k/!A)y!@='           ! really k/sqrt(psi)
!     call fwrite(kp/sqrt(psip),iwidth,3,string(10:))
     call fwrite(kp/sqrt(psip),iwidth,2,string)
!     call jdrwstr(wx2nx(Omegacarr(1)),wy2ny(imag(osolved(ik,1)))+.01,string,1.1)
     ilabel=nint(noc*.75)
     call jdrwstr(wx2nx(Omegacarr(ilabel)),wy2ny(real(osolved(ik,ilabel)))+.01,string,0.)
     call accisflush()
  enddo
  call color(15)
  call jdrwstr(wx2nx(Omegacarr(ilabel)),wy2ny(real(osolved(ik-1,ilabel)))+.04,'k/!a)y!@=',0.)
  call pltend
  open(13,file='fstore.dat',status='unknown')
  do ik=1,nk ! Give legends for value of k/sqrt(psi)
     write(13,'(''legend: '',f5.3,''r'')')karr(ik)
     write(13,'(''legend: '',f5.3,''i'')')karr(ik)
  enddo
  write(13,'(a)')':-x!AW!@/!A)y!@',':-y!Aw!@/!A)y!@'
  write(13,*)noc,nk*2
  do ioc=1,noc
     write(13,*)Omegacarr(ioc)  &
          ,(real(osolved(ik,ioc)),imag(osolved(ik,ioc)),ik=1,nk)
  enddo
  close(13)
end program fomegasolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fomegacont(psip,Omegacp,Typ,kp,lcont,lplot,err,omegap,ormax,oimax)
  use shiftmode
  complex :: omegap
  integer, parameter ::   nor=41,noi=41,nunconv=3
  real :: or(nor),oi(noi),FE,kp
  complex ::  omegacomplex(nor,noi),forcecomplex(nor,noi)
  real, dimension(nor,noi) :: cworka
  integer :: icl
  real :: zclv(20)
  logical :: lcont,lplot,lplot2=.true.,lreadit=.false.
  character*30 string,filename,argument
  real ormax,oimax

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
  read(12,err=101)omegacomplex,forcecomplex
  close(12)
  lreadit=.true.
  write(*,*)'Read forcecomplex from file'
  goto 102
101 write(*,*)'Failed to open or read from file: ', filename
102 continue
  
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
           dioi=-0.04   ! Offset of oi(1) from zero.
           do ioi=1,noi
              oi(ioi)=(ioi-1+dioi)*oimax/(noi-1+dioi)
              omegacomplex(ior,ioi)=complex(or(ior),oi(ioi))
              omega=omegacomplex(ior,ioi)
              call SumHarmonics()
              forcecomplex(ior,ioi)=Fpasstotal+Ftraptotal-FE
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
        write(12)omegacomplex,forcecomplex
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
  endif

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
     Fsum=Fpasstotal+Ftraptotal-FE
     write(*,'(a,i4,5f10.6)')'i,omega,Fsum,err=',i+1,omega,Fsum,err
     if(err.lt..5e-4)exit
  enddo
  if(lplot.and.nzo.lt.2)call polymark(real(omega),imag(omega),1,'@')
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
