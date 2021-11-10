! New version of omegacont.f90 to use for slow hole stability including
! the effects of ion force Fi.
! This version uses only shiftgen.
! Contour the real and imaginary force(s) over a complex omega domain.
! For k, psi and Omegac having some values.
! Solve for the complex omega that makes the complex force zero.


program fomegasolve
  real :: kp,kmid
  logical :: lcont=.true.,lplot=.true.,lerase=.false.
  integer, parameter :: nk=1,noc=1
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
  ormax=0.
  oimax=0.
  call parsefoarguments(psip,vsin,ormax,oimax,lerase,lcont,lplot)
  if(ormax.eq.0)then
     ormax=.3*sqrt(psip)
  endif
  if(oimax.eq.0)then
     oimax=.1*sqrt(psip)
     if(vsin.le.1.2)oimax=.2*sqrt(psip)
     if(vsin.le..9)oimax=.3*sqrt(psip)
  endif
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
             &,ormax,oimax,lerase)
     enddo
  enddo
!  call plotionforce(psip,Typ,vsin,Omegacp)
end program fomegasolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine parsefoarguments(psip,vsin,ormax,oimax,lerase,lcont,lplot)
  character*20 argument
  logical :: lerase,lcont,lplot
  ipfset=3 ! default
  do i=1,iargc()
     call getarg(i,argument)
     if(argument(1:2).eq.'-p')read(argument(3:),*)psip
     if(argument(1:3).eq.'-vs')read(argument(4:),*)vsin
     if(argument(1:3).eq.'-or')read(argument(4:),*)ormax
     if(argument(1:3).eq.'-oi')read(argument(4:),*)oimax
     if(argument(1:3).eq.'-lc')lcont=.not.lcont
     if(argument(1:3).eq.'-lp')lplot=.not.lplot
     if(argument(1:2).eq.'-e')lerase=.not.lerase
     if(argument(1:2).eq.'-c')ipfset=-3
     if(argument(1:2).eq.'-h')goto 1
  enddo
  call pfset(ipfset)
  return
1 write(*,*)'-p psi, -vs vshift, -or -oi real, imag omega,',&
       ' -c no-stopping, -e erase file'
  end subroutine parsefoarguments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fomegacont(psip,Omegacp,Typ,kp,vsin,lcont,lplot,err,ormax,oimax,lerase)
  use shiftgen
  logical :: lerase
  complex :: omegap      ! omegag surrogate maybe read from file
  integer, parameter ::   nor=21,noi=21,nunconv=3
  real :: or(nor),oi(noi),FE,kp
  complex ::  omegacomplex(nor,noi),forcecomplex(nor,noi),Fi
  complex ::  Ftcomplex(nor,noi),Ficomplex(nor,noi)
  real, dimension(nor,noi) :: cworka
  integer :: icl
  real :: zclv(20)
  logical :: lcont,lplot
  logical :: lplot3=.false.,lplot2=.true.
  logical :: lions=.true.,lreadit=.false.
  character*30 string,filename,argument
  real ormax,oimax

  ! Create filename in accordance with passed parameters
  write(filename,'(a,2i2.2,a,i3.3,a,i3.3,a,i3.3,a)')   &
       'F',nor,noi, &
       'Oc',min(999,abs(nint(100*Omegacp))),   &
       'v',min(999,abs(nint(100*vsin))),  &
       'p',min(999,abs(nint(100*psip))),'.arr'
  do i=1,iargc()   ! Check cmdline for filenames.
     call getarg(i,argument)
     if(.not.argument(1:1).eq.'-')then
        filename=argument
        write(*,'(a,$)')'Specified file: '
     endif
  enddo
  ! Try to open the file.
  open(12,file=filename,status='old',form='unformatted',err=101)
  if(.not.lerase)then
     read(12,err=101)norf,noif
     if(norf.ne.nor.or.noif.ne.noi)then
        write(*,*)'Reading from ',filename
        write(*,*)'File array dimensions',norf,nori,' not compatible'
        write(*,*)'Adjust allocation or delete file'
        stop
     endif
     read(12,err=101)or,oi
     read(12,err=101)psip,Omegacp,Typ,kp,omegap,ormax,oimax
     read(12,err=101)omegacomplex,forcecomplex,Ftcomplex,Ficomplex
     read(12,end=100)lions
     close(12)
100  lreadit=.true.
     write(*,*)'Read forcecomplex from file ',filename
  else
     write(*,*)'Overwriting forcecomplex file ',filename
     close(12,status='delete')
  endif
  goto 102
101 write(*,*)'Failed to open or read from file: ', filename
102 continue

  if(lplot)then
     vshift=vsin
     call fvinfplot
  endif
  
  Omegac=Omegacp
  psig=psip
  kg=kp
  Tperpg=Typ
  FE=(psig*kg)**2*128./315.
  isigma=-1

  if(.not.lreadit.and.lcont)then    ! Failed to read from file so calculate
     ! Contruct the forcecomplex matrix
        write(*,'(a)')'ior,  ioi  omegar omegai   Ftotalr  Ftotali'
        Fi=-FE          ! Unless ionforce is called just old k effect.
        do ior=1,nor
           or(ior)=(ior-1)*ormax/(nor-1.)
           dioi=-0.04   ! Offset of oi(1) from zero.
           do ioi=1,noi
              oi(ioi)=(ioi-1+dioi)*oimax/(noi-1+dioi)
              omegacomplex(ior,ioi)=complex(or(ior),oi(ioi))
              if(.false.)then
                 omegag=omegacomplex(ior,ioi)
                 omegaonly=omegag  ! Ignoring kg for now. 
                 Omegacg=Omegac
                 psig=-psig; call SumHarmonicsg(isigma);psig=-psig
                 Ftcomplex(ior,ioi)=Ftotalsumg
              else
                 call electronforce(Ftcomplex(ior,ioi),omegacomplex(ior,ioi)&
                      ,Omegacp,psip,isigma)
              endif
              if(lions)then
                 call ionforce(Fi,omegacomplex(ior,ioi),Omegacp,psip,vsin)
              endif
              Ficomplex(ior,ioi)=Fi
              forcecomplex(ior,ioi)=Ftcomplex(ior,ioi)+Fi
              write(*,'(2i4,2f8.4,2f10.6)')ior,ioi,omegag,forcecomplex(ior,ioi)
           enddo 
        enddo
        write(*,'(a)')'ior,  ioi  omegar omegai   Ftotalr  Ftotali '
        write(*,*)'Omegac,k,psig',Omegac,k,psig
        
        ! Attempt to write but skip if file exists.
        open(12,file=filename,status='new',form='unformatted',err=103)
        write(*,*)'Opened new file: ',filename,' and writing'
        write(12)nor,noi
        write(12)or,oi
        write(12)psip,Omegacp,Typ,kp,omegag,ormax,oimax
        write(12)omegacomplex,forcecomplex,Ftcomplex,Ficomplex
        write(12)lions
        goto 104
103     write(*,*)'New File: ',filename,' cannot be opened; not rewriting.'
104     close(12)
  endif

  write(*,*)'Omegac/omegab=',2*Omegac/sqrt(psig)
  if(lplot)then
     call lplot1(or,oi,nor,noi,vsin,omegac,psig,Ftcomplex)
     call legendline(.1,.9,258,'F!de!d')
     call pltend
     call lplot1(or,oi,nor,noi,vsin,omegac,psig,Ficomplex)
     call legendline(.1,.9,258,'F!di!d')
     call pltend
  
     call lplot1(or,oi,nor,noi,vsin,omegac,psig,forcecomplex)
     call legendline(.1,.9,258,'F!dtotal!d')
  endif
     
! Find root and plot it converging (omegag is set to found omegap implicitly)  
  call iterfindroot(psip,vsin,Omegacp,omegap,isigma,lplot,nunconv)
  write(*,*)'Eigenfrequency=',omegap
  if(lplot)     call pltend
  if(lplot3)then
     call multiframe(2,1,3)
     call ocomplot(or,nor,vsin,omegac,psig,(Ftcomplex(:,1)))
     call legendline(.1,.9,258,'F!de!d at imag(!Aw!@)=0')
!  call orealplot(or,nor,vsin,omegac,psi,real(Ficomplex(:,1)))
     call ocomplot(or,nor,vsin,omegac,psig,(Ficomplex(:,1)/psig**2))
     call legendline(.1,.9,258,'F!di!d')
     call multiframe(0,0,0)
     call pltend()
  endif
  if(lplot2)then
     sqpsi=sqrt(psig)
     tqpsi=psig**0.75
     uqpsi=psig**1.5
     qqpsi=psig**.25
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
     call contourl(real(forcecomplex)/psig**2.5,cworka,nor,nor,noi,zclv,icl, &
          or/tqpsi,oi/uqpsi,icsw)
     call legendline(0.1,-.1,0,'real')
     icl=-1
     zclv(1)=0.
     call dashset(0)
     call contourl(real(forcecomplex)/psig**2.5,cworka,nor,nor,noi,zclv,icl, &
          or/tqpsi,oi/uqpsi,icsw)
     icl=0
     zclv(1)=20
     call color(2)
     call dashset(3)
     call contourl(imag(forcecomplex)/psig**2.5,cworka,nor,nor,noi,zclv,icl, &
          or/tqpsi,oi/uqpsi,icsw)
     call legendline(.65,-.1,0,'imag')
     icl=-1
     zclv(1)=0.
     call dashset(0)
     call contourl(imag(forcecomplex)/psig**2.5,cworka,nor,nor,noi,zclv,icl, &
          or/tqpsi,oi/uqpsi,icsw)
     call color(15)
     call fwrite(kg/qqpsi,iwidth,3,string)
     call legendline(0.05,1.04,258,'k/!Ay!@!u1/4!u='//string)
     call fwrite(omegac,iwidth,2,string)
     call legendline(.45,1.04,258,'!AW!@='//string)
     call fwrite(psig,iwidth,2,string)
     call legendline(.8,1.04,258,'!Ay!@='//string)
     call color(6)
     write(*,*)omegag,.001*sqrt(psip)
     if(imag(omegag).gt.0.0011*sqrt(psip))& 
          call polymark(real(omegag/tqpsi),imag(omegag/uqpsi),1,3)
     call color(15)
     call pltend
  endif
end subroutine fomegacont
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ionforce(Fi,omega,Omegacin,psiin,vsin)
  use shiftgen
  complex :: Fi,omega  !,Ftotalg
  real :: psiin,vsin,Omegacin
  real, parameter :: mime=1836
  omegag=omega*sqrt(mime)
  Omegacg=Omegacin*sqrt(mime)
  omegaonly=omegag
  psig=psiin
  isigma=-1
  vshift=vsin
!  call FgRepelEint(Ftotalg,isigma)  ! Some subtle differences remain.
!  Fi=2.*Ftotalg
  call SumHarmonicsg(isigma)
!  if(abs(Ftotalg-Ftotalsumg).gt.1.e-4)write(*,*)Ftotalg,Ftotalsumg
  Fi=Ftotalsumg
! Undo changes
  omegag=omega
  omegaonly=omegag
  Omegacg=Omegacg/sqrt(mime)
end subroutine ionforce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine electronforce(Felec,omegain,Omegac,psiin,isigma)
  use shiftgen
  complex :: omegain,Felec
  omegag=omegain
  omegaonly=omegag  ! Ignoring kg for now. 
  Omegacg=Omegac
  psig=-psiin; call SumHarmonicsg(isigma);psig=-psig
  Felec=Ftotalsumg
end subroutine electronforce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine lplot1(or,oi,nor,noi,vsin,omegac,psi,forcecomplex)
  real :: or(nor),oi(noi)
  complex ::  forcecomplex(nor,noi)
!  complex ::  Fpcomplex(nor,noi),Ftcomplex(nor,noi),Ficomplex(nor,noi)
  real, dimension(nor,noi) :: cworka
  integer :: icl
  real :: zclv(20)
  character*30 string
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
subroutine iterfindroot(psip,vsin,Omegacp,omegap,isigma,lplot,nunconv)
  real :: psip,vsin,Omegacp
  integer :: isigma,nunconv
  logical lplot
  complex :: omegap,  Fec,Fic,Fsum
     zoif=.001  ! Iteration minimum oi limit factor.
     nzo=0
     omegap=complex(0.7*sqrt(psip)/8.,0.7*sqrt(psip)/16.*1.2)
     call electronforce(Fec,omegap,Omegacp,psip,isigma)
     call ionforce(Fic,omegap,Omegacp,psip,vsin)
     Fsum=Fec+Fic
     do i=1,12
        if(lplot)then
           call color(6)
           call polymark(real(omegap),imag(omegap),1,ichar('0')+i)
        endif
        call complexnewton(Fsum,omegap,err,psip,isigma,vsin,Omegacp)
        if(imag(omegap).lt.zoif*sqrt(psip))then
           nzo=nzo+1
           if(nzo.ge.nunconv)then
              write(*,'(a,i2,a,g10.3)')'Uncoverged after',nzo,' omegai less than',zoif*sqrt(psip)
              err=1.
              goto 1
           endif
           omegap=complex(real(omegap),zoif*sqrt(psip))
        endif
        call electronforce(Fec,omegap,Omegacp,psip,isigma)
        call ionforce(Fic,omegap,Omegacp,psip,vsin)
        Fsum=Fec+Fic
        write(*,'(a,i4,5f10.6)')'i,omegap,Fsum,err=',i,omegap,Fsum,err
        if(err.lt..5e-4)goto 1
     enddo
     write(*,*)'Unconverged after',i-1,' iterations'
1    continue
end subroutine iterfindroot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine complexnewton(Fsum,omegap,err,psip,isigma,vsin,Omegacp)
  ! Take a Newton Step in finding the complex omega root of
  ! complex force by inverting the 2x2 matrix from adjacent evaluations.
  ! Return the Fsum,omegap,err new values.
  complex :: Fsum,omegap
  real :: err,psip,vsin,Omegacp
  real :: eps1=.05,J11,J12,J21,J22,domega1,domega2,det
  complex :: om1,om2,om3,f1,f2,f3,domega,Fec,Fic
  det=1.
  f1=Fsum
  om1=omegap
  ! Calculate the Jacobian's coefficients.
  om2=om1+eps1*abs(real(om1))
  call electronforce(Fec,om2,Omegacp,psip,isigma)
  call ionforce(Fic,om2,Omegacp,psip,vsin)
  f2=Fec+Fic
!  write(*,'(6g12.4)')om1,om2,f2
  if(.not.om2-om1.ne.0)stop 'om2-om1=0'
  J11=real(f2-f1)/real(om2-om1)
  J21=imag(f2-f1)/real(om2-om1)
  om3=om1+complex(0.,eps1*imag(om1))
  call electronforce(Fec,om3,Omegacp,psip,isigma)
  call ionforce(Fic,om3,Omegacp,psip,vsin)
  f3=Fec+Fic
!  write(*,'(6g12.4)')om1,om3,f3
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
  domega1=-(J22*real(f1)-J12*imag(f1))/det
  domega2=-(-J21*real(f1)+J11*imag(f1))/det
!  write(*,*)J11*domega1+J12*domega2,J21*domega1+J22*domega2
  domega=complex(domega1,domega2)
  omegap=om1+domega
  err=abs(domega/omegap)
!  write(*,*)'domega1,domega2,omega',domega1,domega2,omega
end subroutine complexnewton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotionforce(psi,Typ,vsin,Omegacin)
! Repurposed for Ftot and Fe.
  real :: psi,Typ,vsin
  integer, parameter :: nfi=100
  complex, dimension(nfi) :: Ftotarray,Fiarray,Fearray,omegaFi
  omegamax=3.
  write(*,*)'psi=',psi,' vsin=',vsin
  do i=1,nfi
     omegaFi(i)=omegamax*(float(i)/nfi)/sqrt(1836.)+complex(0.,.0001)
     call ionforce(Fiarray(i),omegaFi(i),Omegacin,psi,vsin)
     call electronforce(Fearray(i),omegaFi(i),Omegacin,psi,-1)
     Ftotarray(i)=Fearray(i)+Fiarray(i)
  enddo
  call minmax(Ftotarray,2*nfi,fmin,fmax)
  call pltinit(0.,omegaFi(nfi),fmin,fmax)
  call axis; call axis2
  call axlabels('omega','Ftot, Fe')
  call color(1)
  call polyline(real(omegaFi),real(Ftotarray),nfi)
  call color(2)
  call dashset(1)
  call polyline(real(omegaFi),imag(Ftotarray),nfi)
  if(.true.)then
  call color(3)
  call dashset(2)
  call polyline(real(omegaFi),real(Fearray),nfi)
  call color(4)
  call dashset(3)
  call polyline(real(omegaFi),imag(Fearray),nfi)
  endif
  call dashset(0)
  call color(15)
  call pltend
end subroutine plotionforce
