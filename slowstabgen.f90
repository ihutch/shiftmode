! New version of omegacont.f90 to use for slow hole stability including
! the effects of ion force Fi.
! This version uses only shiftgen.
! Contour the real and imaginary force(s) over a complex omega domain.
! For k, psi and Omegac having some values.
! Solve for the complex omega that makes the complex force zero.
! k-variation is not yet actually implemented. 
! The variation of interest here is more about psip and vshift.

program fomegasolve
  real :: kin,kmid
  logical :: lcont=.true.,lplot=.false.,lerase=.false.,lTiscan=.false.
  integer, parameter :: nk=1,noc=1,npsi=8,nvsin=21,nv0=0
  real :: Omegacarr(noc),karr(nk),psiparray(npsi),vsinarray(nv0:nvsin)
  real :: ormax,oimax
  complex :: omegasolve(nv0:nvsin,npsi),omegap
  character*30 string
  
! Default parameters
  kmid=0.
  Omegacmax=10
  rangek=.5    !fractional k-range
  Typ=1.
  

  isigma=-1
  psip=.25
  vsin=1.25         ! Maxwellian ion component velocity shift
  ormax=0.
  oimax=0.
  call parsefoarguments(psip,vsin,ormax,oimax,kmid,Omegacmax,lerase,lcont,lplot,Ti,lTiscan)
  if(ormax.eq.0)then
     ormax=.3*sqrt(psip)
  endif
  if(oimax.eq.0)then
     oimax=.1*sqrt(psip)
     if(vsin.le.1.2)oimax=.3*sqrt(psip)
     if(vsin.le..9)oimax=.4*sqrt(psip)
  endif
  if(lplot)then
! contouring of psimax and vsmax case.       
     do ik=1,nk
        if(nk.gt.1)then
           kin=kmid*(1.+rangek*(2.*(ik-1.)/max(1.,nk-1.)-1.))
        else
           kin=kmid
        endif
        karr(ik)=kin/sqrt(psip)
!     write(*,*)'kin=',kin,' karr=',karr(ik)
        do ioc=1,noc
           Omegacp=ioc*Omegacmax/noc
           Omegacarr(ioc)=Omegacp/sqrt(psip)
           call fomegacont(psip,Omegacp,Typ,omegap,kin,vsin,lcont,lplot,err&
                &,ormax,oimax,lerase)
        enddo
     enddo
  elseif(lTiscan)then
! omega versus Ti at several vsin
     kin=0.
     Timax=Ti
     Timin=.5
     psip=.05
     omegasolve=0.
     nvs=4
     write(*,*)' psip   vshift   Ti    it     omega'
     do ip=1,nvs ! Really this is vs iterations here
        vs=vsin*(ip-1)/(nvs-1.)
        psiparray(ip)=vs
        do iv=1,nvsin  ! Really this is Ti iterations
           Ti=Timin+(Timax-Timin)*(iv-1.)/(nvsin-.99999)
           if(ip.eq.1)vsinarray(iv)=Ti
           omegap=complex(0.7*sqrt(psip)/8.,.7*sqrt(psip)/8./(1.+vsin))
           call Tset(1.,Ti)
           call iterfindroot(psip,vs,Omegacp,omegap,kin,isigma,lplot,ires)
           if(ires.ne.0)omegasolve(iv,ip)=omegap/sqrt(psip)
           write(*,'(3f8.4,i3,$)')psip,vs,Ti,ires
           write(*,*)omegap
        enddo
     enddo
     call pltinit(vsinarray(1),vsinarray(nvsin),0.,.3+psip)
     call charsize(.019,.019)
     call axis; call axis2
     call axlabels('T!di!d','!Aw!@/!Ay!@!u1/2!u')
     call winset(.true.)
     call legendline(.55,.95,258,'v!ds!d=')
     iline=1
     do ip=1,nvs
        call color(ip)
!        call polyline(vsinarray(1),imag(omegasolve(nvsin,ip))&
!             &*(vsinarray(nvsin)/vsinarray(1:nvsin)),nvsin)
        call fwrite(psiparray(ip),iwidth,2,string)
        call dashset(2*ip-2)
        call polyline(vsinarray(1),real(omegasolve(1:,ip)),nvsin)
        if(real(omegasolve(1,ip)).gt.0.001)then
           call legendline(.5,.95-.05*iline,0,' '//string(1:iwidth)&
                &//' !Aw!@!dr!d')
           iline=iline+1
        endif
        call dashset(2*ip-1)
        call polyline(vsinarray(1),imag(omegasolve(1:,ip)),nvsin)
        call legendline(.5,.95-.05*iline,0,' '//string(1:iwidth)//' !Aw!@!di!d')
        iline=iline+1
        call dashset(0)
     enddo
     call pltend
     
  else
! omega versus vsin plots at various psi.
     kin=0.
     psipmax=psip
     vsmax=vsin
     Omegacp=Omegacmax
     lplot=.false.
     if(nvsin.lt.15)then
        vsmin=0.5
     else
        vsmin=0.1
     endif
     omegasolve=0.
     do ip=1,npsi
        psip=ip*psipmax/npsi
        psiparray(ip)=psip
        do iv=nv0,nvsin
           vsin=vsmin+(vsmax-vsmin)*(iv-1.)/(nvsin-.99999)
           if(iv.eq.0)vsin=0.
           if(ip.eq.1)vsinarray(iv)=vsin
           omegap=complex(0.7*sqrt(psip)/8.,.7*sqrt(psip)/8./(1.+vsin))
           call iterfindroot(psip,vsin,Omegacp,omegap,kin,isigma,lplot,ires)
           if(ires.ne.0)omegasolve(iv,ip)=omegap/sqrt(psip)
           write(*,'(2f8.4,i3,$)')psip,vsin,ires
           write(*,*)omegap
        enddo
     enddo
!  call minmax(omegasolve,2*npsi*nvsin,ommin,ommax)
!  call pltinit(vsinarray(nv0),vsinarray(nvsin),0.,ommax*1.1)
! Now scaled version  
! Uncorrected  call pltinit(vsinarray(nv0),vsinarray(nvsin),0.,.34)
     call pltinit(vsinarray(nv0),vsinarray(nvsin),0.,.3+psipmax)
     call charsize(.019,.019)
     call axis; call axis2
     call axlabels('!Bv!ds!d!@','!Aw!@/!Ay!@!u1/2!u')
     if(Ti.ne.1.)then
        call fwrite(Ti,iwidth,1,string)
        call legendline(.05,.07,258,'T!di!d='//string(1:iwidth))
     endif
     call winset(.true.)
     call legendline(.7,.95,258,'!Ay!@=')
     do ip=1,npsi
        call color(15)
        call fwrite(psiparray(ip),iwidth,3,string)
        call dashset(ip-1)
        call legendline(.65,.95-.05*ip,0,' '//string(1:iwidth))
        call color(1)
        call polyline(vsinarray,real(omegasolve(:,ip)),(nvsin+1-nv0))
        if(ip.eq.1)call legendline(.65,.5,258,'real(!Aw!@)/!Ay!@!u1/2!u')
        call color(4)
        call polyline(vsinarray,imag(omegasolve(:,ip)),(nvsin+1-nv0))
        if(ip.eq.1)call legendline(.1,.9,258,'imag(!Aw!@)/!Ay!@!u1/2!u')
        call dashset(0)
     enddo
     call color(0)
     call pltend
  endif
end program fomegasolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine parsefoarguments(psip,vsin,ormax,oimax,kin,Omegacmax,lerase,lcont,lplot,Ti,lTiscan)
  character*20 argument
  real :: kin
  logical :: lerase,lcont,lplot,lTiscan
  ipfset=3 ! default
  Ti=1.
  do i=1,iargc()
     call getarg(i,argument)
     if(argument(1:2).eq.'-p')read(argument(3:),*)psip
     if(argument(1:2).eq.'-k')read(argument(3:),*)kin
     if(argument(1:3).eq.'-vs')read(argument(4:),*)vsin
     if(argument(1:3).eq.'-or')read(argument(4:),*)ormax
     if(argument(1:3).eq.'-oi')read(argument(4:),*)oimax
     if(argument(1:3).eq.'-oc')read(argument(4:),*)Omegacmax
     if(argument(1:3).eq.'-Ti')then
        read(argument(4:),*)Ti
        write(*,'(a,f8.2)')'Setting Ti to',Ti
        call Tset(1.,Ti)
     endif
     if(argument(1:3).eq.'-lc')lcont=.not.lcont
     if(argument(1:3).eq.'-lp')lplot=.not.lplot
     if(argument(1:3).eq.'-lT')lTiscan=.not.lTiscan
     if(argument(1:2).eq.'-e')lerase=.not.lerase
     if(argument(1:2).eq.'-c')ipfset=-3
     if(argument(1:2).eq.'-h')goto 1
  enddo
  call pfset(ipfset)
  return
1 write(*,*)'-p psi, -vs vshift, -or -oi real, imag omega,',&
       ' -c no-stopping, -e erase file'
  write(*,*)'-Ti ion tempr, -oc Omegac -lp toggle on contours, -lt toggle Tiscan' 
  stop
end subroutine parsefoarguments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fomegacont(psip,Omegacp,Typ,omegap,kin,vsin,lcont,lplot&
     &,err,ormax,oimax,lerase)
  logical :: lerase
  complex :: omegap      ! omegag surrogate maybe read from file
  integer, parameter ::   nor=21,noi=21
  real :: or(nor),oi(noi),kin
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

write(*,*)'fomegacont k=',kin
  ! Create filename in accordance with passed parameters
  write(filename,'(a,2i2.2,a,i3.3,a,i3.3,a,i3.3,a,i3.3,a)')   &
       'F',nor,noi, &
       'Oc',min(999,abs(nint(100*Omegacp))),   &
       'k',min(999,abs(nint(100*kin))),   &
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
     read(12,err=101)psip,Omegacp,Typ,omegap,kin,ormax,oimax
     read(12,err=101)omegacomplex,forcecomplex,Ftcomplex,Ficomplex
     read(12,end=100)lions
     close(12)
100  lreadit=.true.
     write(*,*)'Read forcecomplex from file ',filename
  write(*,*)'kin=',kin,' omegap=',omegap
  else
     write(*,*)'Overwriting forcecomplex file ',filename
     close(12,status='delete')
  endif
  goto 102
101 write(*,*)'Failed to open or read from file: ', filename
102 continue

  if(lplot)call plotfv(vsin)  
  isigma=-1
  FE=kin**2*psip**2*128./315.

  if(.not.lreadit.and.lcont)then    ! Failed to read from file so calculate
     ! Contruct the forcecomplex matrix
        write(*,'(a)')'ior,  ioi  omegar omegai   Ftotalr  Ftotali     k nharm'
        do ior=1,nor
           or(ior)=(ior-1)*ormax/(nor-1.)
           dioi=-0.04   ! Offset of oi(1) from zero.
           do ioi=1,noi
              oi(ioi)=(ioi-1+dioi)*oimax/(noi-1+dioi)
              omegacomplex(ior,ioi)=complex(or(ior),oi(ioi))
              call electronforce(Ftcomplex(ior,ioi),omegacomplex(ior&
                   &,ioi),kin,Omegacp,psip,vsin,isigma)
              if(lions)then
                 call ionforce(Fi,omegacomplex(ior,ioi) ,kin,Omegacp,&
                      & psip,vsin ,isigma)
              endif
              Ficomplex(ior,ioi)=Fi
              forcecomplex(ior,ioi)=Ftcomplex(ior,ioi)+Fi+FE
              write(*,'(2i4,2f8.4,2f10.6,f7.3,i4)')ior,ioi,omegacomplex(ior&
                   &,ioi),forcecomplex(ior,ioi),kin,inharm()
           enddo 
        enddo
        write(*,'(a)')'ior,  ioi  omegar omegai   Ftotalr  Ftotali     k nharm'
        write(*,*)'Omegacp,k,psip',Omegacp,kin,psip
        
        ! Attempt to write but skip if file exists.
        open(12,file=filename,status='new',form='unformatted',err=103)
        write(*,*)'Opened new file: ',filename,' and writing'
        write(12)nor,noi
        write(12)or,oi
        write(12)psip,Omegacp,Typ,omegap,kin,ormax,oimax
        write(12)omegacomplex,forcecomplex,Ftcomplex,Ficomplex
        write(12)lions
        goto 104
103     write(*,*)'New File: ',filename,' cannot be opened; not rewriting.'
104     close(12)
  endif

  write(*,*)'Omegacp/omegab=',2*Omegacp/sqrt(psip)
  if(lplot)then
     call lplot1(or,oi,nor,noi,vsin,omegacp,kin,psip,Ftcomplex/psip**2)
     call legendline(.1,.9,258,'!p!o~!o!qF!de!d/!Ay!@!u2!u')
     call pltend
     call lplot1(or,oi,nor,noi,vsin,omegacp,kin,psip,Ficomplex/psip**2)
     call legendline(.1,.9,258,'!p!o~!o!qF!di!d/!Ay!@!u2!u')
     call pltend

     call lplot1(or,oi,nor,noi,vsin,omegacp,kin,psip,forcecomplex/psip**2)
     if(FE.eq.0.)then
        call legendline(.1,.9,258,'!p!o~!o!qF/!Ay!@!u2!u')
     else
        call legendline(.1,.9,258,'(!p!o~!o!qF+F!dE!d)/!Ay!@!u2!u')
     endif
  endif
     
! Find root and plot it converging (omegag is set to found omegap implicitly)  
!  omegap=complex(0.7*sqrt(psip)/8.,1.*sqrt(psip)/8./(1.+vsin))
  write(*,*)'calling iterfindroot',omegap,kin
  call iterfindroot(psip,vsin,Omegacp,omegap,kin,isigma,lplot,ires)
  write(*,*)'Eigenfrequency=',omegap
  if(lplot)     call pltend
  if(lplot3)then
     call multiframe(2,1,3)
     call ocomplot(or,nor,vsin,omegacp,psip,(Ftcomplex(:,1))/psip**2)
     call legendline(.1,.9,258,'F!de!d at imag(!Aw!@)=0')
!  call orealplot(or,nor,vsin,omegacp,psi,real(Ficomplex(:,1)))
     call ocomplot(or,nor,vsin,omegacp,psip,(Ficomplex(:,1)/psip**2))
     call legendline(.1,.9,258,'F!di!d')
     call multiframe(0,0,0)
     call pltend()
  endif
  if(lplot2)then
     sqpsi=sqrt(psip)
     tqpsi=psip**0.75
     uqpsi=psip**1.5
     qqpsi=psip**.25
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
     call contourl(real(forcecomplex)/psip**2.5,cworka,nor,nor,noi,zclv,icl, &
          or/tqpsi,oi/uqpsi,icsw)
     call legendline(0.1,-.1,0,'real')
     icl=-1
     zclv(1)=0.
     call dashset(0)
     call contourl(real(forcecomplex)/psip**2.5,cworka,nor,nor,noi,zclv,icl, &
          or/tqpsi,oi/uqpsi,icsw)
     icl=0
     zclv(1)=20
     call color(2)
     call dashset(3)
     call contourl(imag(forcecomplex)/psip**2.5,cworka,nor,nor,noi,zclv,icl, &
          or/tqpsi,oi/uqpsi,icsw)
     call legendline(.65,-.1,0,'imag')
     icl=-1
     zclv(1)=0.
     call dashset(0)
     call contourl(imag(forcecomplex)/psip**2.5,cworka,nor,nor,noi,zclv,icl, &
          or/tqpsi,oi/uqpsi,icsw)
     call color(15)
     call fwrite(kin/qqpsi,iwidth,3,string)
     call legendline(0.05,1.04,258,'k/!Ay!@!u1/4!u='//string)
     call fwrite(omegacp,iwidth,2,string)
     call legendline(.45,1.04,258,'!AW!@='//string)
     call fwrite(psip,iwidth,2,string)
     call legendline(.8,1.04,258,'!Ay!@='//string)
     call color(6)
!     write(*,*)omegap,.001*sqrt(psip)
     if(imag(omegap).gt.0.0011*sqrt(psip))& 
          call polymark(real(omegap/tqpsi),imag(omegap/uqpsi),1,3)
     call color(15)
     call pltend
  endif
end subroutine fomegacont
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine iterfindroot(psip,vsin,Omegacp,omegap,kin,isigma,lplot,ires)
  integer, parameter :: niter=12
  real :: psip,vsin,Omegacp,kin
  integer :: isigma,ires,nunconv=3
  logical lplot
  complex :: omegap,  Fec,Fic,Fsum
  real, dimension(0:niter) :: Frit,Fiit
     zoif=.001  ! Iteration minimum oi limit factor.
     nzo=0
     omegap=complex(0.9*sqrt(psip)/8.,.9*sqrt(psip)/8./(1.+vsin))
     Frit(0)=max(real(omegap),-2e-3)
     Fiit(0)=imag(omegap)
     call electronforce(Fec,omegap,kin,Omegacp,psip,vsin,isigma)
     call      ionforce(Fic,omegap,kin,Omegacp,psip,vsin,isigma)
     FE=kin**2*psip**2*128./315.
     Fsum=Fec+Fic+FE
     err=0
     do i=1,niter
        if(lplot)write(*,'(a,2i3,5f10.6)')'i,nharm,omegap,Fsum,err=',i&
             &-1,inharm(),omegap,Fsum,err
        ires=i
        call complexnewton(Fsum,omegap,kin,err,psip,isigma,vsin,Omegacp)
        Frit(i)=max(real(omegap),-2e-3)
        Fiit(i)=imag(omegap)
        if(.not.abs(omegap).lt.1.e6)write(*,*)'Iterfindroot',i,psip,vsin,omegap
        if(imag(omegap).lt.zoif*sqrt(psip))then
           nzo=nzo+1
           if(nzo.ge.nunconv)then
              if(lplot)write(*,'(a,i2,a,g10.3)')'Uncoverged after',nzo,'&
                   & omegai less than',zoif*sqrt(psip)
              err=1.
              ires=0
              goto 1
           endif
           omegap=complex(real(omegap),zoif*sqrt(psip))
        endif
        if(.not.abs(omegap).lt.1.e3)then
           ires=0
           write(*,*)'Iterfind diverging',omegap,err
           goto 1
        endif
        call electronforce(Fec,omegap,kin,Omegacp,psip,vsin,isigma)
        call      ionforce(Fic,omegap,kin,Omegacp,psip,vsin,isigma)
        Fsum=Fec+Fic+FE
!        write(*,*)err
        if(err.lt..5e-4)goto 1
        if(err*abs(omegap).lt.1.e-6)then
           write(*,*)'Apparent convergence at low omegap'
           goto 1
        endif
     enddo
     i=i-1
     if(lplot)write(*,*)'Unconverged after',i,' iterations'
1    continue
     if(lplot)then
        write(*,'(a,i4,5f10.6)')'i,omegap,Fsum,err=',i,omegap,Fsum,err
        if(err.ne.1.and.i.ne.niter)then
           do j=0,i-1
              call color(6)
              call polymark(max(Frit(j),-2e-3),Fiit(j),1,ichar('0')+j)
           enddo
        endif
     endif
end subroutine iterfindroot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine complexnewton(Fsum,omegap,kin,err,psip,isigma,vsin,Omegacp)
  ! Take a Newton Step in finding the complex omega root of
  ! complex force by inverting the 2x2 matrix from adjacent evaluations.
  ! Return the Fsum,omegap,err new values.
  complex :: Fsum,omegap
  real :: err,psip,vsin,Omegacp,kin
  real :: eps1=.05,J11,J12,J21,J22,domega1,domega2,det
  complex :: om1,om2,om3,f1,f2,f3,domega,Fec,Fic
  det=1.
  FE=kin**2*psip**2*128./315.
  f1=Fsum
  om1=omegap
  ! Calculate the Jacobian's coefficients.
  om2=om1+eps1*abs(real(om1))
!write(*,*)'compnewt call, kin=',kin
  call electronforce(Fec,om2,kin,Omegacp,psip,vsin,isigma)
  call      ionforce(Fic,om2,kin,Omegacp,psip,vsin,isigma)
  f2=Fec+Fic+FE
!  write(*,'(6g12.4)')om1,om2,f2
  if(.not.om2-om1.ne.0)stop 'om2-om1=0'
  J11=real(f2-f1)/real(om2-om1)
  J21=imag(f2-f1)/real(om2-om1)
  om3=om1+complex(0.,eps1*abs(imag(om1)))
  call electronforce(Fec,om3,kin,Omegacp,psip,vsin,isigma)
  call      ionforce(Fic,om3,kin,Omegacp,psip,vsin,isigma)
  f3=Fec+Fic+FE
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
!  err=abs(domega/omegap)
  err=abs(domega/omegap)
!  write(*,*)'domega1,domega2,omega',domega1,domega2,omega
end subroutine complexnewton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Plotting routines.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine lplot1(or,oi,nor,noi,vsin,omegacp,kin,psi,forcecomplex)
  real :: or(nor),oi(noi),kin
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
  call fwrite(vsin,iwidth,2,string)
  call legendline(0.02,1.04,258,'v!ds!d='//string)
  call fwrite(omegacp,iwidth,2,string)
  call legendline(.27,1.04,258,'!AW!@='//string)
  call fwrite(kin,iwidth,2,string)
  call legendline(.52,1.04,258,'k='//string)
  call fwrite(psi,iwidth,3,string)
  call legendline(.77,1.04,258,'!Ay!@='//string)
  call dashset(0)
end subroutine lplot1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine orealplot(or,nor,vsin,omegacp,psi,orealforce)
  real ::  or(nor)
  real ::  orealforce(nor)

  call minmax(orealforce,nor,fmin,fmax)
  call pltinit(0.,or(nor),fmin,fmax)
  call axis; call axis2
  call axlabels('real(!Aw!@)','real(Force)')
  call polyline(or,orealforce,nor)
end subroutine orealplot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ocomplot(or,nor,vsin,omegacp,psi,ocomforce)
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
     call ionforce(Fiarray(i),omegaFi(i),kin,Omegacin,psi,vsin,isigma)
     call electronforce(Fearray(i),omegaFi(i),kin,Omegacin,psi,vsin,-1)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
