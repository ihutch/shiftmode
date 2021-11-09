!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine testLofW
    use shiftgen
    use shiftmode
    integer, parameter :: nw=60
    complex, dimension(-ngz:ngz,nw) :: Lgw
    real, dimension(-ngz:ngz,nw) :: zgw,vgw,taugw
    complex :: dForceg(nw),forcet(nw),dFdirect(nw)
    integer :: iwsa(nw)
    real :: Wn(nw)
    real, dimension(-ngz:ngz) :: ones
    logical :: lplotmz=.true.
    ones=1.
    omegag=(5.,0.02)
    omegaonly=omegag
    psig=-.5
    omegad=omegag
    psi=abs(psig)
    isigma=-1
    if(lplotmz)call pltinit(-zm,zm,min(0.,psig),max(1.8*psig,.8*abs(psig)))
    if(lplotmz)call axis
    if(lplotmz)call axlabels('z','W')
    do i=1,nw
       if(psig.gt.0)Wg=1.1*psig*i/nw
       if(psig.lt.0)Wg=psig+1.42*abs(psig)*i/nw
       write(*,*)'Wg=',Wg,' omegag=',omegag
       Wn(i)=Wg
       call LofW(Wg,isigma,dForceg(i))
!       dFdirect(i)=dFordirect
       call Fdirect(Wg,isigma,dFdirect(i))
       if(lplotmz)call color(mod(i-1,15)+1)
       if(lplotmz)call polymark(zg(-ngz:ngz),Wg*ones(-ngz:ngz),2*ngz+1,10)
       iwsa(i)=min(iws,ngz)
       Lgw(-ngz:ngz,i)=Lg(-ngz:ngz) ! Save for plotting.
       zgw(-ngz:ngz,i)=zg(-ngz:ngz)
       vgw(-ngz:ngz,i)=vg(-ngz:ngz)
       if(Wg.gt.0.and.Wg.lt.psig)Lgw(:,i)=Lg+(vg-vg(-ngz))/omegag
       taugw(-ngz:ngz,i)=taug(-ngz:ngz)
!       write(*,*)'dForceg=',dForceg(i),' taug',taug(ngz)
!       write(*,'(10f8.4)')(zg(j),j=-ngz,ngz)
!       write(*,'(10f8.3)')(taug(j),j=-ngz,ngz)
!       write(*,'(10f8.4)')(vg(j),j=-ngz,ngz)
!       write(*,'(10f8.4)')(real(Lg(j)),j=-ngz,ngz)
!       write(*,'(10f8.4)')(imag(Lg(j)),j=-ngz,ngz)
!       write(*,*)'dForceg=',dForceg(i),' taug',taug(ngz)
! Testing against shiftmode. We must use positive velocity in the
! shiftmode calculation because it is not set up to use negative
! integration direction. The sign change is then in the print out. 
       if(Wg.lt.0)then          ! Trapped
          xL=20.
          vpsig=sqrt(2.*(Wg-psig))
          call dFdvpsidvy(vpsig,forcet(i),tb,xlent)
          tdur=tb/2
          xLend=xlent
       elseif(psig.lt.0)then   ! Untrapped.
          vinf=sqrt(2.*Wg)
          xL=zm
          call initialize          
          call dFdvinfdvy(vinf,forcet(i))
          tdur=tau(nx)
          xLend=xL
       endif
       write(*,'(a, 5f10.5)')'End position        ',zg(ngz),xLend
       write(*,'(a, 5f10.4)')'Time duration       ',taug(ngz),tdur
       write(*,'(a, 5f10.6)')'Lg     Lt           ',Lg(ngz),-isigma*Lt(iend+1)
       write(*,'(a, 5f10.6)')'Forceg  real,imag     ',dForceg(i),forcet(i)
       write(*,'(a, 5f10.6)')'Fdirect real,imag     ',dFdirect(i)
          
    enddo
! Put xL back to default
    xL=20

    if(lplotmz)call pltend

    call pfset(3)
    call multiframe(2,1,0)
    if(psig.lt.0)call pltinit(-zm,zm,0.,-isigma*1.5*sqrt(2.*abs(psig)))
    if(psig.ge.0)call pltinit(-zm,zm,-1.5*sqrt(2.*abs(psig)),1.5*sqrt(2.*abs(psig)))
    call axis
    call axlabels('zg','vg')
    do i=1,nw
       call color(mod(i-1,15)+1)
       call polyline(zgw(:,i),vgw(:,i),2*ngz+1)
!       call polyline(zgw(:,i),-isigma*taugw(:,i)/taugw(ngz,i),2*ngz+1)
!       call polymark(zgw(:,i),taugw(:,i)/taugw(ngz,i),2*ngz+1,10)
       
       call polymark(zgw(iwsa(i),i),vgw(iwsa(i),i),1,1)
    enddo
    call color(15)

    call minmax2(taugw,2*ngz+1,2*ngz+1,nw,tmin,tmax)
    call pltinit(-zm,zm,0.,tmax)
    call axis
    call axlabels('zg','tau')
    do i=1,nw
       call color(mod(i-1,15)+1)
       call polyline(zgw(:,i),taugw(:,i),2*ngz+1)
    enddo
    call pltend

    call multiframe(2,1,3)
    call minmax2(real(Lgw),2*ngz+1,2*ngz+1,nw,amin,amax)
    call pltinit(-zm,zm,amin,amax)
    call axis
    call axlabels('z','Real(Lg[+(v-v!d!A;!@!d)/!Aw!@])')
    do i=1,nw
       call color(mod(i-1,15)+1)
       call polyline(zgw(-ngz:ngz,i),real(Lgw(-ngz:ngz,i)),2*ngz+1)
    enddo
    call color(15)
    call minmax2(imag(Lgw),2*ngz+1,2*ngz+1,nw,amin,amax)
    call pltinit(-zm,zm,amin,amax)
    call axis
    call axlabels('z','Imag(Lg[+(v-v!d!A;!@!d)/!Aw!@])')
    do i=1,nw
       call color(mod(i-1,15)+1)
       call polyline(zgw(-ngz:ngz,i),imag(Lgw(-ngz:ngz,i)),2*ngz+1)
    enddo
    call pltend
    call multiframe(0,0,0)

    call minmax(dForceg,2*nw,amin,amax)
    call pltinit(Wn(1),Wn(nw),amin,amax)
    call axis
    call axlabels('W','dForceg,dFdirect')
    call polyline(Wn,real(dForceg),nw)
    call color(4)
    call polyline(Wn,real(dFdirect),nw)
    call dashset(1)
    call polyline(Wn,imag(dFdirect),nw)
    call color(15)
    call polyline(Wn,imag(dForceg),nw)
    call dashset(0)
    call pltend
  end subroutine testLofW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine testFrepel
    use shiftgen
    complex :: Ftotalg
    character*40 annote
    omegag=(.1,0.00000)
    omegaonly=omegag
    psig=.5
    isigma=-1    
    vshift=1.
    write(annote,'(''!Ay!@='',f5.3,'' !Aw!@=('',f5.3'','',f5.3,'')'')')&
         psig,real(omegag),imag(omegag)
    call dcharsize(.018,.018)
    call multiframe(2,1,2)
       call FgRepelEint(Ftotalg,isigma)
       write(*,*)'Repelling Ftotalpg',Ftotalpg
       write(*,*)'Repelling Ftotalrg',Ftotalrg
       call pltinit(vinfarrayr(nge),vinfarrayp(nge),0.,Wgarrayp(nge))
       call axis
       call axlabels('v!d!A;!@!d','W')
       call legendline(0.5,0.9,258,annote(1:lentrim(annote)))
       call winset(.true.)
       call polymark(vinfarrayp,Wgarrayp,nge,1)
       call polyline(vinfarrayp,Wgarrayp,nge)
       call polyline(vinfarrayr,Wgarrayr,nge)
       call polymark(vinfarrayr,Wgarrayr,nge,2)
       call polyline(vinfarrayp,tbp/10.,nge)
       
       call polyline(vinfarrayr,tbr/10.,nge)
       call legendline(.2,.9,0,' t!dorbit!d/10')
!    call pltend
       call minmax(forcegp,2*nge,pmin,pmax)
       call minmax(forcegr,2*nge,rmin,rmax)
       call pltinit(vinfarrayr(nge),vinfarrayp(nge),min(pmin,rmin),max(pmax,rmax))
       call axis
       call axlabels('v!d!A;!@!d','dF/dv!d!a;!@!d')
       call color(1)
       call polyline(vinfarrayr,real(forcegr),nge)
       call polyline(vinfarrayp,real(forcegp),nge)
!    call polymark(vinfarrayr,real(forcegr),nge,1)
!    call polymark(vinfarrayp,real(forcegp),nge,1)
       call legendline(.6,.7,0,' real')
       call color(2)
       call polyline(vinfarrayr,imag(forcegr),nge)
       call polyline(vinfarrayp,imag(forcegp),nge)
!    call polymark(vinfarrayr,imag(forcegr),nge,2)
!    call polymark(vinfarrayp,imag(forcegp),nge,2)
       call legendline(.6,.8,0,' imag')
       call pltend
       call multiframe(0,0,0)
    call fvinfplot
  end subroutine testFrepel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine testAttract
    use shiftgen
    use shiftmode
    real, dimension(nge) :: vpsiarrayp
    complex :: Ftotalg
    character*40 annote
    omegag=(.1,0.001000)
    omegaonly=omegag
    psig=-.5
    isigma=-1    
    vshift=0.
!    write(*,*)'Entered testattract'
    write(annote,'(''!Ay!@='',f5.3,'' !Aw!@=('',f5.3'','',f5.3,'')'')')&
         psig,real(omegag),imag(omegag)
    call dcharsize(.018,.018)
    call multiframe(2,1,2)
!    call FgAttractEint(Ftotalg,isigma)
    call FgEint(Ftotalg,isigma)  ! Generic call is the same.
!    write(*,*)'Return from FgAttractEint. Ftotalg=',Ftotalg
    vpsiarrayp=sqrt(2.*(Wgarrayp(1:nge)-psig))
!    call pltinit(vinfarrayr(nge),vinfarrayp(nge),psig,Wgarrayp(nge))
    call pltinit(vinfarrayr(nge),vpsiarrayp(nge),psig,Wgarrayp(nge))
    call axis
    call axlabels('v!d!Ay!@!d','W')
    call legendline(0.5,0.9,258,annote(1:lentrim(annote)))
    call polymark(vinfarrayr,Wgarrayr,nge,1)
    call polyline(vinfarrayr,Wgarrayr,nge)
    call polyline(vpsiarrayp,Wgarrayp,nge)
    call polymark(vpsiarrayp,Wgarrayp,nge,2)
    
    call minmax(forcegp,2*nge,pmin,pmax)
    call minmax(forcegr,2*nge,rmin,rmax)
!    write(*,*)'pmin,pmax,rmin,rmax',pmin,pmax,rmin,rmax
    call pltinit(vinfarrayr(nge),vpsiarrayp(nge),min(pmin,rmin),max(pmax,rmax))
    call axis
    call axlabels('v!d!Ay!@!d','dF/dv!d!a;!@!d')
    call color(1)
    call polyline(vinfarrayr,real(forcegr),nge)
    call legendline(.6,.7,0,' real trapped')
    call color(2)
    call polyline(vinfarrayr,imag(forcegr),nge)
    call legendline(.6,.8,0,' imag trapped')
    call color(3)
    call polyline(vpsiarrayp,10.*imag(forcegp),nge)
    call legendline(.6,.6,0,' 10xreal passing')
    call color(4)
    call polyline(vpsiarrayp,10.*real(forcegp),nge)
    call legendline(.6,.5,0,' 10ximag passing')
    call pltend
    call multiframe(0,0,0)
!       stop
    call fvinfplot
    psi=-psig                     ! psi is the positive depth
    omega=omegag
    omegad=omega
    omegac=10.
    write(*,*)'testAttract: Calling shiftmode initialize'
    call initialize
    call SumHarmonics
    write(*,*)'Ftraptotal mode',Ftraptotal
    write(*,*)'Ftraptotal gen ',Ftotalrg
    write(*,*)'Fpasstotal mode',Fpasstotal
    write(*,*)'Fpasstotal gen ',Ftotalpg
    write(*,*)'Sum mode       ',Ftraptotal+Fpasstotal
    write(*,*)'Sum gen        ',Ftotalrg+Ftotalpg
!    write(*,*)'Ftotalg        ',Ftotalg
  end subroutine testAttract
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Frepelofomega
    use shiftgen
    integer, parameter :: nor=100
    real, dimension(nor) :: or
    complex, dimension(nor) :: frcomplex,fion
    character*30 string
!    complex :: Ftotalg
    kg=0.
    ormax=10.   !defaults
    oi=0.01
    psig=0.01
    vshift=0.
    isigma=-1
    Fimmobile=(128./315.)/2.       ! Now normalized *psig**2 
! Because frcomplex includes only one velocity direction.
    nvs=1
    call tsparse(ormax,oi,nvs)
    psig=abs(psig)
    vsmax=vshift
    ol=.4
    nl=int(nor*min(ol/ormax,1.))
    do j=1,nvs
       if(nvs.gt.1)vshift=vsmax*(j-1.)/(nvs-1.)
       do i=1,nor
!       or(i)=ormax*(i-1.)/(nor-1.)
          or(i)=ormax*float(i)/(nor)
          omegag=complex(or(i),oi)
          omegaonly=omegag
          call FgRepelEint(frcomplex(i),isigma)
          if(.not.real(frcomplex(i)).lt.1.e20)then
             write(*,*)'Frepelofomega Force Nan?',&
                  i,omegag,omegaonly,frcomplex(i),psig,isigma
             stop
          endif
! Test of ionforce
          call ionforce(fion(i),omegag/sqrt(1836.),omegag/sqrt(1836.),psig,0.,1836.)
          fion(i)=fion(i)/psig**2/2.
          frcomplex(i)=frcomplex(i)/psig**2
          diffmax=max(diffmax,abs(fion(i)-frcomplex(i)))
       enddo
       write(*,*)'diffmax=',diffmax
       if(j.eq.1)then
          write(*,*)'Fimmobile/2=',Fimmobile,' Fdirect=',frcomplex(nor)
          write(*,'(a,f9.5,a,f9.5)')' vshift=',vshift,' psig=',psig
          call pltinit(0.,or(nor),-0.8*Fimmobile,1.2*Fimmobile)
          call axis; call axis2
          call axlabels('real(!Aw!@)/!Aw!@!dpi!d','Force/!Ay!@!u2!u')
          call polymark(ormax,Fimmobile,1,1)
       endif
       call color(mod(j-1,15)+1)
       call polyline(or,real(frcomplex),nor)
       call fwrite(vshift,iwidth,2,string)
       call jdrwstr(wx2nx(ormax*.95),wy2ny(real(frcomplex(nor))),string(1:iwidth),-1.)
       call jdrwstr(wx2nx(ol),wy2ny(imag(frcomplex(nl))),string(1:iwidth),0.)
       if(j.eq.1)call legendline(.5,.1,0,' real')
       call dashset(2)
       call polyline(or,imag(frcomplex),nor)
       if(j.eq.1)call legendline(.5,.15,0,' imag')
       call dashset(0)
    enddo
    call pltend
  end subroutine Frepelofomega
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine testSumHarm
    use shiftgen
    use shiftmode
!    complex :: Ftotalg
    ormax=.1
    psig=-.1
    call tsparse(ormax,oi,nvs)
    if(oi.lt.0.00001)oi=.00001
    omegag=complex(ormax,oi)
    omegaonly=omegag
    write(*,*)'testSumHarm psig                  omegag,              Omegacg'
    write(*,*)psig,omegag,Omegacg
    isigma=-1
!    call FgEint(Ftotalg,isigma)  ! Generic call is the same.
!    write(*,*)'Ftotalg        ',Ftotalg
    call SumHarmonicsg(isigma)
    write(*,*)'FtotalSumg=',Ftotalsumg
    if(psig.lt.0)then
       psi=-psig                     ! psi is the positive depth
       omega=omegag
       omegad=omega
       Omegac=Omegacg
       k=kg
       write(*,*)omega,Omegac,k
       write(*,*)'testSumHarm: Calling shiftmode initialize'
       call initialize
       call SumHarmonics
       write(*,*)'Ftraptotal mode',Ftraptotal
       write(*,*)'Fpasstotal mode',Fpasstotal
       write(*,*)'Sum mode       ',Ftraptotal+Fpasstotal
       write(*,*)'Sum gen        ',FtotalSumg
       write(*,*)'SumM/SumG-1    ',(Ftraptotal+Fpasstotal)/FtotalSumg-1.
    endif
  end subroutine testSumHarm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tsparse(ormax,oi,nvs)
    use shiftgen
    character*30 argument
    do i=1,iargc()
       call getarg(i,argument)
       if(argument(1:2).eq.'-p')read(argument(3:),*)psig
       if(argument(1:2).eq.'-v')read(argument(3:),*)vshift
       if(argument(1:2).eq.'-i')read(argument(3:),*)isigma
       if(argument(1:3).eq.'-zm')read(argument(4:),*)zm
       if(argument(1:3).eq.'-or')read(argument(4:),*)ormax
       if(argument(1:3).eq.'-oi')read(argument(4:),*)oi
       if(argument(1:3).eq.'-oc')read(argument(4:),*)Omegacg
       if(argument(1:3).eq.'-kg')read(argument(4:),*)kg
       if(argument(1:2).eq.'-n')read(argument(3:),*)nvs
       if(argument(1:2).eq.'-h')goto 1
    enddo
    return
1   continue
    write(*,*)' Usage: testshiftgen [-p,-v,-i,-zm,-or,-oi,-oc,-n,-h]'
    write(*,'(a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,i3,a,f8.3,a,f8.3)')&
         'psi=',psig,' zm=',zm,' omax=',ormax,' v=',vshift,' nv=',nvs&
         ,' oc=',Omegacg,' kg=',kg
  end subroutine tsparse
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: nvs=1
  call tsparse(ormax,oi,nvs)
! Some tests might interfere with others.       
  call testLofW
  call testFrepel
  call testAttract
  call testSumHarm
  call Frepelofomega
  call plotionforce(.01,1.,0.)
end program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ionforce(Fi,omega,omegaon,psiin,vsin,mime)
! Calculate ion force for given parameters and ion to electron mass mime.
  use shiftgen
  complex :: Fi,omega,Ftotalg,omegaon
  real :: psiin,vsin,mime
  omegag=omega*sqrt(mime)
  omegaonly=omegaon*sqrt(mime)
  psig=psiin
  isigma=-1
  vshift=vsin
  call FgRepelEint(Ftotalg,isigma)
  Fi=2.*Ftotalg
  if(abs(real(omegag)-nint(real(omegag))).lt.1.e-5)then ! Debugging.
     write(*,*)omegag,Ftotalg
  endif
end subroutine ionforce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotionforce(psi,Typ,vsin)
  real :: psi,Typ,vsin
  integer, parameter :: nfi=100
  complex :: omegaon
  complex, dimension(nfi) :: Fiarray,omegaFi
  omegamax=10
  write(*,*)'psi=',psi,' vsin=',vsin
  do i=1,nfi
     omegaFi(i)=omegamax*(float(i)/nfi)/sqrt(1836.)+complex(0.,.01)/sqrt(1836.)
     omegaon=omegaFi(i)
     call ionforce(Fiarray(i),omegaFi(i),omegaon,psi,vsin,1836.)
     Fiarray(i)=Fiarray(i)/psi
  enddo
  call minmax(Fiarray,2*nfi,fmin,fmax)
  call pltinit(0.,omegamax/sqrt(1836.),fmin,fmax)
  call axis
  call axlabels('omega','Fi')
  call polyline(real(omegaFi),real(Fiarray),nfi)
  call polyline(real(omegaFi),imag(Fiarray),nfi)
  call pltend
end subroutine plotionforce
