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

!    call pfset(3)
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
    character*100 annote,tban
    if(real(omegag).eq.0.)omegag=(.1,0.00000)
    write(*,*)omegag
    omegaonly=omegag
    if(psig.le.0)psig=.5
    isigma=-1
 !   vshift=1.
    write(annote,'(''!Ay!@='',f5.3,'' !Aw!@=('',f5.3'','',f5.3,'')'', '//&
         ' ''  v!ds!d='',f5.3)')psig,real(omegag),imag(omegag),vshift
!    call fvinfplot
    call FgRepelEint(Ftotalg,isigma)
    tbmax=tbr(nge-2)*0.9
    tbfac=5.*nint(tbmax/Wgarrayp(nge)/5.)
!    write(*,*)tbfac,tbmax,int(nge*.9)
    write(*,*)'Repelling Ftotalpg',Ftotalpg
    write(*,*)'Repelling Ftotalrg',Ftotalrg
    call dcharsize(.025,.025)
    call multiframe(2,2,0)
    call pltinit(vinfarrayr(nge),vinfarrayr(1),0.0001,Wgarrayp(nge))
!       call axlabels('v!d!A;!@!d','W')
    call axis
    call axlabels('','W!d!A|!@!d')
    call legendline(0.5,1.06,258,annote(1:lentrim(annote)))
    call winset(.true.)
    call polyline(vinfarrayr,Wgarrayr,nge)
!    call polymark(vinfarrayr,Wgarrayr,nge,ichar('|'))
    call dashset(2)
    call color(5)
    call fwrite(tbfac,iwidth,0,tban)
    call polyline(vinfarrayr,tbr/tbfac,nge)
    call legendline(.5,.3,0,' !Bt!@!dorbit!d/'//tban(1:iwidth))
    call dashset(0)
    call color(15)
    call polymark(vinfarrayr,Wgarrayp(nge)*.98+1.e-6*Wgarrayp,nge,ichar('|'))
!    call pltend
    call minmax(forcegp,2*nge,pmin,pmax)
    call minmax(forcegr,2*nge,rmin,rmax)
    call pltinit(vinfarrayr(nge),vinfarrayr(1),min(pmin,rmin),max(pmax,rmax))
    call axis
    call axlabels('v!d!A;!@!d','dF/dv!d!a;!@!d')
    call legendline(.3,.9,258,'Reflected')
    call color(1)
    call polyline(vinfarrayr,real(forcegr),nge)
    call legendline(.2,.12,0,' real')
    call color(2)
    call dashset(2)
    call polyline(vinfarrayr,imag(forcegr),nge)
    call legendline(.2,.06,0,' imag')
    call dashset(0)
    call color(15)
    call pltinit(vinfarrayp(1),vinfarrayp(nge),0.,Wgarrayp(nge))
    call axis
    call axis2
    call winset(.true.)
    call polymark(vinfarrayp,Wgarrayp(nge)*.98+1.e-6*Wgarrayp,nge,ichar('|'))
    call polyline(vinfarrayp,Wgarrayp,nge)
    call dashset(2)
    call color(5)
    call polyline(vinfarrayp,tbp/tbfac,nge)
    call dashset(0)
    call color(15)
    call pltinit(vinfarrayp(1),vinfarrayp(nge),min(pmin,rmin),max(pmax,rmax))
    call legendline(.3,.9,258,'Passing')
    call axis
    call axis2
    call axlabels('v!d!A;!@!d','')
    call color(1)
    call polyline(vinfarrayp,real(forcegp),nge)
    call color(2)
    call dashset(2)
    call polyline(vinfarrayp,imag(forcegp),nge)
    call dashset(0)
    call pltend
    call multiframe(0,0,0)
  end subroutine testFrepel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine testAttract
    use shiftgen
    use shiftmode
    real, dimension(nge) :: vpsiarrayp
    complex :: Ftotalg
    character*40 annote,ffan
    lioncorrect=.false.
    if(real(omegag).eq.0)omegag=(.1,0.001000)
    omegaonly=omegag
    if(psig.ge.0)psig=-.5
    isigma=-1    
    if(vshift.ne.0.)write(*,*)'WARNING testAttract with vshift=',vshift
!    write(*,*)'Entered testattract'
    write(annote,'(''!Ay!@='',f5.3,'' !Aw!@=('',f5.3'','',f5.3,'')'')')&
         psig,real(omegag),imag(omegag)
    call dcharsize(.02,.02)
!    call FgAttractEint(Ftotalg,isigma)
    call FgEint(Ftotalg,isigma)  ! Generic call is the same.
!    write(*,*)'Return from FgAttractEint. Ftotalg=',Ftotalg
    vpsiarrayp=sqrt(2.*(Wgarrayp(1:nge)-psig))
!    call fvinfplot
    call multiframe(1,2,0)
    call minmax(forcegp,2*nge,pmin,pmax)
    call minmax(forcegr,2*nge,rmin,rmax)
    fpfac=5.*max(int(min(abs(rmax/pmax),abs(rmin/pmin))/5.),1)
    call fwrite(fpfac,iwidth,0,ffan)
    call pltinit(vinfarrayr(nge),vinfarrayr(1)*1.01,min(pmin,rmin),max(pmax,rmax))
    call axis; call axis2
    call axlabels('v!d!Ay!@!d','dF/dv!d!Ay!@!d')
    call legendline(0.5,1.03,258,annote(1:lentrim(annote)))
    call legendline(0.3,.9,258,'Trapped')
    call winset(.true.)
    call polymark(vinfarrayr,(max(pmax,rmax)*.97+vinfarray*1.e-7),nge&
         &,ichar('|'))
    call color(1)
    call polyline(vinfarrayr,real(forcegr),nge)
    call legendline(.05,.1,0,' real')
    call color(2)
    call dashset(2)
    call polyline(vinfarrayr,imag(forcegr),nge)
    call legendline(.05,.05,0,' imag')
    call color(15)
    call dashset(0)
    call pltinit(vpsiarrayp(1),vpsiarrayp(nge),min(pmin,rmin),max(pmax,rmax))
    call axis; call axis2 
    call axlabels('v!d!Ay!@!d','')
    call legendline(0.3,.9,258,'Passing')
    call winset(.true.)
    call polymark(vpsiarrayp,(max(pmax,rmax)*.97+vpsiarray*1.e-7),nge&
         &,ichar('|'))
    call color(3)
    call polyline(vpsiarrayp,fpfac*imag(forcegp),nge)
    call legendline(.4,.1,0,' realx'//ffan(1:iwidth))
    call color(4)
    call dashset(2)
    call polyline(vpsiarrayp,fpfac*real(forcegp),nge)
    call legendline(.4,.05,0,' imagx'//ffan(1:iwidth))
    call dashset(0)
    call pltend
    
    call multiframe(0,0,0)
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
    lioncorrect=.true.
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
    call tsparse(ormax,oi,nvs,isw)
    write(*,*)'vshift,psig,nvs',vshift,psig,nvs
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
! Test of ionforce, which takes electron omega arguments, whereas omega
! here has been specified in ion units.
          call ionforce(fion(i),omegag/sqrt(1836.),omegag/sqrt(1836.),psig,vshift,1836.)
          fion(i)=fion(i)/psig**2/2.
          frcomplex(i)=frcomplex(i)/psig**2
          diffmax=max(diffmax,abs(fion(i)-frcomplex(i)))
       enddo
       write(*,*)'diffmax=',diffmax
       if(j.eq.1)then
          write(*,*)'Fimmobile/2=',Fimmobile,' Fdirect=',frcomplex(nor)
          write(*,'(a,f9.5,a,f9.5)')' vshift=',vshift,' psig=',psig
          call pltinit(0.,or(nor),-0.8*Fimmobile,1.2*Fimmobile)
          call charsize(.02,.02)
          call axis
          call axptset(0.,1.)
          call ticrev
          call altxaxis(1./sqrt(1836.),1./sqrt(1836.))
          call legendline(.35,1.13,258,'real(!Aw!@)/!Aw!@!dpe!d')
          call axptset(1.,0.)
          call ticlabtog
          call altyaxis(1.,1.)
          call ticlabtog
          call axptset(0.,0.)
          call ticrev
          call axlabels('real(!Aw!@)/!Aw!@!dpi!d','!p!o~!o!qF!di!d/!Ay!@!u2!u')
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
    call tsparse(ormax,oi,nvs,isw)
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
  subroutine tsparse(ormax,oi,nvs,isw)
    use shiftgen
    character*30 argument
    isw=0
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
       if(argument(1:2).eq.'-c')call pfset(-3)
       if(argument(1:2).eq.'-s')then
          read(argument(3:),*)j
          isw=isw+j
       endif
       if(argument(1:2).eq.'-h')goto 1
    enddo
    if(isw.eq.0)isw=1
    return
1   continue
    write(*,*)' Usage: testshiftgen [-p,-v,-i,-zm,-or,-oi,-oc,-kg,-s,-n -h]'
    write(*,101) '-p..         set psi         [',psig
    write(*,101) '-v..         set vshift(max) [',vshift
    write(*,101) '-zm..        set zmax        [',zm
    write(*,101) '-or.. -oi..  set omega       [',ormax,oi
    write(*,101) '-oc.. -kg..  set O_c, set k  [',Omegacg,kg
    write(*,102) '-n..         set number of vs[',nvs
    write(*,102) '-s..         set switches    [',isw
    write(*,*)' s=1 Frepel, 2 Fattract, 4 Fr(omega), 8 PlotForce, 16&
         & denem',', LofW, 64 SumHarm'
    write(*,*)' s=128 fvinfplot.'
!    write(*,'(a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,i3,a,f8.3,a,f8.3)') 'psi&
!         &=',psig,' zm=',zm,' omax=',ormax,' v=',vshift,' nv=',nvs ,'&
!         & oc=',Omegacg,' kg=',kg
    101 format(a,6f8.4)
    102 format(a,i5)
  end subroutine tsparse
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine testdenem
  use shiftgen
  integer, parameter :: nvs=100
  real, dimension(nvs) :: vs,denem
  denem=1.
  psig=.1
  vsmax=2.
  do i=1,nvs
     vrshift=i*vsmax/nvs
     vs(i)=vrshift
     call dfefac(denem(i))
  enddo
  call autoplot(vs,denem,nvs)
  call axlabels('vshift','fefac')
  write(*,'(''Factor by which dfe trapped is multiplied for psi='',f8.3)')psig
  call pltend
end subroutine testdenem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use shiftgen
integer :: nvs=1,isw=0
real :: ormax=0.,oi=0.
call tsparse(ormax,oi,nvs,isw)
omegag=complex(ormax,oi)
! Some tests might interfere with others.
  if(isw-2*(isw/2).eq.1) call testFrepel
  isw=isw/2 ! 2
  if(isw-2*(isw/2).eq.1) call testAttract
  isw=isw/2 ! 4
  if(isw-2*(isw/2).eq.1) call Frepelofomega
  isw=isw/2 ! 8
  if(isw-2*(isw/2).eq.1) call plotionforce(.01,1.,0.)
  isw=isw/2 ! 16
  if(isw-2*(isw/2).eq.1) call testdenem
  isw=isw/2 ! 32
  if(isw-2*(isw/2).eq.1) call testLofW
  isw=isw/2 ! 64
  if(isw-2*(isw/2).eq.1) call testSumHarm
  isw=isw/2 ! 128
  if(isw-2*(isw/2).eq.1) call fvinfplot
end program
