! The following are routines that need to use shiftmode module.
! They are needed by shiftgen module routines, but are separated
! so as to make shiftgen formally require no use of shiftmode.
!****************************************************************
! Potentially free-standing routine.
 subroutine Ftrapcompare(i,dfe,dfeperp,obi,resdenom&
       &,resdprev,cdvpsi,Ftotal,dFdvpsig,vpsi)
    use shiftgen
    use shiftmode
    complex :: Ftotal,dFdvpsi,resdenom,resdprev,cdvpsi,dFdvpsig
    if(i.lt.nge)then
    psi=-psig                     ! psi is the positive depth
    omegaonly=omegag
!          write(*,*)'FgA calling shiftmode initialize'
    call initialize
    call dFdvpsidvy(vpsi,dFdvpsi,tbi,xln)
    Fnonres(i)=dFdvpsi*(omegag*dfe-(omegag-omegaonly)*dfeperp)
    Fnonres(i)=Fnonres(i)+sqm1g*real(Fnonres(i)-Fnonres(i-1))  &
         /real(omegabg(i)-omegabg(i-1))*obi
    Ftrap(i)=0.5*(Fnonres(i)/resdenom &
         +Fnonres(i-1)/resdprev)*cdvpsi
    Ftotalmode=Ftotalmode+2.*Ftrap(i)       ! Add to Ftotal integral.
    write(*,'(i3,'' shiftmode'',2f10.5,$)')i,tbi,real(Ftotalmode)
    write(*,*)dFdvpsi
    write(*,'(i3,'' shiftgen '',2f10.5,$)')i,2.*taug(ngz),real(Ftotal)
    write(*,*)dFdvpsig
    else
    write(*,*)'Calling FtEint'
    call FtEint(Ftotalmode,dfperpdWperp,fperp)
    write(*,*)'Fotalmode vs Ftotal (gen)'
    write(*,*)Ftotalmode
    write(*,*)Ftotal
    endif
  end subroutine Ftrapcompare
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! Potentially free-standing routine
  subroutine Fpasscompare(i) ! Do shiftmode version too.
    use shiftgen
    use shiftmode
    complex :: passforce,sumpassforce
    if(i.eq.1)sumpassforce=0.
    psi=-psig                     ! psi is the positive depth
    omega=omegag
    omegad=omega
    call initialize
    call dFdvinfdvy(vinfarrayp(i),passforce)
!This does not work to give agreement.
!          sumpassforce=sumpassforce+passforce*dfdWpar(vinfarrayp(i)&
!               &,fvinf)*dvinf*omegag
    write(*,'(i3,'' shiftmode passing'',f8.4,4e12.5)')i,psig &
         &,passforce!,sumpassforce
    write(*,'(i3,'' shiftgen  passing'',f8.4,4e12.5)')i,psig &
         &,Forcegarray(i)!,Ftp
  end subroutine Fpasscompare
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
