  ! Calculate and plot the magnetized representation in terms of sums
  ! of Bessel functions. Uses  RIBESL(X,ALPHA,NB,IZE,B,NCALC)
  ! in libmodbess  from toms715. 

! The expression of interest is exp(-xi^2)I_m(xi^2) where
! xi=k/\Omega \sqrt(Te/me) = k v_t/\Omega.   v_t=1.
! The required number of orders is nb such that nb*\Omega/k ~ 3v_t.
! So nb ~ 3 xi. Phase velocity of harmonics is n\Omega/k=n/xi.

  integer, parameter :: nbmax=40,nxi=4
  real :: B(nbmax,nxi),vphase(nbmax)
  character*20 string
! Test section for calling the RIBESL routine.
  alpha=0.   ! Integer order.
  ximax=4.
  vmax=3.3
  call pfset(3)
  call pltinit(0.,vmax,0.,.5)
  call charsize(.018,.018)
  call axis
  call axis2
  call axlabels('!Bv!dm!d=m/!Ac!B!dt!d=m!AW!B/kv!dt!d!@', &
       '!Ac!B!dt!d!@ exp(-!Ac!B!dt!d!u2!u)' &
       //'I!dm!d(!Ac!B!dt!d!u2!u)')

  do ixi=1,nxi
     xi=ixi*ximax/nxi
     x=xi**2                  ! Argument
     nb=nint(vmax*xi)+1       ! Number of orders calculated.
     write(*,*) 'For x=       exp(-x)I_n(x): n=0,1,2, ...',nb-1
     nx=20
     xmax=20.
     ize=2                    ! exp(-x)I_n(x) to be calculated.
     do i=1,nb
        vphase(i)=(i-1)/xi
     enddo
     call RIBESL(X,ALPHA,NB,IZE,B(1,ixi),NCALC)
     B(:,ixi)=xi*B(:,ixi)     ! Scale as if an integral.
     if(ncalc.eq.nb)then    ! All orders calculated correctly.
        write(*,'(10f8.4)')x,B(1:nb,ixi)
!        write(*,'(''vphase  '',10f8.4)')vphase(1:nb)
     else
        write(*,*)'Bessel function imprecision',ncalc,' of',nb
     endif
     call color(ixi)
     call polymark(vphase,B(1,ixi),nb,ixi)
     call fwrite(xi,iwidth,2,string)
     call legendline(.7,.9-.05*ixi,ixi,' '//string)
  enddo
  call color(15)
  call legendline(.7,.9,257,'  !Ac!B!dt!d!@')
  call pltend()
  end
