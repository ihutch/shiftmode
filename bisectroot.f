      subroutine bisectroot(x0,x1,vfracfunc,nbi,frac,xval)
! Given two input variables x0,x1, which may be general real arrays, and
! a function interpolation routine vfracfunc that interpolates two such
! variables into a third (xval) and returns the value of the function to
! be solved and the interpolated variable by a call like this:
! fv=vfracfunc(x0,x1,frac,xval), repetitively bisect the interval 0 and
! 1 between x0 and x1, until a root has been found of the equation and
! return the corresponding xval interpolate and interpolation fraction
! frac.

!  On entry:
!    x0,x1 contain the endpoint variables (arrays or scalars).
!    nbi is the maximum number of bisections allowed. So it gives the
!        precision of fraction determination: 1/2^nbi.
!    vfracfunc is the interpolation function
!  On exit:
!    frac if it is between 0 and 1 is the position fraction of the root
!    xval is the corresponding interpolated variable, if 0<=frac<=1.
!    Otherwise no solution has been found.

! The routine can accomodate situations where the x0 and x1 do not yield
! opposite signs of vfracfunc. It searches the region between them (but
! not beyond) in increasingly fine increments for a sign
! crossing. Finding the crossing closest to x0 at the last fineness of
! division. If this search requires more than approximately nbi
! evaluations it fails. If nbi is set at approximately 32 or more, then
! convergence will occur to 32 bit (machine) precision. Roughly nbi is
! the number of significant bits in frac.

!      real x0(*),x1(*),xval(*)
      real frac,fval
      external fval,vfracfunc

      frac0=0.
      fv0=vfracfunc(x0,x1,frac0,xval)
      frac1=1.
      fv1=vfracfunc(x0,x1,frac1,xval)
      ii=0
!      write(*,*)'Starting func values',x0,fv0,x1,fv1
      if(sign(1.,fv0).eq.sign(1.,fv1))then
         write(*,*)'No root',x0,fv0,x1,fv1
         frac=-1.
! Disable sign crossing search.
         return
! Search for sign crossing. 
         do ii=1,1+int(log(float(nbi))/log(2.))
            rjm1=2**ii
            do j=1,2**ii-1,2
               frac=j/rjm1
               fv=vfracfunc(x0,x1,frac,xval)
               if(sign(1.,fv).ne.sign(1.,fv0))then
! Found one. Start the bisection proper.
                  fv1=fv
                  frac0=frac-1/rjm1
                  frac1=frac
                  goto 1
               endif
            enddo
         enddo
! Failed to find a sign crossing.
         frac=-1.
         return
      endif
 1    continue
! Now we enclose a root appropriately. Bisect.
      do i=1,nbi
         frac=(frac0+frac1)*0.5
         fv=vfracfunc(x0,x1,frac,xval)
! Which value to replace.
         if(fv.eq.0)then
            goto 2
         elseif(sign(1.,fv).eq.sign(1.,fv0))then
            if(frac0.eq.frac)goto 2
            frac0=frac
         else
            if(frac1.eq.frac)goto 2
            frac1=frac
         endif
      enddo
 2    continue
      end
!$$$c******************************************************************
!$$$      program bisectiontest
!$$$      real myfracfunc
!$$$      external myfracfunc
!$$$      
!$$$      x0=-2.5
!$$$      x1=200.
!$$$      nbi=33
!$$$      call bisectroot(x0,x1,myfracfunc,nbi,frac,xval)
!$$$      write(*,*)'Finished bisectroot',x0,x1,nbi,frac,xval
!$$$      if(frac.lt.0.)write(*,*)'Failed to find zero crossing.'
!$$$      end
!$$$c******************************************************************
!$$$      real function myfracfunc(x0,x1,frac,xval)
!$$$c Example interpolation function
!$$$      xval=x0+(x1-x0)*frac
!$$$      myfracfunc=xval**2-1.3
!$$$      write(*,*)'frac,xval,myfracfunc',frac,xval,myfracfunc
!$$$      end
!$$$      
