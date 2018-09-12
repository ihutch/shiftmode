! Functionalized version of the hgfunction code.
!****************************************************************
! This is exp(X^2)*erfc(X)
      FUNCTION expERFCC(X)
      Z=ABS(X)      
      T=1./(1.+0.5*Z)
      expERFCC=T*EXP(-1.26551223+T*(1.00002368+T*(.37409196+             &
     &    T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+        &
     &    T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
      IF (X.LT.0.) expERFCC=2.*exp(z**2)-expERFCC
      END
!********************************************************************
      subroutine fhgfunc(psi,xrn,nx,xwidth,gint,hint,gjint,pint)
! Integrate functions over a sech4 hole, using input psi and a total
! range xrn, and nx intervals. chi=sqrt(phi)
! h(chi) is the sum of passing and trapped force
! g(chi) is the total trapped force.
! gj(chi) is the passing force without the boundary term.
! p(chi) is minus the trapped force without the boundary term.
! Thus p=gj-h
! But also, p(chi) is the passing particle density deficit, which is 
! the total number of trapped particles. So it can be compared with
! the trapped number calculated in trapnum.
      real hint,gint,gjint,pint,xrn,psi,xwidth
      integer nx
      tworpi=2./sqrt(3.1415926)
      hint=0.
      gint=0.
      gjint=0.
      pint=0.
      dx=xrn/nx
      do i=1,nx
         x=i*dx
         x=sqrt(psi/cosh(x/xwidth)**4)   ! chi = sqrt(phi)
         epfchi=experfcc(x)
         hint=hint+dx*(-tworpi*x+(2.*x**2-1)*epfchi+1.)
         gint=gint+dx*(tworpi*x-1.+epfchi)
         gjint=gjint+dx*(-tworpi*x+2.*(x**2-1)*epfchi+2.)
         pint=pint+dx*(1.-epfchi)
      enddo
      chi=sqrt(psi)
      epfchi=experfcc(chi)
      hint0=dx*(-tworpi*chi+(2.*chi**2-1)*epfchi+1.)
      gint0=dx*(tworpi*chi-1.+epfchi)
      gjint0=dx*(-tworpi*chi+2.*(chi**2-1)*epfchi+2.)
      pint0=dx*(1.-epfchi)
!      write(*,*)'chimax=',x,'  dx=',dx,'  gint0=',gint0
! Double to represent \pm x-integral, but don't double the zero pos.
      gint=2.*gint+gint0
      hint=2.*hint+hint0
      gjint=2.*gjint+gjint0
      pint=2.*pint+pint0
      end
!************************************************************************
! Evaluate the total number of trapped particles for a hole in which
! the potential is given by the function phifunc(z,psi,phiwidth)
! and the distribution is ffuncw(W,beta), with W in units of Te.
!***********************************************************************
      real function phifunc(zf,psi,phiwidth)
      real zf,psi
      phifunc=psi/cosh(zf/phiwidth)**4
      end
!**********************************************************************
      real function ffuncw(W,beta)
      real W,beta
! Unshifted Schamel distribution function for W in units of T, and beta
! the trapped temperature relative to T.
      if(W.lt.0)then
         ffuncw=exp(-beta*W)
      else
         ffuncw=exp(-W)
      endif
      end
!**********************************************************************
      real function trapnum(psi,beta,phiwidth)
      real psi,beta
      integer nnz,nnv
      parameter (nnz=100,nnv=100)

      zlen=20.
      trapnum=0.
      dz=2.*zlen/nnz
      do i=1,nnz
         z=-zlen+(i-0.5)*dz
         phi=phifunc(z,psi,phiwidth)
         vmax=sqrt(2*phi)
         anz=0.
         dv=vmax/nnv
         do j=1,nnv
            v=(i-0.5)*dv
            w=-phi+v**2/2.
            fv=ffuncw(w,beta)
            anz=anz+fv*dv
         enddo
         trapnum=trapnum+2.*anz*dz
      enddo
      end
!**********************************************************************
      real function flattrap(psi,phiwidth)
      real psi
      parameter (nnz=100)
      zlen=20.
      flattrap=0.
      dz=2.*zlen/nnz
      do i=1,nnz
         z=-zlen+(i-0.5)*dz
         phi=phifunc(z,psi,phiwidth)
         vs=sqrt(2.*phi)
         flattrap=flattrap+2.*vs*dz
      enddo
      end
!*********************************************************************
! Return the low-k growth rate for anisotropic temperature
! normalized to k^2 T_parallel/m_e, for a schamel hole with trapped
! "temperature" beta
      real function growthlk(psi,beta,Ty)
      xrn=10.
      nx=101
      xwidth=4.
      call fhgfunc(psi,xrn,nx,xwidth,gint,hint,gjint,pint)
      txoty=1/Ty
      growthlk=sqrt(Ty*(((1-txoty)*hint+txoty*(1-1/beta)*gint))/hint)
      end
