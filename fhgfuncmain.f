!*********************************************************************
!Main program for testing fhgfunc routines.
      include 'fhgfunc.f'

      program fhgvfuncmain
      parameter (nz=101,npsi=40,nw=100,psimax=.4,zlen=20.,pi=3.1415926)
!      real phi(nz),z(nz)
      real psi(npsi),ptrap(npsi),fptrap(npsi),gint(npsi),hint(npsi)
      real gjint(npsi),pint(npsi)
      real v(nw),fschamel(nw)

      s2pi=sqrt(2.*3.1415926)
      write(*,*)'  psi    beta  flattrap trapnum  ft-tn    gint   pint'
      do j=npsi,1,-1
!      j=1
         psi(j)=psimax*j/npsi
         beta=-1.-(15./16.)*sqrt(pi/psi(j)) ! Shallow hole approx.
!      beta=beta*(.8-.9/abs(beta))    ! Ad hoc correction. Density
!      beta=beta*(1.-.75/abs(beta))    ! Ditto. Force.
         phiwidth=4.
         fptrap(j)=flattrap(psi(j),phiwidth)/s2pi
        call fhgfunc(psi(j),zlen,nz,4.,gint(j),hint(j),gjint(j),pint(j))
         
         wmax=2.*psi(j)
         do i=1,nw
            dw=2.*(psi(j))/nw
            w=-psi(j)+(i-1.)*dw
            v(i)=sqrt(2.*(psi(j)+w))
            fschamel(i)=ffuncw(w,beta)
         enddo
         if(j.eq.npsi)then
            call autoplot(v,fschamel,nw)
            call axlabels('v','fschamel')
         else
            call polyline(v,fschamel,nw)
         endif
         ptrap(j)=trapnum(psi(j),beta,phiwidth)/s2pi
         write(*,'(f8.4,f8.3,$)')psi(j),beta
         write(*,'(6f8.4)')fptrap(j),ptrap(j),fptrap(j)-ptrap(j),gint(j)    &
     &        ,pint(j)
      enddo
      call pltend()
      
      call autoplot(psi,ptrap,npsi)
      call axlabels('psi','trapnum')
      call pltend()

      



      end
