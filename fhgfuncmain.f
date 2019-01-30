!*********************************************************************
!Main program for testing fhgfunc routines.
      include 'fhgfunc.f'

      program fhgvfuncmain
      parameter (nz=101,npsi=80,nw=100,psimax=1.,zlen=20.,pi=3.1415926)
!      real phi(nz),z(nz)
      real psi(npsi),ptrap(npsi),fptrap(npsi),gint(npsi),hint(npsi)
      real gjint(npsi),pint(npsi)
      real v(nw),fschamel(nw)

      call pfset(3)
      s2pi=sqrt(2.*3.1415926)
      write(*,*)'  psi    beta  flattrap trapnum  ft-tn    gint   pint'
     $     ,'     hint   gjint'
      
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
         write(*,'(8f8.4)')fptrap(j),ptrap(j),fptrap(j)-ptrap(j),gint(j)    &
     &        ,pint(j),hint(j),gjint(j)
      enddo
      call pltend()
      
      call autoplot(psi,ptrap,npsi)
      call axlabels('psi','trapnum')
      call pltend()

      call pltinit(0.,psi(npsi),0.,float(nint(gjint(npsi)+0.5)))
      call charsize(.018,.018)
      call polyline(psi,gjint,npsi)
      call polyline(psi,hint,npsi)
      call axis
      call axis2
      call axlabels('!Ay!@','!BH=!AJ!Bhdx ,     J=!AJ!Bjdx')
      call jdrwstr(wx2nx(psi(npsi/2)),wy2ny(hint(npsi/2)),'!BH!@',-4.)
      call jdrwstr(wx2nx(psi(npsi/2)),wy2ny(gjint(npsi/2)),'!BJ!@',-3.)
      call pltend()

      
      write(*,*)'Integrated h            g              gj            p'
     $     ,' at psi=',psimax
      write(*,*)hint(npsi),gint(npsi),gjint(npsi),pint(npsi)

      do i=1,npsi
         ptrap(i)=sqrt(psi(i)/hint(i)*16./3.)
         gint(i)=16.*psi(i)/3.
      enddo
      
      call autoplot(psi,gint,npsi)
      call axlabels('!Ay!@','H        16!Ay!@/3')
      call polyline(psi,hint,npsi)
      call pltend

      call autoplot(psi,ptrap,npsi)
      call axis2
      call axlabels('!Ay!@','(16!Ay!@/3H)!u1/2!u')
      call pltend

      end
