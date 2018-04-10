program tbedoc
  ! Document the bounce time and frequency and the presence of a force
  ! resonance.

  ! We find that the omega_b=2pi/tbe agrees with sqrt(psi)/2.
  ! Also the force resonance occurs at an energy for which the
  ! omega_b = Omegac
  use shiftmode
  implicit none
  complex :: Ftotal
  real :: dfperpdWperp=1.,fperp=1.

  psi=.1

  omega=.01*sqm1
  Omegac=0.1
  omegad=Omegac+omega

!  omegad=omega
  call initialize
  idebug=-2
  call FtEint(Ftotal,dfperpdWperp,fperp)

  write(*,*)'omega_b= 2\pi/tbe'
  write(*,'(10f8.4)')2*3.1415926/tbe
  write(*,*)'sqrt(psi)/2=',sqrt(psi)/2.
  
end program tbedoc
