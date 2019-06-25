subroutine bottom_friction()
  use openarray
  use variables
  use config
  implicit none
  real(kind=8) :: dz_kbm1
  dz_kbm1 = sub(dz, 1,1,kbm1)
  cbc=(0.4/log(0.5*dz_kbm1*dtb/z0b))**2;
  call set(cbc, cbcmin, cbc < cbcmin)
  call set(cbc, cbcmax, cbc > cbcmax)
end subroutine
