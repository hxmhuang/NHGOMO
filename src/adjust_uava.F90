subroutine adjust_uava()
  use openarray
  use variables
  use config
  implicit none
  type(array) zdist
  real(8) realdz
  realdz=sum(dz,3)
  ua=sum(u*dz*AXB(dt)*dum,3)/(realdz*AXB(dt))*dum
  va=sum(v*dz*AYB(dt)*dvm,3)/(realdz*AYB(dt))*dvm
end subroutine
