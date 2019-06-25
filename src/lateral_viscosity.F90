subroutine lateral_viscosity()
  use openarray
  use config
  use variables
  implicit none
  call advct()
  call baropg()
end subroutine
