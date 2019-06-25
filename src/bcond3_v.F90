#include "common.h"
subroutine bcond3_v()
  use openarray
  use variables
  use config
  implicit none
  call set(sub(vf,':',jm,':'), 0.d0)
  call set(sub(vf,':',2,':'), 0.d0)
  call set(sub(vf,':',1,':'), 0.d0)
  call set(sub(vf, 1, ':',':'), 0.d0)
  call set(sub(vf, im,':',':'), 0.d0)
  call set(sub(vf,':',':',kb), 0.d0)
  vf = vf * dvm 
end subroutine
