#include "common.h"
subroutine bcond3_u()
  use openarray
  use variables
  use config
  implicit none
  call set(sub(uf,':',jm,':'), sub(uf,':',jmm1,':') )
  call set(sub(uf,':',1,':'), sub(uf,':',2,':') )
  call set(sub(uf,':',':',kb), 0.d0)
  call set(sub(uf,':',[2,jmm1],':'), sub(uf,':',[2,jmm1],':')*sub(dum,':',[2,jmm1],':') )
end subroutine
