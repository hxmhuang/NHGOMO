#include "common.h"
subroutine bcond6(ff)
  use openarray
  use variables
  use config
  implicit none
  type(array), intent(inout):: ff
  call set(sub(ff,[2,imm1],1,':'),  sub(ff,[2,imm1],2,':'))
  call set(sub(ff,[2,imm1],jm,':'), sub(ff,[2,imm1],jmm1,':'))
  ff  = ff  *fsm + 1.d-10
end subroutine
