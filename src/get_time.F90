#include "common.h"
subroutine get_time(iint)
  use openarray
  use variables
  use config
  implicit none
  integer :: iint
  time_sec=0
  time=dti*iint*1.0/86400.0+time0
  time_sec=time_sec+dti
end subroutine

