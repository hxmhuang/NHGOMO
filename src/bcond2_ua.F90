#include "common.h"
subroutine bcond2_ua()
  use openarray
  use variables
  use config
  implicit none
  call set(sub(uaf, 2   ,[2,jmm1]), sub(tps, 2   ,[2,jmm1]))
  call set(sub(uaf, 1   ,[2,jmm1]), sub(tps, 2   ,[2,jmm1]))
  if(rn==1) then
    call set(sub(uaf, [2,imm1] ,jm), 0.d0)
  else
    call set(sub(uaf, [2,imm1] ,jm), sub(uaf, [2,imm1] ,jmm1))
  endif
  if(rs==1) then
    call set(sub(uaf, [2,imm1] ,1), 0.d0)
  else
    call set(sub(uaf, [2,imm1] ,1), sub(uaf, [2,imm1] , 2))
  endif
  uaf = uaf * dum
end subroutine
