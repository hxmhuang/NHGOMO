#include "common.h"
subroutine bcond2_va()
  use openarray
  use variables
  use config
  implicit none
  type(array):: tmp,tmp1
  tmp=mat_zeros_im_jm_1
  if(rn==1) then
    tmp=sqrt(grav/h)*elf
    call set(sub(vaf, [2,imm1] ,jm), sub(tmp, [2,imm1] ,jm))
  else
    call set(sub(vaf, [2,imm1] ,jm), 0.d0)
  endif
  if(rs==1) then
    call set(sub(vaf, [2,imm1] ,2), -sqrt(grav/sub(h, [2,imm1] ,2))*sub(elf, [2,imm1] ,1))
    call set(sub(vaf, [2,imm1] ,1), sub(vaf, [2,imm1]   ,2))
  else
    call set(sub(vaf, [2,imm1] ,2), 0.d0)
  endif
  vaf = vaf * dvm 
end subroutine
