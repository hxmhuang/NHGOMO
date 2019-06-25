#include "common.h"
subroutine bcond1()
  use openarray
  use variables
  use config
  implicit none
  call set(sub(elf,  1 , ':'), sub(elf, 2   , ':'))
  call set(sub(elf,  im, ':'), sub(elf, imm1, ':'))
  elf=elf*fsm 
  if (rs==1) call set(sub(elf,  [2,imm1] , 1), 0.d0)
  if (rn==1) call set(sub(elf,  [2,imm1] , jm), 0.d0)
end subroutine
