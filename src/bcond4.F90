subroutine bcond4(ff)
  use openarray
  use config
  use variables
  implicit none
  type(array), intent(inout):: ff
  call set(sub(ff, 1,':',':'),  sub(ff, 2,':',':') )
  call set(sub(ff, ':',1,':'),  sub(ff, ':',2,':') )
  call set(sub(ff, im,':',':'), sub(ff, imm1,':',':') )
  call set(sub(ff, ':',jm,':'), sub(ff, ':',jmm1,':') )
  ff = ff * fsm
end subroutine
