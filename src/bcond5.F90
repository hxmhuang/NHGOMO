subroutine bcond5(ff)
  use openarray
  use config
  use variables
  implicit none
  type(array), intent(inout):: ff
  ff=ff*fsm
  call set(sub(ff, ':',1,':'),    sub(ff,   ':',2,':') )
  call set(sub(ff, ':',jm,':'),   sub(ff,   ':',jmm1,':') )
  call set(sub(ff,   im,':',':'),  0.d0 )
  call set(sub(ff,  ':',':',kb), 0.d0)
end subroutine
