subroutine subinv(il,iu,pb,pd,pa,pc)
  use openarray
  implicit none
  integer :: k,lp
  integer, intent(in) :: il,iu
  type(array), intent(in) :: pb,pd,pa
  type(array), intent(inout) :: pc
  type(array) :: r 
  lp=il+1
  do k=lp,iu
    r=sub(pb,':',':',k)/sub(pd,':',':',k-1)
    call set(sub(pd,':',':',k),sub(pd,':',':',k)-r*sub(pa,':',':',k-1))
    call set(sub(pc,':',':',k),sub(pc,':',':',k)-r*sub(pc,':',':',k-1))
  enddo
  call set(sub(pc,':',':',iu),sub(pc,':',':',iu)/sub(pd,':',':',iu))
  do k=iu-1,il,-1
    call set(sub(pc,':',':',k),(sub(pc,':',':',k)-sub(pa,':',':',k)*  &
          sub(pc,':',':',k+1))/sub(pd,':',':',k))
  enddo
end subroutine
