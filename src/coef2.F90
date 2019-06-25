#include "common.h"
subroutine coef2(aa1,bb1,cc1,ff)
  use openarray
  use variables
  use config
  implicit none
  integer k
  type(array),intent(inout)::aa1,bb1,cc1,ff  ! set aa1,bb1,cc1,ff as global var
  type(array) tmp
  do k=1,kb
    call set(sub(aa1,':',':',k), DXB(AXF(etf)) + DXB(AXF(dtf))*sub(z,1,1,k))
    call set(sub(bb1,':',':',k), DYB(AYF(etf)) + DYB(AYF(dtf))*sub(z,1,1,k))
    call set(sub(cc1,':',':',k), (etf-etb)/dti2 * (1.d0 + sub(z,1,1,k))  )
    call set(sub(aa1,1,':',k),   (sub(etf,2,':')-sub(etf,1,':') &
         + (sub(dtf,2,':')-sub(dtf,1,':'))*sub(z,1,1,k) )/sub(dx,1,':') )
    call set(sub(bb1,':',1,k),   (sub(etf,':',2)-sub(etf,':',1) &
         + (sub(dtf,':',2)-sub(dtf,':',1))*sub(z,1,1,k) )/sub(dy,':',1) )
    call set(sub(aa1,im,':',k),   (sub(etf,im,':')-sub(etf,imm1,':') &
         + (sub(dtf,im,':')-sub(dtf,imm1,':'))*sub(z,1,1,k) )/sub(dx,im,':') )
    call set(sub(bb1,':',jm,k),   (sub(etf,':',jm)-sub(etf,':',jmm1) &
         + (sub(dtf,':',jm)-sub(dtf,':',jmm1))*sub(z,1,1,k) )/sub(dy,':',jm) )
  enddo
  ff=AZB(AXF(uf))*aa1+AZB(AYF(vf))*bb1+cc1
  tmp=AXF(uf)*aa1+AYF(vf)*bb1+cc1
  call set(sub(ff,':',':',1), sub(tmp,':',':',1) )
end subroutine
