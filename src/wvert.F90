#include "common.h"
subroutine wvert()
  use openarray
  use config
  use variables
  implicit none
  integer :: k
  type(array) :: u111,u011,u211,pru,uar,fff,rk,rk1
  type(array) :: aaf,bbf,ccf
  aaf=mat_zeros;bbf=mat_zeros
  ccf=mat_zeros;fff=mat_zeros
  u111=mat_zeros;u011=mat_zeros
  u211=mat_zeros;pru=mat_zeros
  uar=dti2/(dtf*dtf)
  call coef2(aaf,bbf,ccf,fff)
  do k=2,kbm1
    rk=KM(k)/(DZZ(k-1)*DZ(k))
    rk1=KM(k+1)/(DZ(k)*DZZ(k))
    call set(sub(u111,':',':',k),1.d0+uar*(rk+rk1))
    call set(sub(u011,':',':',k),-uar*rk)
    call set(sub(u211,':',':',k),-uar*rk1)
    call set(sub(pru,':',':',k),sub(wwf,':',':',k))
  enddo
  call set(sub(u111,':',':',1),1.d0)
  call set(sub(u111,':',':',kb),1.d0)
  call set(sub(u211,':',':',1),0.d0)
  call set(sub(u011,':',':',kb),0.d0)
  call set(sub(pru,':',':',1),sub(fff,':',':',1))
  call set(sub(pru,':',':',kb),sub(fff,':',':',kb))
  call subinv(1,kb,u011,u111,u211,pru)
  wwf=pru
end subroutine
