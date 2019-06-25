#include "common.h"
subroutine profu()
  use openarray
  use config
  use variables
  implicit none
  integer :: k
  type(array) :: tmpdzz
  uf=uf*dum
  c=AXB(km)
  call set(sub(a,':',':',[1,kbm2]),-dti2*( sub(c,':',':',[2,kbm1])+ umol))
  a=a/(dz*dzz*AXB(dtf)*AXB(dtf))
  call set(A(kbm1), 0.d0);call set(A(kb), 0.d0)
  tmpdzz = dzz
  call set(sub(tmpdzz,':',':',[2,kbm1]),sub(dzz,':',':',[1,kbm2]))
  c=-dti2*(c+umol)/(dz*tmpdzz*AXB(dtf)*AXB(dtf))
  call set(C(1) ,0);call set(C(kb) ,0)
  call set(EE(1),A(1)/(A(1)-1.d0) )
  call set(GG(1),(-dti2*wusurf/(-DZ(1)*AXB(dtf))-UF(1))/(A(1)-1.d0))
  do k=2,kbm2
    call set(GG(k), 1.d0 / (A(k) + C(k) * (1.d0 - EE(k-1)) -1.d0))
    call set(EE(k), A(k) * GG(k) * dum)
    call set(GG(k), (C(k) * GG(k-1) - UF(k))* GG(k)*dum)
  enddo
  tps=2*AAM(kbm1)/(dt*DZZ(kbm1)) !need to be compared by wmq
  call set(UF(kbm1) , (C(kbm1)*GG(kbm2)-UF(kbm1))/ &
      (tps*dti2 /(-DZ(kbm1)*AXB(dtf))-&
      1.d0-C(kbm1)*( EE(kbm2)-1.d0)))
  do k=kbm2,1,-1
     call set(UF(k) , (EE(k)*UF(k+1)+GG(k))*dum)
  end do
  wubot=-tps*UF(kbm1)
end subroutine
