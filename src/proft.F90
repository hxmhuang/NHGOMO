#include "common.h"
subroutine proft(f,wfsurf,sward,fsurf,nbc)
  use openarray
  use variables
  use config
  implicit none
  integer,intent(in) :: nbc
  integer :: k
  type(array), intent(inout) :: f
  type(array), intent(in) :: wfsurf,sward,fsurf
  type(array):: tmp
  real*8 tr(5),extc(5)
  tr=(/0.32d0, 0.31d0, 0.29d0, 0.26d0, 0.24d0/)
  extc=(/0.037d0,0.042d0,0.056d0,0.073d0,0.127d0/)
  a=mat_zeros;c=mat_zeros
  tmp=mat_ones*dzz
  call set(sub(a,':',':',[1,kbm2]), -dti2*(sub(kh,':',':',[2,kbm1])+umol))
  a=a/(dz*dzz*dtf*dtf)
  call set(A(kbm1), 0.d0);call set(A(kb), 0.d0)
  call set(sub(c,':',':',[2,kbm1]), sub(tmp,':',':',[1,kbm2]))
  c=(-dti2*(kh+umol))/(dz*c*dtf*dtf)
  call set(C(1), 0.d0);call set(C(kb), 0.d0)
  if(nbc==2 .or. nbc==4) then
    rad=mat_ones*sward*(tr(ntp)*exp(extc(ntp)*z*dtf))
    call set(sub(rad,':',':',kb) , 0.d0)
  endif
  if(nbc==1) then
    call set(EE(1),fsm*A(1)/(A(1)-1.d0))
    call set(GG(1),fsm*(-dti2*wfsurf/(-DZ(1)*dtf)-F(1))/(A(1)-1.d0))
  elseif(nbc==2) then
    call set(EE(1),fsm*A(1)/(A(1)-1.d0))
    call set(GG(1),fsm*(dti2*(wfsurf+RAD(1)-RAD(2))/(DZ(1)*dtf)-F(1))/(A(1)-1.d0))
  else
    call set(EE(1),0.d0)
    call set(GG(1),fsurf)
  endif
  do k=2,kbm2
    call set(GG(k), 1.d0 / (A(k)+C(k)*(1.d0-EE(k-1)) -1.d0))
    call set(EE(k), A(k) * GG(k))
    call set(GG(k),(C(k)*GG(k-1)-F(k)+dti2*(RAD(k)-RAD(k+1))/(dtf*DZ(k)))*GG(k))
  enddo
    call set(F(kbm1), (C(kbm1)*GG(kbm2)-F(kbm1)+dti2*(RAD(kbm1)-RAD(kb)) &
          /(dtf*DZ(kbm1))) / (C(kbm1)*(1.d0-EE(kbm2))-1.d0))
  do k=kbm2,1,-1
    call set(F(k),EE(k)*F(k+1)+GG(k))
  enddo
end subroutine
