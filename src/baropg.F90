#include "common.h"
! Calculates  baroclinic pressure gradient.
subroutine baropg()
  use openarray
  use variables
  use config
  implicit none
  type(array) tmp, tmp1, ttt, axbdt
  tmp=-DZB(zz)* DXB(AZB(rho - rmean)) * AXB(dtb)   &
       + DXB(dtb) * DZB(AXB(rho - rmean)) &
       * AZB(zz) * AZB(dz)! , 3);
  tmp1 = -zz * AXB(dtb) * DXB(rho-rmean)
  call set(sub(tmp, ':',':',1), sub(tmp1, ':',':',1))
  call set(sub(tmp,':',':',kb),0.d0)
  ttt = csum(tmp, 3)
  drhox=ramp * grav * AXB(dt)  * csum(tmp, 3)*dz * dum;
  call set(sub(drhox, ':', ':', kb), 0.d0)
  tmp=-DZB(zz)* DYB(AZB(rho - rmean)) * AYB(dtb)  &
       + DYB(dtb) * DZB(AYB(rho - rmean))  &
       * AZB(zz) * AZB(dz)
  tmp1 = -zz * AYB(dtb) * DYB(rho-rmean)
  call set(sub(tmp, ':',':',1), sub(tmp1, ':',':',1))
  call set(sub(tmp,':',':',kb),0.d0)
  drhoy=ramp * grav * AYB(dt)  * csum(tmp, 3)*dz * dvm;
  call set(sub(drhoy, ':', ':', kb), 0.d0)
end subroutine
