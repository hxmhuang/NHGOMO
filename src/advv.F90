#include "common.h"
subroutine advv()
	use openarray
	use variables
	use config
	implicit none
  type(array)::tmp
  tmp = AYB(w) * AZB(v)
  call set(sub(tmp, ':',':',kb), 0.d0)
  vf =(AYB(dtb)*vb-dti2*( (advy+ drhoy)/dz +AYB(cor*dt*AXF(u))+     &
      grav * AYB(dt)*(DYB(egf+egb) )/2.d0 &
      -DZF( tmp )))/AYB(dtf)
  call set(VF(kb), 0.d0)
end subroutine
