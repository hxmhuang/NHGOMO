#include "common.h"
subroutine advu()
	use openarray
	use variables
	use config
	implicit none
  type(array)::bondu, tmp
  tmp = AXB(w) * AZB(u)
  call set(sub(tmp, ':',':',kb), 0.d0)
	uf= (AXB(dtb)*ub -dti2*( (advx + drhox)/dz - AXB( cor*dt*AYF(v) ) &
        +grav*AXB(dt)* ( DXB(egf+egb) )/2.d0 - &
        DZF( tmp )))/ AXB(dtf)
  call set(UF(kb),0.d0)
  bondu = AXB(w) * AZB(u)
  call set(sub(uf,im,':',':') , sub(bondu,im,':',':'))
end subroutine
