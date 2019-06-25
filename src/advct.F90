subroutine advct()
	use openarray
	use variables
	use config
	implicit none
  type(array)::tmp
  curv = (AYF(v) * DXB(AXF(dy))  &
       - AXF(u) * DYB(AYF(dx))) / dx / dy
  tmp  = (AXF(AXB(dt) * u) * AXF(u) - dt*aam*2.d0*DXF(ub)) * fsm
  advx  = (DXB(tmp) &
       + DYF((AXB(AYB(dt) * v) * AYB(u) - AYB(AXB(dt)) &
       * AYB(AXB(aam))*(DYB(ub) + DXB(vb))) &
       * (2.d0*AYB(dum)*AYB(dum)-AYB(dum*dum))*(2.d0*AXB(dvm)*AXB(dvm)-AXB(dvm*dvm))) &
       - AXB(curv * dt * AYF(v)) ) * dz
  tmp = ( AYF(AYB(dt) * v) * AYF(v) - dt*aam*2.d0*DYF(vb) ) * fsm
  advy  = (DXF((AYB(AXB(dt) * u) * AXB(v) &
       - AYB(AXB(dt))*AYB(AXB(aam))*(DYB(ub) + DXB(vb))) &
       * (2.d0*AYB(dum)*AYB(dum)-AYB(dum*dum))*(2.d0*AXB(dvm)*AXB(dvm)-AXB(dvm*dvm))) &
       + DYB(tmp) + AYB(curv * dt * AXF(u))) * dz
end subroutine
