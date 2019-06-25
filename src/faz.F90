subroutine faz()
  use openarray
  use variables
  use config
  implicit none
  type(array) qvy,qux,quz,qvz,qflwz,qflux,qfluz,qflvy,qflvz
  type(array) ungf,vngf,wng,aaf,bbf,ccf,fff
  aaf=mat_zeros;bbf=mat_zeros
  ccf=mat_zeros;fff=mat_zeros
  call coef2(aaf,bbf,ccf,fff)
  qvy=q*dtf;qux=qvy
  quz = AZB(AXB(q))*AXB(aaf)
  qvz = AZB(AYB(q))*AYB(bbf)
  qflwz = -DZB(q)
  qflux = DXB(qux)
  qfluz = -DZF(quz) 
  qflvy = DYB(qvy)
  qflvz = -DZF(qvz)
  call set(sub(qflwz,':',':',kb),0.d0)
  call set(sub(qflux,2,':',':'), 0.d0)
  call set(sub(qfluz,2,':',':'), 0.d0)
  ungf = (qflux - qfluz)/AXB(dtf)*dti2*(-1.d0) 
  vngf = (qflvy - qflvz)/AYB(dtf)*dti2*(-1.d0) 
  wng = -qflwz/dtf*ramp*dti2
  wq = (wq+ wng)*fsm
  wf = (wq -ccf)*fsm
  uf = (uf + ungf)*dum
  vf = vf + vngf
  vf = vf*dvm
  wwf = (wq + fff - ccf)*fsm
 end subroutine
