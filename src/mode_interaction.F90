subroutine mode_interaction()
  use openarray
  use config
  use variables
  implicit none
  adx2d = sum(advx,  3)
  ady2d = sum(advy,  3)
  drx2d = sum(drhox, 3)
  dry2d = sum(drhoy, 3)
  aam2d = sum(aam*dz*fsm,3)
  call advave()
  adx2d = adx2d-advua
  ady2d = ady2d-advva
  egf=el * ispi
  utf=ua * 2.0 * AXB(d) * isp2i
  vtf=va * 2.0 * AYB(d) * isp2i
end subroutine
