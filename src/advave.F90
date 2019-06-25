subroutine advave()
  use openarray
  use variables
  use config
  implicit none
  type(array) tmp
  tps = AYB(AXB(d)) * AXB(AYB(aam2d)) * (DYB(uab) + DXB(vab))
  tmp = AXF( AXB(d) * ua ) * AXF(ua) &
       - 2.d0 * d * aam2d * DXF(uab)*fsm
  advua = DXB(tmp)  &
       + DYF( (AXB(AYB(d) * va ) * AYB(ua)-tps) &
       * ( 2.d0*AYB(dum)*AYB(dum)- AYB(dum*dum) ) &
       * ( 2.d0*AXB(dvm)*AXB(dvm)- AXB(dvm*dvm) )  )
  tmp = ( AYF( AYB(d) * va ) * AYF(va) &
       - 2.d0 * d * aam2d * DYF(vab) ) * fsm
  advva = DXF( ( AYB( AXB(d) * ua ) * AXB(va) - tps ) &
       *( 2.d0*AYB(dum) * AYB(dum)-AYB(dum*dum) ) &
       *( 2.d0*AXB(dvm) * AXB(dvm)-AXB(dvm*dvm) ) ) &
       +DYB(tmp)
end subroutine
