subroutine external_va(iext)
  use openarray
  use variables
  use config
  implicit none
  integer :: iext
  if(mod(iext,ispadv)==0) then
    advva = DXF( (AYB(AXB(d)*ua)*AXB(va) &
          -AYB(AXB(d))*AXB(AYB(aam2d))*(DYB(uab) + DXB(vab))) &
          *(2*AYB(dum )*AYB(dum )-AYB(dum *dum ))  &
          *(2*AXB(dvm )*AXB(dvm )-AXB(dvm *dvm ))) &
          +DYB( (AYF(AYB(d)*va)* AYF(va)-2.0*d * aam2d * DYF(vab))*fsm )
    call set(sub(advva, 1, ':',1), 0.d0)
    call set(sub(advva, im, ':',1), 0.d0)
    call set(sub(advva, ':', 1,1), 0.d0)
    call set(sub(advva, ':', jm,1), 0.d0)
  endif
  vaf= (AYB(h+elb)* vab &
       -2.0* dte*(ady2d + advva + AYB(cor* d* AXF(ua)) &
       + grav*AYB(d)*((1.0-2.0*alpha)* DYB(el) &
       + alpha*(DYB(elb)+ DYB(elf))) &
       + dry2d +(wvsurf-wvbot))) / AYB(h+elf)
  call bcond2_va()
end subroutine
