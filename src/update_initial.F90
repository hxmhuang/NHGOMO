#include "common.h"
subroutine update_initial()
  use openarray
  use variables
  use config
  implicit none
  integer::k
  ua=uab;va=vab
  el=elb;et=etb;etf=et
  d=h+el;dt=h+et
  dtb=h+etb;dtf=h+etf
  q2b = 1.d-8;q2lb = 1.d-8
  kh = 1.d-7;km = 2.d-6
  l=q2lb/(q2b+small)
  kq = 0.2d0*l*sqrt(q2b)
  aam = aam_init
  q2 = q2b;q2l = q2lb
  t = tb;s = sb
  u = ub;v = vb
  call dens(rho, t, s)
  call baropg()
  drx2d = sum(drhox, 3)
  dry2d = sum(drhoy, 3)
end subroutine
