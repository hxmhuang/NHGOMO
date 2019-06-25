subroutine smoth_update_ts(ff,f,fb)
  use openarray
  use variables
  use config
  implicit none
  type(array), intent(inout) :: ff, f, fb
  fb = f+0.5d0*smoth*(dtf*ff+dtb*fb-2*dt*f)/dt
  f = ff
end subroutine
subroutine smoth_update_uvw()
  use config
  use variables
  implicit none
  u=u+0.5d0*smoth*(uf+ub-2.d0*u-sum((uf+ub-2.d0*u)*dz*dum, 3))
  call set(sub(u,':',':',kb),0.d0)
  v=v+0.5d0*smoth*(vf+vb-2.d0*v-sum((vf+vb-2.d0*v)*dz*dvm, 3))
  call set(sub(v,':',':',kb),0.d0)
  ub = u; u = uf; vb = v; v = vf
  wwb= ww; ww=wwf; wb=w; w=wf
end subroutine
