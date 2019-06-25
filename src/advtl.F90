subroutine advtl(ff,f,fb)
	use openarray
	use variables
	use config
	implicit none
  integer i,j,k
  type(array),intent(inout)::ff,f,fb
  type(array)::xflux,yflux,zflux,zflux1,fdtd1,fdtd2
  type(array)::delfi,delfiavg,fimin,fimax,test1,test2
  type(array)::tmp,tmp1,fdt,tmpu,tmpflux
  call set(sub(f,':',':',kb),  sub(f,':',':',kbm1) )
  call set(sub(fb,':',':',kb), sub(fb,':',':',kbm1) )
	fdt=f*dt
  tmpu=mat_zeros;tmpflux=mat_zeros
  xflux=mat_zeros;yflux=mat_zeros
  zflux=mat_zeros;zflux1=mat_zeros
!---------------- 	ADVECTIVE FLUXES -----------
!  ----------x-direction:-----------
	fimin=fdt;fimax=fdt
  fdtd1=mat_zeros;fdtd2=mat_zeros
  call set(sub(fdtd2,[2,imm1],':',':'), sub(fdt,[3,im],  ':',':') )
  call set(sub(fdtd1,[2,imm1],':',':'), sub(fdt,[1,imm2],':',':') )
	delfiavg=0.5*(fdtd2 - fdtd1)
	call set(fimin, fdtd1, fdtd1<fimin)
	call set(fimin, fdtd2, fdtd2<fimin)
	call set(fimax, fdtd1, fdtd1>fimax)
	call set(fimax, fdtd2, fdtd2>fimax)
  tmp=mat_zeros;tmp1=mat_zeros
	call set(tmp,2*(fdt-fimin),fdt>fimin)
	call set(tmp,abs(delfiavg),abs(delfiavg)<tmp)
	call set(tmp1,2*(fimax-fdt),fimax>fdt)
	call set(tmp,tmp1,tmp>tmp1)
	delfi=1.d0*(delfiavg>=0.d0)*tmp-1.d0*(delfiavg<0.d0)*tmp
	call set(sub(delfi,1,':',':'),  sub(delfi,2,':',':') )
	call set(sub(delfi,im,':',':'), sub(delfi,imm1,':',':') )
	tmp=1.d0*(u>=0.d0)
  call set(sub(tmpu,[1,imm1],':',':'), sub(u,[2,im],':',':') )
  tmp1 = 0.5d0* delfi * (1.d0- dti* tmpu/dx  )  
  tmpflux = tmp1 + 1.d0*fdt 
  call set( sub(xflux,[2,im],':',':'), sub(u,[2,im],':',':') * &
            sub(tmp,[2,im],':',':') * sub(tmpflux,[1,imm1],':',':') )
	xflux=xflux+u*(fdt-0.5*delfi*(1+u*dti/dx))*(1-tmp)
!  ----------y-direction:-----------
	fimin=fdt;fimax=fdt
	fdtd1=shift(fdt, 0, 1,0);fdtd2=shift(fdt, 0,-1,0)
	delfiavg=0.5*(fdtd2 - fdtd1)
	call set(fimin, fdtd1, fdtd1<fimin)
	call set(fimin, fdtd2, fdtd2<fimin)
	call set(fimax, fdtd1, fdtd1>fimax)
	call set(fimax, fdtd2, fdtd2>fimax)
  tmp=mat_zeros;tmp1=mat_zeros
	call set(tmp,2*(fdt-fimin),fdt>fimin)
	call set(tmp,abs(delfiavg),abs(delfiavg)<tmp)
	call set(tmp1,2*(fimax-fdt),fimax>fdt)
	call set(tmp,tmp1,tmp>tmp1)
	delfi=1.d0*(delfiavg>=0.d0)*tmp-1.d0*(delfiavg<0.d0)*tmp
	call set(sub(delfi,':',1,':'),  sub(delfi,':',2,':') )
	call set(sub(delfi,':',jm,':'), sub(delfi,':',jmm1,':') )
	tmp= 1.d0*(v>=0.d0)
  tmpu=mat_zeros
  call set(sub(tmpu,':',[1,jmm1],':'), sub(v,':',[2,jm],':') )
  tmpflux = fdt + 0.5* delfi * (1- tmpu * dti /dy  )  
  call set( sub(yflux,':',[2,jm],':'), sub(v,':',[2,jm],':') * &
            sub(tmp,  ':',[2,jm],':')* sub(tmpflux,':',[1,jmm1],':') )
	yflux=yflux+v*(fdt-0.5*delfi*(1+v*dti/dy))*(1-tmp)
!--------------ADD HORIZONTAL DifFUSION--------
	xflux=xflux-AXB(aam)*AXB(h)*DXB(fb)*dum
	yflux=yflux-AYB(aam)*AYB(h)*DYB(fb)*dvm
!  ----------vertical advection:----
	fimin=f;fimax=f
  call set(sub(fdtd2,':',':',[2,kbm1]), sub(f, ':',':',[3,kb]) )
  call set(sub(fdtd1,':',':',[2,kbm1]), sub(f, ':',':',[1,kbm2]) )
	delfiavg=0.5d0*(fdtd2 - fdtd1)
	call set(fimin, fdtd1, fdtd1<fimin)
	call set(fimin, fdtd2, fdtd2<fimin)
	call set(fimax, fdtd1, fdtd1>fimax)
	call set(fimax, fdtd2, fdtd2>fimax)
  tmp=mat_zeros;tmp1=mat_zeros
	call set(tmp,2.d0*(f-fimin),f>fimin)
	call set(tmp,abs(delfiavg),abs(delfiavg)<tmp)
	call set(tmp1,2.d0*(fimax-f),fimax>f)
	call set(tmp,tmp1,tmp>tmp1)
	delfi=1.d0*(delfiavg>=0.d0)*tmp-1.d0*(delfiavg<0.d0)*tmp
	call set(sub(delfi,':',':',1),  sub(delfi,':',':',2) )
	call set(sub(delfi,':',':',kb), sub(delfi,':',':',kbm1) )
  tmp= 1.d0*(w<0.d0)
  tmpu=mat_zeros
  tmpu=shift( w,0,0,-1)
  tmpflux = -0.5d0* (1.d0 + tmpu/dt * dti /dz )* delfi + f
  tmpu=shift(tmpflux, 0, 0,1)
  zflux=w*tmp*tmpu
  zflux1= w*(f+0.5d0*delfi* (1.d0 -w/dt/dz*dti) )*(1-tmp)+zflux
  call set(sub(zflux,':',':',[1,kbm1]),  sub(zflux1,':',':',[1,kbm1]) )
  call set(sub(zflux,':',':',1), sub(zflux,':',':',1)*(1.d0-sub(tmp,':',':',1)))
  call set(sub(zflux,':',':',kb),sub(zflux,':',':',kb)*sub(tmp,':',':',kb))
  ff=(fdt-dti*(DXF(xflux)+DYF(yflux)-1.d0*DZF(zflux)))/dtf
end subroutine
