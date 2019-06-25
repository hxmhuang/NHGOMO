module config
  implicit none
  integer :: im,jm,kb,kbm1,kbm2,imm1,jmm1,imm2,jmm2
  integer :: mode, isplit,iskp,jskp,ntp,nbct,nbcs,ispadv
  real(kind=8) :: dte,tend,time0,prtd1,prtd2,swtch,iswtch
  real(kind=8) :: tprni,corfac,alpha,tbias,sbias,rhoref
  real(kind=8) :: small,pi,grav,z0b,cbcmin,cbcmax,horcon,umol
  real(kind=8) :: smoth, aam_init, ramp, vmaxl
  character(len=60) :: in_path,problem
  namelist /setting/ im,jm,kb,dte,isplit,tend,mode, problem
  namelist /setting/ prtd1,prtd2,time0,iswtch
  namelist /setting/ swtch,iskp,jskp,tprni,ntp,corfac
  namelist /setting/ nbct,nbcs,ispadv,alpha,tbias,sbias
  namelist /setting/ small,pi,grav,z0b,rhoref
  namelist /setting/ cbcmin,cbcmax,horcon,umol
  namelist /setting/ smoth,aam_init,ramp,vmaxl,in_path
contains
 subroutine LoadConfig()
   open(1, file='config.txt')
   read(1, setting)
   imm1 = im - 1;  imm2 = im - 2
   jmm1 = jm - 1;  jmm2 = jm - 2
   kbm1 = kb - 1;  kbm2 = kb - 2
 end subroutine
end module config
