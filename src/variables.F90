#include "common.h"
module variables
  use openarray
  type(array) :: dt, dtb,dtf
  integer rn,  rs
  type(array) aam2d, advua, advva, adx2d, ady2d
  type(array) art, aru, arv, cbc, cor, d, drx2d
  type(array) dry2d   , dum     , dvm
  type(array) dx      , dy      
  type(array) egb
  type(array) egf     , el      , elb     , elf
  type(array) et      , etb     , etf     , fluxua
  type(array) fluxva  , fsm     , h       
  type(array) psi
  type(array) ssurf   , swrad   , vfluxb
  type(array) tps     , tsurf   , ua      , vfluxf
  type(array) uab     , uaf     , utb     , utf
  type(array) va      , vab     , vaf
  type(array) vtb     , vtf     , wssurf  , wtsurf
  type(array) wubot   , wusurf  , wvbot   , wvsurf
  type(array) wubot1  , wvbot1
  type(array) aam  , advx , advy , advw , a
  type(array) c    , drhox, drhoy
  type(array) ee   , gg   , kh   , km
  type(array) kq   , l    , q2b  , q2, q2f, q2lf
  type(array) q2lb , q2l  , rho  , rmean
  type(array) sb   , sclim, s    , tb
  type(array) tclim, t    , ub   , uf
  type(array) u    , vb   , vf   , v    , q    , qb  , qf
  type(array) w    , wb   , wf   , ww   , wwb  , wwf , wq
  type(array) z, zz, dz, dzz
  type(array) curv, curv2d, curv_x, curv_y
  type(array) z_3d, dt_3d
  type(array) numnelt
  type(array)   :: mat_zeros, mat_ones, axf_u,ayf_v
  type(array)  ::  mat_zeros_im_jm_1 
  type(array)  ::  mat_zeros_im_1_1  
  type(array)  ::  mat_zeros_im_1_kb 
  type(array)  ::  mat_zeros_1_jm_kb 
  type(array)  ::  mat_zeros_1_jm_1  
  type(array)  ::  mat_zeros_1_1_kb 
  type(array)  ::  mat_ones_1_jm_kb
  type(array)  ::  mat_ones_im_1_kb
  type(array)  ::  swrad0
  real(kind=8) :: dti, dte2, dti2, ispi, isp2i
  real(kind=8) :: period, time,time_sec
  integer :: iprint, iend
  type(array) :: tf, sf, rad
  real(kind=8),allocatable :: z1(:),zz1(:),dzz1(:),dz1(:)
  real(kind=8),allocatable :: apr(:),rr(:),qq(:)
  integer,allocatable :: ia(:),ja(:),iaa(:)
  type(array) :: dxb_axf_dy, dyb_ayf_dx
  integer :: total_nums,numnz,cnz,cnn,nums_s
contains
  subroutine init_variables()
    use config
    implicit none
    character(len=1000) :: fnc
    call LoadConfig()
    mat_zeros         = zeros(im, jm, kb, sw=1, dt=2)
    mat_ones          = ones(im, jm, kb, sw=1, dt=2)
    mat_zeros_im_jm_1 = sub(mat_zeros, ':', ':', 1)
    mat_zeros_im_1_1  = sub(mat_zeros, ':', 1, 1)
    mat_zeros_im_1_kb = sub(mat_zeros, ':', 1, ':')
    mat_zeros_1_jm_kb = sub(mat_zeros, 1, ':', ':')
    mat_zeros_1_jm_1  = sub(mat_zeros, 1, ':', 1)
    mat_zeros_1_1_kb  = sub(mat_zeros, 1, 1,   ':')
    mat_ones_1_jm_kb  = sub(mat_ones,1,':',':')
    mat_ones_im_1_kb  = sub(mat_ones,':',1,':')
    numnelt= zeros(im,jm,kb,sw=1,dt=0) 
    aam=mat_zeros  ;advx=mat_zeros ;advy=mat_zeros ;a=mat_zeros    ;
    c=mat_zeros    ;drhox=mat_zeros;drhoy=mat_zeros;advw=mat_zeros ;
    ee=mat_zeros   ;gg=mat_zeros   ;kh=mat_zeros   ;km=mat_zeros   ;
    kq=mat_zeros   ;l=mat_zeros    ;q2b=mat_zeros  ;q2=mat_zeros   ;
    q2lb=mat_zeros ;q2l=mat_zeros  ;rho=mat_zeros  ;rmean=mat_zeros;
    sb=mat_zeros   ;sclim=mat_zeros;s=mat_zeros    ;tb=mat_zeros   ;
    tclim=mat_zeros;t=mat_zeros    ;ub=mat_zeros   ;uf=mat_zeros   ;
    u=mat_zeros    ;vb=mat_zeros   ;vf=mat_zeros   ;v=mat_zeros    ;
    w=mat_zeros    ;q2f=mat_zeros  ;q2lf=mat_zeros ;q=mat_zeros
    ww=mat_zeros   ;wwb=mat_zeros  ;wwf=mat_zeros  ;qb=mat_zeros   
    tf = mat_zeros; sf = mat_zeros ;wq=mat_zeros   ;qf=mat_zeros
    wb = mat_zeros; wf = mat_zeros;
    aam2d=mat_zeros_im_jm_1   ;advua=mat_zeros_im_jm_1   ;
    ady2d=mat_zeros_im_jm_1   ;art=mat_zeros_im_jm_1     ;
    cbc=mat_zeros_im_jm_1     ;
    cor=mat_zeros_im_jm_1     ;
    d=mat_zeros_im_jm_1       ;drx2d=mat_zeros_im_jm_1   ;
    dry2d=mat_zeros_im_jm_1   ;dt=mat_zeros_im_jm_1      ;
    dx=mat_zeros_im_jm_1      ;dy=mat_zeros_im_jm_1      ;
    egf=mat_zeros_im_jm_1     ;el=mat_zeros_im_jm_1      ;
    et=mat_zeros_im_jm_1      ;etb=mat_zeros_im_jm_1     ;
    fluxva=mat_zeros_im_jm_1  ;fsm=mat_zeros_im_jm_1     ;
    ssurf=mat_zeros_im_jm_1   ;
    tps=mat_zeros_im_jm_1     ;tsurf=mat_zeros_im_jm_1   ;
    uab=mat_zeros_im_jm_1     ;uaf=mat_zeros_im_jm_1     ;
    va=mat_zeros_im_jm_1      ;vab=mat_zeros_im_jm_1     ;
    vtb=mat_zeros_im_jm_1     ;vtf=mat_zeros_im_jm_1     ;
    wubot=mat_zeros_im_jm_1   ;wusurf=mat_zeros_im_jm_1  ;
    wubot1=mat_zeros_im_jm_1  ;wvbot1=mat_zeros_im_jm_1  ;
    advva=mat_zeros_im_jm_1   ;adx2d=mat_zeros_im_jm_1   ;
    aru=mat_zeros_im_jm_1     ;arv=mat_zeros_im_jm_1     ;
    dum=mat_zeros_im_jm_1     ;dvm=mat_zeros_im_jm_1     ;
    egb=mat_zeros_im_jm_1     ;
    elb=mat_zeros_im_jm_1     ;elf=mat_zeros_im_jm_1     ;
    etf=mat_zeros_im_jm_1     ;fluxua=mat_zeros_im_jm_1  ;
    h=mat_zeros_im_jm_1       ;psi=mat_zeros_im_jm_1     ;    
    swrad=mat_zeros_im_jm_1   ;vfluxb=mat_zeros_im_jm_1  ;
    swrad0 = mat_zeros_im_jm_1
    ua=mat_zeros_im_jm_1      ;vfluxf=mat_zeros_im_jm_1  ;
    utb=mat_zeros_im_jm_1     ;utf=mat_zeros_im_jm_1     ;
    vaf=mat_zeros_im_jm_1     ;                         
    wssurf=mat_zeros_im_jm_1  ;wtsurf=mat_zeros_im_jm_1  ;
    wvbot=mat_zeros_im_jm_1   ;wvsurf=mat_zeros_im_jm_1  ;
    d   = mat_zeros_im_jm_1;
    dtf = mat_zeros_im_jm_1;
    dtb = mat_zeros_im_jm_1;
    rad = mat_zeros           
    curv = mat_zeros
    curv_x = mat_zeros;  curv_y = mat_zeros
    curv2d = mat_zeros_im_jm_1
    fnc = trim(in_path)//trim(problem)//".nc"
    if(get_rank()==0) then
       print*, "trim(fnc)=", trim(fnc)
    endif
    call read_init(trim(fnc))
    dti=dte*isplit
    dte2=dte*2.d0
    dti2=dti*2.d0
    iend=int(tend/dti)
    iprint=int(prtd1*24.0*3600.0/dti+0.5)
    ispi=1.d0/isplit
    isp2i=1.d0/(2.d0*isplit)
  end subroutine
end module
