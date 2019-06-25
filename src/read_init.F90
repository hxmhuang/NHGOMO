subroutine read_init(fnc)
  use openarray
  use variables
  use config
  use openarray
  use mpi
  implicit none
  character(len=*) :: fnc
  integer :: ii, ierr,i,j,k
  real*8, allocatable :: tmp_z1(:,:,:), &
       tmp_zz1(:,:,:), tmp_dz1(:,:,:), &
       tmp_dzz1(:,:,:)
  if(get_rank() == 0) &
       print*, "start reading init variables..."
  z       = load(fnc,'z');
  zz      = load(fnc,'zz');
  dz      = load(fnc,'dz');
  dzz     = load(fnc,'dzz');
  allocate(z1(kb), zz1(kb), dz1(kb), dzz1(kb))
  allocate(tmp_z1(1,1,kb), tmp_zz1(1,1,kb), &
       tmp_dz1(1,1,kb), tmp_dzz1(1,1,kb))
  tmp_z1   = z;tmp_zz1  = zz
  tmp_dz1  = dz;tmp_dzz1 = dzz
  z1   = tmp_z1(1,1,:)
  zz1  = tmp_zz1(1,1,:)
  dz1  = tmp_dz1(1,1,:)
  dzz1 = tmp_dzz1(1,1,:)
  dx      = load(fnc,'dx')       ; !call disp(dx, 'dx = ')
  dy      = load(fnc,'dy')       ; !call disp(dy, 'dy = ')
  cor     = load(fnc,'cor')      ; !call disp(cor, 'cor = ')
  h       = load(fnc,'h')        ;
  fsm     = load(fnc,'fsm')      ;
  dum     = load(fnc,'dum')      ;
  dvm     = load(fnc,'dvm')      ;
  art     = load(fnc,'art')      ;
  aru     = load(fnc,'aru')      ; !call disp(aru, 'aru = ')
  arv     = load(fnc,'arv')      ; !call disp(arv, 'arv = ')
  rn     = load(fnc,'rn')      ;!call disp(rfn, 'rfn = ')
  rs     = load(fnc,'rs')      ;!call disp(rfs, 'rfs = ')   
  tb      = load(fnc,'tb')       ;
  sb      = load(fnc,'sb')       ;
  tclim   = load(fnc,'tclim')    ;
  sclim   = load(fnc,'sclim')    ;
  ub      = load(fnc,'ub');
  vb      = load(fnc,'vb');
  uab     = load(fnc,'uab');
  vab     = load(fnc,'vab');
  elb     = load(fnc,'elb');
  etb     = load(fnc,'etb');
  ssurf   = load(fnc,'ssurf')    ;
  tsurf   = load(fnc,'tsurf')    ;
  rho     = load(fnc,'rho')    ;
  rmean   = load(fnc,'rmean')    ;
end subroutine
