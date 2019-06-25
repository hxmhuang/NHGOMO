#include "common.h"
#define STENCIL_BOX 1
program nh_gomo
  use openarray
  use config
  use variables
  implicit none
  integer :: ierr
  integer :: myrank, mysize
  real(kind=8) :: vtot, atot, taver, saver, eaver, tsalt
  integer :: iext, iint, imax, jmax
  real(kind=8) :: vamax
  integer :: max_step
  character(len=80) :: filename,string

  call oa_init(MPI_COMM_WORLD, [-1,-1,1])
  
  call set_stencil(STENCIL_BOX, 1)

  call oa_get_option(max_step, "s", -1)

  call init_variables()

  call grid_init('C', dx, dy, dz)

  call init_fields()
  
  call fsm_compute()
 
  call update_initial()

  call bottom_friction()

  if(get_rank()==0) then
    open(unit=60,file='conservation7200.txt')
    write(60,'(2A23,3A12,A31,A10)')'vtot','atot','eaver','taver','saver','tsalt','vamax'
  endif

!***********************************************************************
!                                                                      *
!                  begin numerical integration                         *
!                                                                      *
!***********************************************************************
  do iint = 1, iend
    if(get_rank()==0)then
    print*,"iint=",iint
    endif
    call get_time(iint)

    call surface_forcing()

    call adjust_uava()
    
    call lateral_viscosity()
    
    call mode_interaction()

    !********************** begin external mode *************************
    do iext = 1,isplit
    
      call external_el()

      call external_ua(iext)

      call external_va(iext)

      call external_update(iext, vamax, imax, jmax)
 
    enddo
    !*****************end external (2-d) mode calculation*****************
 
    !***********continue with internal (3-d) mode calculation*************
    if(vamax <= vmaxl) then

      if(((iint/=1).or.(time0/=0.e0)).and.(mode/=2)) then

        !----------------------------------------------------------------
        !               compute sf and tf
        !----------------------------------------------------------------
!        dtf=h+etf 
        call advtl(sf,s,sb)

        call advtl(tf,t,tb)
       
        call proft(sf,wssurf,swrad,ssurf,nbcs)
 
        call proft(tf,wtsurf,swrad,tsurf,nbct) 
        
        call bcond4(tf)

        call bcond4(sf)
 
        call smoth_update_ts(sf,s,sb)

        call smoth_update_ts(tf,t,tb)

        call dens(rho,t,s)

        call bcond4(rho)

        !----------------------------------------------------------------
        !               compute uf and vf
        !----------------------------------------------------------------
        
        call advu()
   
        call advv()

        call profu()

        call profv()

        call bcond3_u()

        call bcond3_v()

        !----------------------------------------------------------------
        !               compute wwf
        !----------------------------------------------------------------

        call advctw()

        call advew()

        call wvert()

        call bcond5(wwf)

        !----------------------------------------------------------------
        !               compute non-hydrastatic pressure q
        !----------------------------------------------------------------
       
        call construct_matrix(iint)

        call faz()

        call bcond3_u()

        call bcond3_v()

        call bcond5(w)

        call smoth_update_uvw()

      endif

      call internal_update()

    endif
    if(mod(iint,500)==0)then
      write(string,*) iint 
      filename="output/s_"//trim(adjustl(string))//".nc"
      call save(sf,filename,'salinity')
    endif
    call print_section(iint, iext, vamax,imax,jmax, &
           vtot, atot, taver, saver, eaver, tsalt)

    if(get_rank() == 0) then
      write(60,'(F25.4,F22.7,F12.7,F12.7,F16.10,F29.7,F10.5)') &
            vtot, atot, eaver,taver,saver ,tsalt,vamax
    endif

  enddo

  call oa_finalize()

end program nh_gomo
