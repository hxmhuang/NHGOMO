#include "common.h"
subroutine construct_matrix(iint)
  use openarray
  use variables
  use config
  use pvariable
  use petscsys
  implicit none

  integer :: i,j,i_global,j_global,k,m,n,ierr
  integer :: iint
  type(array) :: aa1,aa2,aa3,aa4,bb1,bb2,bb3,bb4,&
    ga1,ga2,gb1,gb2,gc1,gc2,gen
  type(array) :: r,r1,r2
  type(array) :: aaf,bbf,ccf,fff
  type(array) :: ddx,ddy,dq
  type(array) :: fsm_tmp
  integer :: ind(3),ind_global(6),xs,ys
  real(kind=8) :: tmp0,tmp1,tmp2,tmp3,tmp4
  real(kind=8),pointer :: paa1(:,:,:),paa2(:,:,:),paa3(:,:,:),&
                    paa4(:,:,:),pbb1(:,:,:),pbb2(:,:,:),&
                    pbb3(:,:,:),pbb4(:,:,:),pga1(:,:,:),&
                    pga2(:,:,:),pgb1(:,:,:),pgb2(:,:,:),&
                    pgc1(:,:,:),pgc2(:,:,:),pgen(:,:,:),&
                    paaf(:,:,:),pbbf(:,:,:),pdx(:,:,:), &
                    pddx(:,:,:),pddy(:,:,:),pdy(:,:,:), &
                    pdq(:,:,:) ,pfsm(:,:,:),pr(:,:,:),pq(:,:,:)
  integer,pointer:: pnumnelt(:,:,:)
  aa1=mat_zeros;gen=mat_zeros;aa2=mat_zeros
  aa3=mat_zeros;aa4=mat_zeros;bb1=mat_zeros
  bb2=mat_zeros;bb3=mat_zeros;bb4=mat_zeros
  ga1=mat_zeros;ga2=mat_zeros;gb1=mat_zeros
  gc1=mat_zeros;gc2=mat_zeros;aaf=mat_zeros
  gb2=mat_zeros;gc1=mat_zeros;gen=mat_zeros
  bbf=mat_zeros;ccf=mat_zeros;fff=mat_zeros
  tol=0.05d0 
  ddx=AXF(dx); ddy=AXF(dy)
  fsm_tmp=fsm; dq=dtf
  call get_local_buffer(paa1,aa1);call get_local_buffer(paa2,aa2)
  call get_local_buffer(paa3,aa3);call get_local_buffer(paa4,aa4)
  call get_local_buffer(pbb1,bb1);call get_local_buffer(pbb2,bb2)
  call get_local_buffer(pbb3,bb3);call get_local_buffer(pbb4,bb4)
  call get_local_buffer(pga1,ga1);call get_local_buffer(pga2,ga2)
  call get_local_buffer(pgb1,gb1);call get_local_buffer(pgb2,gb2)
  call get_local_buffer(pgc1,gc1);call get_local_buffer(pgc2,gc2)
  call get_local_buffer(pgen,gen);call get_local_buffer(paaf,aaf)
  call get_local_buffer(pbbf,bbf);call get_local_buffer(pddx,ddx)
  call get_local_buffer(pddy,ddy);call get_local_buffer(pdq,dq)
  call get_local_buffer(pfsm,fsm);call get_local_buffer(pdx,dx)
  call get_local_buffer(pdy,dy)  ;call get_local_buffer(pnumnelt,numnelt)
  
  call update_ghost(ddx)
  call update_ghost(ddy)
  call set(sub(fsm,1,':'),0.d0)
  call set(sub(fsm,im,':'),0.d0)
  call set(sub(fsm,':',1),0.d0)
  call set(sub(fsm,':',jm),0.d0)
  call set(sub(dq,1,':'),sub(dtf,2,':'))
  call set(sub(dq,im,':'),sub(dtf,imm1,':'))
  call set(sub(dq,':',1),sub(dtf,':',2))
  call set(sub(dq,':',jm),sub(dtf,':',jmm1))
  
  call update_ghost(dq)
  call update_ghost(fsm)
  call coef1(aaf,bbf,ccf,fff)       !weitiao
  call update_ghost(aaf)
  call update_ghost(bbf)
  call update_ghost(numnelt)
  ia=0;ja=0;apr=0
!  call disp(numnelt,"numnelt=")

!petsc init--------------------------------------------------------
  PETSC_COMM_WORLD=MPI_COMM_WORLD
  p_N=total_nums
  if(iint==2)then
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  call VecCreate(MPI_COMM_WORLD,p_b,ierr)
  call VecSetSizes(p_b,cnn,p_N,ierr)
  call VecSetFromOptions(p_b,ierr)
  call VecDuplicate(p_b,p_x,ierr)
  call KSPCreate(MPI_COMM_WORLD,p_ksp,ierr)
  call KSPSetOperators(p_ksp,p_A,p_A,ierr)
  allocate(ix(0:cnn-1))
  do p_i=0,cnn-1
     ix(p_i)=p_i-1+nums_s
  enddo
  endif
!------------------------------------------------------------------
  ind = shape(pfsm)
  ind_global = get_box_corners(fsm)
  xs=ind_global(1);ys=ind_global(3)
  m=0;n=0
  do k=2,kb
    do j=2,ind(2)-1
      j_global=j+ys
      do i=2,ind(1)-1
        i_global=i+xs
        if(pfsm(i,j,2)>0.9d0)then
          n=n+1
          ia(n)=m
          paa1(i+1,j,k+1)=0.25d0*paaf(i+1,j,k+1)/pdx(i,j,2)
          paa2(i-1,j,k+1)=-0.25d0*paaf(i-1,j,k+1)/pdx(i,j,2)
          pbb1(i,j+1,k+1)=0.25d0*pbbf(i,j+1,k+1)/pdy(i,j,2)
          pbb2(i,j-1,k+1)=-0.25d0*pbbf(i,j-1,k+1)/pdy(i,j,2)
          pga1(i+1,j,k)=dz1(k-1)/pdx(i,j,2)/pddx(i,j,2)*pdq(i+1,j,2)
          pga2(i-1,j,k)=dz1(k-1)/pdx(i,j,2)/pddx(i-1,j,2)*pdq(i-1,j,2)
          pgb1(i,j+1,k)=dz1(k-1)/pdy(i,j,2)/pddy(i,j,2)*pdq(i,j+1,2)
          pgb2(i,j-1,k)=dz1(k-1)/pdy(i,j,2)/pddy(i,j-1,2)*pdq(i,j-1,2)
          pgc1(i,j,k+1)=1.e0/(dzz1(k-1)*pdq(i,j,2))*ramp     !weitiao
          if(k>2)then
            paa3(i+1,j,k-1)=-0.25d0*paaf(i+1,j,k-1)/pdx(i,j,2)
            paa4(i-1,j,k-1)=0.25d0*paaf(i-1,j,k-1)/pdx(i,j,2)
            pbb3(i,j+1,k-1)=-0.25d0*pbbf(i,j+1,k-1)/pdy(i,j,2)
            pbb4(i,j-1,k-1)=0.25d0*pbbf(i,j-1,k-1)/pdy(i,j,2)
            pgc2(i,j,k-1)=1.d0/(dzz1(k-2)*pdq(i,j,2))*ramp
            pgen(i,j,k) = - pdq(i,j,2)*dz1(k-1)*((1.d0/pddx(i,j,2)+1.d0/pddx(i-1,j,2))/ &
                  pdx(i,j,2)+(1.d0/pddy(i,j,2)+1.d0/pddy(i,j-1,2))/pdy(i,j,2))- &
                  (1.d0/dzz1(k-2)+1.d0/dzz1(k-1))/(pdq(i,j,2))*ramp
          else
            pgen(i,j,k) = - pdq(i,j,2)*dz1(k-1)*((1.d0/pddx(i,j,2)+1.d0/pddx(i-1,j,2))/ &
                  pdx(i,j,2)+(1.d0/pddy(i,j,2)+1.d0/pddy(i,j-1,2))/pdy(i,j,2))-&
                  (1.d0/dzz1(k-1)+1.d0/dzz1(k-1))/(pdq(i,j,2))*ramp
          endif
          if(k==kb)then
            pga1(i+1,j,k)=pga1(i+1,j,k)+paa1(i+1,j,k+1)
            pga2(i-1,j,k)=pga2(i-1,j,k)+paa2(i-1,j,k+1)
            pgb1(i,j+1,k)=pgb1(i,j+1,k)+pbb1(i,j+1,k+1)
            pgb2(i,j-1,k)=pgb2(i,j-1,k)+pbb2(i,j-1,k+1)
            pgen(i,j,k)=pgen(i,j,k)+pgc1(i,j,k+1)
          endif
          if(k==3)then

            paa3(i+1,j,k-1)=0.d0
            paa4(i-1,j,k-1)=0.d0
            pbb3(i,j+1,k-1)=0.d0
            pbb4(i,j-1,k-1)=0.d0
            pgc2(i,j,k-1)=0.d0

            !paa2(i+1,j,k-1)=0.d0
            !pbb2(i,j+1,k-1)=0.d0 !need to debug
            !pbb4(i,j-1,k-1)=0.d0
            !pgc2(i,j,k-1)=0.d0
          endif

          if(i_global==3)then
            if(k>2)then
              pgc2(i,j,k-1)=pgc2(i,j,k-1)+paa4(i-1,j,k-1)
            endif
            pgc1(i,j,k+1)=pgc1(i,j,k+1)+paa2(i-1,j,k+1)
            pgen(i,j,k)=pgen(i,j,k)+pga2(i-1,j,k)
          else
            if(pfsm(i-1,j,2)<0.9d0)then
              pgen(i,j,k)=pgen(i,j,k)+pga2(i-1,j,k)
            endif
          endif
     
          if(i_global==im)then
            if(k>2)then 
              pgc2(i,j,k-1)=pgc2(i,j,k-1)+paa3(i+1,j,k-1)
            endif
            pgc1(i,j,k+1)=pgc1(i,j,k+1)+paa1(i+1,j,k+1)
            pgen(i,j,k)=pgen(i,j,k)+pga1(i+1,j,k)
          else
            if(pfsm(i+1,j,2)<0.9d0)then
              pgen(i,j,k)=pgen(i,j,k)+pga1(i+1,j,k)
            endif
          endif

          if(j_global==3)then
            if(k>2)then
              pgc2(i,j,k-1)=pgc2(i,j,k-1)+pbb4(i,j-1,k-1)
            endif
            pgc1(i,j,k+1)=pgc1(i,j,k+1)+pbb2(i,j-1,k+1)
            pgen(i,j,k)=pgen(i,j,k)+pgb2(i,j-1,k)
          else
            if(pfsm(i,j-1,2)<0.9d0)then
              pgen(i,j,k)=pgen(i,j,k)+pgb2(i,j-1,k)
            endif
          endif

          if(j_global==jm)then
            if(k>2)then
              pgc2(i,j,k-1)=pgc2(i,j,k-1)+pbb3(i,j+1,k-1)
            endif
            pgc1(i,j,k+1)=pgc1(i,j,k+1)+pbb1(i,j+1,k+1)
            pgen(i,j,k)=pgen(i,j,k)+pgb1(i,j+1,k)
          else
            if(pfsm(i,j+1,2)<0.9d0)then
              pgen(i,j,k)=pgen(i,j,k)+pgb1(i,j+1,k)
            endif
          endif
        if(k<kb)then 
          if(pfsm(i+1,j,2)>0.9d0)then
            m=m+1
            apr(m)=paa1(i+1,j,k+1)
            ja(m)=pnumnelt(i+1,j,k+1)-1
          endif

          if(pfsm(i,j+1,2)>0.9d0)then
            m=m+1
            apr(m)=pbb1(i,j+1,k+1)
            ja(m)=pnumnelt(i,j+1,k+1)-1
          endif

          if(pfsm(i,j,2)>0.9d0)then
            m=m+1
            apr(m)=pgc1(i,j,k+1)
            ja(m)=pnumnelt(i,j,k+1)-1
          endif

          if(pfsm(i,j-1,2)>0.9d0)then
            m=m+1
            apr(m)=pbb2(i,j-1,k+1)
            ja(m)=pnumnelt(i,j-1,k+1)-1
          endif

          if(pfsm(i-1,j,2)>0.9d0)then
            m=m+1
            apr(m)=paa2(i-1,j,k+1)
            ja(m)=pnumnelt(i-1,j,k+1)-1
          endif
        endif

          if(pfsm(i+1,j,2)>0.9d0)then
            m=m+1
            apr(m)=pga1(i+1,j,k)
            ja(m)=pnumnelt(i+1,j,k)-1
          endif

          if(pfsm(i,j+1,2)>0.9d0)then
            m=m+1
            apr(m)=pgb1(i,j+1,k)
            ja(m)=pnumnelt(i,j+1,k)-1
          endif
 
          m=m+1
          apr(m)=pgen(i,j,k)
          ja(m)=pnumnelt(i,j,k)-1

          if(pfsm(i,j-1,2)>0.9d0)then
            m=m+1
            apr(m)=pgb2(i,j-1,k)
            ja(m)=pnumnelt(i,j-1,k)-1
          endif

          if(pfsm(i-1,j,2)>0.9d0)then
            m=m+1
            apr(m)=pga2(i-1,j,k)
            ja(m)=pnumnelt(i-1,j,k)-1
          endif

          if(k>=3)then
            if(pfsm(i+1,j,2)>0.9d0)then
              m=m+1
              apr(m)=paa3(i+1,j,k-1)
              ja(m)=pnumnelt(i+1,j,k-1)-1
            endif

            if(pfsm(i,j+1,2)>0.9d0)then
              m=m+1
              apr(m)=pbb3(i,j+1,k-1)
              ja(m)=pnumnelt(i,j+1,k-1)-1
            endif

            if(pfsm(i,j,2)>0.9d0)then
              m=m+1
              apr(m)=pgc2(i,j,k-1)
              ja(m)=pnumnelt(i,j,k-1)-1
            endif

            if(pfsm(i,j-1,2)>0.9d0)then
              m=m+1
              apr(m)=pbb4(i,j-1,k-1)
              ja(m)=pnumnelt(i,j-1,k-1)-1
           endif

           if(pfsm(i-1,j,2)>0.9d0)then
             m=m+1
             apr(m)=paa4(i-1,j,k-1)
             ja(m)=pnumnelt(i-1,j,k-1)-1
           endif
         endif
       endif
     enddo
   enddo
 enddo
 ia(cnn+1)=cnz
 call MatCreateMPIAIJWithArrays(MPI_COMM_WORLD,cnn,cnn,p_N,p_N,ia,ja,apr,p_A,ierr)
! call MatView(p_A,PETSC_VIEWER_STDOUT_WORLD,ierr)
  r1=uf*AXB(dq)*dz
  r2=vf*AYB(dq)*dz
  call coef2(aaf,bbf,ccf,fff) 
  wq=-fff+wwf+ccf
  r=1.d0/dti2*(DXF(r1)+DYF(r2)-DZF(wq)*dz)
  call get_local_buffer(pr,r)    
  call get_local_buffer(pq,q) 
  m=0
  do k=2,kb
    do j=2,ind(2)-1
      do i=2,ind(1)-1
        if(pfsm(i,j,2)>0.9d0)then
          m=m+1
          iaa(m)= pnumnelt(i,j,k)-1
          rr(m) = pr(i,j,k)
        endif
      enddo
    enddo
  enddo 
  call VecSetValues(p_b,cnn,iaa,rr,INSERT_VALUES,ierr)
  call VecAssemblyBegin(p_b,ierr)
  call VecAssemblyEnd(p_b,ierr)
!  call VecView(p_b,PETSC_VIEWER_STDOUT_WORLD,ierr)

  !ksp solve-------------------:-----------------------------------
  call KSPSetOperators(p_ksp,p_A,p_A,ierr)
  call KSPSetInitialGuessNonzero(p_ksp,PETSC_TRUE,ierr)
  call KSPSetTolerances(p_ksp,tol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,200,ierr)
  call KSPSetFromOptions(p_ksp,ierr)
  call KSPSolve(p_ksp,p_b,p_x,ierr)
!  call KSPView(p_ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
  call VecGetValues(p_x,cnn,ix,qq,ierr)
  call MatDestroy(p_A,ierr)
  !---------------------------------------------------------------
  m=0
  do k=2,kb
    do j=2,ind(2)-1
      do i=2,ind(1)-1
        if(pfsm(i,j,2)>0.9d0)then
          m=m+1
          pq(i,j,k)=qq(m)
        endif
      enddo
    enddo
  enddo

  !ksp destroy----------------------------------------------------
  if(iint==iend)then
  deallocate(ix)
  call VecDestroy(p_b,ierr)
  call VecDestroy(p_x,ierr)
  call KSPDestroy(p_ksp,ierr)
  endif
  !---------------------------------------------------------------

  call set(sub(q,1,':',':'),sub(q,2,':',':'))
  call set(sub(q,im,':',':'),sub(q,imm1,':',':'))
  call set(sub(q,':',1,':'),sub(q,':',2,':'))
  call set(sub(q,':',jm,':'),sub(q,':',jmm1,':'))

  call set(sub(q,':',':',kb),sub(q,':',':',kbm1))
  call set(sub(q,':',':',1),0.d0)

  fsm=fsm_tmp
end subroutine
