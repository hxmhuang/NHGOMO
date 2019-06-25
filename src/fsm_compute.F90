subroutine fsm_compute()
  use openarray
  use config
  use variables
  use mpi
  implicit none
  type(array) :: fsm_tmp
  integer :: i,j,k,m,size,rank,IEEOR,ind(6),local_ind(3)
  real(kind=8) :: tmp0,tmp1,tmp2,tmp3,tmp4
  integer,allocatable :: nums(:)
  real(kind=8),pointer :: pfsm(:,:,:)
  integer,pointer :: pnumnelt(:,:,:)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,IEEOR)
  rank=get_rank()
  allocate(nums(size))
  fsm_tmp=fsm
  call set(sub(fsm_tmp,1,':'),0.d0)
  call set(sub(fsm_tmp,im,':'),0.d0)
  call set(sub(fsm_tmp,':',1),0.d0)
  call set(sub(fsm_tmp,':',jm),0.d0)
  call update_ghost(fsm_tmp)
  ind=get_box_corners(fsm_tmp)
  call get_local_buffer(pfsm,fsm_tmp)
  call get_local_buffer(pnumnelt,numnelt)
  local_ind=shape(pfsm)
  cnn=0;cnz=0
  do j=ind(3)+1,ind(4)
  do i=ind(1)+1,ind(2)
    tmp0 = local_sub(fsm_tmp,i,j,1)
    if(tmp0>0.9d0)then
      cnn=cnn+kbm1
      cnz=cnz+3*kbm2+1
      tmp1= local_sub(fsm_tmp,i-1,j,1)
      tmp2= local_sub(fsm_tmp,i+1,j,1)
      tmp3= local_sub(fsm_tmp,i,j-1,1)
      tmp4= local_sub(fsm_tmp,i,j+1,1)
      if(tmp1>0.9d0) cnz=cnz+3*kbm2+1
      if(tmp2>0.9d0) cnz=cnz+3*kbm2+1
      if(tmp3>0.9d0) cnz=cnz+3*kbm2+1
      if(tmp4>0.9d0) cnz=cnz+3*kbm2+1
    endif
  enddo
  enddo
  CALL MPI_ALLGATHER(cnn,1,MPI_INTEGER,nums,1,MPI_INTEGER,MPI_COMM_WORLD,IEEOR)
!  CALL MPI_ALLREDUCE(cnz,numnz,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IEEOR)
  total_nums=sum(nums)
  nums_s=sum(nums(1:rank))+1
  m=nums_s
  allocate(apr(cnz),ia(cnn+1),ja(cnz),rr(cnn),qq(cnn),iaa(cnn))
  do k=2,kb
  do j=2,local_ind(2)-1
  do i=2,local_ind(1)-1
    if(pfsm(i,j,2)>0.9d0)then
      pnumnelt(i,j,k)=m
      m=m+1
    endif
  enddo
  enddo
  enddo
end
