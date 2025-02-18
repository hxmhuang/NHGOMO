module pvariable
  !use petscvec
  use petscsys
  implicit none
#include "petsc.h"
  PetscInt :: aaa1
  PetscInt  :: p_size,p_rank,p_N,p_numperrow,p_i,p_j,p_nn
  PetscScalar :: p_aprA,p_aprb
  Mat       :: p_A 
  Vec       :: p_b,p_x,p_x_all
  KSP       :: P_ksp
  PC        :: pc
  VecScatter:: p_ctx
  IS        :: p_from,p_to
  PetscReal :: tol
  PetscInt,dimension(:),allocatable :: ix

  !call PCSetType(pc,PCJACOBI,ierr)
  !call KSPSetType(p_ksp,KSPCG,ierr)
end module
