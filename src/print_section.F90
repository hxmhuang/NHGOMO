subroutine print_section(iint,iext,vamax,imax,jmax,vtot,atot,taver,saver,eaver,tsalt)
  use openarray
  use config
  use variables
  implicit none
  real(kind=8),intent(out) :: eaver, taver, saver, tsalt,vtot, atot
  type(array)::darea
  real(kind=8) :: realdz
  real(kind=8),intent(in) :: vamax
  integer,intent(in) :: imax, jmax, iint, iext
  integer :: ierr
  if (iint >= iswtch) then
     iprint=int(prtd2*24.0*3600.0/dti + 0.5)
  endif
     realdz=sum(dz,3)
     darea=dx*dy*fsm
     vtot =sum(sum(realdz * darea * dt,2),1)
     atot =sum(sum(darea,2),1)
     taver=sum(sum(sum(dz*tb*darea*dt,3),2),1)/vtot
     saver=sum(sum(sum(dz*sb*darea*dt,3),2),1)/vtot
     tsalt=(saver+sbias)*vtot
     eaver=sum(sum(et*darea,2),1)/atot
  if (mod(iint,int(10))==0 .or. vamax >= vmaxl) then
     if(get_rank() == 0) then 
        write(*,*) '******************************************************************************'
        write(*,'(A,F16.7,A,I8,A,i8,A,I8)') &
             " time=",time,    ";   iint=",iint, &
             ";   iext=",iext, ";   iprint=",iprint
        write(*,*) '******************************************************************************'
        write(*,'(A,F25.12,A,F22.12,A,F30.12/A,F25.12,		&
             A,F22.12,A,F30.12)') &
             '   vtot=',vtot,  ';   atot=',atot, &
             ';  eaver=',eaver, '  taver=',taver, &
             ';  saver=',saver, ';  tsalt=',tsalt
     endif
     if (vamax > vmaxl) then
        if(get_rank() == 0) then
           write(*,*) '*********************************************************'
           write(*,'(A,F16.7,A,I8,A,I8,A,I8)')"time=",time,";   iint=",iint, &
               ";   iext=",iext,";   iprint=",iprint
           write(*,*) '************************************************'
           write(*,*) '************ abnormal job end ******************' 
           write(*,*) '************* user terminated ******************'
           write(*,*) '************************************************'
           write(*,'(A,E12.3, A,I5,A,I5)') &
                'vamax=', vamax, ';    imax=', imax, ';    jmax=', jmax
        endif
        stop
     endif
  endif
end subroutine print_section
