
# PETSC_INC=-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
# PETSC_LIB=-L${PETSC_DIR}/${PETSC_ARCH}/lib
#export OPEN_ARRAY=/home/siofive/OpenArray_CXX/build
#export OPEN_ARRAY=/GPFS/cess/wangmq/OpenArray_CXX/build
#export OPEN_ARRAY=/GPFS/cess/wangmq/openarray_cxx_cess/build
#export OPEN_ARRAY=/GPFS/cess/wangmq/build.4.23/

export OPEN_ARRAY=/GPFS/cess/wangmq/OpenArray/build.4.23/

SRCDIR           = ./src/
OBJDIR           = ./lib/
BINDIR           = ./bin/

EXT_LIB = -L${EXT_PATH}/lib64/
JIT_LIB = ${EXT_PATH}/lib64/

#petsc----------------
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
#---------------------

EXE              = nh_gomo
#FC	         = mpif90
#FLINKER          = mpif90 
FC	         = mpiifort -O2 -g -DBOOST_LOG_DYN_LINK -w
FLINKER          = mpiifort -O2 -g -w
CFLAGS	         =
FFLAGS	         =-Wno-tabs -I ${EXT_PATH}/include \
		  -J ${OBJDIR} -I ${OPEN_ARRAY} -g \
		 -fbacktrace -ffree-line-length-0 
CPPFLAGS         =
FPPFLAGS         =
CLEANFILES       = nh_gomo *.o *.mod *.nc
NP               = 1

OOO              = \
		print_section.o read_var.o \
		dens.o baropg.o bottom_friction.o lateral_bc.o \
		get_time.o advct.o advt2.o lateral_viscosity.o \
		advave.o bcond1.o smoth_update.o smol_adif.o\
		bcond2_ua.o bcond2_va.o external_el.o \
		external_ua.o external_va.o external_update.o \
		bcond3_u.o bcond3_v.o bcond4.o bcond6.o \
		internal_q.o adjust_uv.o adjust_ufvf.o \
		mode_interaction.o internal_t.o internal_u.o \
		internal_v.o internal_w.o \
		internal_update.o \
		surface_forcing.o update_initial.o \
		oa_kernels.o


OBJ = ${addprefix ${OBJDIR}/, config.o \
        variables.o dens.o read_init.o  \
        init_fields.o update_initial.o \
        bottom_friction.o get_time.o \
        fsm_compute.o surface_forcing.o \
        advct.o baropg.o \
        lateral_viscosity.o advave.o \
        mode_interaction.o bcond1.o \
        external_el.o bcond2_ua.o \
        external_ua.o external_va.o \
        bcond2_va.o external_update.o \
        advtl.o proft.o smoth_update.o \
        advu.o advv.o profu.o profv.o pvariable.o \
        advctw.o advew.o coef1.o coef2.o subinv.o wvert.o \
        construct_matrix.o faz.o adjust_uava.o \
        internal_update.o bcond5.o  \
        bcond6.o \
        bcond4.o  \
        bcond3_u.o bcond3_v.o \
        print_section.o \
        nh_gomo.o}



OBJMAIN	= ${OBJ}

.DEFAULT_GOAL := all

all :
	@./fypp1 -p  -m re -m string -m io -m os --create-parents \
	src/oa_kernels.fypp src/oa_kernels.F90
	@make main


${OBJDIR}/%.o: ${SRCDIR}/%.F90
	-${FC} ${FFLAGS} -c -wd1572 -g -I/GPFS/cess/liuc/Libaraies/petsc-3.6.4/include -I/GPFS/cess/liuc/Libaraies/petsc-3.6.4/linux-intel/include -o $@ $<


%.o: %.mod

main : ${OBJ} chkopts
	-${FLINKER} -o ${BINDIR}/nh_gomo ${OBJ} ${PETSC_KSP_LIB} \
	${EXT_LIB} -I ${OPEN_ARRAY} -L${OPEN_ARRAY} \
	-lopenarray -lm -ldl -lstdc++ \
	-lboost_program_options -lboost_system -lboost_log -lboost_log_setup -lboost_thread -ljit -lpnetcdf # /home/siofive/pnetcdf.intel/lib/libpnetcdf.a

clean ::
	-rm lib/*.o
	-rm *.mod
	-rm bin/${EXE}

