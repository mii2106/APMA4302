include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CFLAGS += -pedantic -std=c99

tri: tri.o
	-${CLINKER} -o bvp bvp.o  ${PETSC_LIB}
	${RM} bvp.o


# testing

runtri_1:
	-@../testit.sh bvp "-a_mat_view ::ascii_dense" 1 1

runtri_2:
	-@../testit.sh tri "-bvp_m 1000 -ksp_rtol 1.0e-4 -ksp_type cg -pc_type bjacobi -sub_pc_type jacobi -ksp_converged_reason" 2 2

runloadsolve_1:
	-@./bvp -ksp_view_mat binary:A.dat -ksp_view_rhs binary:b.dat > /dev/null
	-@../testit.sh loadsolve "-verbose -fA A.dat -fb b.dat -ksp_view_mat -ksp_view_rhs -ksp_view_solution" 1 1

test_tri: runtri_1 runtri_2


test:  test_tri 

# etc

.PHONY: clean distclean runtri_1 runtri_2 test  test_tri 

distclean: clean

clean::
	@rm -f *~ tri  *tmp
	@rm -f *.dat *.dat.info

