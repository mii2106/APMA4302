include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CFLAGS += -pedantic -std=c99


reaction: reaction2d.o
	-${CLINKER} -o reaction2d reaction2d.o ${PETSC_LIB}
	${RM} reaction.o

# testing

runreaction2d_1:
	-@../testit.sh reaction2d "-snes_converged_reason -da_refine 1" 1 1

runreaction2d_2:
	-@../testit.sh reaction2d "-da_grid_x 7 -mat_is_symmetric" 1 2

runreaction2d_3:
	-@../testit.sh reaction2d "-snes_fd_color -da_refine 1 -snes_view" 4 3

runreaction2d_4:
	-@../testit.sh reaction2d "-snes_converged_reason -snes_linesearch_type cp" 1 4

runreaction2d_5:
	-@../testit.sh reaction2d "-snes_test_jacobian" 1 5

runreaction2d_6:
	-@../testit.sh reaction2d "-snes_converged_reason -rct_noRinJ" 1 6


test_reaction2d: runreaction2d_1 runreaction2d_2 runreaction2d_3 runreaction2d_4 runreaction2d_5 runreaction2d_6

test: test_reaction2d

.PHONY: runreaction2d_1 runreaction2d_2 runreaction2d_3 runreaction2d_4 runreaction2d_5 runreaction2d_6 test test_reaction2d

distclean: clean

clean::
	@rm -f *~ reaction2d *tmp

