include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CFLAGS += -pedantic -std=c99

expx: HW1.o
	-${CLINKER} -o HW1 HW1.o  ${PETSC_LIB}
	${RM} HW1.o
# etc

.PHONY: clean

clean::
	@rm -f *~ expx balanced

