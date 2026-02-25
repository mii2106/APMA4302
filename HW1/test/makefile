include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CFLAGS += -pedantic -std=c99

expx: expx.o
	-${CLINKER} -o expx expx.o  ${PETSC_LIB}
	${RM} expx.o
# etc

.PHONY: clean

clean::
	@rm -f *~ expx balanced

