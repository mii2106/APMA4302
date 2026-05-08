include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CFLAGS += -pedantic -std=c99

biharm: biharm.o poissonfunctions.o
	-${FLINKER} -o biharm biharm.o poissonfunctions.o ${PETSC_LIB}
	${RM} biharm.o poissonfunctions.o



.PHONY: clean distclean 

distclean: clean

clean::
	@rm -f *~  biharm *tmp
	@rm -f *.pdf *.txt
