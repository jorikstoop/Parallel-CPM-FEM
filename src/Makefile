OBJECTS = CPM_dH.o \
	FE_assembly.o \
	FE_local.o \
	FE_nodes2dofs.o \
	FE_solver.o \
	cellforces.o \
	cellmoves.o \
	init.o \
	mt.o \
	mylib.o \
	cpmfem.o \
	FE_csr.o \
	write.o

# consider using higher compiler omptimization level (-O2 or -O3)
CFLAGS = -O 
OMPFLAGS = -fopenmp

cpmfem: $(OBJECTS)
	mpicc $(CFLAGS) -g -o cpmfem $(OBJECTS) $(OMPFLAGS) -lm

clean:
	rm -f cpmfem
	rm -f *.o

.c.o:
	mpicc $(CFLAGS) $(OMPFLAGS) -c $< -g -o $@
