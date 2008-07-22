include ../config.in

all:  modules \
      coicsr.o \
      construct_matrix.o \
      locate_irn_jcn.o \
      reduce.o \
      scale_global_matrix.o \
      sort_I_mrgrnk.o \
      global_matrix_structure.o
           
modules : mod_global_distributed_matrix.o
 
clean :
	rm *.o *.mod
