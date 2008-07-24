include config.in

DIRS = datatypes models/$(MODEL) models communication elements grids matrix solvers plots diagnostics

MAIN = jorek_$(MODEL)

.PHONY : modules sources $(MAIN)

all :   modules sources $(MAIN)

modules :
	for dir in $(DIRS); do     \
          ($(MAKE) -C $$dir modules) \
        done

sources :	
	for dir in $(DIRS); do \
          ($(MAKE) -C $$dir all) \
        done
	
clean :	
	rm $(MAIN) ; \
	for dir in $(DIRS); do   \
          ($(MAKE) -C $$dir clean) \
        done
	
	
$(MAIN) : jorek2_main.f90
	$(FC) $(FFLAGS)   \
	jorek2_main.f90   \
	datatypes/*.o     \
	$(MODEL_DIR)/*.o  \
	models/*.o        \
	communication/*.o \
	elements/*.o      \
	grids/*.o         \
	matrix/*.o        \
	solvers/*.o       \
	plots/*.o         \
	diagnostics/*.o   \
	 -o $(MAIN) $(INCLUDES) $(LIBS)
