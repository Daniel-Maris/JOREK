include config.in

DIRS = datatypes models/$(MODEL) models communication elements grids matrix solvers plots diagnostics

MAIN = jorek_$(MODEL)

.PHONY : modules sources $(MAIN) jorek2_diagno jorek2vtk jorek2vtk_3d

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
	$(FC) $(FFLAGS)      \
	jorek2_main.f90      \
	datatypes/*.o        \
	models/$(MODEL)/*.o  \
	models/*.o           \
	communication/*.o    \
	elements/*.o         \
	grids/*.o            \
	matrix/*.o           \
	solvers/*.o          \
	plots/*.o            \
	diagnostics/*.o      \
	 -o $(MAIN) $(INCLUDES) $(LIBS)

jorek2vtk : modules sources
	$(FC) $(FFLAGS)                   \
	diagnostics/jorek2vtk.f90         \
	datatypes/mod_parameters.o        \
	datatypes/mod_data_structure.o    \
	elements/mod_gauss.o              \
	models/$(MODEL)/mod_phys_module.o \
	elements/mod_basis_at_gaussian.o  \
	communication/import_restart.o    \
	elements/initialise_basis.o       \
	elements/basis_functions.o        \
	elements/basis_functions2.o       \
        elements/interp.o                 \
	elements/interp_RZ.o              \
	 -o $(JOREK_DIR)/jorek2vtk $(INCLUDES) $(LIBS)
	 
jorek2vtk_3d : modules sources
	$(FC) $(FFLAGS)                   \
	diagnostics/jorek2vtk_3d.f90      \
	datatypes/mod_parameters.o        \
	datatypes/mod_data_structure.o    \
	elements/mod_gauss.o              \
	models/$(MODEL)/mod_phys_module.o \
	elements/mod_basis_at_gaussian.o  \
	communication/import_restart.o    \
	elements/initialise_basis.o       \
	elements/basis_functions.o        \
	elements/basis_functions2.o       \
        elements/interp.o                 \
	elements/interp_RZ.o              \
	 -o $(JOREK_DIR)/jorek2vtk_3d $(INCLUDES) $(LIBS)

jorek2_diagno : modules sources
	$(FC) $(FFLAGS)                   \
	diagnostics/jorek2_diagno.f90     \
	datatypes/mod_parameters.o        \
	datatypes/mod_data_structure.o    \
	communication/export_helena.o     \
	plots/plot_velocity_profile.o     \
	diagnostics/find_flux_surfaces.o  \
	diagnostics/flux_surface_add_line.o  \
	diagnostics/flux_surface_add_point.o  \
	diagnostics/find_axis.o           \
	diagnostics/find_xpoint.o         \
	diagnostics/RZ_minmax.o           \
	diagnostics/psi_minmax.o          \
	diagnostics/integrals.o           \
	solvers/mnewtax.o                 \
	solvers/solvP3.o                  \
	solvers/solve_M2.o                \
	solvers/root.o                    \
	elements/*.o                      \
	models/$(MODEL)/mod_phys_module.o \
	communication/import_restart.o    \
	grids/find_RZ.o                   \
	grids/meshac.o                    \
	grids/fgauss.o                    \
	 -o $(JOREK_DIR)/jorek2_diagno $(INCLUDES) $(LIBS)
