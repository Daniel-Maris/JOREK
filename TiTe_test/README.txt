In all cases, JOREK is first run with model 083, and the intear_init namelist is given as input. The hard-coded parameters are as follows:
    n_tor = 5
    n_coord_tor = 21
    n_period = 5
    n_coord_period = 5
    n_plane = 40
    l_pol_domm = 9
This computes the initial condidtions for the reduced MHD variables. The simulation is then restarted with model 183, using the namelist intear_prerun, with the binary being compiled with the same hard-coded parameters as above, except for the model.
Finally, after the pre-run is finished the simulation is restarted with the same binary as the previous step and the namelist intear_production is used for this run.

Note that the initialization phase (gvec2jorek.dat and intear_init) as well as the Dommaschk coefficients (dc_W7A_unst_10kPa) are the same for both cases.

Note that all parameters in the namelist having to do with initial density, temperature and FF' profiles, grid construction and boundary are ignored/overwritten by the GVEC import, except for n_flux and n_tht, which must match Ns and Ntheta in the gvec2jorek.dat file, respectively. Further, resistivity, viscosity and other diffusive parameters play no role in the initialization.