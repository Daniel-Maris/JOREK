# --- General settings
jorekmodel="600"
description="kinetic_main with rep coupling scheme using model$jorekmodel simulating a circular cross section plasma, tearing mode onset."
options="with_vpar=.false. with_TiTe=.false. with_neutrals=.false. with_impurities=.false. with_refluid=.false."
mpitasks=8
num_threads=8
binaries="kinetic_main"
extra_restart="part_restart.h5"
requiredfiles="input part_restart.h5"
extra_remote_files="part_restart.h5"
particle_example="kinetic_main"
particle_example_dir="particles/examples"

# NOTE: This case starts from late in a particle run, therefore requiring both a fluid and particle restart file
#       The initial_run - which would normally generate fluid restart - is therefore not used
#       Instead, inital_run just copies the begin.h5 file to fluid restart, rather than re-generating it from scratch


# --- Compile the code for the test case
function compile_jorek () {
  ./util/config.sh model=$jorekmodel n_tor=3 n_coord_tor=1 l_pol_domm=0 n_plane=4 n_period=1 n_coord_period=1   $options  || exit 1
  make $compilopt $debugoptions kinetic_main                                                                              || exit 1
}


# --- Initial run only required when preparing or updating the test case
function initial_run () {
  cp ${codedir}/reg_tests/testcases/particle_rep_600/begin.h5 ./jorek_restart.h5                                          || exit 1
  echo "Initial run not compatible with this test case, as both a fluid and particle restart are required."               || exit 1
}


# --- Carry out the test case
function restart_run () {
  ${codedir}/util/setinput.sh input restart_particles=.t. restart=.t. nstep_n=2 tstep_n=1                                 || exit 1
  export OMP_NUM_THREADS=$num_threads                                                                                     || exit 1
  echo "setting OMP_NUM_THREADS=$num_threads, due to the requirements of the test"                                        || exit 1
  $MPIRUN $mpitasks ./kinetic_main < input | tee logfile                                                                  || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  compare_results_generic 1.e-10                                                                                           || exit 1
}
