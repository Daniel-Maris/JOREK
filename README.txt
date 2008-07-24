*******************************************************
*                    JOREK2                           *
*                                                     *
* copyright : Guido Huysmans (Association Euratom/CEA *
*******************************************************

to build jorek2:

- edit the config.in file to match your environment
  (some examples are available in the directory configs)
  
  choose a model : model199 (reduced MHD, no v_parallel)
                   model300 (reduced MHD, with V_parallel)

- edit the model dependent mod_parameters.f90 file to define
  n_tor    : the number of toroidal harmnonics (sin+cos) 
  n_plane  : the number of toroidal planes (for FFT)
  n_period : the periodicity of the torus
  
- make all (in top directory)

to execute jorek2:

mpirun -np n_cpu jorek2_model199 < namelist/model199/intear

this is a simple testcase (example output is available in test_output)

the number of cpu's should be a multiple of the number of harmonics ((n_tor-1)/2+1)
