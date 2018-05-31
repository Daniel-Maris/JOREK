#!/bin/bash
eval `tclsh /work/imas/opt/modules-tcl/modulecmd.tcl $(basename $SHELL) autoinit`

module purge

module use /work/imas/opt/EasyBuild/modules/all
module use /work/imas/etc/modules/all

module load cURL/7.28.1-goolf-1.5.16
module load intel/17.0.4
module load mpich2/3.1.3-intel
module load scotch/5.1.12b
module load metis/5.1.0
module load mumps/4.10.0
module load pastix/5.2.2.16
module load FFTW/3.3.4-gompi-1.5.16
module load hdf5/1.10.2-intel-12.0.2

module list

export LANG=C
export JOREK_HOST=iter-hpc
export compilethreads=4
export MAKEFLAGS="-j$compilethreads"
export PRERUN="export OMP_NUM_THREADS=4"
export MPIRUN="mpirun -np "
export BATCHCOMMAND="qsub"

export http_proxy=${JOREK_HTTP_PROXY}
export https_proxy=${JOREK_HTTP_PROXY}
