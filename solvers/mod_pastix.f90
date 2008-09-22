module pastix_module             ! PastiX specific variables
  integer(kind=8)       :: pastix_data
  integer               :: pastix_iparm(64)
  real*8                :: pastix_dparm(64)
  integer,allocatable   :: pastix_perm_vars(:), pastix_iperm_vars(:)
  integer,allocatable   :: sparskit_work(:)
  integer, allocatable  :: ihwb(:),iwk(:)
  logical               :: use_pastix, pastix_initialised, pastix_analysed, pastix_smp_only
  integer               :: pastix_facto
  parameter (pastix_facto=2)
  integer               :: pastix_rhs
  parameter (pastix_rhs=0)
  integer               :: pastix_sym
  parameter (pastix_sym=1)
  integer               :: pastix_nthrd=1
 ! parameter (pastix_nthrd=1)
  integer               :: pastix_iter
  parameter (pastix_iter=250)
  integer               :: pastix_ricar
  parameter (pastix_ricar=0)
  integer               :: pastix_iluk
  parameter (pastix_iluk=3)
  integer               :: pastix_amalg
  parameter (pastix_amalg=5)
  integer               :: pastix_endsolve
  parameter (pastix_endsolve=5)
  real*8                :: pastix_epsilon
  parameter (pastix_epsilon=1.d-12)
  real*8                :: pastix_pivot
  parameter (pastix_pivot=1.d-32)
endmodule pastix_module
