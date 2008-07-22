module pastix_module             ! PastiX specific variables
  integer(kind=8)       :: pastix_data
  integer               :: pastix_iparm(64)
  real*8                :: pastix_dparm(64)
  integer,allocatable   :: pastix_perm_vars(:), pastix_iperm_vars(:)
  integer,allocatable   :: sparskit_work(:)
  integer, allocatable  :: ihwb(:),iwk(:)
  logical               :: use_pastix, pastix_initialised, pastix_analysed
endmodule pastix_module