module vacuum_response_module
!-----------------------------------------------------------------------
! module contains the vacuum response matrix for all toroidal harmonics
!-----------------------------------------------------------------------

integer             :: n_dof_bnd                  ! degrees of freedom on the boundary
real*8, allocatable :: vacuum_response(:,:,:)     ! the vacuum response matrix (idrive,ireponse,itor)

end module
