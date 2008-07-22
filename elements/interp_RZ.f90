subroutine interp_RZ(node_list,element_list,i_elm,s,t,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
!-----------------------------------------------------------------------------------------
! subroutine calculates the interpolation within one element (i_elm) for a given position
! (s,t) in the local coordinates
!-----------------------------------------------------------------------------------------
use data_structure
implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list

real*8  :: G(4,4), G_s(4,4), G_t(4,4), G_st(4,4), G_ss(4,4), G_tt(4,4)
real*8  :: s,t, R, R_s, R_t, R_st, R_ss, R_tt, Z,Z_s,Z_t, Z_st, Z_ss, Z_tt
integer :: kv, iv, kf, i_elm

call basisfunctions2(s,t,G(1:4,1:4),G_s(1:4,1:4),G_t(1:4,1:4),G_st(1:4,1:4),G_ss(1:4,1:4),G_tt(1:4,1:4))

R = 0.d0; R_s = 0.d0; R_t = 0.d0; R_st = 0.d0; R_ss = 0.d0; R_tt = 0.d0
Z = 0.d0; Z_s = 0.d0; Z_t = 0.d0; Z_st = 0.d0; Z_ss = 0.d0; Z_tt = 0.d0

do kv = 1,n_vertex_max  ! 4 vertices

  iv = element_list%element(i_elm)%vertex(kv)  ! the node number

  do kf = 1, n_order+1       ! 4 basis functions

    R    = R    + node_list%node(iv)%x(kf,1) * element_list%element(i_elm)%size(kv,kf) * G(kv,kf)
    R_s  = R_s  + node_list%node(iv)%x(kf,1) * element_list%element(i_elm)%size(kv,kf) * G_s(kv,kf)
    R_t  = R_t  + node_list%node(iv)%x(kf,1) * element_list%element(i_elm)%size(kv,kf) * G_t(kv,kf)
    R_st = R_st + node_list%node(iv)%x(kf,1) * element_list%element(i_elm)%size(kv,kf) * G_st(kv,kf)
    R_ss = R_ss + node_list%node(iv)%x(kf,1) * element_list%element(i_elm)%size(kv,kf) * G_ss(kv,kf)
    R_tt = R_tt + node_list%node(iv)%x(kf,1) * element_list%element(i_elm)%size(kv,kf) * G_tt(kv,kf)

    Z    = Z    + node_list%node(iv)%x(kf,2) * element_list%element(i_elm)%size(kv,kf) * G(kv,kf)
    Z_s  = Z_s  + node_list%node(iv)%x(kf,2) * element_list%element(i_elm)%size(kv,kf) * G_s(kv,kf)
    Z_t  = Z_t  + node_list%node(iv)%x(kf,2) * element_list%element(i_elm)%size(kv,kf) * G_t(kv,kf)
    Z_st = Z_st + node_list%node(iv)%x(kf,2) * element_list%element(i_elm)%size(kv,kf) * G_st(kv,kf)
    Z_ss = Z_ss + node_list%node(iv)%x(kf,2) * element_list%element(i_elm)%size(kv,kf) * G_ss(kv,kf)
    Z_tt = Z_tt + node_list%node(iv)%x(kf,2) * element_list%element(i_elm)%size(kv,kf) * G_tt(kv,kf)

  enddo

enddo

return
end
