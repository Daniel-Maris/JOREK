subroutine energy(node_list,element_list,W_mag,W_kin)
!---------------------------------------------------------------
!
!---------------------------------------------------------------
use data_structure
use Gauss
use basis_at_gaussian
use phys_module

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

real*8     :: x_g(n_gauss,n_gauss),        x_s(n_gauss,n_gauss),        x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss),        y_s(n_gauss,n_gauss),        y_t(n_gauss,n_gauss)
real*8     :: eq_g(n_var,n_gauss,n_gauss), eq_s(n_var,n_gauss,n_gauss), eq_t(n_var,n_gauss,n_gauss)
real*8     :: density_eq(n_gauss,n_gauss)

integer    :: i, j, k, in, ms, mt, iv, inode, ife, n_elements
real*8     :: W_kin(n_tor), W_mag(n_tor), xjac, BigR, wst
real*8     :: ps0_x, ps0_y, u0_x, u0_y

W_mag = 0.d0
W_kin = 0.d0

do ife =1,  element_list%n_elements

  element = element_list%element(ife)

  do iv = 1, n_vertex_max
    inode     = element%vertex(iv)
    nodes(iv) = node_list%node(inode)
  enddo

  x_g(:,:) = 0.d0;    x_s(:,:) = 0.d0;    x_t(:,:) = 0.d0;
  y_g(:,:) = 0.d0;    y_s(:,:) = 0.;      y_t(:,:) = 0.d0;
  eq_g(:,:,:) = 0.d0; eq_s(:,:,:) = 0.d0; eq_t(:,:,:) = 0.d0;

  do i=1,n_vertex_max
    do j=1,n_order+1
      do ms=1, n_gauss
        do mt=1, n_gauss

          x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H(i,j,ms,mt)
          y_g(ms,mt) = y_g(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H(i,j,ms,mt)

          x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_t(i,j,ms,mt)
          y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_t(i,j,ms,mt)

        enddo
      enddo
    enddo
  enddo

  do in=1,n_tor

    eq_g(:,:,:) = 0.d0; eq_s(:,:,:) = 0.d0; eq_t(:,:,:) = 0.d0;

    do i=1,n_vertex_max
      do j=1,n_order+1
        do ms=1, n_gauss
          do mt=1, n_gauss

            do k=1,n_var
              eq_g(k,ms,mt)  = eq_g(k,ms,mt)  + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)
              eq_s(k,ms,mt)  = eq_s(k,ms,mt)  + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)
              eq_t(k,ms,mt)  = eq_t(k,ms,mt)  + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)
            enddo
        
	    if (in .eq. 1) then
              density_eq(ms,mt) = abs(eq_g(5,ms,mt))
            endif

          enddo
        enddo
      enddo
    enddo

!--------------------------------------------------- sum over the Gaussian integration points

    do ms=1, n_gauss

      do mt=1, n_gauss

        wst = wgauss(ms)*wgauss(mt)

        xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        BigR = x_g(ms,mt)

        ps0_x = (   y_t(ms,mt) * eq_s(1,ms,mt) - y_s(ms,mt) * eq_t(1,ms,mt) ) / xjac
        ps0_y = ( - x_t(ms,mt) * eq_s(1,ms,mt) + x_s(ms,mt) * eq_t(1,ms,mt) ) / xjac
        u0_x  = (   y_t(ms,mt) * eq_s(2,ms,mt) - y_s(ms,mt) * eq_t(2,ms,mt) ) / xjac
        u0_y  = ( - x_t(ms,mt) * eq_s(2,ms,mt) + x_s(ms,mt) * eq_t(2,ms,mt) ) / xjac

        W_mag(in) = W_mag(in) +                     (ps0_x*ps0_x + ps0_y*ps0_y ) / BigR    * xjac * wst
        W_kin(in) = W_kin(in) + density_eq(ms,mt) * (u0_x*u0_x   + u0_y*u0_y)    * BigR**3 * xjac * wst

      enddo
    enddo
  enddo

enddo

do in=1,n_tor
  if (mode(in) .ne. 0) then
    W_mag(in) = 0.5d0 * W_mag(in)
    W_kin(in) = 0.5d0 * W_kin(in)
  endif
enddo

return
end
