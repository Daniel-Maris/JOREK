module mod_chi
  use mod_semianalytical
  use mod_parameters
  use phys_module, only: domm, dcoef, F0, R_geo
  implicit none
  
  integer, parameter :: m_tor = (n_coord_tor - 1)/2
  
  type type_Cfunc
    real*8,  dimension(:), allocatable :: coef, lcoef
    integer, dimension(:), allocatable :: pwr, lpwr
  end type type_Cfunc
  
  type type_Bfunc
    real*8,  dimension(:), allocatable :: coef
    integer, dimension(:), allocatable :: zpwr
    type(type_Cfunc), dimension(:), allocatable :: rfunc
  end type type_Bfunc
  
  type(type_Cfunc), dimension(0:l_pol_domm,0:m_tor), private :: CD, CN
  type(type_Bfunc), dimension(0:n_order-1,0:n_order-1,0:l_pol_domm,0:m_tor), private :: D, N
  
  contains
  
  subroutine init_chi_basis()
    implicit none
    integer :: i, m, l, k, j, i_ord, j_ord, clsz, csz
    real*8,  dimension(:), allocatable :: tmpcoef
    integer, dimension(:), allocatable :: tmppwr
    
    ! Construct the C^D_{m,l}(R) and C^N_{m,l}(R) functions
    do i=0,m_tor
      m = i*n_coord_period
      do l=0,l_pol_domm
        allocate(CD(l,i)%coef(2*(l+1))); allocate(CD(l,i)%lcoef(l+1))
        allocate(CD(l,i)%pwr(2*(l+1)));  allocate(CD(l,i)%lpwr(l+1))
        allocate(CN(l,i)%coef(2*(l+1))); allocate(CN(l,i)%lcoef(l+1))
        allocate(CN(l,i)%pwr(2*(l+1)));  allocate(CN(l,i)%lpwr(l+1))
        do k=0,l
          CD(l,i)%coef(2*k+1) = -(alpha(k)*(gamma_st(l-m-k) - alpha(l-m-k)) - gamma(k)*alpha_st(l-m-k) + alpha(k)*beta_st(l-k))
          CD(l,i)%pwr(2*k+1)  = 2*k + m
          CD(l,i)%coef(2*k+2) = beta(k)*alpha_st(l-k);     CD(l,i)%pwr(2*k+2) = 2*k - m
          CD(l,i)%lcoef(k+1)  = -alpha(k)*alpha_st(l-m-k); CD(l,i)%lpwr(k+1)  = 2*k + m
          
          CN(l,i)%coef(2*k+1) = alpha(k)*gamma(l-m-k) - gamma(k)*alpha(l-m-k) + alpha(k)*beta(l-k); CN(l,i)%pwr(2*k+1)  = 2*k + m
          CN(l,i)%coef(2*k+2) = -beta(k)*alpha(l-k);                                                CN(l,i)%pwr(2*k+2)  = 2*k - m
          CN(l,i)%lcoef(k+1)  = alpha(k)*alpha(l-m-k);                                              CN(l,i)%lpwr(k+1)   = 2*k + m
          
!          CD(i,l) = CD(i,l) - (alpha(k)*(alpha_st(l-m-k)*lnR + gamma_st(l-m-k) - alpha(l-m-k)) - gamma(k)*alpha_st(l-m-k) &
!                  + alpha(k)*beta_st(l-k))*R**(2*k+m) + beta(k)*alpha_st(l-k)*R**(2*k-m)
!          CN(i,l) = CN(i,l) + (alpha(k)*(alpha(l-m-k)*lnR + gamma(l-m-k)) - gamma(k)*alpha(l-m-k) + alpha(k)*beta(l-k))*R**(2*k+m) &
!                  - beta(k)*alpha(l-k)*R**(2*k-m)
        end do
      end do
    end do
    
    ! Construct the D_{m,l}(R,z) and N_{m,l}(R,z) functions
    do i=0,m_tor
      do l=0,l_pol_domm
        allocate(D(0,0,l,i)%coef(l/2 + 1));  allocate(N(0,0,l,i)%coef(l/2 + 1))
        allocate(D(0,0,l,i)%zpwr(l/2 + 1));  allocate(N(0,0,l,i)%zpwr(l/2 + 1))
        allocate(D(0,0,l,i)%rfunc(l/2 + 1)); allocate(N(0,0,l,i)%rfunc(l/2 + 1))
        do k=0,l/2
          D(0,0,l,i)%coef(k+1)  = 1.d0/fact(l - 2*k); N(0,0,l,i)%coef(k+1)  = 1.d0/fact(l - 2*k)
          D(0,0,l,i)%zpwr(k+1)  = l - 2*k;            N(0,0,l,i)%zpwr(k+1)  = l - 2*k
          D(0,0,l,i)%rfunc(k+1) = CD(k,i);            N(0,0,l,i)%rfunc(k+1) = CN(k,i)
        end do
      end do
    end do
    
    ! Differentiate D and N with respect to R
    do i=0,m_tor
      do l=0,l_pol_domm
        do i_ord=1,n_order-1
          D(i_ord,0,l,i) = D(i_ord-1,0,l,i); N(i_ord,0,l,i) = N(i_ord-1,0,l,i)
          do k=0,l/2
            clsz = size(D(i_ord,0,l,i)%rfunc(k+1)%lcoef)
            csz  = size(D(i_ord,0,l,i)%rfunc(k+1)%coef)
            allocate(tmpcoef(clsz+csz)); allocate(tmppwr(clsz+csz))
            tmpcoef(clsz+1:clsz+csz) = D(i_ord,0,l,i)%rfunc(k+1)%coef; tmppwr(clsz+1:clsz+csz) = D(i_ord,0,l,i)%rfunc(k+1)%pwr
            call move_alloc(tmpcoef,D(i_ord,0,l,i)%rfunc(k+1)%coef); call move_alloc(tmppwr,D(i_ord,0,l,i)%rfunc(k+1)%pwr)
            do j=1,clsz
              D(i_ord,0,l,i)%rfunc(k+1)%coef(j) = D(i_ord,0,l,i)%rfunc(k+1)%lcoef(j)
              D(i_ord,0,l,i)%rfunc(k+1)%pwr(j)  = D(i_ord,0,l,i)%rfunc(k+1)%lpwr(j) - 1
              if (D(i_ord,0,l,i)%rfunc(k+1)%lpwr(j) .eq. 0) then
                D(i_ord,0,l,i)%rfunc(k+1)%lcoef(j) = 0.d0
              else
                D(i_ord,0,l,i)%rfunc(k+1)%lcoef(j) = D(i_ord,0,l,i)%rfunc(k+1)%lpwr(j) * D(i_ord,0,l,i)%rfunc(k+1)%lcoef(j)
                D(i_ord,0,l,i)%rfunc(k+1)%lpwr(j)  = D(i_ord,0,l,i)%rfunc(k+1)%lpwr(j) - 1
              end if
            end do
            do j=clsz+1,clsz+csz
              if (D(i_ord,0,l,i)%rfunc(k+1)%pwr(j) .eq. 0) then
                D(i_ord,0,l,i)%rfunc(k+1)%coef(j) = 0.d0
              else
                D(i_ord,0,l,i)%rfunc(k+1)%coef(j) = D(i_ord,0,l,i)%rfunc(k+1)%pwr(j) * D(i_ord,0,l,i)%rfunc(k+1)%coef(j)
                D(i_ord,0,l,i)%rfunc(k+1)%pwr(j)  = D(i_ord,0,l,i)%rfunc(k+1)%pwr(j) - 1
              end if
            end do
            
            clsz = size(N(i_ord,0,l,i)%rfunc(k+1)%lcoef)
            csz  = size(N(i_ord,0,l,i)%rfunc(k+1)%coef)
            allocate(tmpcoef(clsz+csz)); allocate(tmppwr(clsz+csz))
            tmpcoef(clsz+1:clsz+csz) = N(i_ord,0,l,i)%rfunc(k+1)%coef; tmppwr(clsz+1:clsz+csz) = N(i_ord,0,l,i)%rfunc(k+1)%pwr
            call move_alloc(tmpcoef,N(i_ord,0,l,i)%rfunc(k+1)%coef); call move_alloc(tmppwr,N(i_ord,0,l,i)%rfunc(k+1)%pwr)
            do j=1,clsz
              N(i_ord,0,l,i)%rfunc(k+1)%coef(j) = N(i_ord,0,l,i)%rfunc(k+1)%lcoef(j)
              N(i_ord,0,l,i)%rfunc(k+1)%pwr(j)  = N(i_ord,0,l,i)%rfunc(k+1)%lpwr(j) - 1
              if (N(i_ord,0,l,i)%rfunc(k+1)%lpwr(j) .eq. 0) then
                N(i_ord,0,l,i)%rfunc(k+1)%lcoef(j) = 0.d0
              else
                N(i_ord,0,l,i)%rfunc(k+1)%lcoef(j) = N(i_ord,0,l,i)%rfunc(k+1)%lpwr(j) * N(i_ord,0,l,i)%rfunc(k+1)%lcoef(j)
                N(i_ord,0,l,i)%rfunc(k+1)%lpwr(j)  = N(i_ord,0,l,i)%rfunc(k+1)%lpwr(j) - 1
              end if
            end do
            do j=clsz+1,clsz+csz
              if (N(i_ord,0,l,i)%rfunc(k+1)%pwr(j) .eq. 0) then
                N(i_ord,0,l,i)%rfunc(k+1)%coef(j) = 0.d0
              else
                N(i_ord,0,l,i)%rfunc(k+1)%coef(j) = N(i_ord,0,l,i)%rfunc(k+1)%pwr(j) * N(i_ord,0,l,i)%rfunc(k+1)%coef(j)
                N(i_ord,0,l,i)%rfunc(k+1)%pwr(j)  = N(i_ord,0,l,i)%rfunc(k+1)%pwr(j) - 1
              end if
            end do
          end do
        end do
      end do
    end do
    
    ! Differentiate D and N with respect to z
    do i=0,m_tor
      do l=0,l_pol_domm
        do j_ord=1,n_order-1
          do i_ord=0,n_order-1
            D(i_ord,j_ord,l,i) = D(i_ord,j_ord-1,l,i); N(i_ord,j_ord,l,i) = N(i_ord,j_ord-1,l,i)
            do k=0,l/2
              if (D(i_ord,j_ord,l,i)%zpwr(k+1) .eq. 0) then
                D(i_ord,j_ord,l,i)%coef(k+1) = 0.d0
              else
                D(i_ord,j_ord,l,i)%coef(k+1) = D(i_ord,j_ord,l,i)%zpwr(k+1) * D(i_ord,j_ord,l,i)%coef(k+1)
                D(i_ord,j_ord,l,i)%zpwr(k+1) = D(i_ord,j_ord,l,i)%zpwr(k+1) - 1
              end if
              if (N(i_ord,j_ord,l,i)%zpwr(k+1) .eq. 0) then
                N(i_ord,j_ord,l,i)%coef(k+1) = 0.d0
              else
                N(i_ord,j_ord,l,i)%coef(k+1) = N(i_ord,j_ord,l,i)%zpwr(k+1) * N(i_ord,j_ord,l,i)%coef(k+1)
                N(i_ord,j_ord,l,i)%zpwr(k+1) = N(i_ord,j_ord,l,i)%zpwr(k+1) - 1
              end if
            end do
          end do
        end do
      end do
    end do  
    
    contains
    
    pure real*8 function alpha(n)
      implicit none
      integer, intent(in) :: n
      
      if (n .lt. 0) then
        alpha = 0.0
      else
        alpha = (-1)**n/(fact(n+m)*fact(n)*2**(2*n+m))
      end if
    end function alpha
    
    pure real*8 function alpha_st(n)
      implicit none
      integer, intent(in) :: n
      
      alpha_st = (2*n + m)*alpha(n)
    end function alpha_st
    
    pure real*8 function beta(n)
      implicit none
      integer, intent(in) :: n
      integer :: pwr
      
      if (n .lt. 0 .or. n .ge. m) then
        beta = 0.0
      else
        pwr = 2*n-m+1
        beta = fact(m-n-1)/(fact(n)*(2.d0**pwr))
      end if
    end function beta
    
    pure real*8 function beta_st(n)
      implicit none
      integer, intent(in) :: n
      
      beta_st = (2*n - m)*beta(n)
    end function beta_st
    
    pure real*8 function gamma(n)
      implicit none
      integer, intent(in) :: n
      integer             :: i
      
      gamma = 0.0
      do i=1,n
        gamma = gamma + 1.0/i + 1.0/(m+i)
      end do
      gamma = gamma*alpha(n)/2.0
    end function gamma
    
    pure real*8 function gamma_st(n)
      implicit none
      integer, intent(in) :: n
      
      gamma_st = (2*n + m)*gamma(n)
    end function gamma_st
  end subroutine init_chi_basis
  
  function get_chi(R,z,phi)
    implicit none
    real*8,  intent(in) :: R, z, phi
    real*8, dimension(0:n_order-1,0:n_order-1,0:n_order-1) :: get_chi
    real*8, dimension(0:n_order-1) :: dksinmp, dkcosmp, V_ml
    real*8  :: Rn, zn, cval, D_ml, N_ml_1
    integer :: i, j, k, m, l, i_ord, j_ord, k_ord
    
    get_chi = 0.d0
    get_chi(0,0,0) = phi; get_chi(0,0,1) = 1.d0 ! Include the phi term
    
    if (domm) then
      Rn = R/R_geo
      zn = z/R_geo
      do i=0,m_tor
        m = i*n_coord_period
        do k_ord=0,n_order-1
          dksinmp(k_ord) = m**k_ord*(((1+(-1)**k_ord)/2)*sin(m*phi) + ((1-(-1)**k_ord)/2)*cos(m*phi))
          dkcosmp(k_ord) = m**k_ord*(((1+(-1)**k_ord)/2)*cos(m*phi) + ((1-(-1)**k_ord)/2)*sin(m*phi))
        end do
        do l=0,l_pol_domm
          do j_ord=0,n_order-1
            do i_ord=0,n_order-1
              D_ml = 0.d0; N_ml_1 = 0.d0
              do k=0,l/2
                cval = 0.d0
                do j=1,size(D(i_ord,j_ord,l,i)%rfunc(k+1)%lcoef)
                  cval = cval + D(i_ord,j_ord,l,i)%rfunc(k+1)%lcoef(j)*log(Rn)*Rn**D(i_ord,j_ord,l,i)%rfunc(k+1)%lpwr(j)
                end do
                do j=1,size(D(i_ord,j_ord,l,i)%rfunc(k+1)%coef)
                  cval = cval + D(i_ord,j_ord,l,i)%rfunc(k+1)%coef(j)*Rn**D(i_ord,j_ord,l,i)%rfunc(k+1)%pwr(j)
                end do
                D_ml = D_ml + D(i_ord,j_ord,l,i)%coef(k+1)*cval*zn**D(i_ord,j_ord,l,i)%zpwr(k+1)
              end do
              do k=0,(l-1)/2
                cval = 0.d0
                do j=1,size(N(i_ord,j_ord,abs(l-1),i)%rfunc(k+1)%lcoef)
                  cval = cval + N(i_ord,j_ord,abs(l-1),i)%rfunc(k+1)%lcoef(j)*log(Rn)*Rn**N(i_ord,j_ord,abs(l-1),i)%rfunc(k+1)%lpwr(j)
                end do
                do j=1,size(N(i_ord,j_ord,abs(l-1),i)%rfunc(k+1)%coef)
                  cval = cval + N(i_ord,j_ord,abs(l-1),i)%rfunc(k+1)%coef(j)*Rn**N(i_ord,j_ord,abs(l-1),i)%rfunc(k+1)%pwr(j)
                end do
                N_ml_1 = N_ml_1 + N(i_ord,j_ord,abs(l-1),i)%coef(k+1)*cval*zn**N(i_ord,j_ord,abs(l-1),i)%zpwr(k+1)
              end do
              do k_ord=0,n_order-1
                V_ml(k_ord) = (dcoef(1,l,i)*dkcosmp(k_ord) + dcoef(2,l,i)*dksinmp(k_ord))*D_ml &
                            + (dcoef(3,l,i)*dkcosmp(k_ord) + dcoef(4,l,i)*dksinmp(k_ord))*N_ml_1
              end do
              get_chi(i_ord,j_ord,:) = get_chi(i_ord,j_ord,:) + V_ml/(R_geo**(i_ord+j_ord))
            end do
          end do
        end do
      end do
    end if
    
    get_chi = F0*get_chi
  end function get_chi
  
  pure real*8 function fact(n)
    implicit none
    integer, intent(in) :: n
    integer             :: i
    
    fact = 1.0
    do i=2,n
      fact = fact*i
    end do
  end function fact
  
end module mod_chi
