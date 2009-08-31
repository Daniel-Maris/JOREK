subroutine basisfunctions1(s,H,H_s,H_ss)
!--------------------------------------------------------------
! subroutine which defines the basis functions in one dimension
!
! derived from a mixed Bezier/Cubic finite element representation
!   index 1 : counts the vertex
!   index 2 : counts the variables (p,u,v,w)
!
! the functions are defined on the interval [0,1]
!----------------------------------------------------------------
implicit none
real*8 :: s, H(2,2), H_s(2,2), H_ss(2,2)


!---------------------------------------------------------- vertex (1)
H(1,1)   =(-1.d0 + s)**2*(1.d0 + 2.d0*s)
H_s(1,1) =6.d0*(-1.d0 + s)*s
H_ss(1,1)=6.d0*(-1.d0 + 2.d0*s)

H(1,2)   =3.d0*(-1.d0 + s)**2*s
H_s(1,2) =3.d0*(-1.d0 + s)*(-1.d0 + 3.d0*s)
H_ss(1,2)=6.d0*(-2.d0 + 3.d0*s)

!---------------------------------------------------------- vertex (2)
H(2,1)   =-s**2*(-3.d0 + 2.d0*s)
H_s(2,1) =-6.d0*(-1.d0 + s)*s
H_ss(2,1)=-6.d0*(-1.d0 + 2.d0*s)

H(2,2)   =-3.d0*(-1.d0 + s)*s**2
H_s(2,2) =-3.d0*s*(-2.d0 + 3.d0*s)
H_ss(2,2)=-6.d0*(-1.d0 + 3.d0*s)

return
end
