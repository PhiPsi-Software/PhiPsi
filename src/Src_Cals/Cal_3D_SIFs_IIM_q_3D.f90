!-----------------------------------------------------------
! Brief: Evaluate the 3D weight function q for IIM domain.
!
! Parameters:
!   Input:  r     - radial distance from the crack front
!           s     - coordinate along the crack front
!           R_in  - inner radius where q=1
!           R_out - outer radius where q=0
!           L_f   - half-length of front core region
!   Output: q_val - product of radial and front weight factors
!
! Notes:   q is a separable piecewise-linear function of r and s.
!-----------------------------------------------------------

subroutine Cal_3D_SIFs_IIM_q_3D(r, s, R_in, R_out, L_f, q_val)
! Compute 3D q-function.
! 2026-01-28.

use Global_Float_Type

implicit none
real(kind=FT), intent(in) :: r, s, R_in, R_out, L_f
real(kind=FT), intent(out) :: q_val

real(kind=FT) :: q_radial, q_front

if (r <= R_in) then
    q_radial = ONE
else if (r >= R_out) then
    q_radial = ZR
else
    q_radial = (R_out - r) / (R_out - R_in)
endif

if (abs(s) <= L_f) then
    q_front = ONE
else if (abs(s) >= TWO*L_f) then
    q_front = ZR
else
    q_front = (TWO*L_f - abs(s)) / L_f
endif

q_val = q_radial * q_front
    
end subroutine Cal_3D_SIFs_IIM_q_3D


