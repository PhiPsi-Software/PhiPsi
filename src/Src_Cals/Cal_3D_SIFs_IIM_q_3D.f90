

subroutine Cal_3D_SIFs_IIM_q_3D(r, s, R_in, R_out, L_f, q_val)

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


