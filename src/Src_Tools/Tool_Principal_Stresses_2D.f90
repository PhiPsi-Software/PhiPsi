 
subroutine Tool_Principal_Stresses_2D(S_xx,S_yy,S_xy,S_1,S_3,theta_stress)
use Global_Float_Type
use Global_Common

implicit none
real(kind=FT),intent(in)::S_xx,S_yy,S_xy
real(kind=FT),intent(out)::S_1,S_3,theta_stress
real(kind=FT) Sxx_minus_Syy

S_1=(S_xx+S_yy)/TWO+sqrt(((S_xx-S_yy)/TWO)**2+S_xy**2)
S_3=(S_xx+S_yy)/TWO-sqrt(((S_xx-S_yy)/TWO)**2+S_xy**2)

if(S_xx >= S_yy)then
    if(S_xx==S_yy) then
        Sxx_minus_Syy = Tol_15
    else
        Sxx_minus_Syy  = S_xx - S_yy 
    endif
    theta_stress = atan(TWO*S_xy/Sxx_minus_Syy)/TWO+pi/TWO
elseif(S_xx < S_yy)then
    theta_stress = atan(TWO*S_xy/(S_xx-S_yy))/TWO
end if

RETURN
end subroutine Tool_Principal_Stresses_2D