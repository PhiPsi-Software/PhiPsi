!-----------------------------------------------------------
! Brief: Evaluate the smooth-tip enrichment function R and its
!        spatial derivatives R_x, R_y for a given parent coord.
!
! Parameters:
!   Input:  lambda   - transition width parameter
!           Le       - element characteristic length
!           cos_omega, sin_omega - crack direction cosines
!           kesi     - parent coordinate along crack
!   Output: R        - value of smooth enrichment function
!           R_x, R_y - spatial derivatives of R
!
! Notes:   Uses a quintic smooth-step blending; returns early
!          on degenerate inputs (lambda, kesi, Le out of range).
!-----------------------------------------------------------

subroutine Cal_R_and_R_x_R_y(lambda,Le,cos_omega,sin_omega,kesi,R,R_x,R_y)
!2026-05-26.

use Global_Float_Type
use Global_Common

implicit none
real(kind=FT),intent(in) ::lambda,kesi,Le,cos_omega,sin_omega
real(kind=FT),intent(out)::R,R_x,R_y
real(kind=FT) yita,S,kesi_x,kesi_y,d_S
      
R   = ONE
R_x = ZR
R_y = ZR

if (lambda < ZR) then
    print *, '    ERROR-2026052601 :: lambda should be larger than 0!'
    print *, '                        ERROR in Cal_R_and_R_x_R_y.f90!'
    call Warning_Message('S',Keywords_Blank)
endif

if (kesi < ZR) then
    R = lambda
    return
endif

if (abs(Le) <= Tol_10) then
    print *, '    ERROR-2026052603 :: The value of Le is too small!'
    print *, '                        ERROR in Cal_R_and_R_x_R_y.f90!'
    print *, '                        Le: ',Le
    call Warning_Message('S',Keywords_Blank)
endif

if (Le < ZR) then
    print *, '    ERROR-2026052604 :: Le should be larger than 0!'
    print *, '                        ERROR in Cal_R_and_R_x_R_y.f90!'
    print *, '                        Le: ',Le
    call Warning_Message('S',Keywords_Blank)
endif

if (abs(lambda) <= Tol_10) then
    R   = ZR
    R_x = ZR
    R_y = ZR
    return
endif

if (kesi > lambda) then
    R   = ZR
    R_x = ZR
    R_y = ZR
    return
endif

yita = (lambda - kesi)/lambda
S    = 6.0D0*yita**5 - 15.0D0*yita**4 + 10.0D0*yita**3
d_S  = 30.0D0*yita**2*(1.0D0 - yita)**2
kesi_x = cos_omega/Le
kesi_y = sin_omega/Le

R   = lambda + (ONE - lambda)*S
R_x = -((ONE - lambda)/lambda)*d_S*kesi_x
R_y = -((ONE - lambda)/lambda)*d_S*kesi_y




end subroutine Cal_R_and_R_x_R_y