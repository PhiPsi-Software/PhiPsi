

subroutine Cal_3D_SIFs_IIM_Auxiliary_Fields_Part1(Mode, r, theta, G, nu, kappa, Aux_Sig, Aux_Eps, Aux_Grad_u)

use Global_Float_Type
implicit none

integer, intent(in) :: Mode
real(kind=FT), intent(in) :: r, theta, G, nu, kappa
real(kind=FT), intent(out) :: Aux_Sig(3,3), Aux_Eps(3,3), Aux_Grad_u(3,3)
real(kind=FT) :: sqrt_r, sqrt_2pir, K_aux
real(kind=FT) :: cos_t, sin_t, cos_t2, sin_t2, cos_3t2, sin_3t2
real(kind=FT) :: fac_stress, fac_strain
real(kind=FT) :: drdx1, drdx2, dtdx1, dtdx2
real(kind=FT) :: du1dr, du1dt, du2dr, du2dt, du3dr, du3dt

K_aux = ONE
sqrt_r = sqrt(r)
sqrt_2pir = sqrt(TWO * pi * r)
cos_t = cos(theta)
sin_t = sin(theta)
cos_t2 = cos(theta / TWO)
sin_t2 = sin(theta / TWO)
cos_3t2 = cos(THR * theta / TWO)
sin_3t2 = sin(THR * theta / TWO)

fac_stress = K_aux / sqrt_2pir
fac_strain = K_aux / (TWO * G) * sqrt(r / (TWO * pi))

drdx1 = cos_t
drdx2 = sin_t
dtdx1 = -sin_t / r
dtdx2 = cos_t / r

Aux_Sig = ZR
Aux_Eps = ZR
Aux_Grad_u = ZR

if(Mode == 1) then
    Aux_Sig(1,1) = fac_stress * cos_t2 * (ONE - sin_t2 * sin_3t2)
    Aux_Sig(2,2) = fac_stress * cos_t2 * (ONE + sin_t2 * sin_3t2)
    Aux_Sig(1,2) = fac_stress * sin_t2 * cos_t2 * cos_3t2
    Aux_Sig(2,1) = Aux_Sig(1,2)
    Aux_Sig(3,3) = nu * (Aux_Sig(1,1) + Aux_Sig(2,2))
    du1dr = (ONE/(FOU*G))*(ONE/sqrt_2pir)*cos_t2*(kappa - cos_t)
    du1dt = fac_strain * (-HLF * sin_t2 * (kappa - cos_t) + cos_t2 * sin_t)
    du2dr = (ONE/(FOU*G))*(ONE/sqrt_2pir)*sin_t2*(kappa - cos_t)
    du2dt = fac_strain * (HLF * cos_t2 * (kappa - cos_t) + sin_t2 * sin_t)
    
    Aux_Grad_u(1,1) = du1dr * drdx1 + du1dt * dtdx1
    Aux_Grad_u(1,2) = du1dr * drdx2 + du1dt * dtdx2
    Aux_Grad_u(2,1) = du2dr * drdx1 + du2dt * dtdx1
    Aux_Grad_u(2,2) = du2dr * drdx2 + du2dt * dtdx2
    
    Aux_Eps(1,1) = Aux_Grad_u(1,1)
    Aux_Eps(2,2) = Aux_Grad_u(2,2)
    Aux_Eps(1,2) = HLF * (Aux_Grad_u(1,2) + Aux_Grad_u(2,1))
    Aux_Eps(2,1) = Aux_Eps(1,2)
    Aux_Eps(3,3) = ZR
    
elseif(Mode == 2) then
    Aux_Sig(1,1) = -fac_stress * sin_t2 * (TWO + cos_t2 * cos_3t2)
    Aux_Sig(2,2) = fac_stress * sin_t2 * cos_t2 * cos_3t2
    Aux_Sig(1,2) = fac_stress * cos_t2 * (ONE - sin_t2 * sin_3t2)
    Aux_Sig(2,1) = Aux_Sig(1,2)
    Aux_Sig(3,3) = nu * (Aux_Sig(1,1) + Aux_Sig(2,2))
    
    du1dr = (ONE/(FOU*G))*(ONE/sqrt_2pir)*sin_t2*(kappa + TWO + cos_t)
    du1dt = fac_strain * (HLF * cos_t2 * (kappa + TWO + cos_t) - sin_t2 * sin_t)
    du2dr = -(ONE/(FOU*G))*(ONE/sqrt_2pir)*cos_t2*(kappa - TWO + cos_t)
    du2dt = fac_strain * (HLF * sin_t2 * (kappa - TWO + cos_t) + cos_t2 * sin_t)
    
    Aux_Grad_u(1,1) = du1dr * drdx1 + du1dt * dtdx1
    Aux_Grad_u(1,2) = du1dr * drdx2 + du1dt * dtdx2
    Aux_Grad_u(2,1) = du2dr * drdx1 + du2dt * dtdx1
    Aux_Grad_u(2,2) = du2dr * drdx2 + du2dt * dtdx2
    
    Aux_Eps(1,1) = Aux_Grad_u(1,1)
    Aux_Eps(2,2) = Aux_Grad_u(2,2)
    Aux_Eps(1,2) = HLF * (Aux_Grad_u(1,2) + Aux_Grad_u(2,1))
    Aux_Eps(2,1) = Aux_Eps(1,2)
    Aux_Eps(3,3) = ZR
else
    Aux_Sig(1,3) = -fac_stress * sin_t2
    Aux_Sig(3,1) = Aux_Sig(1,3)
    Aux_Sig(2,3) = fac_stress * cos_t2
    Aux_Sig(3,2) = Aux_Sig(2,3)
    
    du3dr = (ONE/G) *sqrt(TWO/pi)* HLF / sqrt_r * sin_t2
    du3dt = (ONE/G) *sqrt(TWO*r/pi)* HLF * cos_t2
    
    Aux_Grad_u(3,1) = du3dr * drdx1 + du3dt * dtdx1
    Aux_Grad_u(3,2) = du3dr * drdx2 + du3dt * dtdx2
    
    Aux_Eps(1,3) = HLF * Aux_Grad_u(3,1)
    Aux_Eps(3,1) = Aux_Eps(1,3)
    Aux_Eps(2,3) = HLF * Aux_Grad_u(3,2)
    Aux_Eps(3,2) = Aux_Eps(2,3)
endif
    
end subroutine Cal_3D_SIFs_IIM_Auxiliary_Fields_Part1



