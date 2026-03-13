!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !

subroutine Cal_3D_SIFs_IIM_Auxiliary_Fields_Part2(Mode, r_in, theta, G, nu, kappa, &
                     Aux_dSig11_dX1, Aux_dSig12_dX2, Aux_dSig13_dX3,    &
                     Aux_dSig12_dX1, Aux_dSig22_dX2, Aux_dSig23_dX3,    &
                     Aux_dSig13_dX1, Aux_dSig23_dX2, Aux_dSig33_dX3,    &
                     Aux_dEps11_dX1, Aux_dEps12_dX1, Aux_dEps13_dX1,    &
                     Aux_dEps22_dX1, Aux_dEps23_dX1, Aux_dEps33_dX1,    &
                     Aux_ddu1_dX1dX1, Aux_ddu2_dX1dX1, Aux_ddu3_dX1dX1, &
                     Aux_ddu1_dX2dX1, Aux_ddu2_dX2dX1, Aux_ddu3_dX2dX1, & 
                     Aux_ddu1_dX3dX1, Aux_ddu2_dX3dX1, Aux_ddu3_dX3dX1) 
!2026-01-28.
use Global_Float_Type
implicit none
integer, intent(in) :: Mode
real(kind=FT), intent(in) :: r_in, theta, G, nu, kappa
real(kind=FT), intent(out) :: Aux_dSig11_dX1, Aux_dSig12_dX2, Aux_dSig13_dX3
real(kind=FT), intent(out) :: Aux_dSig12_dX1, Aux_dSig22_dX2, Aux_dSig23_dX3
real(kind=FT), intent(out) :: Aux_dSig13_dX1, Aux_dSig23_dX2, Aux_dSig33_dX3
real(kind=FT), intent(out) :: Aux_dEps11_dX1, Aux_dEps12_dX1, Aux_dEps13_dX1
real(kind=FT), intent(out) :: Aux_dEps22_dX1, Aux_dEps23_dX1, Aux_dEps33_dX1
real(kind=FT), intent(out) :: Aux_ddu1_dX1dX1, Aux_ddu2_dX1dX1, Aux_ddu3_dX1dX1
real(kind=FT), intent(out) :: Aux_ddu1_dX2dX1, Aux_ddu2_dX2dX1, Aux_ddu3_dX2dX1 
real(kind=FT), intent(out) :: Aux_ddu1_dX3dX1, Aux_ddu2_dX3dX1, Aux_ddu3_dX3dX1
real(kind=FT) :: sqrt_r, sqrt_2pir, K_aux
real(kind=FT) :: cos_t, sin_t, cos_t2, sin_t2, cos_3t2, sin_3t2
real(kind=FT) :: fac_stress, fac_strain
real(kind=FT) :: drdx1, drdx2, drdx3, dtdx1, dtdx2, dtdx3
real(kind=FT) :: du1dr, du1dt, du2dr, du2dt, du3dr, du3dt
real(kind=FT) :: dsigma11dr,dsigma12dr,dsigma13dr
real(kind=FT) :: dsigma22dr,dsigma23dr
real(kind=FT) :: dsigma11dt,dsigma12dt,dsigma13dt
real(kind=FT) :: dsigma22dt,dsigma23dt
real(kind=FT) :: Aux_dSig11_dX3,Aux_dSig22_dX3
real(kind=FT) :: factorA,factor32,factor34,factorC,factor14,r32,r12
real(kind=FT) :: Aux_Sig(3,3)
real(kind=FT) :: ddrdx1dx1,ddtdx1dx1,ddrdx1dx2,ddtdx1dx2,ddrdx1dx3,ddtdx1dx3
real(kind=FT) :: ddu1drdr,ddu1drdt,ddu1dtdt,ddu2drdr,ddu2drdt,ddu2dtdt,ddu3drdr,ddu3drdt,ddu3dtdt
real(kind=FT) :: ddrdx3dx1,ddtdx3dx1,ddrdx2dx1,ddtdx2dx1
real(kind=FT) :: r

Aux_dSig11_dX1 = ZR; Aux_dSig12_dX2 = ZR; Aux_dSig13_dX3 = ZR
Aux_dSig12_dX1 = ZR; Aux_dSig22_dX2 = ZR; Aux_dSig23_dX3 = ZR
Aux_dSig13_dX1 = ZR; Aux_dSig23_dX2 = ZR; Aux_dSig33_dX3 = ZR
Aux_dEps11_dX1 = ZR; Aux_dEps12_dX1 = ZR; Aux_dEps13_dX1 = ZR
Aux_dEps22_dX1 = ZR; Aux_dEps23_dX1 = ZR; Aux_dEps33_dX1 = ZR
Aux_ddu1_dX1dX1 = ZR; Aux_ddu2_dX1dX1 = ZR; Aux_ddu3_dX1dX1 = ZR
Aux_ddu1_dX2dX1 = ZR; Aux_ddu2_dX2dX1 = ZR; Aux_ddu3_dX2dX1 = ZR
Aux_ddu1_dX3dX1 = ZR; Aux_ddu2_dX3dX1 = ZR; Aux_ddu3_dX3dX1 = ZR

r = r_in
if (r <=1.0D-12) then
    r = 1.0D-12
endif

K_aux = ONE

sqrt_r    = sqrt(r)
sqrt_2pir = sqrt(TWO * pi * r)
factorA   = ONE/sqrt_2pir
factor32  = THR/TWO  
factor34  = THR/FOU 
factor14  = ONE/FOU 
factorC   = ONE/(TWO*G*sqrt(TWO*pi))
r32       = r**(-THR/TWO)
r12       = r**(-ONE/TWO)

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
drdx3 = ZR
dtdx1 = -sin_t / r
dtdx2 = cos_t / r
dtdx3 = ZR

ddrdx1dx1 =  sin_t**2/r
ddtdx1dx1 =  sin(TWO*theta)/(r**2)
ddrdx1dx2 = -sin_t*cos_t/r
ddtdx1dx2 = -cos(TWO*theta)/(r**2)
ddrdx2dx1 =  ddrdx1dx2
ddtdx2dx1 =  ddtdx1dx2
ddrdx1dx3 =  ZR
ddtdx1dx3 =  ZR
ddrdx3dx1 =  ddrdx1dx3
ddtdx3dx1 =  ddtdx1dx3
    
Aux_Sig = ZR


if(Mode == 1) then

    
    Aux_Sig(1,1) = fac_stress * cos_t2 * (ONE - sin_t2 * sin_3t2)
    Aux_Sig(2,2) = fac_stress * cos_t2 * (ONE + sin_t2 * sin_3t2)
    Aux_Sig(1,2) = fac_stress * sin_t2 * cos_t2 * cos_3t2
    Aux_Sig(2,1) = Aux_Sig(1,2)
    Aux_Sig(3,3) = nu * (Aux_Sig(1,1) + Aux_Sig(2,2))
    
    dsigma11dr = -HLF*Aux_Sig(1,1)/r
    dsigma12dr = -HLF*Aux_Sig(1,2)/r
    dsigma13dr =  ZR
    dsigma22dr = -HLF*Aux_Sig(2,2)/r
    dsigma23dr =  ZR
    dsigma11dt = factorA*(-HLF*sin_t2*(ONE-sin_t2*sin_3t2) - cos_t2*(HLF*cos_t2*sin_3t2 + factor32*sin_t2*cos_3t2))
    dsigma12dt = factorA*(HLF*cos_t*cos_3t2 - factor34*sin_t*sin_3t2)
    dsigma13dt =  ZR
    dsigma22dt = factorA*(-HLF*sin_t2*(ONE+sin_t2*sin_3t2) + cos_t2*(HLF*cos_t2*sin_3t2 + factor32*sin_t2*cos_3t2))
    dsigma23dt =  ZR
    
    Aux_dSig11_dX1 =  dsigma11dr * drdx1 + dsigma11dt * dtdx1
    Aux_dSig12_dX2 =  dsigma12dr * drdx2 + dsigma12dt * dtdx2
    Aux_dSig13_dX3 =  dsigma13dr * drdx3 + dsigma13dt * dtdx3
    Aux_dSig12_dX1 =  dsigma12dr * drdx1 + dsigma12dt * dtdx1
    Aux_dSig22_dX2 =  dsigma22dr * drdx2 + dsigma22dt * dtdx2
    Aux_dSig23_dX3 =  dsigma23dr * drdx3 + dsigma23dt * dtdx3
    Aux_dSig13_dX1 =  dsigma13dr * drdx1 + dsigma13dt * dtdx1
    Aux_dSig23_dX2 =  dsigma23dr * drdx2 + dsigma23dt * dtdx2
    Aux_dSig11_dX3 =  dsigma11dr * drdx3 + dsigma11dt * dtdx3
    Aux_dSig22_dX3 =  dsigma22dr * drdx3 + dsigma22dt * dtdx3
    Aux_dSig33_dX3 =  nu*(Aux_dSig11_dX3 + Aux_dSig22_dX3)
    
    du1dr = (ONE/(FOU*G))*(ONE/sqrt_2pir)*cos_t2*(kappa - cos_t)
    du1dt = fac_strain * (-HLF * sin_t2 * (kappa - cos_t) + cos_t2 * sin_t)
    du2dr = (ONE/(FOU*G))*(ONE/sqrt_2pir)*sin_t2*(kappa - cos_t)
    du2dt = fac_strain * (HLF * cos_t2 * (kappa - cos_t) + sin_t2 * sin_t)

    ddu1drdr = -factorC/FOU*r32*cos_t2*(kappa - cos_t) 
    ddu1drdt =  factorC/TWO*r12*(-HLF*sin_t2*(kappa - cos_t) + cos_t2*sin_t) 
    ddu1dtdt =  factorC*sqrt(r)*(-factor14*cos_t2*(kappa - cos_t) - sin_t2*sin_t + cos_t2*cos_t)  
    ddu2drdr = -factorC/FOU*r32*sin_t2*(kappa - cos_t) 
    ddu2drdt =  factorC/TWO*r12*(HLF*cos_t2*(kappa - cos_t) + sin_t2*sin_t) 
    ddu2dtdt =  factorC*sqrt(r)*(-factor14*sin_t2*(kappa - cos_t) + cos_t2*sin_t +sin_t2*cos_t) 
    ddu3drdr =  ZR
    ddu3drdt =  ZR
    ddu3dtdt =  ZR
    
elseif(Mode == 2) then
    Aux_Sig(1,1) = -fac_stress * sin_t2 * (TWO + cos_t2 * cos_3t2)
    Aux_Sig(2,2) = fac_stress * sin_t2 * cos_t2 * cos_3t2
    Aux_Sig(1,2) = fac_stress * cos_t2 * (ONE - sin_t2 * sin_3t2)
    Aux_Sig(2,1) = Aux_Sig(1,2)
    Aux_Sig(3,3) = nu * (Aux_Sig(1,1) + Aux_Sig(2,2))
    
    dsigma11dr = -HLF*Aux_Sig(1,1)/r
    dsigma12dr = -HLF*Aux_Sig(1,2)/r
    dsigma13dr =  ZR
    dsigma22dr = -HLF*Aux_Sig(2,2)/r
    dsigma23dr =  ZR
    dsigma11dt = factorA*(-HLF*cos_t2*(TWO+cos_t2*cos_3t2) - sin_t2*(-HLF*sin_t2*cos_3t2 - factor32*cos_t2*sin_3t2))
    dsigma12dt = factorA*(-HLF*sin_t2*(ONE-sin_t2*sin_3t2) + cos_t2*(-HLF*cos_t2*sin_3t2 - factor32*sin_t2*cos_3t2))
    dsigma13dt =  ZR
    dsigma22dt = factorA*(HLF*cos_t*cos_3t2 - factor32*sin_t2*cos_t2*sin_3t2)
    dsigma23dt =  ZR
    
    Aux_dSig11_dX1 =  dsigma11dr * drdx1 + dsigma11dt * dtdx1
    Aux_dSig12_dX2 =  dsigma12dr * drdx2 + dsigma12dt * dtdx2
    Aux_dSig13_dX3 =  dsigma13dr * drdx3 + dsigma13dt * dtdx3
    Aux_dSig12_dX1 =  dsigma12dr * drdx1 + dsigma12dt * dtdx1
    Aux_dSig22_dX2 =  dsigma22dr * drdx2 + dsigma22dt * dtdx2
    Aux_dSig23_dX3 =  dsigma23dr * drdx3 + dsigma23dt * dtdx3
    Aux_dSig13_dX1 =  dsigma13dr * drdx1 + dsigma13dt * dtdx1
    Aux_dSig23_dX2 =  dsigma23dr * drdx2 + dsigma23dt * dtdx2
    Aux_dSig11_dX3 =  dsigma11dr * drdx3 + dsigma11dt * dtdx3
    Aux_dSig22_dX3 =  dsigma22dr * drdx3 + dsigma22dt * dtdx3
    Aux_dSig33_dX3 =  nu*(Aux_dSig11_dX3 + Aux_dSig22_dX3)
    
    du1dr = (ONE/(FOU*G))*(ONE/sqrt_2pir)*sin_t2*(kappa + TWO + cos_t)
    du1dt = fac_strain * (HLF * cos_t2 * (kappa + TWO + cos_t) - sin_t2 * sin_t)
    du2dr = -(ONE/(FOU*G))*(ONE/sqrt_2pir)*cos_t2*(kappa - TWO + cos_t)
    du2dt = fac_strain * (HLF * sin_t2 * (kappa - TWO + cos_t) + cos_t2 * sin_t)
    
    ddu1drdr = -factorC/FOU*r32*sin_t2*(kappa + TWO + cos_t) 
    ddu1drdt =  factorC/TWO*r12*(HLF*cos_t2*(kappa + TWO + cos_t)  - sin_t2*sin_t)  
    ddu1dtdt =  factorC*sqrt(r)*(-factor14*sin_t2*(kappa + TWO + cos_t) - cos_t2*sin_t - sin_t2*cos_t) 
    ddu2drdr =  factorC/FOU*r32*cos_t2*(kappa - TWO + cos_t) 
    ddu2drdt =  factorC/TWO*r12*(HLF*sin_t2*(kappa -TWO + cos_t) + cos_t2*sin_t) 
    ddu2dtdt =  factorC*sqrt(r)*(factor14*cos_t2*(kappa - TWO + cos_t) - sin_t2*sin_t + cos_t2*cos_t) 
    ddu3drdr =  ZR
    ddu3drdt =  ZR
    ddu3dtdt =  ZR
    
else
    Aux_Sig(1,3) = -fac_stress * sin_t2
    Aux_Sig(3,1) = Aux_Sig(1,3)
    Aux_Sig(2,3) = fac_stress * cos_t2
    Aux_Sig(3,2) = Aux_Sig(2,3)
    
    
    Aux_dSig11_dX1 =  ZR
    Aux_dSig12_dX2 =  ZR
    Aux_dSig13_dX3 =  ZR
    Aux_dSig12_dX1 =  ZR
    Aux_dSig22_dX2 =  ZR
    Aux_dSig23_dX3 =  ZR
    Aux_dSig13_dX1 = -HLF*Aux_Sig(1,3)/r*cos_t + factorA/(TWO*r)*cos_t2*sin_t
    Aux_dSig23_dX2 = -HLF*Aux_Sig(2,3)/r*sin_t - factorA/(TWO*r)*sin_t2*cos_t
    Aux_dSig11_dX3 =  ZR
    Aux_dSig22_dX3 =  ZR
    Aux_dSig33_dX3 =  nu*(Aux_dSig11_dX3 + Aux_dSig22_dX3)
    
    
    du1dr = ZR
    du1dt = ZR
    du2dr = ZR
    du2dt = ZR
    du3dr = (ONE/G) *sqrt(TWO/pi)* HLF / sqrt_r * sin_t2
    du3dt = (ONE/G) *sqrt(TWO*r/pi)* HLF * cos_t2
    
    ddu1drdr =  ZR
    ddu1drdt =  ZR
    ddu1dtdt =  ZR
    ddu2drdr =  ZR
    ddu2drdt =  ZR
    ddu2dtdt =  ZR
    ddu3drdr =  -ONE/(FOU*G)*sqrt(TWO/pi)*r32*sin_t2
    ddu3drdt =   ONE/(FOU*G)*sqrt(TWO/pi)*r12*cos_t2
    ddu3dtdt =  -ONE/(FOU*G)*sqrt(TWO/pi)*sqrt(r)*sin_t2
endif

Aux_ddu1_dX1dX1 = ddu1drdr*(drdx1)**2 + TWO*ddu1drdt*drdx1*dtdx1 + ddu1dtdt*(dtdx1)**2 + &
                  du1dr*ddrdx1dx1 + du1dt*ddtdx1dx1 
Aux_ddu2_dX1dX1 = ddu2drdr*(drdx1)**2 + TWO*ddu2drdt*drdx1*dtdx1 + ddu2dtdt*(dtdx1)**2 + &
                  du2dr*ddrdx1dx1 + du2dt*ddtdx1dx1
Aux_ddu3_dX1dX1 = ddu3drdr*(drdx1)**2 + TWO*ddu3drdt*drdx1*dtdx1 + ddu3dtdt*(dtdx1)**2 + &
                  du3dr*ddrdx1dx1 + du3dt*ddtdx1dx1
Aux_ddu1_dX2dX1= ddu1drdr*drdx2*drdx1 + ddu1drdt*(drdx2*dtdx1+dtdx2*drdx1) + &
                 ddu1dtdt*dtdx2*dtdx1 + du1dr*ddrdx2dx1 + du1dt*ddtdx2dx1 
Aux_ddu2_dX2dX1= ddu2drdr*drdx2*drdx1 + ddu2drdt*(drdx2*dtdx1+dtdx2*drdx1) + &
                 ddu2dtdt*dtdx2*dtdx1 + du2dr*ddrdx2dx1 + du2dt*ddtdx2dx1 
Aux_ddu3_dX2dX1= ddu3drdr*drdx2*drdx1 + ddu3drdt*(drdx2*dtdx1+dtdx2*drdx1) + &
                 ddu3dtdt*dtdx2*dtdx1 + du3dr*ddrdx2dx1 + du3dt*ddtdx2dx1 
Aux_ddu1_dX3dX1= ddu1drdr*drdx3*drdx1 + ddu1drdt*(drdx3*dtdx1+dtdx3*drdx1) + &
                 ddu1dtdt*dtdx3*dtdx1 + du1dr*ddrdx3dx1 + du1dt*ddtdx3dx1 
Aux_ddu2_dX3dX1= ddu2drdr*drdx3*drdx1 + ddu2drdt*(drdx3*dtdx1+dtdx3*drdx1) + &
                 ddu2dtdt*dtdx3*dtdx1 + du2dr*ddrdx3dx1 + du2dt*ddtdx3dx1 
Aux_ddu3_dX3dX1= ddu3drdr*drdx3*drdx1 + ddu3drdt*(drdx3*dtdx1+dtdx3*drdx1) + &
                 ddu3dtdt*dtdx3*dtdx1 + du3dr*ddrdx3dx1 + du3dt*ddtdx3dx1 

Aux_dEps11_dX1 = Aux_ddu1_dX1dX1
Aux_dEps12_dX1 = HLF*(Aux_ddu1_dX2dX1 + Aux_ddu2_dX1dX1)
Aux_dEps13_dX1 = HLF*(Aux_ddu1_dX3dX1 + Aux_ddu3_dX1dX1)
Aux_dEps22_dX1 = Aux_ddu2_dX2dX1
Aux_dEps23_dX1 = HLF*(Aux_ddu2_dX3dX1 + Aux_ddu3_dX2dX1)
Aux_dEps33_dX1 = Aux_ddu3_dX3dX1

    
end subroutine Cal_3D_SIFs_IIM_Auxiliary_Fields_Part2



