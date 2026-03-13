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
 
subroutine Cal_Shape_Functions_Hex8(xi, eta, zeta, N, dN_dxi, dN_deta, dN_dzeta)
!2026-01-28.

use Global_Float_Type      
implicit none
real(kind=FT), intent(in) :: xi, eta, zeta
real(kind=FT), intent(out) :: N(8), dN_dxi(8), dN_deta(8), dN_dzeta(8)

real(kind=FT) :: xip, xim, etap, etam, zetap, zetam

xip = ONE + xi
xim = ONE - xi
etap = ONE + eta
etam = ONE - eta
zetap = ONE + zeta
zetam = ONE - zeta

N(1) = ZP125 * xim * etam * zetam
N(2) = ZP125 * xip * etam * zetam
N(3) = ZP125 * xip * etap * zetam
N(4) = ZP125 * xim * etap * zetam
N(5) = ZP125 * xim * etam * zetap
N(6) = ZP125 * xip * etam * zetap
N(7) = ZP125 * xip * etap * zetap
N(8) = ZP125 * xim * etap * zetap

dN_dxi(1) = -ZP125 * etam * zetam
dN_dxi(2) =  ZP125 * etam * zetam
dN_dxi(3) =  ZP125 * etap * zetam
dN_dxi(4) = -ZP125 * etap * zetam
dN_dxi(5) = -ZP125 * etam * zetap
dN_dxi(6) =  ZP125 * etam * zetap
dN_dxi(7) =  ZP125 * etap * zetap
dN_dxi(8) = -ZP125 * etap * zetap

dN_deta(1) = -ZP125 * xim * zetam
dN_deta(2) = -ZP125 * xip * zetam
dN_deta(3) =  ZP125 * xip * zetam
dN_deta(4) =  ZP125 * xim * zetam
dN_deta(5) = -ZP125 * xim * zetap
dN_deta(6) = -ZP125 * xip * zetap
dN_deta(7) =  ZP125 * xip * zetap
dN_deta(8) =  ZP125 * xim * zetap

dN_dzeta(1) = -ZP125 * xim * etam
dN_dzeta(2) = -ZP125 * xip * etam
dN_dzeta(3) = -ZP125 * xip * etap
dN_dzeta(4) = -ZP125 * xim * etap
dN_dzeta(5) =  ZP125 * xim * etam
dN_dzeta(6) =  ZP125 * xip * etam
dN_dzeta(7) =  ZP125 * xip * etap
dN_dzeta(8) =  ZP125 * xim * etap
    
end subroutine Cal_Shape_Functions_Hex8
