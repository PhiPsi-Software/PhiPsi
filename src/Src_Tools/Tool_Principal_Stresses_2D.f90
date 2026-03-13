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
 
subroutine Tool_Principal_Stresses_2D(S_xx,S_yy,S_xy,S_1,S_3,theta_stress)
! This function calculates the principal stress.
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