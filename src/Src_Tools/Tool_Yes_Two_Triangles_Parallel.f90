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
 
subroutine Tool_Yes_Two_Triangles_Parallel(Tri_1,Tri_2,Yes_Parallel)   
! Determine whether two 3D triangles are parallel by checking the direction of their outward
! normals.
!2022-08-05.

!......................
! Variable Declaration
!......................
use Global_Float_Type
use Global_Common
implicit none
real(kind=FT),intent(in)::Tri_1(3,3),Tri_2(3,3)
logical,intent(out)::Yes_Parallel            
real(kind=FT) Tri_1_n_Vector(3), Tri_2_n_Vector(3)
real(kind=FT) Tri_1_a_Vector(3), Tri_1_b_Vector(3)
real(kind=FT) Tri_2_a_Vector(3), Tri_2_b_Vector(3)
real(kind=FT) c_angle
Yes_Parallel  = .False.
     
Tri_1_a_Vector(1:3) = Tri_1(2,1:3) - Tri_1(1,1:3)
Tri_1_b_Vector(1:3) = Tri_1(3,1:3) - Tri_1(1,1:3)
call Vector_Cross_Product_3(Tri_1_a_Vector,Tri_1_b_Vector,Tri_1_n_Vector)   

Tri_2_a_Vector(1:3) = Tri_2(2,1:3) - Tri_2(1,1:3)
Tri_2_b_Vector(1:3) = Tri_2(3,1:3) - Tri_2(1,1:3)
call Vector_Cross_Product_3(Tri_2_a_Vector,Tri_2_b_Vector,Tri_2_n_Vector)   

call Tool_Angle_of_Vectors_a_and_b_3D(Tri_1_n_Vector,Tri_2_n_Vector,c_angle,2)

if(abs(c_angle)<=Tol_8) then
    Yes_Parallel  = .True.
    return
endif

if(abs(c_angle-pi)<=Tol_8) then
    Yes_Parallel  = .True.
    return
endif

return 
end SUBROUTINE Tool_Yes_Two_Triangles_Parallel       
