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
 
subroutine D3_Get_Point_Domain_Number(point,Domain_Number)
! Get the partition number of the point.
!2022-11-24. NEWFTU2022112401.


!***************
! Public Module
!***************
use Global_Float_Type
use Global_Model
use Global_Common

!**********************
! Variable Declaration
!**********************
implicit none
real(kind=FT),intent(in)::point(3)
integer,intent(out):: Domain_Number
real(kind=FT) c_Tol
real(kind=FT) x,y,z


!*****************************************
! Get the partition number of the element
!*****************************************
x = point(1)
y = point(2)
z = point(3)

c_Tol    = Tol_7
if( ((x - Model_Center_x) > -c_Tol) .and. ((y - Model_Center_y) > -c_Tol) .and. ((z - Model_Center_z) > -c_Tol)) then
    Domain_Number  = 1
    return
elseif( ((x - Model_Center_x) <  c_Tol) .and. ((y - Model_Center_y) >-c_Tol) .and. ((z - Model_Center_z) > -c_Tol)) then
    Domain_Number  = 2
    return
elseif( ((x - Model_Center_x) <  c_Tol) .and. ((y - Model_Center_y) < c_Tol) .and. ((z - Model_Center_z) > -c_Tol)) then
    Domain_Number  = 3
    return
elseif( ((x - Model_Center_x) > -c_Tol) .and. ((y - Model_Center_y) < c_Tol) .and. ((z - Model_Center_z) > -c_Tol)) then
    Domain_Number  = 4
    return
elseif( ((x - Model_Center_x) > -c_Tol) .and. ((y - Model_Center_y) >-c_Tol) .and. ((z - Model_Center_z) <  c_Tol)) then
    Domain_Number  = 5
    return
elseif( ((x - Model_Center_x) <  c_Tol) .and. ((y - Model_Center_y) >-c_Tol) .and. ((z - Model_Center_z) <  c_Tol)) then
    Domain_Number  = 6
    return
elseif( ((x - Model_Center_x) <  c_Tol) .and. ((y - Model_Center_y) < c_Tol) .and. ((z - Model_Center_z) <  c_Tol)) then
    Domain_Number  = 7
    return
elseif( ((x - Model_Center_x) > -c_Tol) .and. ((y - Model_Center_y) < c_Tol) .and. ((z - Model_Center_z) <  c_Tol)) then
    Domain_Number  = 8
    return
else
    print *, '    Error :: cannot get domain number in D3_Get_Point_Domain_Number.f90!'
    call Warning_Message('S',Keywords_Blank)
endif

    
RETURN
END subroutine D3_Get_Point_Domain_Number