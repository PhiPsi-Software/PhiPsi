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
 
SUBROUTINE D3_Offset_Initial_Crack_to_Ele_Center(i_C,c_Crack_Type)
! Adjust the centroid of the crack surface to the centroid of the corresponding element.
! 2022-08-02.
!
!
! c_Crack_Type = 1  ! Rectangular crack surface or polygonal crack surface
! c_Crack_Type = 2  !Circular crack surface
! c_Crack_Type = 3  ! Elliptical crack surface (To be done)

!#############################
! Read public variable module
!#############################
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Crack_Common
use Global_Crack_3D
use Global_Filename
use Global_Cal_Ele_Num_by_Coors_3D  

!###########################
! Variable Type Declaration
!###########################
implicit none
integer,intent(in)::i_C,c_Crack_Type
integer num_Cr_Poi
real(kind=FT) CR_Center_Old(3)
real(kind=FT) c_Elem_Centroid(3)
real(kind=FT) Center_Offset_Vector(3)
integer OUT_Elem,i_Cr_Poi
integer Ele_Num_Cache

!############################################
! If it is a rectangular or polygonal crack.
!############################################
Ele_Num_Cache = 1
if (c_Crack_Type == 1) then
    num_Cr_Poi = Each_Cr_Poi_Num(i_C)
    CR_Center_Old(1) =sum(Crack3D_Coor(i_C,1:num_Cr_Poi,1))/num_Cr_Poi
    CR_Center_Old(2) =sum(Crack3D_Coor(i_C,1:num_Cr_Poi,2))/num_Cr_Poi
    CR_Center_Old(3) =sum(Crack3D_Coor(i_C,1:num_Cr_Poi,3))/num_Cr_Poi
    ! Calculate the element number based on the coordinates
    call Cal_Ele_Num_by_Coors_3D(CR_Center_Old(1),CR_Center_Old(2),CR_Center_Old(3),Ele_Num_Cache,OUT_Elem)
    ! element centroid
    c_Elem_Centroid(1:3) = Elem_Centroid(OUT_Elem,1:3)
    ! Obtain offset vector
    Center_Offset_Vector(1:3) = c_Elem_Centroid(1:3) - CR_Center_Old(1:3)
    ! Update rectangle or polygon coordinates
    do i_Cr_Poi=1,num_Cr_Poi
        Crack3D_Coor(i_C,i_Cr_Poi,1:3) = Crack3D_Coor(i_C,i_Cr_Poi,1:3) + Center_Offset_Vector(1:3)
    enddo
    !To be done.
endif
      
!###########################################################
! If it is a circular crack: just update the coordinates of 
! the center, and leave the other variables unchanged.
!###########################################################
if (c_Crack_Type == 2) then
    ! Original crack center.
    CR_Center_Old(1:3) = Crack3D_Cir_Coor(i_C,1:3)
    
    ! Calculate the element number based on the coordinates
    call Cal_Ele_Num_by_Coors_3D(CR_Center_Old(1),CR_Center_Old(2),CR_Center_Old(3),Ele_Num_Cache,OUT_Elem)
    
    ! element centroid
    c_Elem_Centroid(1:3) = Elem_Centroid(OUT_Elem,1:3)
    
    ! Update center coordinates
    Crack3D_Cir_Coor(i_C,1:3) = c_Elem_Centroid(1:3)
    
    ! Normal direction.
    !n_Vector(1:3) = Crack3D_Cir_Coor(i_C,4:6)

    ! Update center coordinates
    !Crack3D_Cir_Coor(i_C,1:3) = CR_Center_New(1:3) 
endif


!############################################################
! If it is an elliptical crack, it is not supported for now.
!############################################################
if (c_Crack_Type == 3) then
    print *,'     ERROR(2022080201) :: Elliptical cracks not supported yet!'
    call Warning_Message('S',Keywords_Blank) 
endif
        

      
RETURN
END SUBROUTINE D3_Offset_Initial_Crack_to_Ele_Center
