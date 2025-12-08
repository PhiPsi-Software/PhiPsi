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
 
SUBROUTINE D3_Allocate_Ele_Memory(i_E,Value_Type)
! Allocate memory for element-related variables. 2022-09-04.



!-----------------------------
! Read Public Variable Module
!-----------------------------
use Global_Float_Type
use Global_Crack_3D

implicit none
integer,intent(in):: i_E,Value_Type

!--------------------------------------------    
! Allocate memory for different jagged array 
! variables according to Value_Type.
!--------------------------------------------    
select case(Value_Type)      

!////////////////////////////////////////////////////////////////////
! Value_Type=1, related to crack-tip enhancement elements and nodes.
!////////////////////////////////////////////////////////////////////
case(1)
    !IMPROV2022120301.
    if(.not. allocated(Solid_El_Vertex_Num(i_E)%row))then
        allocate(Solid_El_Vertex_Num(i_E)%row(Solid_El_Max_num_Crs))
        Solid_El_Vertex_Num(i_E)%row(1:Solid_El_Max_num_Crs) = 0
    endif
    
    if(.not. allocated(Solid_El_Vertex_Coor(i_E)%row))then
        allocate(Solid_El_Vertex_Coor(i_E)%row(Solid_El_Max_num_Crs,5,3))
        Solid_El_Vertex_Coor(i_E)%row(1:Solid_El_Max_num_Crs,1:5,1:3) = ZR
    endif
    if(.not. allocated(Solid_El_Vertex_Nor_Vec(i_E)%row))then
        allocate(Solid_El_Vertex_Nor_Vec(i_E)%row(Solid_El_Max_num_Crs,5,3))
        Solid_El_Vertex_Nor_Vec(i_E)%row(1:Solid_El_Max_num_Crs,1:5,1:3) = ZR
    endif
    if(.not. allocated(Solid_El_Vertex_x_Vec(i_E)%row))then
        allocate(Solid_El_Vertex_x_Vec(i_E)%row(Solid_El_Max_num_Crs,5,3))
        Solid_El_Vertex_x_Vec(i_E)%row(1:Solid_El_Max_num_Crs,1:5,1:3) = ZR
    endif
    if(.not. allocated(Solid_El_Vertex_y_Vec(i_E)%row))then
        allocate(Solid_El_Vertex_y_Vec(i_E)%row(Solid_El_Max_num_Crs,5,3))
        Solid_El_Vertex_y_Vec(i_E)%row(1:Solid_El_Max_num_Crs,1:5,1:3) = ZR
    endif
    if(.not. allocated(Solid_El_Vertex_z_Vec(i_E)%row))then
        allocate(Solid_El_Vertex_z_Vec(i_E)%row(Solid_El_Max_num_Crs,5,3))
        Solid_El_Vertex_z_Vec(i_E)%row(1:Solid_El_Max_num_Crs,1:5,1:3) = ZR
    endif
    if(.not. allocated(Solid_El_Pre_Vertex_Coor(i_E)%row))then
        allocate(Solid_El_Pre_Vertex_Coor(i_E)%row(Solid_El_Max_num_Crs,5,3))
        Solid_El_Pre_Vertex_Coor(i_E)%row(1:Solid_El_Max_num_Crs,1:5,1:3) = ZR
    endif

    if(.not. allocated(Solid_El_Pre_Vertex_Nor_Vec(i_E)%row))then
        allocate(Solid_El_Pre_Vertex_Nor_Vec(i_E)%row(Solid_El_Max_num_Crs,5,3))
        Solid_El_Pre_Vertex_Nor_Vec(i_E)%row(1:Solid_El_Max_num_Crs,1:5,1:3) = ZR
    endif
    if(.not. allocated(Solid_El_Pre_Vertex_x_Vec(i_E)%row))then
        allocate(Solid_El_Pre_Vertex_x_Vec(i_E)%row(Solid_El_Max_num_Crs,5,3))
        Solid_El_Pre_Vertex_x_Vec(i_E)%row(1:Solid_El_Max_num_Crs,1:5,1:3) = ZR
    endif
    if(.not. allocated(Solid_El_Pre_Vertex_y_Vec(i_E)%row))then
        allocate(Solid_El_Pre_Vertex_y_Vec(i_E)%row(Solid_El_Max_num_Crs,5,3))
        Solid_El_Pre_Vertex_y_Vec(i_E)%row(1:Solid_El_Max_num_Crs,1:5,1:3) = ZR
    endif
    if(.not. allocated(Solid_El_Pre_Vertex_z_Vec(i_E)%row))then
        allocate(Solid_El_Pre_Vertex_z_Vec(i_E)%row(Solid_El_Max_num_Crs,5,3))
        Solid_El_Pre_Vertex_z_Vec(i_E)%row(1:Solid_El_Max_num_Crs,1:5,1:3) = ZR
    endif
    if(.not. allocated(Solid_El_Nex_Vertex_Coor(i_E)%row))then
        allocate(Solid_El_Nex_Vertex_Coor(i_E)%row(Solid_El_Max_num_Crs,5,3))
        Solid_El_Nex_Vertex_Coor(i_E)%row(1:Solid_El_Max_num_Crs,1:5,1:3) = ZR
    endif
    
    if(.not. allocated(Solid_El_Nex_Vertex_Nor_Vec(i_E)%row))then
        allocate(Solid_El_Nex_Vertex_Nor_Vec(i_E)%row(Solid_El_Max_num_Crs,5,3))
        Solid_El_Nex_Vertex_Nor_Vec(i_E)%row(1:Solid_El_Max_num_Crs,1:5,1:3) = ZR
    endif
    if(.not. allocated(Solid_El_Nex_Vertex_x_Vec(i_E)%row))then
        allocate(Solid_El_Nex_Vertex_x_Vec(i_E)%row(Solid_El_Max_num_Crs,5,3))
        Solid_El_Nex_Vertex_x_Vec(i_E)%row(1:Solid_El_Max_num_Crs,1:5,1:3) = ZR
    endif
    if(.not. allocated(Solid_El_Nex_Vertex_y_Vec(i_E)%row))then
        allocate(Solid_El_Nex_Vertex_y_Vec(i_E)%row(Solid_El_Max_num_Crs,5,3))
        Solid_El_Nex_Vertex_y_Vec(i_E)%row(1:Solid_El_Max_num_Crs,1:5,1:3) = ZR
    endif
    if(.not. allocated(Solid_El_Nex_Vertex_z_Vec(i_E)%row))then
        allocate(Solid_El_Nex_Vertex_z_Vec(i_E)%row(Solid_El_Max_num_Crs,5,3))
        Solid_El_Nex_Vertex_z_Vec(i_E)%row(1:Solid_El_Max_num_Crs,1:5,1:3) = ZR
    endif
    
    if(.not. allocated(Solid_El_Tip_BaseLine(i_E)%row))then
        allocate(Solid_El_Tip_BaseLine(i_E)%row(Solid_El_Max_num_Crs,2,3))
        Solid_El_Tip_BaseLine(i_E)%row(1:Solid_El_Max_num_Crs,1:2,1:3) = ZR
    endif

    if(.not. allocated(Solid_El_Tip_BaseLine_Nor_Vec(i_E)%row))then
        allocate(Solid_El_Tip_BaseLine_Nor_Vec(i_E)%row(Solid_El_Max_num_Crs,3))
        Solid_El_Tip_BaseLine_Nor_Vec(i_E)%row(1:Solid_El_Max_num_Crs,1:3) = ZR
    endif
    if(.not. allocated(Solid_El_Tip_BaseLine_x_Vec(i_E)%row))then
        allocate(Solid_El_Tip_BaseLine_x_Vec(i_E)%row(Solid_El_Max_num_Crs,3))
        Solid_El_Tip_BaseLine_x_Vec(i_E)%row(1:Solid_El_Max_num_Crs,1:3) = ZR
    endif
    if(.not. allocated(Solid_El_Tip_BaseLine_y_Vec(i_E)%row))then
        allocate(Solid_El_Tip_BaseLine_y_Vec(i_E)%row(Solid_El_Max_num_Crs,3))
        Solid_El_Tip_BaseLine_y_Vec(i_E)%row(1:Solid_El_Max_num_Crs,1:3) =ZR
    endif
    if(.not. allocated(Solid_El_Tip_BaseLine_z_Vec(i_E)%row))then
        allocate(Solid_El_Tip_BaseLine_z_Vec(i_E)%row(Solid_El_Max_num_Crs,3))
        Solid_El_Tip_BaseLine_z_Vec(i_E)%row(1:Solid_El_Max_num_Crs,1:3) = ZR
    endif
    
    if(.not. allocated(Solid_El_Tip_BaseLine_T_theta(i_E)%row))then
        allocate(Solid_El_Tip_BaseLine_T_theta(i_E)%row(Solid_El_Max_num_Crs,3))
        Solid_El_Tip_BaseLine_T_theta(i_E)%row(1:Solid_El_Max_num_Crs,1:3) = ZR
    endif
 
    if(.not. allocated(Solid_El_Tip_BaseLine_T_Matrix(i_E)%row))then
        allocate(Solid_El_Tip_BaseLine_T_Matrix(i_E)%row(Solid_El_Max_num_Crs,3,3))
        Solid_El_Tip_BaseLine_T_Matrix(i_E)%row(1:Solid_El_Max_num_Crs,1:3,1:3) = ZR
    endif
!////////////////////////////////////////////////////
! Value_Type=2, related to fluid elements and nodes.
!////////////////////////////////////////////////////
case(2)

!//////////////
!Value_Type=3, 
!//////////////
case(3)


!//////////////
!Value_Type=4, 
!//////////////
case(4)

end select
      


RETURN
END SUBROUTINE D3_Allocate_Ele_Memory
