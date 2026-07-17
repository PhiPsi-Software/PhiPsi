!-----------------------------------------------------------
! Brief: Allocate ragged arrays for a given 3D solid element.
!
! Parameters:
!   Input:  i_E        - element index (1-based)
!           Value_Type - storage class (1 crack-tip enriched data,
!                        2 fluid data, 3 junction data, 4 other)
!
! Notes:   Allocates Solid_El_* ragged arrays (vertex coords, normal,
!          tangent vectors, base/tip baseline) sized by Solid_El_Max_num_Crs.
!-----------------------------------------------------------

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
