!-----------------------------------------------------------
! Brief: Compute the x/y/z coordinate range of a natural fracture.
!
! Parameters:
!   Input:  i_C        - natural-fracture index (1-based)
!   Output: Coor_Ranges- min/max for each axis: (axis, [min,max])
!
! Notes:   Only supports rectangular or polygonal natural fractures.
!          Aborts with an error if vertex count is non-positive.
!-----------------------------------------------------------

SUBROUTINE D3_Get_NF_Coor_Ranges(i_C,Coor_Ranges)
! Obtain the coordinate range of natural fractures. Only rectangular or polygonal natural fractures
! are supported.
!2023-01-09.

!-----------------------------
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Crack_Common
use Global_Crack_3D   

!---------------------------
! Variable Type Declaration
!---------------------------
implicit none
integer i_C
real(kind=FT),intent(out):: Coor_Ranges(3,2)
integer num_Vertex

!-----------------------------------------------------------
! Only supports rectangular or polygonal natural fractures.
!-----------------------------------------------------------
num_Vertex = Each_NaCr3D_Poi_Num(i_C)  
if(num_Vertex <=0)then
    print *,'    ERROR-2023010901 :: Not rectangle or polygon NF!'
    print *,'                        In D3_Get_NF_Coor_Ranges.f90!'
    call Warning_Message('S',Keywords_Blank)      
endif

!--------------------------
! Extract coordinate range
!--------------------------
! x range
Coor_Ranges(1,1) =  minval(Na_Crack3D_Coor(i_C,1:num_Vertex,1))
Coor_Ranges(1,2) =  maxval(Na_Crack3D_Coor(i_C,1:num_Vertex,1))     
! y-range
Coor_Ranges(2,1) =  minval(Na_Crack3D_Coor(i_C,1:num_Vertex,2))
Coor_Ranges(2,2) =  maxval(Na_Crack3D_Coor(i_C,1:num_Vertex,2))   
! z range
Coor_Ranges(3,1) =  minval(Na_Crack3D_Coor(i_C,1:num_Vertex,3))
Coor_Ranges(3,2) =  maxval(Na_Crack3D_Coor(i_C,1:num_Vertex,3))           

RETURN
END SUBROUTINE D3_Get_NF_Coor_Ranges
