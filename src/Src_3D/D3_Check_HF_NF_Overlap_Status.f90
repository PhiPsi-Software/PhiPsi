!-----------------------------------------------------------
! Brief: Test whether HF and NF bounding boxes overlap in 3D.
!
! Parameters:
!   Input:  HF_Coor_Ranges - HF bounding range (3 x [min,max])
!           NF_Coor_Ranges - NF bounding range (3 x [min,max])
!   Output: Yes_Overlap    - true if HF and NF overlap on all three axes
!
! Notes:   Uses Tool_Yes_Two_Ranges_Overlapped_Double on x, y, and z.
!          Supports rectangular/polygonal natural fractures only.
!-----------------------------------------------------------

SUBROUTINE D3_Check_HF_NF_Overlap_Status(HF_Coor_Ranges,NF_Coor_Ranges,Yes_Overlap)
! Used to detect the coordinate coincidence of hydraulic fractures and natural fractures.
! Only supports rectangular or polygonal natural fractures.
! 2023-01-09.

!*****************************
! Read public variable module
!*****************************
use Global_Float_Type

!**********************
! Variable Declaration
!**********************
implicit none
logical,intent(out)::Yes_Overlap
real(kind=FT),intent(in)::HF_Coor_Ranges(3,2),NF_Coor_Ranges(3,2)
logical x_Logical_Yes,y_Logical_Yes,z_Logical_Yes

!-----------
! Detection
!-----------
Yes_Overlap = .False.
call Tool_Yes_Two_Ranges_Overlapped_Double(HF_Coor_Ranges(1,1:2),NF_Coor_Ranges(1,1:2),x_Logical_Yes)  
if (x_Logical_Yes .eqv. .False.) then
    return
endif
call Tool_Yes_Two_Ranges_Overlapped_Double(HF_Coor_Ranges(2,1:2),NF_Coor_Ranges(2,1:2),y_Logical_Yes)  
if (y_Logical_Yes .eqv. .False.) then
    return
endif
call Tool_Yes_Two_Ranges_Overlapped_Double(HF_Coor_Ranges(3,1:2),NF_Coor_Ranges(3,1:2),z_Logical_Yes)  
if (z_Logical_Yes .eqv. .False.) then
    return
endif  
Yes_Overlap = .True.


RETURN
END SUBROUTINE D3_Check_HF_NF_Overlap_Status
