!-----------------------------------------------------------
! Brief: Compute the volume of a quadrilateral pyramid with apex A and base BCDE
!
! Parameters:
!   Input:  A(3)     - pyramid apex coordinates
!   Input:  B(3),C(3),D(3),E(3) - four base vertex coordinates
!   Output: Volume   - resulting pyramid volume
!
! Notes:   Splits base into two triangles and uses perpendicular distance from apex.
!-----------------------------------------------------------

subroutine Tool_Volume_Pyramid(A,B,C,D,E,Volume)

! Calculate the volume of a quadrilateral pyramid
! Theory, see my notes V3-P138
use Global_Float_Type
IMPLICIT NONE

real(kind=FT),intent(in)::A(3),B(3),C(3),D(3),E(3)
real(kind=FT),intent(out)::Volume

real(kind=FT) Area_BCDE,Area_BCD,Area_CDE
real(kind=FT) Distance

call Tool_Area_Tri_3D(B,C,D,Area_BCD)
call Tool_Area_Tri_3D(C,D,E,Area_CDE)

Area_BCDE = Area_BCD + Area_CDE



if(Area_BCD>=Tol_20) then
    call Tool_Dis_Point_to_3D_Tri_only_Dis(A,B,C,D,Distance)
else
    call Tool_Dis_Point_to_3D_Tri_only_Dis(A,C,D,E,Distance)
endif


Volume = ONE/THR*Area_BCDE*abs(Distance)

RETURN
END subroutine Tool_Volume_Pyramid
