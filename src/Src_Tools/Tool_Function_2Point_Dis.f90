!-----------------------------------------------------------
! Brief: Euclidean distance between two 2D points
!
! Parameters:
!   Input:  A(2) - first point coordinates
!   Input:  B(2) - second point coordinates
!   Returns: real - Euclidean distance
!-----------------------------------------------------------

function Tool_Function_2Point_Dis(A,B)
!This function cals distance between point A and B.

use Global_Float_Type
implicit none
real(kind=FT),intent(in)::A(2),B(2)
real(kind=FT) :: Tool_Function_2Point_Dis

Tool_Function_2Point_Dis = sqrt((A(1)-B(1))**2+(A(2)-B(2))**2)

end function Tool_Function_2Point_Dis
