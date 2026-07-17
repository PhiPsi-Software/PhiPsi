!-----------------------------------------------------------
! Brief: Euclidean distance between two 3D points
!
! Parameters:
!   Input:  A(3) - first point coordinates
!   Input:  B(3) - second point coordinates
!   Returns: real - Euclidean distance
!-----------------------------------------------------------

function Tool_Function_2Point_Dis_3D(A,B)
!This function cals distance between point A and B.

use Global_Float_Type
implicit none
real(kind=FT),intent(in)::A(3),B(3)
real(kind=FT) :: Tool_Function_2Point_Dis_3D

Tool_Function_2Point_Dis_3D = sqrt((A(1)-B(1))**2 +(A(2)-B(2))**2 +(A(3)-B(3))**2)

end function Tool_Function_2Point_Dis_3D
