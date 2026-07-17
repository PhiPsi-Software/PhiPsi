!-----------------------------------------------------------
! Brief: Returns the centroid of a quadrilateral with vertices ordered counter-clockwise.
!
! Parameters:
!   Input:  x, y     - arrays of 4 vertex coordinates
!   Output: Centroid - 2D centroid coordinates
!
! Notes:   Uses the simple mean of the four vertex coordinates; a precise
!          analytic formula is present in the source but commented out.
!-----------------------------------------------------------

subroutine Tool_Centroid_Quad(x,y,Centroid)

! Calculate the centroid of a quadrilateral (with the four points arranged counterclockwise)
! Reference: http://www.doc88.com/p-3761246197349.html
! The formula is only applicable to four points arranged counterclockwise.

use Global_Float_Type
IMPLICIT NONE
real(kind=FT),intent(in):: x(4)
real(kind=FT),intent(in):: y(4)
real(kind=FT),intent(out)::Centroid(2)
Centroid(1) = sum(x(1:4))/FOU
Centroid(2) = sum(y(1:4))/FOU


RETURN
END subroutine Tool_Centroid_Quad
