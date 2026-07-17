!-----------------------------------------------------------
! Brief: Estimate the element average length at a given point.
!
! Parameters:
!   Input:  x       - x-coordinate of the point
!           y       - y-coordinate of the point
!   Output: Ave_Length - characteristic length (sqrt of area)
!
! Notes:   Locates the element containing (x,y) and returns the
!          square root of its polygon area as a length measure.
!-----------------------------------------------------------

subroutine Cal_Length_of_ELe_by_Coors_Ave(x,y,Ave_Length)
!     Calculate the average length of the element which the 
!     point(x,y) is located in. 
!     Written by Fang Shi on May 14, 2017

use Global_Float_Type
use Global_Model

implicit none

real(kind=FT) x,y
real(kind=FT),intent(out)::Ave_Length
integer :: OUT_Elem
integer :: NN(4)
real(kind=FT) :: area

call Cal_Ele_Num_by_Coors(x,y,OUT_Elem) 

NN  = Elem_Node(OUT_Elem,1:4)

call Tool_Area_Polygon(Coor(NN,1),Coor(NN,2),4,area)

Ave_Length = sqrt(area)

return 
end subroutine Cal_Length_of_ELe_by_Coors_Ave                       
