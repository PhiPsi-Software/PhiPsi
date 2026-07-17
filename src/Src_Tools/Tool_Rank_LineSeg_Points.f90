!-----------------------------------------------------------
! Brief: Reorder unordered segment points end-to-end along AB.
!
! Parameters:
!   In/Out: Line_AB   - Point set of the segment; reordered in place
!   Input:  num_Points - Number of points defining the segment
!
! Notes:   First and last rows are kept as endpoints A and B;
!   intermediate points are sorted by ascending distance to A so
!   the resulting sequence is geometrically continuous.
!-----------------------------------------------------------

subroutine Tool_Rank_LineSeg_Points(Line_AB,num_Points)
!     Reorder the unordered points that make up line segment AB from start to finish so that they connect end to end.
!     Sort according to the distance from each point to endpoint 1
use Global_Float_Type      
implicit none
real(kind=FT),intent(inout):: Line_AB(num_Points,2)
integer,intent(in):: num_Points

real(kind=FT) new_Line_AB(num_Points,2),L(num_Points), A(2),B(2)
integer i_P,Location_Min

A(1) = Line_AB(1,1)
A(2) = Line_AB(1,2)
B(1) = Line_AB(num_Points,1)
B(2) = Line_AB(num_Points,2) 

new_Line_AB(1,1:2) = A
new_Line_AB(num_Points,1:2) = B


L(1:num_Points) =  1.0D8
do i_P=2,num_Points-1
    L(i_P) = sqrt((Line_AB(i_P,2)-A(2))**2+ (Line_AB(i_P,1)-A(1))**2)
end do

do i_P=2,num_Points-1
    Location_Min = minLoc(L(1:num_Points),1) 
    new_Line_AB(i_P,1:2) = Line_AB(Location_Min,1:2)
    L(Location_Min) = 1.0D8 
end do


Line_AB = new_Line_AB

return 
end subroutine Tool_Rank_LineSeg_Points              
