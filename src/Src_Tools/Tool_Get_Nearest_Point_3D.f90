!-----------------------------------------------------------
! Brief: Find the index of the 3D point nearest to a given query point.
!
! Parameters:
!   Input:  num_Point - number of candidate points
!   Input:  Points - array of 3D point coordinates
!   Input:  A(3) - query point
!   Output: Point_Num - index of the nearest point in Points
!
! Notes:   Computes all distances explicitly and selects the
!   minimum using MINLOC; inlined Euclidean distance.
!-----------------------------------------------------------

subroutine Tool_Get_Nearest_Point_3D(num_Point,Points,A, Point_Num)
!     Find the point closest to A from a set of 3D points.
!     2022-04-30.
!     NEWFTU2022043007.

use Global_Float_Type
use Global_Common

implicit none
integer,intent(in)::num_Point
real(kind=FT),intent(in)::Points(num_Point,3),A(3)       
integer,intent(out)::Point_Num
real(kind=FT) :: Dis_Vector(num_Point)
real(kind=FT) :: Tool_Function_2Point_Dis_3D
integer :: i_P

if(num_Point<=1)then
    print *, '    Error :: num_Point<=1!'
    print *, '             in Tool_Get_Nearest_Point_3D.f!'
    call Warning_Message('S',Keywords_Blank)      
endif


do i_P=1,num_Point
    Dis_Vector(i_P)=sqrt((A(1)-Points(i_P,1))*(A(1)-Points(i_P,1)) + (A(2)-Points(i_P,2))*(A(2)-Points(i_P,2)) &
    + (A(3)-Points(i_P,3))*(A(3)-Points(i_P,3)) )
enddo



Point_Num = MINLOC(Dis_Vector(1:num_Point),1)

return 
end subroutine Tool_Get_Nearest_Point_3D                     
