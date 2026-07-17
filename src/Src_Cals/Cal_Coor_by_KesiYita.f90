!-----------------------------------------------------------
! Brief: Map a 2D point from local to global coordinates for a quad element.
!
! Parameters:
!   Input:  kesi,yita   - local (parametric) coordinates
!           X_NODES     - x-coordinates of the 4 element nodes
!           Y_NODES     - y-coordinates of the 4 element nodes
!   Output: Out_x,Out_y - global Cartesian coordinates of the point
!
! Notes:   Standard Q4 isoparametric mapping using bilinear shape functions.
!-----------------------------------------------------------

subroutine Cal_Coor_by_KesiYita(kesi,yita,X_NODES,Y_NODES, Out_x,Out_y)
!     Calculate the coordinates in the global coordinate system based on the coordinates in the local coordinate system
use Global_Float_Type
implicit none
real(kind=FT),intent(in)::kesi,yita,X_NODES(4),Y_NODES(4)
real(kind=FT),intent(out)::Out_x,Out_y
real(kind=FT) :: N(4)

N  = 0.25D0 * [(1-kesi)*(1-yita),(1+kesi)*(1-yita), (1+kesi)*(1+yita),(1-kesi)*(1+yita)]

Out_x  = N(1)*X_NODES(1)+N(2)*X_NODES(2)+ N(3)*X_NODES(3)+N(4)*X_NODES(4)
Out_y  = N(1)*Y_NODES(1)+N(2)*Y_NODES(2)+ N(3)*Y_NODES(3)+N(4)*Y_NODES(4)

return 
end subroutine Cal_Coor_by_KesiYita                          
