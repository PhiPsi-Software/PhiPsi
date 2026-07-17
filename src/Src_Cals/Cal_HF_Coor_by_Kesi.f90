!-----------------------------------------------------------
! Brief: Map a 1D local coordinate (kesi) on a hydraulic
!        fracture element to its 2D global coordinates.
!
! Parameters:
!   Input:  kesi  - local coordinate in [-1, 1]
!           x1,y1 - coordinates of the first element node
!           x2,y2 - coordinates of the second element node
!   Output: Out_x - global x of the mapped point
!           Out_y - global y of the mapped point
!-----------------------------------------------------------

subroutine Cal_HF_Coor_by_Kesi(kesi,x1,y1,x2,y2, Out_x,Out_y)
!     Calculate the global coordinates of the hydraulic fracture Gauss point based on the local coordinates kesi
!     x1 and y1 are the coordinates of the first node of the hydraulic fracture element. The node is also called the calculation point.
!     x2 and y2 are the coordinates of the first node of the hydraulic fracture element.
use Global_Float_Type
implicit none
real(kind=FT),intent(in)::kesi,x1,y1,x2,y2
real(kind=FT),intent(out)::Out_x,Out_y
real(kind=FT) x0,y0
real(kind=FT) angle,length

x0     = HLF*(x1+x2)
y0     = HLF*(y1+y2)
angle  = atan2(y2-y1,x2-x1)
length = sqrt((x2-x1)**2 + (y2-y1)**2)

Out_x  = x0  + HLF*length*kesi*cos(angle)
Out_y  = y0  + HLF*length*kesi*sin(angle)

return 
end subroutine Cal_HF_Coor_by_Kesi                         
