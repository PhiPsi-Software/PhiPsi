!-----------------------------------------------------------
! Brief: Test whether a 2D point lies on the boundary of a polygon
!
! Parameters:
!   Input:  x,y       - point coordinates
!   Input:  xpol,ypol - polygon vertex coordinates
!   Input:  npol      - number of vertices
!   Output: Yes_ONOUT - true if the point lies on any polygon edge
!
! Notes:   Iterates edges and delegates to Tool_Yes_On_Line.
!-----------------------------------------------------------

subroutine Tool_Yes_On_Poly(x,y,xpol,ypol,npol,Yes_ONOUT)
!     Determine whether a point is on the boundary of a polygon
use Global_Float_Type
integer i,npol
real(kind=FT) x,y,xpol(npol),ypol(npol)
logical,intent(out)::Yes_ONOUT

logical Yes_ON


Yes_ONOUT = .False.

do i=1,npol-1 

    call Tool_Yes_On_Line(x,y, [xpol(i),ypol(i)], [xpol(i+1),ypol(i+1)],Yes_ON)
    if (Yes_ON.eqv..True.) then
        Yes_ONOUT = .True.
        exit
    endif

end do 

return 
end SUBROUTINE Tool_Yes_On_Poly                          
