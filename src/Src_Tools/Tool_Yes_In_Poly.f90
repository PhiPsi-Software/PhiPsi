!-----------------------------------------------------------
! Brief: Test whether a 2D point lies inside a closed polygon (ray-casting test)
!
! Parameters:
!   Input:  x,y       - point coordinates
!   Input:  xpol,ypol - polygon vertex coordinates
!   Input:  npol      - number of vertices
!   Output: Yes_INOUT - true if the point is inside
!
! Notes:   Performs a quick bounding-box reject before the crossing test.
!-----------------------------------------------------------

subroutine Tool_Yes_In_Poly(x,y,xpol,ypol,npol,Yes_INOUT)
!     Determine whether a point is inside a polygon. Note that the polygon must be closed.
!     npol equals the number of edges plus one or the number of vertices plus one
use Global_Float_Type
implicit none
integer i,j,npol
real(kind=FT) x,y,xpol(npol),ypol(npol)
logical Yes_INOUT

Yes_INOUT = .False.
if((x.gt.maxval(xpol)) .or. (x.lt.minval(xpol)) .or. (y.gt.maxval(ypol)) .or. (y.lt.minval(ypol))) then
    Yes_INOUT = .False.
    goto 100

else
    j = npol-1 
    do i=1,npol-1
        if ( ((ypol(i).gt.y).neqv. (ypol(j).gt.y)) .and. (x .lt. (xpol(j)-xpol(i)) * (y-ypol(i)) &
        / (ypol(j)-ypol(i)) + xpol(i)) ) then
        Yes_INOUT = .NOT. Yes_INOUT
    end if
    j = i
end do
end if

100 continue

return 
end SUBROUTINE Tool_Yes_In_Poly                          
