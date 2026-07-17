!-----------------------------------------------------------
! Brief: Compute 2D distance from point c to segment ab
!
! Parameters:
!   Input:  a(2),b(2)  - segment endpoints
!   Input:  c(2)       - query point coordinates
!   Output: Dis        - minimum normal distance to the segment
!
! Notes:   Uses normal projection onto the unit edge vector.
!-----------------------------------------------------------

subroutine Tool_Dis_Point_to_Seg(a,b,c,Dis)
    ! Calculate the distance from point c to the line segment ab (normal distance, the minimum of the vertex distances)
    ! http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment

    use Global_Float_Type
    implicit none
    real(kind=FT), intent(in) :: a(2), b(2), c(2)
    real(kind=FT), intent(out) :: Dis
    real(kind=FT):: t(2), n(2), ac(2)
    real(kind=FT):: dd
    t   = b - a
    dd  = sqrt(t(1)**2+t(2)**2)
    t   = t/dd
    n   = (/-t(2), t(1)/)
    ac  = c - a
    Dis = abs(ac(1)*n(1)+ac(2)*n(2))

    return
end SUBROUTINE Tool_Dis_Point_to_Seg
