!-----------------------------------------------------------
! Brief: Locates the two points on a circle at a given arc length from a reference point on the circle.
!
! Parameters:
!   Input:  x, y, r       - circle center coordinates and radius
!   Input:  Point_I       - reference point on the circle
!   Input:  Length_of_Arc - signed arc length from Point_I
!   Output: Point_Status  - 1 if Point_I lies on the circle, -1 otherwise
!   Output: P1_P2         - 2x2 array of the two resulting points
!
! Notes:   Uses circle-circle intersection of the original circle with a
!          concentric circle of chord-radius 2*r*sin(theta/2) to find candidates.
!-----------------------------------------------------------

subroutine Tool_2_Points_of_Cir_give_Point_L( x,y,r, Point_I,Length_of_Arc,Point_Status, P1_P2)
    ! Given a circle and a point I on the circle, calculate the coordinates of the two points on the circle where the arc length from I is Length_of_Arc.
    ! Usage: See the diagram at the lower left corner of note V5-58_P117.
    ! Algorithm principle: The essence of this problem is to calculate the intersection points of two circles.

    use Global_Float_Type
    use Global_Elem_Area_Vol

    implicit none
    real(kind=FT),intent(in)  :: x,y,r,Point_I(2),Length_of_Arc
    real(kind=FT),intent(out) :: P1_P2(2,2)
    integer,intent(out)       :: Point_Status
    real(kind=FT) :: Tool_Function_2Point_Dis
    real(kind=FT) :: c_r
    real(kind=FT) :: x1,y1,r1
    real(kind=FT) :: Radian_Angle
    real(kind=FT) :: c_Inters(2,2)
    integer :: c_Circles_Status,c_num_Inters


    P1_P2(1:2,1:2) = ZR
    c_r = Tool_Function_2Point_Dis(Point_I,[x,y])
    if(abs(c_r - r)<=Tol_11)then
        Point_Status = 1
        Radian_Angle = Length_of_Arc/(TWO*pi*r)*(TWO*pi)
        x1 = Point_I(1)
        y1 = Point_I(2)
        r1 = TWO*(r*sin(Radian_Angle/TWO))
call Tool_Intersection_Circle_and_Circle( x,y,r, x1,y1,r1,c_Circles_Status, c_num_Inters,c_Inters)
        P1_P2 = c_Inters
    else
        Point_Status =-1
        P1_P2(1:2,1:2) = ZR
    endif


    return
    end SUBROUTINE Tool_2_Points_of_Cir_give_Point_L
