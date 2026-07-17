!-----------------------------------------------------------
! Brief: Finds all intersection points of a 3D line segment AB with a brick (8-node hex) element.
!
! Parameters:
!   Input:  A, B          - endpoints of the line segment
!   Input:  c_El          - brick element index
!   Output: Yes_Inter     - true if at least one intersection exists
!   Output: num_Inter     - number of distinct intersections found
!   Output: InterSection_P - up to 2 intersection points (3D)
!
! Notes:   Decomposes the brick into 12 triangles and tests each via
!          Tool_Intersection_of_AB_and_Triangle_3D, deduplicating hits.
!-----------------------------------------------------------

subroutine Tool_All_Intersections_of_AB_and_Brick_Ele(A, B, c_El, Yes_Inter, num_Inter, InterSection_P)
    use Global_Float_Type
    use Global_Elem_Area_Vol
    use Global_Model

    implicit none

    real(kind=FT), intent(in)  :: A(3), B(3)
    integer,       intent(in)  :: c_El
    logical,       intent(out) :: Yes_Inter
    integer,       intent(out) :: num_Inter
    real(kind=FT), intent(out) :: InterSection_P(2,3)

    real(kind=FT) :: c_X_NODES(8), c_Y_NODES(8), c_Z_NODES(8)
    real(kind=FT) :: Point1(3), Point2(3), Point3(3)
    real(kind=FT) :: c_InterSection_P(3)
    logical       :: c_Yes_Inter, Yes_Exist
    integer       :: i_Tri

    integer, save :: Triangle(12,3)
    logical, save :: Triangle_Initialized = .false.

    Yes_Inter      = .false.
    num_Inter      = 0
    InterSection_P = ZR

    if (.not. Triangle_Initialized) then
        Triangle( 1,1:3) = (/ 1, 2, 5 /)
        Triangle( 2,1:3) = (/ 2, 6, 5 /)
        Triangle( 3,1:3) = (/ 2, 3, 6 /)
        Triangle( 4,1:3) = (/ 3, 7, 6 /)
        Triangle( 5,1:3) = (/ 3, 7, 8 /)
        Triangle( 6,1:3) = (/ 3, 8, 4 /)
        Triangle( 7,1:3) = (/ 4, 8, 5 /)
        Triangle( 8,1:3) = (/ 4, 5, 1 /)
        Triangle( 9,1:3) = (/ 5, 6, 8 /)
        Triangle(10,1:3) = (/ 6, 7, 8 /)
        Triangle(11,1:3) = (/ 1, 2, 4 /)
        Triangle(12,1:3) = (/ 2, 3, 4 /)
        Triangle_Initialized = .true.
    end if

    c_X_NODES = G_X_NODES(1:8, c_El)
    c_Y_NODES = G_Y_NODES(1:8, c_El)
    c_Z_NODES = G_Z_NODES(1:8, c_El)

    do i_Tri = 1, 12

        Point1(1) = c_X_NODES(Triangle(i_Tri,1))
        Point1(2) = c_Y_NODES(Triangle(i_Tri,1))
        Point1(3) = c_Z_NODES(Triangle(i_Tri,1))

        Point2(1) = c_X_NODES(Triangle(i_Tri,2))
        Point2(2) = c_Y_NODES(Triangle(i_Tri,2))
        Point2(3) = c_Z_NODES(Triangle(i_Tri,2))

        Point3(1) = c_X_NODES(Triangle(i_Tri,3))
        Point3(2) = c_Y_NODES(Triangle(i_Tri,3))
        Point3(3) = c_Z_NODES(Triangle(i_Tri,3))

        c_Yes_Inter = .false.
call Tool_Intersection_of_AB_and_Triangle_3D(A, B, Point1, Point2, Point3, c_Yes_Inter, c_InterSection_P)

        if (c_Yes_Inter) then
            Yes_Inter = .true.
            Yes_Exist = .false.

            if (num_Inter >= 1) then
                if (all(abs(InterSection_P(1,1:3) - c_InterSection_P(1:3)) < Tol_11)) then
                    Yes_Exist = .true.
                end if
            end if

            if ((.not. Yes_Exist) .and. num_Inter >= 2) then
                if (all(abs(InterSection_P(2,1:3) - c_InterSection_P(1:3)) < Tol_11)) then
                    Yes_Exist = .true.
                end if
            end if

            if (.not. Yes_Exist) then
                if (num_Inter < 2) then
                    num_Inter = num_Inter + 1
                    InterSection_P(num_Inter,1:3) = c_InterSection_P(1:3)
                end if
            end if

        end if

    end do

    return
end subroutine Tool_All_Intersections_of_AB_and_Brick_Ele