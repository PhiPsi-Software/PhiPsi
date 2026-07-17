!-----------------------------------------------------------
! Brief: Determine which of 8 octants contains a given 3D point.
!
! Parameters:
!   Input:  point         - 3D point coordinates (x,y,z)
!   Output: Domain_Number - octant index in 1..8
!
! Notes:   The model center Model_Center_* splits the model into
!          2x2x2=8 partitions; the closest octant to the point
!          is returned. Aborts on points lying outside all octants.
!-----------------------------------------------------------

subroutine D3_Get_Point_Domain_Number(point,Domain_Number)
! Get the partition number of the point.
!2022-11-24. NEWFTU2022112401.


!***************
! Public Module
!***************
use Global_Float_Type
use Global_Model
use Global_Common

!**********************
! Variable Declaration
!**********************
implicit none
real(kind=FT),intent(in)::point(3)
integer,intent(out):: Domain_Number
real(kind=FT) c_Tol
real(kind=FT) x,y,z


!*****************************************
! Get the partition number of the element
!*****************************************
x = point(1)
y = point(2)
z = point(3)

c_Tol    = Tol_7
if( ((x - Model_Center_x) > -c_Tol) .and. ((y - Model_Center_y) > -c_Tol) .and. ((z - Model_Center_z) > -c_Tol)) then
    Domain_Number  = 1
    return
elseif( ((x - Model_Center_x) <  c_Tol) .and. ((y - Model_Center_y) >-c_Tol) .and. ((z - Model_Center_z) > -c_Tol)) then
    Domain_Number  = 2
    return
elseif( ((x - Model_Center_x) <  c_Tol) .and. ((y - Model_Center_y) < c_Tol) .and. ((z - Model_Center_z) > -c_Tol)) then
    Domain_Number  = 3
    return
elseif( ((x - Model_Center_x) > -c_Tol) .and. ((y - Model_Center_y) < c_Tol) .and. ((z - Model_Center_z) > -c_Tol)) then
    Domain_Number  = 4
    return
elseif( ((x - Model_Center_x) > -c_Tol) .and. ((y - Model_Center_y) >-c_Tol) .and. ((z - Model_Center_z) <  c_Tol)) then
    Domain_Number  = 5
    return
elseif( ((x - Model_Center_x) <  c_Tol) .and. ((y - Model_Center_y) >-c_Tol) .and. ((z - Model_Center_z) <  c_Tol)) then
    Domain_Number  = 6
    return
elseif( ((x - Model_Center_x) <  c_Tol) .and. ((y - Model_Center_y) < c_Tol) .and. ((z - Model_Center_z) <  c_Tol)) then
    Domain_Number  = 7
    return
elseif( ((x - Model_Center_x) > -c_Tol) .and. ((y - Model_Center_y) < c_Tol) .and. ((z - Model_Center_z) <  c_Tol)) then
    Domain_Number  = 8
    return
else
    print *, '    Error :: cannot get domain number in D3_Get_Point_Domain_Number.f90!'
    call Warning_Message('S',Keywords_Blank)
endif

    
RETURN
END subroutine D3_Get_Point_Domain_Number