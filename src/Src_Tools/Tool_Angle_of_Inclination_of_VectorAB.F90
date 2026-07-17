!-----------------------------------------------------------
! Brief: Computes the inclination angle (in degrees, 0-360) of the 2D vector from A to B.
!
! Parameters:
!   Input:  A, B  - 2D point coordinates
!   Output: Angle - inclination angle in degrees, range 0 to 360
!
! Notes:   Branches on slope sign and quadrant to convert the principal arctangent
!          result into the full 0-360 degree range, and guards against NaN.
!-----------------------------------------------------------

subroutine Tool_Angle_of_Inclination_of_VectorAB(A,B,Angle)
! Calculate the inclination of line segment AB, the inclination of the vector pointing from
! A to B (in degrees), with the range of the inclination being 0-360 degrees.
! The arctangent function directly gives a tilt angle range of (-90, 90) degrees.
! https://www.easycalculation.com/analytical/slope-length-equation.php
! Formula:
! m = (Y2-Y1) / (X2-X1)
! Distance = Square Root ( (X2-X1)2 + (Y2-Y1)2 )
! Angle = arctan ( m )
! Line of Equation is y = mx + b
! Where, m = slope

use Global_Float_Type

implicit none
real(kind=FT),intent(in)  :: A(2),B(2)
real(kind=FT),intent(out) :: Angle
real(kind=FT) :: X1,Y1,X2,Y2,m
real(kind=FT) :: c_Angle

X1 = A(1)
Y1 = A(2)
X2 = B(1)
Y2 = B(2)

if (abs(X2-X1)<=Tol_11)then
    m = (Y2-Y1) / Tol_11
else
    m = (Y2-Y1) / (X2-X1)
endif
c_Angle = atan(m) *Con_180/pi

if (abs(X2-X1)<=Tol_11)then
    if(Y1<=Y2)then
        Angle = Con_90
    else
        Angle = Con_270
    endif
elseif(abs(X2-X1)>Tol_11)then
    if(X1 <= X2) then
        if(Y1<=Y2)then
            Angle = c_Angle
        else
            Angle = Con_360 - abs(c_Angle)
        endif
    elseif(X1 > X2) then
        if(Y1<=Y2)then
            Angle = Con_180 - abs(c_Angle)
        else
            Angle = Con_180 + abs(c_Angle)
        endif
    endif
endif

if (isnan(Angle)) then
    print *, '    Error :: Tool_Angle_of_Inclination_of_VectorAB.f'
endif

return
end SUBROUTINE Tool_Angle_of_Inclination_of_VectorAB
