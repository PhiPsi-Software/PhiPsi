function angle_rad_3d ( p1, p2, p3 )

!*****************************************************************************80
!
!! ANGLE_RAD_3D returns the angle in radians between two rays in 3D.
!
!  Discussion:
!
!    The routine always computes the SMALLER of the two angles between
!    two rays.  Thus, if the rays make an (exterior) angle of
!    1.5 pi radians, the (interior) angle of 0.5 pi radians will be reported.
!
!    X dot Y = Norm(X) * Norm(Y) * Cos ( Angle(X,Y) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = FT ) P1(3), P2(3), P3(3), points defining an angle.
!    The rays are P1 - P2 and P3 - P2.
!
!    Output, real ( kind = FT ) ANGLE_RAD_3D, the angle between the two rays,
!    in radians.  This value will always be between 0 and PI.  If either ray has
!    zero length, then the angle is returned as zero.
use Global_Float_Type
implicit none

integer, parameter :: dim_num = 3

real ( kind = FT ) angle_rad_3d
real ( kind = FT ) dot
real ( kind = FT ) p1(dim_num)
real ( kind = FT ) p2(dim_num)
real ( kind = FT ) p3(dim_num)
real ( kind = FT ) r8_acos
real ( kind = FT ) v1norm
real ( kind = FT ) v2norm

v1norm = sqrt ( sum ( ( p1(1:dim_num) - p2(1:dim_num) )**2 ) )

if ( v1norm == 0.0D+00 ) then
    angle_rad_3d = 0.0D+00
    return
end if

v2norm = sqrt ( sum ( ( p3(1:dim_num) - p2(1:dim_num) )**2 ) )

if ( v2norm == 0.0D+00 ) then
    angle_rad_3d = 0.0D+00
    return
end if

dot = sum ( ( p1(1:dim_num) - p2(1:dim_num) ) * ( p3(1:dim_num) - p2(1:dim_num) ) )

angle_rad_3d = r8_acos ( dot / ( v1norm * v2norm ) )

return
end function angle_rad_3d 


subroutine line_exp_point_near_3d ( p1, p2, p, pn, dist, t )

use Global_Float_Type
implicit none

integer, parameter :: dim_num = 3

real ( kind = FT ) bot
real ( kind = FT ) dist
logical line_exp_is_degenerate_nd
real ( kind = FT ) p(dim_num)
real ( kind = FT ) p1(dim_num)
real ( kind = FT ) p2(dim_num)
real ( kind = FT ) pn(dim_num)
real ( kind = FT ) t

if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_EXP_POINT_NEAR_3D - Fatal error!'
    write ( *, '(a)' ) '  The line is degenerate.'
    stop 1
end if
bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )


t = sum ( ( p(1:dim_num) - p1(1:dim_num) ) * ( p2(1:dim_num) - p1(1:dim_num) ) ) / bot
pn(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )
dist = sqrt ( sum ( ( pn(1:dim_num) - p(1:dim_num) )**2 ) )

return
end subroutine line_exp_point_near_3d



function line_exp_is_degenerate_nd ( dim_num, p1, p2 )

use Global_Float_Type
implicit none

integer dim_num

logical line_exp_is_degenerate_nd
real ( kind = FT ) p1(dim_num)
real ( kind = FT ) p2(dim_num)

line_exp_is_degenerate_nd = ( all ( p1(1:dim_num) == p2(1:dim_num) ) )

return
end function line_exp_is_degenerate_nd 



function r8_acos ( c )

use Global_Float_Type
implicit none

real ( kind = FT ) c
real ( kind = FT ) c2
real ( kind = FT ) r8_acos

c2 = c
c2 = max ( c2, -1.0D+00 )
c2 = min ( c2, +1.0D+00 )

r8_acos = acos ( c2 )

return
end function r8_acos
