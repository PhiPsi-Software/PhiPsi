!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !
 
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

  dot = sum ( ( p1(1:dim_num) - p2(1:dim_num) ) &
            * ( p3(1:dim_num) - p2(1:dim_num) ) )

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
  

  t = sum ( ( p(1:dim_num) - p1(1:dim_num) ) &
          * ( p2(1:dim_num) - p1(1:dim_num) ) ) / bot
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