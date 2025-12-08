 
subroutine Tool_chrpak_r8_to_s_left ( r8, s )

  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r8
  character ( len = * ) s
  integer s_length
  character ( len = 14 ) s2

  s_length = len ( s )

  if ( s_length < 14 ) then

    do i = 1, s_length
      s(i:i) = '*'
    end do

  else if ( r8 == 0.0D+00 ) then
    s(1:14) = '     0.0      '
  else
    write ( s2, '(g14.6)' ) r8
    s(1:14) = s2
  end if
  s = adjustl ( s )

  return
end subroutine Tool_chrpak_r8_to_s_left

