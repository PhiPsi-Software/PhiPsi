 
subroutine Tool_chrpak_s_is_r ( s, r, lval )

  implicit none

  integer, parameter :: rk = kind ( 1.0E+00 )

  integer ierror
  integer length
  logical lval
  real ( kind = rk ) r
  character ( len = * ) s
  integer s_length

  s_length = len_trim ( s )

  call Tool_chrpak_s_to_r4 ( s, r, ierror, length )

  if ( ierror == 0 .and. s_length <= length ) then
    lval = .true.
  else
    lval = .false.
    r = 0.0E+00
  end if

  return
end subroutine Tool_chrpak_s_is_r