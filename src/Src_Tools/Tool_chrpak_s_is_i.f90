 
function Tool_chrpak_s_is_i ( s, i )

  implicit none

  integer i
  integer ierror
  integer length
  character ( len = * ) s
  logical Tool_chrpak_s_is_i
  integer s_length

  s_length = len_trim ( s )

  call Tool_chrpak_s_to_i4 ( s, i, ierror, length )

  if ( ierror == 0 .and. s_length <= length ) then
    Tool_chrpak_s_is_i = .true.
  else
    Tool_chrpak_s_is_i = .false.
    i = 0
  end if

  return
end function Tool_chrpak_s_is_i