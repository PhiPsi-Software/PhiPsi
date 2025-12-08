 
function Tool_chrpak_s_digits_count ( s )


  implicit none

  character c
  logical Tool_chrpak_ch_is_digit
  integer n
  character ( len = * ) s
  integer Tool_chrpak_s_digits_count
  integer s_len
  integer s_pos

  s_len = len_trim ( s )

  s_pos = 0
  n = 0

  do while ( s_pos < s_len )

    s_pos = s_pos + 1

    c = s(s_pos:s_pos)

    if ( Tool_chrpak_ch_is_digit ( c ) ) then
      n = n + 1
    end if

  end do

  Tool_chrpak_s_digits_count = n

  return

end function Tool_chrpak_s_digits_count