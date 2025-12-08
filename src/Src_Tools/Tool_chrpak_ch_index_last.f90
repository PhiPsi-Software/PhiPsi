 
function Tool_chrpak_ch_index_last( s, ch )

  implicit none

  character ch
  integer Tool_chrpak_ch_index_last
  integer i
  character ( len = * ) s
  integer s_length

  Tool_chrpak_ch_index_last = -1
  s_length = len_trim ( s )

  do i = s_length, 1, -1

    if ( s(i:i) == ch ) then
      Tool_chrpak_ch_index_last = i
      return
    end if

  end do

  return
end function Tool_chrpak_ch_index_last