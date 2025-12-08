 
function Tool_chrpak_ch_index_first ( s, ch )

  implicit none

  character ch
  integer Tool_chrpak_ch_index_first
  integer i
  character ( len = * ) s
  integer s_length

  Tool_chrpak_ch_index_first = - 1
  s_length = len_trim ( s )

  do i = 1, s_length

    if ( s(i:i) == ch ) then
      Tool_chrpak_ch_index_first = i
      return
    end if

  end do

  return
end function Tool_chrpak_ch_index_first