 
function Tool_chrpak_s_indexi( s, sub )

  implicit none

  integer i
  integer llen2
  character ( len = * ) s
  logical Tool_chrpak_s_eqi
  integer Tool_chrpak_s_indexi
  integer s_length
  character ( len = * ) sub

  Tool_chrpak_s_indexi = 0

  s_length = len_trim ( s )
  llen2 = len_trim ( sub )
  if ( s_length == 0 ) then
    s_length = len ( s )
  end if

  if ( llen2 == 0 ) then
    llen2 = len ( sub )
  end if

  if ( s_length < llen2 ) then
    return
  end if

  do i = 1, s_length + 1 - llen2

    if ( Tool_chrpak_s_eqi ( s(i:i+llen2-1), sub ) ) then
      Tool_chrpak_s_indexi = i
      return
    end if

  end do

  return
end function Tool_chrpak_s_indexi

