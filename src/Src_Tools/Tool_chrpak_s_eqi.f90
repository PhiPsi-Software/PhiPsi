 
function Tool_chrpak_s_eqi( s1, s2 )

  implicit none

  character c1
  character c2
  integer i
  integer lenc
  logical Tool_chrpak_s_eqi
  character ( len = *  ) s1
  integer s1_length
  character ( len = *  ) s2
  integer s2_length

  s1_length = len ( s1 )
  s2_length = len ( s2 )
  lenc = min ( s1_length, s2_length )

  Tool_chrpak_s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call Tool_chrpak_ch_cap ( c1 )
    call Tool_chrpak_ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc + 1, s1_length
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, s2_length
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  Tool_chrpak_s_eqi = .true.

  return
end function Tool_chrpak_s_eqi