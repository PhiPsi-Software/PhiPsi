 
function Tool_chrpak_ch_eqi_chrpak ( c1, c2 )

  implicit none

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical Tool_chrpak_ch_eqi_chrpak

  c1_cap = c1
  c2_cap = c2

  call Tool_chrpak_ch_cap_chrpak ( c1_cap )
  call Tool_chrpak_ch_cap_chrpak ( c2_cap )

  if ( c1_cap == c2_cap ) then
    Tool_chrpak_ch_eqi_chrpak = .true.
  else
    Tool_chrpak_ch_eqi_chrpak = .false.
  end if

  return
end function Tool_chrpak_ch_eqi_chrpak

