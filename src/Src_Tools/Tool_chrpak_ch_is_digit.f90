 
function Tool_chrpak_ch_is_digit ( ch )


  implicit none

  character ch
  logical Tool_chrpak_ch_is_digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then
    Tool_chrpak_ch_is_digit = .true.
  else
    Tool_chrpak_ch_is_digit = .false.
  end if

  return
end function Tool_chrpak_ch_is_digit