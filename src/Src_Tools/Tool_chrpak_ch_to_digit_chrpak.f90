 
subroutine Tool_chrpak_ch_to_digit_chrpak ( ch, digit )

  implicit none

  character ch
  integer digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then

    digit = iachar ( ch ) - 48

  else

    digit = - 1

  end if

  return
end subroutine Tool_chrpak_ch_to_digit_chrpak
