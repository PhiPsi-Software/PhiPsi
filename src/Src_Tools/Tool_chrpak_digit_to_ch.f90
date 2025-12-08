 
subroutine Tool_chrpak_digit_to_ch(digit, ch)

  implicit none

  character ch
  integer digit

  if ( 0 <= digit .and. digit <= 9 ) then

    ch = achar ( digit + 48 )

  else

    ch = '*'

  end if

  return
end subroutine Tool_chrpak_digit_to_ch

