 
subroutine Tool_chrpak_ch_cap ( ch )

  implicit none

  character ch
  integer itemp

  itemp = iachar ( ch )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    ch = achar ( itemp - 32 )
  end if

  return
  
end subroutine Tool_chrpak_ch_cap

