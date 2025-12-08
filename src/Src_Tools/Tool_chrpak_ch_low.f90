 
subroutine Tool_chrpak_ch_low ( ch )


implicit none

character ch
integer i

i = iachar ( ch )

if ( 65 <= i .and. i <= 90 ) then
    ch = achar ( i + 32 )
end if

return

end subroutine Tool_chrpak_ch_low