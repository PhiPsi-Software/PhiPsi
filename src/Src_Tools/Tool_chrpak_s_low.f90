 
subroutine Tool_chrpak_s_low ( s )



implicit none

integer i
character ( len = * ) s
integer s_length

s_length = len_trim ( s )

do i = 1, s_length
    call Tool_chrpak_ch_low ( s(i:i) )
end do

return

end subroutine Tool_chrpak_s_low