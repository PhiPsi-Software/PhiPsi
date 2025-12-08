 
subroutine Tool_chrpak_s_blank_delete ( s )

  implicit none

  character ch
  integer get
  integer put
  character ( len = * ) s
  integer s_length
  character, parameter :: tab = achar ( 9 )

  put = 0
  s_length = len_trim ( s )

  do get = 1, s_length

    ch = s(get:get)

    if ( ch /= ' ' .and. ch /= tab ) then
      put = put + 1
      s(put:put) = ch
    end if

  end do

  s(put+1:s_length) = ' '

  return
end subroutine Tool_chrpak_s_blank_delete

