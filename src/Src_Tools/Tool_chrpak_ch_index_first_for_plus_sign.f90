 
function Tool_chrpak_ch_index_first_for_plus_sign(s)

  implicit none

  character(len=1) ch
  integer Tool_chrpak_ch_index_first_for_plus_sign
  integer i
  character ( len = * ) s
  integer s_length
  
  ch = '+'
  
  Tool_chrpak_ch_index_first_for_plus_sign = - 1
  s_length = len_trim ( s )

  do i = 1, s_length

    if ( s(i:i) == ch ) then
      
      if(i==1) then
      else
          if(s(i-1:i-1)=='E' .or. s(i-1:i-1)=='D' .or. s(i-1:i-1)=='e' .or. s(i-1:i-1)=='d' )then
          else
              Tool_chrpak_ch_index_first_for_plus_sign = i
              return
          endif
      endif
      
      
    end if

  end do

  return
end function Tool_chrpak_ch_index_first_for_plus_sign