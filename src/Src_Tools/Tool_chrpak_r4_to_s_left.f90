 
subroutine Tool_chrpak_r4_to_s_left ( r4, s )

  implicit none

  integer, parameter :: rk = kind ( 1.0E+00 )

  real ( kind = rk ) r4
  character ( len = * ) s
  character ( len = 14 ) s2
  integer Sign_Location
  integer Tool_chrpak_ch_index_first

  if ( real ( int ( r4 ), kind = rk ) == r4 ) then
    write ( s2, '(i14)' ) int ( r4 )
  else if ( abs ( r4 ) < 999999.5E+00 ) then
    write ( s2, '(f14.6)' ) r4
  else
    write ( s2, '(g14.6)' ) r4
  end if

  s = adjustl ( s2 )
  
  Sign_Location = Tool_chrpak_ch_index_first(s,'E')
  if(Sign_Location>=2) then
      if(s(Sign_Location+1:Sign_Location+1)=="+") then
          s(Sign_Location+1:Sign_Location+1)=""
          call Tool_chrpak_s_blank_delete(s)
      endif
  endif
  
  Sign_Location = Tool_chrpak_ch_index_first(s,'e')
  if(Sign_Location>=2) then
      if(s(Sign_Location+1:Sign_Location+1)=="+") then
          s(Sign_Location+1:Sign_Location+1)=""
          call Tool_chrpak_s_blank_delete(s)
      endif
  endif  
  
  return
end subroutine Tool_chrpak_r4_to_s_left

