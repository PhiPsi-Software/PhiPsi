 
subroutine Tool_chrpak_s_replace_i ( s, sub1, sub2 )

  implicit none

  integer ilo
  integer len1
  integer len2
  integer lens
  integer Tool_chrpak_s_indexi
  character ( len = * ) s
  character ( len = * ) sub1
  character ( len = * ) sub2

  lens = len ( s )

  if ( lens <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_REPLACE_I - Serious error!'
    write ( *, '(a)' ) '  Null string not allowed!'
    return
  end if

  len1 = len ( sub1 )

  if ( len1 <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_REPLACE_I - Serious error!'
    write ( *, '(a)' ) '  Null SUB1 not allowed!'
    return
  end if

  len2 = len ( sub2 )

  ilo = Tool_chrpak_s_indexi ( s, sub1 )
  if ( ilo /= 0 ) then
    s(ilo+len2:lens+len2-len1) = s(ilo+len1:lens)
    s(ilo:ilo+len2-1) = sub2(1:len2)
  end if

  return
end subroutine Tool_chrpak_s_replace_i


