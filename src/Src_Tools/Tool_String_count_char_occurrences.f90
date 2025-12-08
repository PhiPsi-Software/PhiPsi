 
recursive subroutine Tool_String_count_char_occurrences(input_string, input_chr, num_char)

character(len=*), intent(in) :: input_string
character(len=1), intent(in) :: input_chr
integer, intent(out) :: num_char

integer :: i, len_input
character(len=1) :: current_char

len_input = len_trim(input_string)
num_char = 0

do i = 1, len_input
    current_char = input_string(i:i)
    if (current_char == input_chr) then
        num_char = num_char + 1
    endif
end do


return 
end SUBROUTINE Tool_String_count_char_occurrences                  
