 
recursive subroutine Tool_String_extract_content_in_parentheses(Input_String,Output_String,left_index,right_index,found) 

implicit none  

character(len=*), intent(in) :: input_string
character(len=*), intent(out) :: output_string
integer, intent(out) :: left_index, right_index
logical, intent(out) :: found

integer :: i, len_input, stack_count
character(len=1) :: current_char
character(len=256) :: stack
character(len=256) :: Temp_string,Temp_output_string
integer Temp_left_index,Temp_right_index
logical Temp_found

len_input = len_trim(input_string)
stack = ''
stack_count = 0
found = .false.

do i = 1, len_input
    current_char = input_string(i:i)
    if (current_char == '(') then
        if (stack_count == 0) then
            left_index = i
        endif
        stack = trim(stack) // '('
        stack_count = stack_count + 1
    elseif (current_char == ')') then
        if (stack_count > 0) then
            stack_count = stack_count - 1
            if (stack_count == 0) then
                stack = adjustl(stack)
                right_index = i
                found = .true.
                exit
            else
                stack = trim(stack) // ')'
            endif
        endif
    elseif (stack_count > 0) then
        stack = trim(stack) // current_char
    endif
end do

if (found) then
    output_string = input_string(left_index + 1:right_index - 1)
    Temp_string(1:len_trim(output_string)) = output_string(1:len_trim(output_string))
    call Tool_String_extract_content_in_parentheses(Temp_string,Temp_Output_String,Temp_left_index,Temp_right_index,Temp_found) 
    if(Temp_found) then
        left_index  = left_index  + Temp_left_index
        right_index = left_index + (Temp_right_index-Temp_left_index)
        found = Temp_found
        output_string = input_string(left_index + 1:right_index - 1)
    endif
    
else
    output_string = ''
endif


return 
end SUBROUTINE Tool_String_extract_content_in_parentheses                   
