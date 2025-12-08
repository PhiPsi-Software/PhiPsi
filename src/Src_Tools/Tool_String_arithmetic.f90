 
recursive subroutine Tool_String_arithmetic(Instruction)

use Global_Float_Type
use Global_Elem_Area_Vol
use Global_Common

implicit none
character(len=256),intent(inout)::Instruction
integer Sign_Location
integer ierror
integer Read_Value_length
integer, parameter :: rk = kind ( 1.0E+00 )
real(kind=rk) Read_Value_After,Read_Value_Before,Value_Try_1,Value_Try_2,Result_Value
integer i_Try
logical Yes_Value_1,Yes_Value_2
integer Index_Start
integer Value_Before_Start_Index
integer Value_After_End_Index
logical Yes_Found_Value_Before,Yes_Found_Value_After
character(len=256) Result_String
integer i_Calculate
character(len=1) first_Charc
integer string_length,num_added_char
character(len=256) Tem_String
integer Sign_Location_plus,Sign_Location_minus
logical logical_plus,logical_minus
character(len=256) Input_Instruction
integer Tool_chrpak_ch_index_first
integer Tool_chrpak_ch_index_first_for_plus_sign,Tool_chrpak_ch_index_first_for_minus_sign
character(len=256) Output_String_for_parentheses,parentheses_result
logical parentheses_found
integer left_index, right_index
integer I_Try_Parentheses
integer num_char_left,num_char_right

Input_Instruction = Instruction


call Tool_String_count_char_occurrences(Input_Instruction,'(', num_char_left)
call Tool_String_count_char_occurrences(Input_Instruction,')', num_char_right)
if(num_char_left/=num_char_right) then
    print *, '    Error :: number of "C" differs from number of  ")" in the *kpp file!'
    print *, "             Input_Instruction:",Input_Instruction
    print *, "             num_char_left:    ",num_char_left
    print *, "             num_char_right:   ",num_char_right
    call Warning_Message('S',Keywords_Blank)
endif
do I_Try_Parentheses =1,50
    call Tool_String_extract_content_in_parentheses(Input_Instruction,Output_String_for_parentheses,&
                                                    left_index, right_index,parentheses_found) 
    if(parentheses_found) then
        call Tool_String_arithmetic(Output_String_for_parentheses)
        parentheses_result = Output_String_for_parentheses
        
        Tem_String = Input_Instruction
        Input_Instruction(left_index:left_index) = ''
        Input_Instruction(left_index+1:left_index+len_trim(parentheses_result)) = parentheses_result(1:len_trim(parentheses_result))
        Input_Instruction(left_index+len_trim(parentheses_result)+1:left_index+len_trim(parentheses_result)+1) = ''
        Input_Instruction(left_index+len_trim(parentheses_result)+2:) = Tem_String(right_index+1:)
        
        call Tool_chrpak_s_blank_delete(Input_Instruction) 
    endif
enddo

do i_Calculate= 1,100
    Yes_Found_Value_After = .false.
    Sign_Location = Tool_chrpak_ch_index_first(Input_Instruction,'*')
    if(Sign_Location==1)then
        exit
    elseif(Sign_Location<=0) then
        exit
    endif

    if(Sign_Location>=2) then
        call Tool_chrpak_s_to_r4(Input_Instruction(Sign_Location+1:),Read_Value_After, ierror, Read_Value_length)
        if(Read_Value_length<=0)then
            print *, '    Error :: Cannot find value after sign "*"!'
            call Warning_Message('S',Keywords_Blank)
        else
            Yes_Found_Value_After = .True.
            Value_After_End_Index   = Sign_Location+1+Read_Value_length-1
            if(Input_Instruction(Value_After_End_Index:Value_After_End_Index)==",")then
                Value_After_End_Index =  Value_After_End_Index -1
            endif
        endif

        Yes_Found_Value_Before = .false.
        do i_Try=1,100
            Index_Start = Sign_Location-(i_Try)
            if (Index_Start <=0) exit
            call Tool_chrpak_s_is_r(Input_Instruction(Index_Start:Sign_Location-1), Value_Try_1,Yes_Value_1)
            Value_Before_Start_Index = Index_Start
            if(Yes_Value_1) then
                if(Index_Start==1)then
                    Read_Value_Before = Value_Try_1
                    Yes_Found_Value_Before = .True.
                    exit
                else
                    first_Charc = Input_Instruction(Index_Start-1:Index_Start-1)
                    if((first_Charc =='+') .or. (first_Charc =='-')) then
                        Read_Value_Before = Value_Try_1
                        Yes_Found_Value_Before = .True.
                        exit
                    endif
                    call Tool_chrpak_s_is_r(Input_Instruction(Index_Start-1:Sign_Location-1), Value_Try_2,Yes_Value_2)
                    if(Yes_Value_2 .eqv. .false.) then
                        Read_Value_Before = Value_Try_1
                        Yes_Found_Value_Before = .True.
                        exit
                    endif
                endif
            endif
        enddo

        if(Yes_Found_Value_Before .eqv. .false.)then
            print *, '    Error :: Cannot find value before sign "*"!'
            call Warning_Message('S',Keywords_Blank)
        endif

        if(Yes_Found_Value_Before .and. Yes_Found_Value_After) then
            Result_Value = Read_Value_Before*Read_Value_After
            call Tool_chrpak_r4_to_s_left (Result_Value,Result_String)
            Result_String = trim(adjustl(Result_String))
            string_length = Value_After_End_Index- Value_Before_Start_Index +1

            if(len_trim(Result_String) <string_length) then
                Input_Instruction(Value_Before_Start_Index:Value_After_End_Index) = Result_String
            else
                num_added_char = len_trim(Result_String) - string_length
                Tem_String = Input_Instruction
                Input_Instruction(Value_Before_Start_Index:Value_After_End_Index+num_added_char) = Result_String
                Input_Instruction(Value_After_End_Index+num_added_char+1:) = Tem_String(Value_After_End_Index+1:)
            endif

            call Tool_chrpak_s_blank_delete(Input_Instruction)
        endif
    endif
enddo

do i_Calculate= 1,100
    Yes_Found_Value_After = .false.
    Sign_Location = Tool_chrpak_ch_index_first(Input_Instruction,'/')
    if(Sign_Location==1)then
        exit
    elseif(Sign_Location<=0) then
        exit
    endif

    if(Sign_Location>=2) then
        call Tool_chrpak_s_to_r4 (Input_Instruction(Sign_Location+1:),Read_Value_After, ierror, Read_Value_length)
        if(Read_Value_length<=0)then
            print *, '    Error :: Cannot find value after sign "/"!'
            call Warning_Message('S',Keywords_Blank)
        else
            Yes_Found_Value_After = .True.
            Value_After_End_Index   = Sign_Location+1+Read_Value_length-1
            if(Input_Instruction(Value_After_End_Index:Value_After_End_Index)==",")then
                Value_After_End_Index =  Value_After_End_Index -1
            endif
        endif

        Yes_Found_Value_Before = .false.
        do i_Try=1,100
            Index_Start = Sign_Location-(i_Try)
            if (Index_Start <=0) exit
            call Tool_chrpak_s_is_r(Input_Instruction(Index_Start:Sign_Location-1), Value_Try_1,Yes_Value_1)
            Value_Before_Start_Index = Index_Start
            if(Yes_Value_1) then
                if(Index_Start==1)then
                    Read_Value_Before = Value_Try_1
                    Yes_Found_Value_Before = .True.
                    exit
                else
                    first_Charc = Input_Instruction(Index_Start-1:Index_Start-1)
                    if((first_Charc =='+') .or. (first_Charc =='-')) then
                        Read_Value_Before = Value_Try_1
                        Yes_Found_Value_Before = .True.
                        exit
                    endif
                    call Tool_chrpak_s_is_r(Input_Instruction(Index_Start-1:Sign_Location-1), Value_Try_2,Yes_Value_2)
                    if(Yes_Value_2 .eqv. .false.) then
                        Read_Value_Before = Value_Try_1
                        Yes_Found_Value_Before = .True.
                        exit
                    endif
                endif
            endif
        enddo

        if(Yes_Found_Value_Before .eqv. .false.)then
            print *, '    Error :: Cannot find value before sign "/"!'
            call Warning_Message('S',Keywords_Blank)
        endif

        if(Yes_Found_Value_Before .and. Yes_Found_Value_After) then
            Result_Value = Read_Value_Before/Read_Value_After
            call Tool_chrpak_r4_to_s_left (Result_Value,Result_String)
            Result_String = trim(adjustl(Result_String))

            string_length = Value_After_End_Index- Value_Before_Start_Index +1
            if(len_trim(Result_String) <string_length) then
                Input_Instruction(Value_Before_Start_Index:Value_After_End_Index) = Result_String
            else
                num_added_char = len_trim(Result_String) - string_length
                Tem_String = Input_Instruction
                Input_Instruction(Value_Before_Start_Index:Value_After_End_Index+num_added_char) = Result_String
                Input_Instruction(Value_After_End_Index+num_added_char+1:) = Tem_String(Value_After_End_Index+1:)
            endif

            call Tool_chrpak_s_blank_delete(Input_Instruction)

        endif
    endif
enddo

do i_Calculate= 1,100
    Yes_Found_Value_After = .false.
    Sign_Location_plus = Tool_chrpak_ch_index_first_for_plus_sign(Input_Instruction)
    
    Sign_Location_minus=  Tool_chrpak_ch_index_first_for_minus_sign(Input_Instruction)
    

    logical_plus  = .False.
    logical_minus = .False.

    if(Sign_Location_plus==Sign_Location_minus) then
        exit
    endif

    if(Sign_Location_plus<=0 .and. Sign_Location_minus >0)then
        logical_minus = .True.
        Sign_Location = Sign_Location_minus
    endif

    if(Sign_Location_plus > 0 .and. Sign_Location_minus<=0)then
        logical_plus  = .True.
        Sign_Location = Sign_Location_plus
    endif

    if(Sign_Location_plus > 0 .and. Sign_Location_minus>0) then
        if(Sign_Location_plus < Sign_Location_minus) then
            logical_plus  = .True.
            Sign_Location = Sign_Location_plus
        elseif(Sign_Location_plus > Sign_Location_minus) then
            logical_minus = .True.
            Sign_Location = Sign_Location_minus
        endif
    endif
    
    if(Sign_Location_plus > 0 .and. Sign_Location_minus>0) then
        if(Sign_Location_minus==1)then
            logical_plus  = .True.
            Sign_Location = Sign_Location_plus
        endif
    endif
    
    
    if(Sign_Location_minus==1) then
        logical_minus  = .True.
    endif


    if(Sign_Location==1)then

    elseif(Sign_Location<=0) then
        exit
    endif

    
    if(Sign_Location>=2) then

        call Tool_chrpak_s_to_r4(Input_Instruction(Sign_Location+1:),Read_Value_After, ierror, Read_Value_length)

        if(Read_Value_length<=0)then
            if(logical_plus) then
                print *, '    Error :: Cannot find value after sign "+"!'
            elseif(logical_minus)then
                print *, '    Error :: Cannot find value after sign "-"!'
            endif
            call Warning_Message('S',Keywords_Blank)
        else
            Yes_Found_Value_After = .True.
            Value_After_End_Index   = Sign_Location+1+Read_Value_length-1
            if(Input_Instruction(Value_After_End_Index:Value_After_End_Index)==",")then
                Value_After_End_Index =  Value_After_End_Index -1
            endif
        endif

        Yes_Found_Value_Before = .false.
        do i_Try=1,100
            Index_Start = Sign_Location-(i_Try)
            if (Index_Start <=0) exit
            call Tool_chrpak_s_is_r(Input_Instruction(Index_Start:Sign_Location-1), Value_Try_1,Yes_Value_1)
            Value_Before_Start_Index = Index_Start
            if(Yes_Value_1) then
                if(Index_Start==1)then
                    Read_Value_Before = Value_Try_1
                    Yes_Found_Value_Before = .True.
                    exit
                else
                    call Tool_chrpak_s_is_r(Input_Instruction(Index_Start-1:Sign_Location-1), Value_Try_2,Yes_Value_2)
                    if(Yes_Value_2 .eqv. .false.) then
                        Read_Value_Before = Value_Try_1
                        Yes_Found_Value_Before = .True.
                        exit
                    endif
                endif
            endif
        enddo

        if(Yes_Found_Value_Before .eqv. .false.)then
            if(logical_plus) then
                print *, '    Error :: Cannot find value before sign "+"!'
                print *, "             Input_Instruction:",Input_Instruction(1:len_trim(Input_Instruction))
            elseif(logical_minus)then
                print *, '    Error :: Cannot find value before sign "-"!'
                print *, "             Input_Instruction:",Input_Instruction(1:len_trim(Input_Instruction))
            endif
            call Warning_Message('S',Keywords_Blank)
        endif

        if(Yes_Found_Value_Before .and. Yes_Found_Value_After) then
            if(logical_plus) then
                Result_Value = Read_Value_Before + Read_Value_After
            elseif(logical_minus)then
                Result_Value = Read_Value_Before - Read_Value_After
            endif
            call Tool_chrpak_r4_to_s_left (Result_Value,Result_String)
            Result_String = trim(adjustl(Result_String))


            string_length = Value_After_End_Index- Value_Before_Start_Index +1

            if(len_trim(Result_String) <string_length) then
                Input_Instruction(Value_Before_Start_Index:Value_After_End_Index) = Result_String
            else
                num_added_char = len_trim(Result_String) - string_length
                Tem_String = Input_Instruction
                Input_Instruction(Value_Before_Start_Index:Value_After_End_Index+num_added_char) = Result_String
                Input_Instruction(Value_After_End_Index+num_added_char+1:) = Tem_String(Value_After_End_Index+1:)
            endif

            call Tool_chrpak_s_blank_delete(Input_Instruction)

        endif
    endif
enddo

Instruction = Input_Instruction

return 
end SUBROUTINE Tool_String_arithmetic                       
