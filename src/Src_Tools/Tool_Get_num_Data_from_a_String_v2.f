 
      subroutine Tool_Get_num_Data_from_a_String_v2(c_String,num_Char,
     &                                           num_Data,Yes_Even)
       
      use Global_Float_Type
      IMPLICIT NONE
      integer,intent(in)::num_Char
      character(*),intent(in)::c_String
      integer,intent(out)::num_Data
      logical,intent(out)::Yes_Even
      integer i_Char,num_Space
      
      
      
      
      num_Space = 0
      do i_Char = 1,num_Char-1
          if(c_String(i_Char:i_Char).eq.' ' .and. 
     &       c_String(i_Char+1:i_Char+1).ne.' ' )then
              num_Space = num_Space +1
          endif
      enddo
      num_Data = num_Space + 1
      
      Yes_Even = .False.
      if(mod(num_Data,2)==0) Yes_Even = .True.
      
      RETURN
      end subroutine Tool_Get_num_Data_from_a_String_v2 