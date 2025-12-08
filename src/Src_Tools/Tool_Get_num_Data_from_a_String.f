 
      subroutine Tool_Get_num_Data_from_a_String(c_String,num_Char,
     &                                           num_Data,Yes_Even)
       
      use Global_Float_Type
      IMPLICIT NONE
      integer,intent(in)::num_Char
      character(len=num_Char),intent(in)::c_String
      integer,intent(out)::num_Data
      logical,intent(out)::Yes_Even
      integer i_Char,num_Comma
      
      num_Comma = 0
      do i_Char = 1,num_Char
          if(c_String(i_Char:i_Char).eq.',')then
              num_Comma = num_Comma +1
          endif
      enddo
      num_Data = num_Comma + 1
      
      Yes_Even = .False.
      if(mod(num_Data,2)==0) Yes_Even = .True.
      
      RETURN
      end subroutine Tool_Get_num_Data_from_a_String