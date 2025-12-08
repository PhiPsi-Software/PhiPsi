 
      subroutine Vector_Range_Number_of_Value_Dou(value_to_be_find,n,
     &               values,Yes_In,
     &               Range_Num)   
    
      use Global_Float_Type
      implicit none

      integer,intent(in)::n
      real(kind=FT),intent(in)::value_to_be_find,values(n)
      integer,intent(out)::Range_Num
      logical,intent(out)::Yes_In
      integer i               
      
      Yes_In = .False.

      do i=1,n-1
          if(value_to_be_find >= values(i)   .and. 
     &       value_to_be_find <= values(i+1)  )then
              Range_Num = i
              Yes_In = .True.
              return
          endif
      end do
      
      return 
      end subroutine Vector_Range_Number_of_Value_Dou
    


