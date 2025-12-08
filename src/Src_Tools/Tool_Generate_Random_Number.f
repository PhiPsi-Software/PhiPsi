 
      subroutine Tool_Generate_Random_Number(out_number)
      use Global_Float_Type 
      use Global_Common
      
      
 
 
 
      
      implicit none
      
      real(kind=FT),intent(out):: out_number
      integer k
      
      
      k = Seed / 127773

      Seed = 16807 * ( Seed - k * 127773 ) - k * 2836

      if (Seed .lt. 0) then
        Seed = Seed+ 2147483647
      end if
      
      out_number = dble ( Seed ) * 4.656612875D-10
      
      return 
      end SUBROUTINE Tool_Generate_Random_Number                        
