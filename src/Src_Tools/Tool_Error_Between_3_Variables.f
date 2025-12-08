 
      subroutine Tool_Error_Between_3_Variables(A,B,C,Error)

      use Global_Float_Type      
      implicit none
      real(kind=FT),intent(in):: A,B,C
      real(kind=FT),intent(out):: Error
      real(kind=FT) Error_AB,Error_BC,Error_AC
      
      if (max(A,B) ==ZR) then
          Error_AB = ZR
      else
          Error_AB = abs(A-B) / max(A,B)
      endif
      
      if (max(B,C) ==ZR) then
          Error_BC = ZR
      else
          Error_BC = abs(B-C) / max(B,C)
      endif
      
      if (max(A,C) ==ZR) then
          Error_AC= ZR
      else
          Error_AC = abs(A-C) / max(A,C)
      endif      
      
      Error = max(Error_AB,Error_BC,Error_AC)
      
      return 
      end SUBROUTINE Tool_Error_Between_3_Variables                         
