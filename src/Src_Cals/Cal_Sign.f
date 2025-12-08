 
      subroutine Cal_Sign(Variable,Sign_O)
      
      use Global_Float_Type     
      implicit none
      
      real(kind=FT),intent(in)::Variable
      real(kind=FT),intent(out)::Sign_O

      if(Variable.eq.ZR) then
          Sign_O = ZR
      elseif (Variable > ZR) then
          Sign_O = ONE
      elseif (Variable < ZR) then
          Sign_O = -ONE      
      end if
      
      return 
      end SUBROUTINE Cal_Sign                
