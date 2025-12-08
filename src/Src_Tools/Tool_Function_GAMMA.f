 
      subroutine Tool_Function_GAMMA(X,GA)
      use Global_Float_Type
      real(kind=FT),intent(in)::X
      real(kind=FT),intent(out)::GA

      GA = GAMMA(X)


      return
      end SUBROUTINE Tool_Function_GAMMA
