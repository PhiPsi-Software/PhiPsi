 
      subroutine Tool_Function_BETA(P,Q,BT)

      use Global_Float_Type
      IMPLICIT real(kind=FT) (A-H,O-Z)
      
      CALL Tool_Function_GAMMA(P,GP)
      CALL Tool_Function_GAMMA(Q,GQ)
      PPQ=P+Q
      CALL Tool_Function_GAMMA(PPQ,GPQ)
      BT=GP*GQ/GPQ
      
      return 
      end SUBROUTINE Tool_Function_BETA                         
