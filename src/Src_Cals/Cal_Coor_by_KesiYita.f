 
      subroutine Cal_Coor_by_KesiYita(kesi,yita,X_NODES,Y_NODES,
     &                                Out_x,Out_y)
      use Global_Float_Type
      implicit none
      real(kind=FT),intent(in)::kesi,yita,X_NODES(4),Y_NODES(4)
      real(kind=FT),intent(out)::Out_x,Out_y
      real(kind=FT) N(4)
      
      N  = 0.25D0 * [(1-kesi)*(1-yita),(1+kesi)*(1-yita),
     &               (1+kesi)*(1+yita),(1-kesi)*(1+yita)]
      
      Out_x  = N(1)*X_NODES(1)+N(2)*X_NODES(2)+
     &         N(3)*X_NODES(3)+N(4)*X_NODES(4)
      Out_y  = N(1)*Y_NODES(1)+N(2)*Y_NODES(2)+
     &         N(3)*Y_NODES(3)+N(4)*Y_NODES(4)
     
      return 
      end SUBROUTINE Cal_Coor_by_KesiYita                          
