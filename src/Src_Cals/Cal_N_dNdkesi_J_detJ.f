 
      subroutine Cal_N_dNdkesi_J_detJ(kesi,yita,
     &                                X_NODES,Y_NODES,
     &                                detJ,J,N,dNdkesi)
     
      use Global_Float_Type      
      implicit none
      
      real(kind=FT),intent(in)::kesi,yita,X_NODES(4),Y_NODES(4)
      real(kind=FT),intent(out)::detJ,dNdkesi(4,2)
      real(kind=FT),intent(out):: N(2,8)
      real(kind=FT) N1,N2,N3,N4
      real(kind=FT) o
      
      real(kind=FT) temp(2,4),Coor(2,4),J(2,2) 
      
      o   = ZR
      
      Coor(1,:) = [X_NODES(1),X_NODES(2),X_NODES(3),X_NODES(4)]
      Coor(2,:) = [Y_NODES(1),Y_NODES(2),Y_NODES(3),Y_NODES(4)]
      
      N1 = (ONE-kesi)*(ONE-yita)/FOU
      N2 = (ONE+kesi)*(ONE-yita)/FOU
      N3 = (ONE+kesi)*(ONE+yita)/FOU
      N4 = (ONE-kesi)*(ONE+yita)/FOU
      N(1,:) = [N1,o,N2,o,N3,o,N4,o]
      N(2,:) = [o,N1,o,N2,o,N3,o,N4]

      temp(1,1:4) =[(yita-ONE)/FOU,(-yita+ONE)/FOU,
     &                    (yita+ONE)/FOU,(-yita-ONE)/FOU]
      temp(2,1:4) =[(kesi-ONE)/FOU,(-kesi-ONE)/FOU,
     &                    (kesi+ONE)/FOU,(-kesi+ONE)/FOU]
              
      dNdkesi = transpose(temp)

      J = MATMUL(Coor,dNdkesi)
      
      call Matrix_Det_2x2(J,detJ)  
      
      return 
      end SUBROUTINE Cal_N_dNdkesi_J_detJ                     
