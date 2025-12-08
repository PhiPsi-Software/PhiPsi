 
      subroutine Cal_detJ(kesi,yita,X_NODES,Y_NODES,detJ)
     
      use Global_Float_Type
      implicit none
      
      real(kind=FT),intent(in)::kesi,yita,X_NODES(4),Y_NODES(4)
      real(kind=FT),intent(out)::detJ       
      real(kind=FT) temp(2,4),Coor(2,4),J(2,2),dNdkesi(4,2)
      
      
      Coor(1,1:4) = [X_NODES(1),X_NODES(2),X_NODES(3),X_NODES(4)]
      Coor(2,1:4) = [Y_NODES(1),Y_NODES(2),Y_NODES(3),Y_NODES(4)]

      temp(1,:) =[(yita-ONE)/FOU,(-yita+ONE)/FOU,
     &            (yita+ONE)/FOU,(-yita-ONE)/FOU]
      temp(2,:) =[(kesi-ONE)/FOU,(-kesi-ONE)/FOU,
     &            (kesi+ONE)/FOU,(-kesi+ONE)/FOU]

      dNdkesi = transpose(temp)

      J(1:2,1:2) = MATMUL(Coor(1:2,1:4),dNdkesi(1:4,1:2))  

      call Matrix_Det_2x2(J,detJ)       

      return 
      end SUBROUTINE Cal_detJ                  
