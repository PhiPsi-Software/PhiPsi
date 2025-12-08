 
      subroutine Cal_detJ_3D(kesi,yita,zeta,
     &                    X_NODES,Y_NODES,Z_NODES,
     &                    detJ)
     

      use Global_Float_Type
      implicit none
      
      real(kind=FT),intent(in)::kesi,yita,zeta
      real(kind=FT),intent(in)::X_NODES(8),Y_NODES(8),Z_NODES(8)
      real(kind=FT),intent(out)::detJ       
      real(kind=FT) temp(3,8),Coor(3,8),J(3,3),dNdkesi(8,3)
      
      Coor(1,:) = [X_NODES(1),X_NODES(2),X_NODES(3),X_NODES(4),
     &             X_NODES(5),X_NODES(6),X_NODES(7),X_NODES(8)]
      Coor(2,:) = [Y_NODES(1),Y_NODES(2),Y_NODES(3),Y_NODES(4),
     &             Y_NODES(5),Y_NODES(6),Y_NODES(7),Y_NODES(8)]
      Coor(3,:) = [Z_NODES(1),Z_NODES(2),Z_NODES(3),Z_NODES(4),
     &             Z_NODES(5),Z_NODES(6),Z_NODES(7),Z_NODES(8)]
      
      temp(1,1:8)=[ -(ONE - yita)*(ONE - zeta)/EIG,
     &               (ONE - yita)*(ONE - zeta)/EIG,
     &               (ONE + yita)*(ONE - zeta)/EIG,
     &              -(ONE + yita)*(ONE - zeta)/EIG,
     &              -(ONE - yita)*(ONE + zeta)/EIG,
     &               (ONE - yita)*(ONE + zeta)/EIG,
     &               (ONE + yita)*(ONE + zeta)/EIG,
     &              -(ONE + yita)*(ONE + zeta)/EIG]  
     
      temp(2,1:8)=[ -(ONE - kesi)*(ONE - zeta)/EIG,
     &              -(ONE + kesi)*(ONE - zeta)/EIG,
     &               (ONE + kesi)*(ONE - zeta)/EIG,
     &               (ONE - kesi)*(ONE - zeta)/EIG,
     &              -(ONE - kesi)*(ONE + zeta)/EIG,
     &              -(ONE + kesi)*(ONE + zeta)/EIG,
     &               (ONE + kesi)*(ONE + zeta)/EIG,
     &               (ONE - kesi)*(ONE + zeta)/EIG] 
     
      temp(3,1:8)=[ -(ONE - kesi)*(ONE - yita)/EIG,
     &              -(ONE + kesi)*(ONE - yita)/EIG,
     &              -(ONE + kesi)*(ONE + yita)/EIG,
     &              -(ONE - kesi)*(ONE + yita)/EIG,
     &               (ONE - kesi)*(ONE - yita)/EIG,
     &               (ONE + kesi)*(ONE - yita)/EIG,
     &               (ONE + kesi)*(ONE + yita)/EIG,
     &               (ONE - kesi)*(ONE + yita)/EIG]   
              
      dNdkesi = transpose(temp)

      J = MATMUL(Coor,dNdkesi)     

      call Matrix_Det_3x3(J,detJ) 
      
      
      return 
      end SUBROUTINE Cal_detJ_3D                  
