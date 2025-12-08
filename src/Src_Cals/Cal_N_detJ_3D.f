 
      subroutine Cal_N_detJ_3D(kesi,yita,zeta,
     &                         X_NODES,Y_NODES,Z_NODES,
     &                         detJ,N)
     
      use Global_Float_Type      
      implicit none
      
      real(kind=FT),intent(in)::kesi,yita,zeta,
     &                         X_NODES(8),Y_NODES(8),Z_NODES(8)
      real(kind=FT),intent(out)::detJ
      real(kind=FT),intent(out):: N(3,24)
      real(kind=FT) dNdkesi(8,3),J(3,3) 
      real(kind=FT) N1,N2,N3,N4,N5,N6,N7,N8
      real(kind=FT) o
      
      real(kind=FT) temp(3,8),Coor(3,8)
      
      o      = ZR
      dNdkesi(1:8,1:3) = ZR
      Coor(1,:)=[X_NODES(1),X_NODES(2),X_NODES(3),X_NODES(4),
     &           X_NODES(5),X_NODES(6),X_NODES(7),X_NODES(8)]
      Coor(2,:)=[Y_NODES(1),Y_NODES(2),Y_NODES(3),Y_NODES(4),
     &           Y_NODES(5),Y_NODES(6),Y_NODES(7),Y_NODES(8)]
      Coor(3,:)=[Z_NODES(1),Z_NODES(2),Z_NODES(3),Z_NODES(4),
     &           Z_NODES(5),Z_NODES(6),Z_NODES(7),Z_NODES(8)]
      
        N1 = (ONE-kesi)*(ONE-yita)*(ONE-zeta)/EIG
        N2 = (ONE+kesi)*(ONE-yita)*(ONE-zeta)/EIG
        N3 = (ONE+kesi)*(ONE+yita)*(ONE-zeta)/EIG
        N4 = (ONE-kesi)*(ONE+yita)*(ONE-zeta)/EIG
      N5 = (ONE-kesi)*(ONE-yita)*(ONE+zeta)/EIG
      N6 = (ONE+kesi)*(ONE-yita)*(ONE+zeta)/EIG
      N7 = (ONE+kesi)*(ONE+yita)*(ONE+zeta)/EIG
      N8 = (ONE-kesi)*(ONE+yita)*(ONE+zeta)/EIG
      
        N(1,1:24) =
     &        [N1,o,o,N2,o,o,N3,o,o,N4,o,o,N5,o,o,N6,o,o,N7,o,o,N8,o,o]
        N(2,1:24) =
     &        [o,N1,o,o,N2,o,o,N3,o,o,N4,o,o,N5,o,o,N6,o,o,N7,o,o,N8,o]
      N(3,1:24) =
     &        [o,o,N1,o,o,N2,o,o,N3,o,o,N4,o,o,N5,o,o,N6,o,o,N7,o,o,N8]

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
      end SUBROUTINE Cal_N_detJ_3D            
