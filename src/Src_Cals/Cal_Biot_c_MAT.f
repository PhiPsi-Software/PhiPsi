 
      subroutine Cal_Biot_c_MAT
      
      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Crack
      use Global_Material
      use Global_Field_Problem

      implicit none
      integer i_E
      real(kind=FT) c_thick,c_B(3,8)
      real(kind=FT) c_X_NODES(4),c_Y_NODES(4)
      integer c_NN(4) 
      real(kind=FT) kesi(4),yita(4),weight(4)  
      integer i_G
      real(kind=FT) c_kesi,c_yita
      integer local_U(8),local_P(4),G_Counter
      real(kind=FT) tem_PNPxy(8),c_temp(8,4),detJ,N_4(4)
      
      Biot_c_MAT(1:Num_Node*2,1:Num_Node) = ZR
      
      if (Key_Integral_Sol  == 2)then
          call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,
     &                               kesi,yita,weight)
      elseif (Key_Integral_Sol  == 3)then
          call Cal_Gauss_Points_QUAD_for_SUBQUAD(Num_Sub_Quads,
     &                           kesi,yita,weight)
          Num_Gauss_P_FEM = Num_Sub_Quads*4
      endif
      G_Counter = 0
      do i_E = 1,Num_Elem
          c_thick = thick(Elem_Mat(i_E))   
          c_NN    = G_NN(:,i_E)
          c_X_NODES = G_X_NODES(:,i_E)
          c_Y_NODES = G_Y_NODES(:,i_E)           
          local_U=[c_NN(1)*2-1,c_NN(1)*2,c_NN(2)*2-1,c_NN(2)*2,
     &             c_NN(3)*2-1,c_NN(3)*2,c_NN(4)*2-1,c_NN(4)*2]     
          local_P=[c_NN(1),c_NN(2),c_NN(3),c_NN(4)] 
          call Cal_detJ(c_kesi,c_yita,c_X_NODES,c_Y_NODES,detJ)  
          do i_G = 1,Num_Gauss_P_FEM  
              G_Counter = G_Counter +1
              c_kesi = kesi(i_G); c_yita = yita(i_G)
              call Cal_Ele_B_N4(i_E,c_X_NODES,c_Y_NODES,
     &                          c_kesi,c_yita,c_B) 
              tem_PNPxy(1) = c_B(1,1)
              tem_PNPxy(2) = c_B(2,2)
              tem_PNPxy(3) = c_B(1,3)
              tem_PNPxy(4) = c_B(2,4)
              tem_PNPxy(5) = c_B(1,5)
              tem_PNPxy(6) = c_B(2,6)
              tem_PNPxy(7) = c_B(1,7)
              tem_PNPxy(8) = c_B(2,8)
              N_4(1)  = (ONE-c_kesi)*(ONE-c_yita)/FOU
              N_4(2)  = (ONE+c_kesi)*(ONE-c_yita)/FOU
              N_4(3)  = (ONE+c_kesi)*(ONE+c_yita)/FOU
              N_4(4)  = (ONE-c_kesi)*(ONE+c_yita)/FOU
              call Vectors_Multi(tem_PNPxy,8,N_4,4,c_temp)   
              Biot_c_MAT(local_U,local_P) = 
     &                Biot_c_MAT(local_U,local_P)  + 
     &                c_temp(1:8,1:4)*c_thick*weight(i_G)*detJ
          end do
      enddo
          
      return 
      end SUBROUTINE Cal_Biot_c_MAT        
