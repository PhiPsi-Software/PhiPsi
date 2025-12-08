!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !
 
      SUBROUTINE D3_Get_InSitu_Node_Reactions
      ! Obtain the constraint reactions for all nodes and store them in Reactions_Nodes_3D(num_Node, 3).
      ! All nodes of the model need to be constrained.
      ! 2022-06-03.
      ! Abandon this plan.


c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Stress
      use Global_Material
      
c     --------------------------
C     Variable Type Declaration
c     --------------------------
      implicit none
      integer i_E
      real(kind=FT) c_D(6,6),U(24)
      real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8),
     &                 c_kesi,c_yita,c_zeta,c_Stress(6)
      integer c_NN(8)
      integer G_Counter,i_G
      real(kind=FT) kesi(Num_Gauss_P_FEM_3D),
     &              yita(Num_Gauss_P_FEM_3D),
     &              zeta(Num_Gauss_P_FEM_3D),
     &              weight(Num_Gauss_P_FEM_3D)
      real(kind=FT) c_temp(24),detJ
      integer local(24)      
      real(kind=FT) c_B(6,24)
      real(kind=FT) J(3,3),dNdkesi(8,3),N(3,24),Inverse_J(3,3),  
     &              dNdx(8,3)   
      real(kind=FT) c_Node_Forces(num_Node*3)
      c_Node_Forces(1:num_Node*3) = ZR
      
      
      
      call Cal_Gauss_Points_3D_8nodes(Num_Gauss_P_FEM_3D,
     &                                kesi,yita,zeta,weight)          
      G_Counter = 0
        
          
      do i_E = 1,Num_Elem   
          c_D     = D(Elem_Mat(i_E),1:6,1:6)     
          c_NN    = G_NN(1:8,i_E)
          c_X_NODES = G_X_NODES(:,i_E)
          c_Y_NODES = G_Y_NODES(:,i_E)    
          c_Z_NODES = G_Z_NODES(:,i_E) 
          !Traditional index locations
          local=[c_NN(1)*3-2,c_NN(1)*3-1,c_NN(1)*3,
     &           c_NN(2)*3-2,c_NN(2)*3-1,c_NN(2)*3,
     &           c_NN(3)*3-2,c_NN(3)*3-1,c_NN(3)*3,
     &           c_NN(4)*3-2,c_NN(4)*3-1,c_NN(4)*3,
     &           c_NN(5)*3-2,c_NN(5)*3-1,c_NN(5)*3,
     &           c_NN(6)*3-2,c_NN(6)*3-1,c_NN(6)*3,
     &           c_NN(7)*3-2,c_NN(7)*3-1,c_NN(7)*3,
     &           c_NN(8)*3-2,c_NN(8)*3-1,c_NN(8)*3]
          ! Gauss Point Loop
          do i_G = 1,Num_Gauss_P_FEM_3D 
              G_Counter = G_Counter +1
              c_kesi = kesi(i_G)
              c_yita = yita(i_G)
              c_zeta = zeta(i_G)
              ! Calculates N, dNdkesi, J and the determinant of Jacobian matrix.   
              call Cal_N_dNdkesi_J_detJ_3D(c_kesi,c_yita,c_zeta,
     &                              c_X_NODES,c_Y_NODES,c_Z_NODES,
     &                              detJ,J,N,dNdkesi)    
              !The inverse matrix of J.
              call Matrix_Inverse_3x3(J,Inverse_J)
              !Calculate dNdx.
              dNdx = MATMUL(dNdkesi,Inverse_J)
              !B matrix for conventional element. 
              c_B(1:6,1:24)   =  ZR
              c_B(1,1:24:3)   =  dNdx(1:8,1)
              c_B(2,2:24:3)   =  dNdx(1:8,2)
              c_B(3,3:24:3)   =  dNdx(1:8,3)
              c_B(4,1:24:3)   =  dNdx(1:8,2)
              c_B(4,2:24:3)   =  dNdx(1:8,1)
              c_B(5,2:24:3)   =  dNdx(1:8,3)
              c_B(5,3:24:3)   =  dNdx(1:8,2)
              c_B(6,1:24:3)   =  dNdx(1:8,3)
              c_B(6,3:24:3)   =  dNdx(1:8,1)    
              !Calculate stress.
              c_Stress = [InSitu_Strs_Gaus_xx(i_E,i_G),
     &                    InSitu_Strs_Gaus_yy(i_E,i_G),
     &                    InSitu_Strs_Gaus_zz(i_E,i_G),
     &                    InSitu_Strs_Gaus_xy(i_E,i_G),
     &                    InSitu_Strs_Gaus_yz(i_E,i_G),
     &                    InSitu_Strs_Gaus_xz(i_E,i_G)]
              !Calculate internal force.
              c_temp(1:24) =  MATMUL(transpose(c_B(1:6,1:24)),c_Stress)
              !Update F.
              c_Node_Forces(local) = c_Node_Forces(local) - 
     &                               c_temp(1:24)*weight(i_G)*detJ     
          end do
      enddo
      
      RETURN
      END SUBROUTINE D3_Get_InSitu_Node_Reactions
