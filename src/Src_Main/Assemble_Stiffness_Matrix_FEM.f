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
 
      SUBROUTINE Assemble_Stiffness_Matrix_FEM(isub,globalK,
     &                              T_Freedom,Total_Num_G_P)
c     Assemble the stiffness matrix.
      use Global_Float_Type
      use Global_Model
      use Global_Filename
      use Global_Common
      use Global_Material
      use omp_lib
      
      implicit none
      
      integer,intent(in)::isub,T_Freedom
      real(kind=FT),intent(out):: globalK(T_Freedom,T_Freedom)
      integer,intent(out)::Total_Num_G_P
      integer i_E
      real(kind=FT) c_thick,c_D(3,3)
      real(kind=FT) c_X_NODES(4),c_Y_NODES(4)
      integer c_NN(4) 
      real(kind=FT) kesi(Num_Gauss_P_FEM),yita(Num_Gauss_P_FEM),
     &              weight(Num_Gauss_P_FEM)             
      real(kind=FT) localK(8,8)
      integer local(8),i_row,i_col,nIndex
      real(kind=FT) penalty_k
      integer i_CP_node,c_CP_node,first_node
      integer first_node_DOF,c_node_DOF
      integer i_Boux_nonzero,i_Bouy_nonzero,c_node
      integer c_i
      integer num_Ele_Killed
      
      print *,'    Assembling the global stiffness matrix......'  
      Total_Num_G_P = 0
      globalK(1:T_Freedom,1:T_Freedom) = ZR
      
      call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,kesi,yita,weight)
      nIndex = 0
      ! Life-and-Death element preparation
      if(Key_EKILL==1)then
          num_Ele_Killed = count(Ele_Killed_Each_Load_Step(1:isub,:)>0)
      endif
      
      !..................
      ! Inter-unit cycle
      !..................
      do i_E = 1,Num_Elem
          ! Starting number of Gauss points for each element
          Ele_GP_Start_Num(i_E) = Total_Num_G_P + 1
          Total_Num_G_P = Total_Num_G_P +Num_Gauss_P_FEM
      enddo

      !.............................
      ! OpenMP multi-core computing
      !.............................
!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_E,i_row,i_col,localK,
!$OMP&              c_D,c_NN,c_X_NODES,c_Y_NODES,c_thick,local
!$OMP&              ) 
      do i_E = 1,Num_Elem
          ! Starting number of Gauss points for each element
          !Ele_GP_Start_Num(i_E) = Total_Num_G_P + 1
          !Total_Num_G_P = Total_Num_G_P + 4
          c_thick = thick(Elem_Mat(i_E))
          c_D     = D(Elem_Mat(i_E),1:3,1:3)   
          
          
          ! Elastic modulus of the Weibull distribution. 2024-06-25. NEWFTU2024062402.
          if(Flag_Weibull_E)then
              if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
                  c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
              endif
          endif
          
          ! Life-and-Death element processing: weakening D matrix (elastic modulus)
          if(Key_EKILL==1 .and. num_Ele_Killed>=1) then
            if(any(Ele_Killed_Each_Load_Step(1:isub,1:num_Ele_Killed)
     &                                                      ==i_E))then
                c_D= c_D*EKILL_Weaken_Factor
            endif
          endif  
          
          ! Damage unit processing: weakening D matrix (elastic modulus), 2020-02-28
          if(Key_Element_Break==1 .and. Elem_Break(i_E))then
              !c_D= c_D*EKILL_Weaken_Factor
              c_D= c_D*1.0D-3
          endif   
          
          c_NN    = G_NN(:,i_E)
          c_X_NODES = G_X_NODES(:,i_E)
          c_Y_NODES = G_Y_NODES(:,i_E)           
            !Traditional index locations
          local=[c_NN(1)*2-1,c_NN(1)*2,c_NN(2)*2-1,c_NN(2)*2,
     &           c_NN(3)*2-1,c_NN(3)*2,c_NN(4)*2-1,c_NN(4)*2]
          !Get the element stiffness matrix of the current element	
            call Cal_Ele_Stiffness_Matrix_N4(c_X_NODES,c_Y_NODES,
     &                           c_thick,c_D,kesi,yita,weight,
     &                           localK)
          do i_row = 1,8
              do i_col = 1,8
                  globalK(local(i_row),local(i_col)) = 
     &                   globalK(local(i_row),local(i_col)) +
     &                   localK(i_row,i_col)
              end do
          end do
      end do
!$omp end parallel do
      
      ! Save the number of Gauss points for each element to a global variable
      num_GP_Elem(1:num_Elem) = Num_Gauss_P_FEM
      
      ! Filter the calculation results
      !option 1
C     do i=1,T_Freedom
C         do j=1,T_Freedom
C             if (abs(globalK(i,j))< = 1.0D-8) then
C                  globalK(i,j) = ZR 
C             endif
C         enddo
C     end do

      ! Option 2, speed greatly increased
      where (abs(globalK)< = 1.0D-8)
          globalK = ZR 
      end where

      !.................................................................
      ! Handling of non-zero displacement boundary conditions:
      ! Finite Element Methods in Engineering, Zeng Pan, Fourth Edition
      !.................................................................
      if(Num_Boux_nonzero >0)then
          !penalty_k =  1.0D4*maxval(globalK)
          penalty_k =  penalty_k_bou_nonzero 
          do i_Boux_nonzero =1,Num_Boux_nonzero
              c_node = int(Bou_x_nonzero(i_Boux_nonzero,1))
              c_node_DOF = 2*c_node-1
              globalK(c_node_DOF,c_node_DOF) =
     &                globalK(c_node_DOF,c_node_DOF)+penalty_k              
          enddo
      endif
      if(Num_Bouy_nonzero >0)then
          !penalty_k =  1.0D4*maxval(globalK)
          penalty_k =  penalty_k_bou_nonzero 
          do i_Bouy_nonzero =1,Num_Bouy_nonzero
              c_node = int(Bou_y_nonzero(i_Bouy_nonzero,1))
              c_node_DOF = 2*c_node
              globalK(c_node_DOF,c_node_DOF) =
     &                globalK(c_node_DOF,c_node_DOF)+penalty_k              
          enddo
      endif
      
            
      !...................................................................................
      ! Penalty stiffness method for node coupling, theoretical basis:
      !http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch09.d/IFEM.Ch09.pdf
      !...................................................................................
      ! Used for option 1 - directly defining coupled nodes through a keyword file (in this method, only
      ! one set of coupling can be defined per direction)
      if (num_CP_x_nodes>0)then
          penalty_k =  1.0D4*maxval(globalK)
          first_node = CP_x_nodes(1)
          do i_CP_node =2,num_CP_x_nodes
              c_CP_node = CP_x_nodes(i_CP_node)
              first_node_DOF = 2*first_node-1
              c_node_DOF = 2*c_CP_node-1
              globalK(first_node_DOF,first_node_DOF) =
     &                globalK(first_node_DOF,first_node_DOF)+penalty_k
              globalK(c_node_DOF,c_node_DOF) =
     &                globalK(c_node_DOF,c_node_DOF)+penalty_k
              globalK(first_node_DOF,c_node_DOF) =
     &                globalK(first_node_DOF,c_node_DOF)-penalty_k
              globalK(c_node_DOF,first_node_DOF) =
     &                globalK(c_node_DOF,first_node_DOF)-penalty_k
          enddo
      endif
      if (num_CP_y_nodes > 0)then
          penalty_k =  1.0D4*maxval(globalK)
          first_node = CP_y_nodes(1)
          do i_CP_node =2,num_CP_y_nodes
              c_CP_node = CP_y_nodes(i_CP_node)
              first_node_DOF = 2*first_node
              c_node_DOF = 2*c_CP_node
              globalK(first_node_DOF,first_node_DOF) =
     &                globalK(first_node_DOF,first_node_DOF)+penalty_k
              globalK(c_node_DOF,c_node_DOF) =
     &                globalK(c_node_DOF,c_node_DOF)+penalty_k
              globalK(first_node_DOF,c_node_DOF) =
     &                globalK(first_node_DOF,c_node_DOF)-penalty_k
              globalK(c_node_DOF,first_node_DOF) =
     &                globalK(c_node_DOF,first_node_DOF)-penalty_k
          enddo
      endif
      ! Used for option 2 – defining coupled nodes through dofx, dofy, dofz files (this method can only
      ! define multiple sets of coupling for each direction)
      if (num_CP_set_x>0)then
          penalty_k =  1.0D4*maxval(globalK)
          do c_i = 1,num_CP_set_x
              first_node = CP_nodes_x(c_i,1)
              do i_CP_node =2, num_nodes_CP_set_x(c_i)
                  c_CP_node = CP_nodes_x(c_i,i_CP_node)
                  first_node_DOF = 2*first_node-1
                  c_node_DOF = 2*c_CP_node-1
                  globalK(first_node_DOF,first_node_DOF) =
     &                globalK(first_node_DOF,first_node_DOF)+penalty_k
                  globalK(c_node_DOF,c_node_DOF) =
     &                globalK(c_node_DOF,c_node_DOF)+penalty_k
                  globalK(first_node_DOF,c_node_DOF) =
     &                globalK(first_node_DOF,c_node_DOF)-penalty_k
                  globalK(c_node_DOF,first_node_DOF) =
     &                globalK(c_node_DOF,first_node_DOF)-penalty_k
              enddo
          enddo
      endif      
      if (num_CP_set_y>0)then
          penalty_k =  1.0D4*maxval(globalK)
          do c_i = 1,num_CP_set_y
              first_node = CP_nodes_y(c_i,1)
              do i_CP_node =2, num_nodes_CP_set_y(c_i)
                  c_CP_node = CP_nodes_y(c_i,i_CP_node)
                  first_node_DOF = 2*first_node
                  c_node_DOF = 2*c_CP_node
                  globalK(first_node_DOF,first_node_DOF) =
     &                globalK(first_node_DOF,first_node_DOF)+penalty_k
                  globalK(c_node_DOF,c_node_DOF) =
     &                globalK(c_node_DOF,c_node_DOF)+penalty_k
                  globalK(first_node_DOF,c_node_DOF) =
     &                globalK(first_node_DOF,c_node_DOF)-penalty_k
                  globalK(c_node_DOF,first_node_DOF) =
     &                globalK(c_node_DOF,first_node_DOF)-penalty_k
              enddo
          enddo
      endif        
      
      RETURN
      END SUBROUTINE Assemble_Stiffness_Matrix_FEM
