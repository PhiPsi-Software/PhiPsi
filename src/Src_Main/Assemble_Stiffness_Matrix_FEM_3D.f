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
 
      SUBROUTINE Assemble_Stiffness_Matrix_FEM_3D(isub,c_globalK,
     &                                 c_Total_Freedom,Total_Num_G_P)
c     Assemble the stiffness matrix (3D problem).
      use Global_Float_Type
      use Global_Model
      use Global_Filename
      use Global_Common
      use Global_Material
      use omp_lib
      ! Program Interface Testing. 2022-09-06.
      use Global_Cal_Ele_Stiffness_Matrix_3D_8nodes
      
      implicit none
      !include 'omp_lib.h'
      integer,intent(in)::isub,c_Total_Freedom
      integer,intent(out)::Total_Num_G_P
      real(kind=FT) ,intent(out)::c_globalK(c_Total_Freedom,
     &                                      c_Total_Freedom)
      integer i_E
      integer R_num,C_num
      real(kind=FT) c_D(6,6)
      real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
      integer c_NN(8)  
      real(kind=FT) kesi(Num_Gauss_P_FEM_3D),yita(Num_Gauss_P_FEM_3D),
     &               zeta(Num_Gauss_P_FEM_3D),weight(Num_Gauss_P_FEM_3D)       
      real(kind=FT) localK(24,24)
      integer local(24),i_row,i_col,nIndex 
      integer mat_num
      real(kind=FT) Rot_c_D_Comp(6,6),c_D_Comp(6,6),Volume_Ratio
      real(kind=FT) T_Matrix(6,6),TT_Matrix(6,6)
      
      c_globalK(1:Total_FD,1:Total_FD) = ZR
      
      
      call Cal_Gauss_Points_3D_8nodes(Num_Gauss_P_FEM_3D,
     &                                kesi,yita,zeta,weight)
      nIndex = 0
      Total_Num_G_P = 0
      
      !.............................
      ! OpenMP multi-core computing
      !.............................
!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_E,i_row,i_col,localK,
!$OMP&              c_D,c_NN,c_X_NODES,c_Y_NODES,c_Z_NODES,local,
!$OMP&              Volume_Ratio,c_D_Comp,T_Matrix,TT_Matrix,
!$OMP&              Rot_c_D_Comp,mat_num)       
!$OMP&            SCHEDULE(static)     

      do i_E = 1,Num_Elem
          mat_num = Elem_Mat(i_E)
          c_D(1:6,1:6)     = D(Elem_Mat(i_E),1:6,1:6)    
          !If the mat is composite material.
          if (Material_Type(mat_num)==5)then
              Volume_Ratio = Material_Para_Added(mat_num,10)
              c_D_Comp = D_Comp(mat_num,1:6,1:6)
              T_Matrix = Ele_ComMat_RotMatrix(i_E,1:6,1:6)
              TT_Matrix= TRANSPOSE(T_Matrix)
              Rot_c_D_Comp = MATMUL(TT_Matrix,c_D_Comp)
              Rot_c_D_Comp = MATMUL(Rot_c_D_Comp,T_Matrix)
              c_D =(ONE-Volume_Ratio)*c_D + Volume_Ratio*Rot_c_D_Comp
          endif
          c_NN    = G_NN(1:8,i_E)
          c_X_NODES = G_X_NODES(1:8,i_E)
          c_Y_NODES = G_Y_NODES(1:8,i_E)    
          c_Z_NODES = G_Z_NODES(1:8,i_E) 
          !$omp critical      
          Total_Num_G_P = Total_Num_G_P + Num_Gauss_P_FEM_3D
          !$omp end critical  
          !Traditional index locations
          local=[c_NN(1)*3-2,c_NN(1)*3-1,c_NN(1)*3,
     &           c_NN(2)*3-2,c_NN(2)*3-1,c_NN(2)*3,
     &           c_NN(3)*3-2,c_NN(3)*3-1,c_NN(3)*3,
     &           c_NN(4)*3-2,c_NN(4)*3-1,c_NN(4)*3,
     &           c_NN(5)*3-2,c_NN(5)*3-1,c_NN(5)*3,
     &           c_NN(6)*3-2,c_NN(6)*3-1,c_NN(6)*3,
     &           c_NN(7)*3-2,c_NN(7)*3-1,c_NN(7)*3,
     &           c_NN(8)*3-2,c_NN(8)*3-1,c_NN(8)*3]    
          ! Calculate the stiffness matrix of the current element
          call Cal_Ele_Stiffness_Matrix_3D_8nodes(i_E,
     &         Num_Gauss_P_FEM_3D,c_X_NODES,c_Y_NODES,c_Z_NODES,
     &                           c_D,kesi,yita,zeta,weight,
     &                           localK)  
!$omp critical          
          ! Global stiffness matrix of the assembly
          do i_row = 1,24
              do i_col = 1,24
                  c_globalK(local(i_row),local(i_col)) = 
     &                   c_globalK(local(i_row),local(i_col)) +
     &                   localK(i_row,i_col)
              end do
          end do  
!$omp end critical          
      end do
!$omp end parallel do 
      
    
      
      RETURN
      END SUBROUTINE Assemble_Stiffness_Matrix_FEM_3D
