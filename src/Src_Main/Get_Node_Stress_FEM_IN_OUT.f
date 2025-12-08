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
 
      SUBROUTINE Get_Node_Stress_FEM_IN_OUT(Yes_Add_Insitu,isub,c_DISP,
     &        Stress_xx_N,Stress_yy_N,Stress_xy_N,Stress_vm_N)
c     Computational Node Stress
C     Store into local variables: Stress_xx_N, Stress_yy_N, Stress_xy_N, Stress_vm_N

c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Dynamic
      use Global_Material
      use Global_Stress
      use Global_Disp
      
c     --------------------------
C     Variable Type Declaration
c     --------------------------
      implicit none
      integer, intent(in)::isub
      logical, intent(in)::Yes_Add_Insitu
      real(kind=FT), intent(in)::c_DISP(num_Node*2)
      real(kind=FT), intent(out)::Stress_xx_N(num_Node),
     &                               Stress_yy_N(num_Node),
     &                               Stress_xy_N(num_Node),
     &                               Stress_vm_N(num_Node)
      integer i_E
      real(kind=FT) c_thick,c_D(3,3),U(8)
      real(kind=FT) c_v   
      real(kind=FT) c_X_NODES(4),c_Y_NODES(4),T_Kesi(4),T_Yita(4),
     &                 c_kesi,c_yita,c_Stress(3)
      integer c_NN(4),i_N,C_Node,c_Count(num_Node) 
      real(kind=FT) c_v_all(Num_Node)
            
      print *,'    Calculating stress of all nodes...'
      
c     -----------------
C     Initialized to 0
c     -----------------
      Stress_xx_N(1:num_Node) = ZR
      Stress_yy_N(1:num_Node) = ZR
      Stress_xy_N(1:num_Node) = ZR
      Stress_vm_N(1:num_Node) = ZR
      c_Count(1:num_Node)     = 0
      
c     -----------------
C     Inter-unit cycle
c     -----------------
      T_Kesi = [-ONE,  ONE,  ONE, -ONE]
      T_Yita = [-ONE, -ONE,  ONE,  ONE]
      do i_E = 1,Num_Elem
          c_thick = thick(Elem_Mat(i_E))
          c_D     = D(Elem_Mat(i_E),:,:)
          
          ! Elastic modulus of the Weibull distribution. 2024-06-25. NEWFTU2024062402.
          if(Flag_Weibull_E)then
              if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
                  c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
              endif
          endif
          
          c_v     = Material_Para(Elem_Mat(i_E),2)
          c_NN    = G_NN(:,i_E)
          c_X_NODES = G_X_NODES(:,i_E)
          c_Y_NODES = G_Y_NODES(:,i_E)  
          U = [c_DISP(c_NN(1)*2-1),c_DISP(c_NN(1)*2),
     &         c_DISP(c_NN(2)*2-1),c_DISP(c_NN(2)*2),
     &         c_DISP(c_NN(3)*2-1),c_DISP(c_NN(3)*2),
     &         c_DISP(c_NN(4)*2-1),c_DISP(c_NN(4)*2)]
            ! Node Loop
          do i_N = 1,4
                C_Node = c_NN(i_N)
              c_v_all(C_Node) = c_v
            c_kesi = T_Kesi(i_N)                                               
            c_yita = T_Yita(i_N)
            call Cal_Ele_Stress_N4(i_E,i_N,c_X_NODES,c_Y_NODES,
     &                                 c_D,c_kesi,c_yita,U,
     &                                 c_Stress)      
            c_Count(C_Node) = c_Count(C_Node) + 1
              Stress_xx_N(C_Node) = Stress_xx_N(C_Node) + 
     &                                 c_Stress(1)
            Stress_yy_N(C_Node) = Stress_yy_N(C_Node) + 
     &                                 c_Stress(2)
            Stress_xy_N(C_Node) = Stress_xy_N(C_Node) + 
     &                                 c_Stress(3)    

     
              ! Add initial stress field
              if(Key_InSitu_Strategy==2 .and. Yes_Add_Insitu)then
                  Stress_xx_Node(C_Node) = Stress_xx_Node(C_Node) +
     &                                     Str_xx_InSitu(C_Node)
                  Stress_yy_Node(C_Node) = Stress_yy_Node(C_Node) +
     &                                     Str_yy_InSitu(C_Node)
                  Stress_xy_Node(C_Node) = Stress_xy_Node(C_Node) +
     &                                     Str_xy_InSitu(C_Node)
              endif
            end do
          
          
          
      end do  
c     --------------------
C     Node Stress Average
c     --------------------
      do i_N=1,num_Node   
          Stress_xx_N(i_N) = Stress_xx_N(i_N)/c_Count(i_N)
            Stress_yy_N(i_N) = Stress_yy_N(i_N)/c_Count(i_N)
            Stress_xy_N(i_N) = Stress_xy_N(i_N)/c_Count(i_N)
          call Tool_von_Mises(Stress_xx_N(i_N),
     &                            Stress_yy_N(i_N),
     &                            Stress_xy_N(i_N),
     &                            c_v_all(i_N),Stress_vm_N(i_N))
      end do
     
      RETURN
      END SUBROUTINE Get_Node_Stress_FEM_IN_OUT
