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
 
      SUBROUTINE Get_Node_Stress_FEM(isub)
c     Computational Node Stress
C     Store into global variables: Stress_xx_Node, Stress_yy_Node, Stress_xy_Node, Stress_vm_Node

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
      use omp_lib
      
c     --------------------------
C     Variable Type Declaration
c     --------------------------
      implicit none
      !include 'omp_lib.h'
      integer, intent(in)::isub
      integer i_E
      real(kind=FT) c_thick,c_D(3,3),U(8)
      real(kind=FT) c_v   
      real(kind=FT) c_X_NODES(4),c_Y_NODES(4),T_Kesi(4),T_Yita(4),
     &                 c_kesi,c_yita,c_Stress(3)
      integer c_NN(4),i_N,C_Node,c_Count(num_Node) 
      real(kind=FT) c_v_all(Num_Node)
      integer num_Ele_Killed,Yes_Killed_Ele      
      real(kind=FT) c_T_Alpha,c_TStress(3)
      
      print *,'    Calculating stress of all nodes...'
c     -----------------
C     Initialized to 0
c     -----------------
      Stress_xx_Node(1:num_Node) = ZR
      Stress_yy_Node(1:num_Node) = ZR
      Stress_xy_Node(1:num_Node) = ZR
      Stress_vm_Node(1:num_Node) = ZR
      ! Thermal stress
      if(Key_Thermal_Stress==1)then
          TStress_xx_Node(1:num_Node) = ZR
          TStress_yy_Node(1:num_Node) = ZR
          TStress_xy_Node(1:num_Node) = ZR
          TStress_vm_Node(1:num_Node) = ZR
      endif
      c_Count(1:num_Node)        = 0
      ! Life-and-Death element preparation
      if(Key_EKILL==1)then
          num_Ele_Killed = count(Ele_Killed_Each_Load_Step(1:isub,:)>0)
      endif  
      
c     -----------------
C     Inter-unit cycle
c     -----------------
      T_Kesi = [-ONE,  ONE,  ONE, -ONE]
      T_Yita = [-ONE, -ONE,  ONE,  ONE]
      !.............................
      ! OpenMP multi-core computing
      !.............................
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,i_N,c_Stress,
!$OMP&              c_thick,c_D,c_NN,c_X_NODES,c_Y_NODES ,U,
!$OMP&              C_Node,c_kesi,c_yita) 
      do i_E = 1,Num_Elem
          !thread_id = omp_get_thread_num()
          c_thick = thick(Elem_Mat(i_E))
          c_D     = D(Elem_Mat(i_E),:,:)   
          
          ! Elastic modulus of the Weibull distribution. 2024-06-25. NEWFTU2024062402.
          if(Flag_Weibull_E)then
              if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
                  c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
              endif
          endif
          
          c_v     = Material_Para(Elem_Mat(i_E),2)
          c_T_Alpha = T_Alpha(Elem_Mat(i_E))
          c_NN    = G_NN(:,i_E)
          c_X_NODES = G_X_NODES(:,i_E)
          c_Y_NODES = G_Y_NODES(:,i_E)             
          U = [DISP(c_NN(1)*2-1),DISP(c_NN(1)*2),
     &         DISP(c_NN(2)*2-1),DISP(c_NN(2)*2),
     &         DISP(c_NN(3)*2-1),DISP(c_NN(3)*2),
     &         DISP(c_NN(4)*2-1),DISP(c_NN(4)*2)]
          ! Life and Death element Processing: Gaussian Stress is 0
          Yes_Killed_Ele = 0
          if(Key_EKILL==1 .and. num_Ele_Killed>=1) then
            if(any(Ele_Killed_Each_Load_Step(1:isub,1:num_Ele_Killed)
     &                                                      ==i_E))then
                Yes_Killed_Ele = 1
            endif
          endif      
          
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
            Stress_xx_Node(C_Node) = Stress_xx_Node(C_Node) + 
     &                                 c_Stress(1)
            Stress_yy_Node(C_Node) = Stress_yy_Node(C_Node) + 
     &                                 c_Stress(2)
            Stress_xy_Node(C_Node) = Stress_xy_Node(C_Node) + 
     &                                 c_Stress(3)  
     
     
              ! Dealing with thermal stress
              if(Key_Thermal_Stress==1)then
                  c_TStress(1:3) =ZR
                  if (Key_Type_2D==1) then
                    c_TStress = c_T_Alpha*
     &                         Thermal_Str_Temper(Elem_Mat(i_E))*
     &                         MATMUL(c_D,[ONE,ONE,ZR])
                  elseif (Key_Type_2D==2) then
                    c_TStress = (ONE+c_v)*c_T_Alpha*
     &                         Thermal_Str_Temper(Elem_Mat(i_E))*
     &                         MATMUL(c_D,[ONE,ONE,ZR])
                  endif
                  TStress_xx_Node(C_Node) = TStress_xx_Node(C_Node) + 
     &                                      c_TStress(1)
                  TStress_yy_Node(C_Node) = TStress_yy_Node(C_Node) + 
     &                                      c_TStress(2)
                  TStress_xy_Node(C_Node) = TStress_xy_Node(C_Node) + 
     &                                      c_TStress(3)   
                  Stress_xx_Node(C_Node) = Stress_xx_Node(C_Node) -
     &                                     c_TStress(1)
                  Stress_yy_Node(C_Node) = Stress_yy_Node(C_Node) -
     &                                     c_TStress(2)
                  Stress_xy_Node(C_Node) = Stress_xy_Node(C_Node) -
     &                                     c_TStress(3)        
              endif
     
              ! Add the initial stress field
              if(Key_InSitu_Strategy==2)then
                  Stress_xx_Node(C_Node) = Stress_xx_Node(C_Node) +
     &                                     Str_xx_InSitu(C_Node)
                  Stress_yy_Node(C_Node) = Stress_yy_Node(C_Node) +
     &                                     Str_yy_InSitu(C_Node)
                  Stress_xy_Node(C_Node) = Stress_xy_Node(C_Node) +
     &                                     Str_xy_InSitu(C_Node)
              endif
              ! If the element is killed, the Gaussian stress is 0.
              if(Yes_Killed_Ele==1) then
                  Stress_xx_Node(C_Node) = ZR
                  Stress_yy_Node(C_Node) = ZR
                  Stress_xy_Node(C_Node) = ZR
              endif               
            end do
      end do  
!$omp end parallel do   

      ! Node Stress Average
!$omp parallel do DEFAULT(SHARED) PRIVATE(i_N) 
      do i_N=1,num_Node   
           Stress_xx_Node(i_N) = Stress_xx_Node(i_N)/c_Count(i_N)
           Stress_yy_Node(i_N) = Stress_yy_Node(i_N)/c_Count(i_N)
           Stress_xy_Node(i_N) = Stress_xy_Node(i_N)/c_Count(i_N)
           call Tool_von_Mises(Stress_xx_Node(i_N),
     &                        Stress_yy_Node(i_N),
     &                        Stress_xy_Node(i_N),
     &                        c_v_all(i_N),Stress_vm_Node(i_N))
           if(Key_Thermal_Stress==1)then
              TStress_xx_Node(i_N) = TStress_xx_Node(i_N)/c_Count(i_N)
              TStress_yy_Node(i_N) = TStress_yy_Node(i_N)/c_Count(i_N)
              TStress_xy_Node(i_N) = TStress_xy_Node(i_N)/c_Count(i_N)
              call Tool_von_Mises(TStress_xx_Node(i_N),
     &                              TStress_yy_Node(i_N),
     &                              TStress_xy_Node(i_N),
     &                              c_v_all(i_N),TStress_vm_Node(i_N))
          endif
      end do
!$omp end parallel do
  
      RETURN
      END SUBROUTINE Get_Node_Stress_FEM
