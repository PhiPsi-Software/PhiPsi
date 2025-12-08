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
 
      SUBROUTINE Get_Gauss_Stress_FEM(isub)
c     Computational Node Stress
C     Store in global variables: Stress_xx_Gauss, Stress_yy_Gauss, Stress_xy_Gauss, Stress_vm_Gauss

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
      integer, intent(in)::isub
      !include 'omp_lib.h'
      integer i_E
      real(kind=FT) c_D(3,3),U(8)
      real(kind=FT) c_X_NODES(4),c_Y_NODES(4),
     &                 c_kesi,c_yita,c_Stress(3)
      integer c_NN(4)
      real(kind=FT) kesi_N_Enr(Num_Gauss_P_FEM),
     &              yita_N_Enr(Num_Gauss_P_FEM),
     &              weight_N_Enr(Num_Gauss_P_FEM)    
      integer G_Counter,i_G
      real(kind=FT) kesi(Num_Gauss_P_FEM),yita(Num_Gauss_P_FEM)
      real(kind=FT) c_v
      integer num_Ele_Killed,Yes_Killed_Ele
      real(kind=FT) c_T_Alpha,c_TStress(3)
      real(kind=FT) c_S_1,c_S_3,c_theta_stress,c_Ele_Gauss(3)
      
      print *,'    Calculating stress of all Gauss points...'
      
      call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,
     &                           kesi_N_Enr,yita_N_Enr,weight_N_Enr)
     
      G_Counter = 0
      ! Life-and-Death element preparation
      if(Key_EKILL==1)then
          num_Ele_Killed = count(Ele_Killed_Each_Load_Step(1:isub,:)>0)
      endif  
      
c     -----------------
C     Inter-unit cycle
c     -----------------
      !.............................
      ! OpenMP multi-core computing
      !.............................
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,i_G,c_Stress,
!$OMP&              c_D,c_NN,c_X_NODES,c_Y_NODES ,U,
!$OMP&              c_kesi,c_yita,G_Counter) 
      do i_E = 1,Num_Elem
          c_D     = D(Elem_Mat(i_E),:,:)   
          
          ! Elastic modulus of the Weibull distribution. 2024-06-25. NEWFTU2024062402.
          if(Flag_Weibull_E)then
              if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
                  c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
              endif
          endif
          
          c_v     = Material_Para(Elem_Mat(i_E),2)
          c_T_Alpha  = T_Alpha(Elem_Mat(i_E))  
          c_NN    = G_NN(:,i_E)
          c_X_NODES = G_X_NODES(:,i_E)
          c_Y_NODES = G_Y_NODES(:,i_E)          
          U = [DISP(c_NN(1)*2-1),DISP(c_NN(1)*2),
     &         DISP(c_NN(2)*2-1),DISP(c_NN(2)*2),
     &         DISP(c_NN(3)*2-1),DISP(c_NN(3)*2),
     &         DISP(c_NN(4)*2-1),DISP(c_NN(4)*2)]
     
          c_Ele_Gauss(1:3) = ZR
          
          ! Life and Death element Processing: Gaussian Stress is 0
          Yes_Killed_Ele = 0
          if(Key_EKILL==1 .and. num_Ele_Killed>=1) then
            if(any(Ele_Killed_Each_Load_Step(1:isub,1:num_Ele_Killed)
     &                                                      ==i_E))then
                Yes_Killed_Ele = 1
            endif
          endif 
          kesi(1:Num_Gauss_P_FEM)   = kesi_N_Enr
          yita(1:Num_Gauss_P_FEM)   = yita_N_Enr
          ! Loop over each Gauss point.
          do i_G = 1,Num_Gauss_P_FEM
              !G_Counter =  G_Counter +1 
              G_Counter =  (i_E-1)*Num_Gauss_P_FEM +i_G 
            c_kesi = kesi(i_G)                                          
            c_yita = yita(i_G)
            call Cal_Ele_Stress_N4(i_E,i_G,c_X_NODES,c_Y_NODES,
     &                                 c_D,c_kesi,c_yita,U,
     &                                 c_Stress)   
              Stress_xx_Gauss(G_Counter) = c_Stress(1)
            Stress_yy_Gauss(G_Counter) = c_Stress(2)
            Stress_xy_Gauss(G_Counter) = c_Stress(3)
              ! Accumulation of stresses at each Gauss point of the current element (used to calculate the average
              ! Gauss stress of each element)
              c_Ele_Gauss(1) = c_Ele_Gauss(1)+c_Stress(1)
              c_Ele_Gauss(2) = c_Ele_Gauss(2)+c_Stress(1)
              c_Ele_Gauss(3) = c_Ele_Gauss(3)+c_Stress(1)
              
              ! Thermal stress, Theory: Equation 15.1.98 from the Finite Element Method Basic Tutorial (5th
              ! Edition), 2019-09-24
              if(Key_Thermal_Stress==1)then
                  c_TStress(1:3) =ZR
                  if(Key_Type_2D==1)then
                    c_TStress = c_T_Alpha*
     &                         Thermal_Str_Temper(Elem_Mat(i_E))*
     &                         MATMUL(c_D,[ONE,ONE,ZR])
                  elseif(Key_Type_2D==2)then
                    c_TStress = (ONE+c_v)*c_T_Alpha*
     &                         Thermal_Str_Temper(Elem_Mat(i_E))*
     &                         MATMUL(c_D,[ONE,ONE,ZR])
                  endif
                  Stress_xx_Gauss(G_Counter) = c_Stress(1)-c_TStress(1)
                Stress_yy_Gauss(G_Counter) = c_Stress(2)-c_TStress(2)
                Stress_xy_Gauss(G_Counter) = c_Stress(3)-c_TStress(3)     
                  ! Accumulation of stresses at each Gauss point of the current element (used to calculate the average
                  ! Gauss stress of each element)
                  c_Ele_Gauss(1) = c_Ele_Gauss(1)-c_TStress(1)
                  c_Ele_Gauss(2) = c_Ele_Gauss(2)-c_TStress(2)
                  c_Ele_Gauss(3) = c_Ele_Gauss(3)-c_TStress(3)                
              endif

              
              ! Add initial stress field
              if(Key_InSitu_Strategy==2)then
                  Stress_xx_Gauss(G_Counter)=Stress_xx_Gauss(G_Counter)
     &                   + InSitu_Strs_Gaus_xx(i_E,i_G)
                  Stress_yy_Gauss(G_Counter)=Stress_yy_Gauss(G_Counter)
     &                   + InSitu_Strs_Gaus_yy(i_E,i_G) 
                  Stress_xy_Gauss(G_Counter)=Stress_xy_Gauss(G_Counter)
     &                   + InSitu_Strs_Gaus_xy(i_E,i_G) 
                  ! Accumulation of stresses at each Gauss point of the current element (used to calculate the average
                  ! Gauss stress of each element)
                  c_Ele_Gauss(1) = c_Ele_Gauss(1) + 
     &                                      InSitu_Strs_Gaus_xx(i_E,i_G)
                  c_Ele_Gauss(2) = c_Ele_Gauss(2) + 
     &                                      InSitu_Strs_Gaus_yy(i_E,i_G) 
                  c_Ele_Gauss(3) = c_Ele_Gauss(3) + 
     &                                      InSitu_Strs_Gaus_xy(i_E,i_G)       
              endif
              ! If the element is killed or destroyed, the Gauss stress is 0.
              if(Yes_Killed_Ele==1  .or. Elem_Break(i_E)) then
                  Stress_xx_Gauss(G_Counter) = ZR
                  Stress_yy_Gauss(G_Counter) = ZR
                  Stress_xy_Gauss(G_Counter) = ZR
                  c_Ele_Gauss(1) = ZR
                  c_Ele_Gauss(2) = ZR
                  c_Ele_Gauss(3) = ZR
              endif              
              call Tool_von_Mises(Stress_xx_Gauss(G_Counter),
     &                                Stress_yy_Gauss(G_Counter),
     &                                Stress_xy_Gauss(G_Counter),
     &                                c_v,Stress_vm_Gauss(G_Counter)) 

          end do  
          ! Average Gauss stress of each element
          Elem_Ave_Gauss_Stress(i_E,1:3) = c_Ele_Gauss(1:3)/
     &                                     Num_Gauss_P_FEM
          ! Calculate the average principal stress for each element
          call Tool_Principal_Stresses_2D(Elem_Ave_Gauss_Stress(i_E,1),
     &                                    Elem_Ave_Gauss_Stress(i_E,2),
     &                                    Elem_Ave_Gauss_Stress(i_E,3),
     &                                    c_S_1,c_S_3,c_theta_stress)
          Elem_Ave_Gauss_S1(i_E) =  c_S_1          
      end do  
!$omp end parallel do        

      RETURN
      END SUBROUTINE Get_Gauss_Stress_FEM
