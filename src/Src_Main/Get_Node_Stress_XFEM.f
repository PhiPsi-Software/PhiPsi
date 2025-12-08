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
 
      SUBROUTINE Get_Node_Stress_XFEM(isub,c_DISP)
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
      use Global_Crack
      use Global_Crack_Common
      use Global_Inclusion
      use Global_Cross
c     --------------------------
C     Variable Type Declaration
c     --------------------------
      implicit none
      integer, intent(in)::isub
      real(kind=FT),intent(in)::c_DISP(Total_FD)
      integer i_E,i_C
      real(kind=FT) c_thick,c_D(3,3),U(80)
      real(kind=FT) c_v   
      real(kind=FT) c_X_NODES(4),c_Y_NODES(4),T_Kesi(4),T_Yita(4),
     &                 c_kesi,c_yita,c_Stress(3),B(3,80),
     &                 tem_B(3,80)
      integer c_NN(4),i_N,C_Node,c_Count(num_Node),num_B,num_tem_B
      integer:: Location_ESM(MDOF_2D)
      integer Location_ESM_C_Cr_NoFEM(60)
      integer num_Loc_ESM_C_Cr_NoFEM
      integer num_Loc_ESM
      integer::Location_ESM_C_Crack(80)
      integer num_Loc_ESM_C_Crack,i_H
      integer i_Incl
      logical Yes_Gauss_in_Incl
      integer c_Incl_Num,c_MatNum
      real(kind=FT) c_T_Alpha,c_TStress(3)
      integer i_Cross
      integer num_Ele_Killed,Yes_Killed_Ele
      real(kind=FT) c_v_all(Num_Node)
      
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
      EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gauss_Points)= .false.
      do i_E = 1,Num_Elem
          c_thick = thick(Elem_Mat(i_E))
          c_D     = D(Elem_Mat(i_E),:,:)  
          
          ! Elastic modulus of the Weibull distribution. 2024-06-25. NEWFTU2024062402.
          if(Flag_Weibull_E)then
              if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
                  c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
              endif
          endif
          
          c_v     = v(Elem_Mat(i_E),1)
          c_T_Alpha  = T_Alpha(Elem_Mat(i_E)) 
          c_NN    = G_NN(:,i_E)
          c_X_NODES = G_X_NODES(:,i_E)
          c_Y_NODES = G_Y_NODES(:,i_E)           
          ! Life and Death element Processing: Gaussian Stress is 0
          Yes_Killed_Ele = 0
          if(Key_EKILL==1 .and. num_Ele_Killed>=1) then
            if(any(Ele_Killed_Each_Load_Step(1:isub,1:num_Ele_Killed)
     &                                                      ==i_E))then
                Yes_Killed_Ele = 1
            endif
          endif           
            !Decide the location of each element stiffness matrix in the global stiffness matrix.
          Location_ESM(1:MDOF_2D)  = 0
          num_Loc_ESM = 0       
          do i_C =1,num_Crack 
              call Location_Element_Stiff_Matrix(i_E,i_C,
     &                                        c_POS(:,i_C),
     &                                        Location_ESM_C_Crack,
     &                                        num_Loc_ESM_C_Crack,
     &                                        Location_ESM_C_Cr_NoFEM,
     &                                        num_Loc_ESM_C_Cr_NoFEM)
              Location_ESM(num_Loc_ESM+1:
     &                       num_Loc_ESM+num_Loc_ESM_C_Crack) = 
     &                       Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
              num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
          end do
          
          !Decide the location of each element stiffness matrix in the global stiffness matrix for hole.   
          if(num_Hole/=0)then
            do i_H =1,num_Hole 
              call Location_Element_Stiff_Matrix_Hl(i_E,i_H,
     &                                    c_POS_Hl(:,i_H),
     &                                    Location_ESM_C_Crack,
     &                                    num_Loc_ESM_C_Crack,
     &                                    Location_ESM_C_Cr_NoFEM,
     &                                    num_Loc_ESM_C_Cr_NoFEM)
              ! Includes FEM degrees of freedom
              Location_ESM(num_Loc_ESM+1:
     &                  num_Loc_ESM+num_Loc_ESM_C_Crack) = 
     &                  Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
              num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
            end do
          endif
          !Decide the location of each element stiffness matrix in the global stiffness matrix for cross.   
          if(num_Cross/=0)then
            do i_Cross =1,num_Cross 
              call Location_Element_Stiff_Matrix_Cross(i_E,i_Cross,
     &                                    c_POS_Cross(:,i_Cross),
     &                                    Location_ESM_C_Crack,
     &                                    num_Loc_ESM_C_Crack,
     &                                    Location_ESM_C_Cr_NoFEM,
     &                                    num_Loc_ESM_C_Cr_NoFEM)
              ! Includes FEM degrees of freedom
              Location_ESM(num_Loc_ESM+1:
     &                  num_Loc_ESM+num_Loc_ESM_C_Crack) = 
     &                  Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
              num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
            end do
          endif
          !Decide the location of each element stiffness matrix in the global stiffness matrix for inclusion.   
          if(num_Inclusion/=0)then
            do i_Incl =1,num_Inclusion
              call Location_Element_Stiff_Matrix_Incl(i_E,i_Incl,
     &                                    c_POS_Incl(:,i_Incl),
     &                                    Location_ESM_C_Crack,
     &                                    num_Loc_ESM_C_Crack,
     &                                    Location_ESM_C_Cr_NoFEM,
     &                                    num_Loc_ESM_C_Cr_NoFEM)              
              ! Includes FEM degrees of freedom
              Location_ESM(num_Loc_ESM+1:
     &                  num_Loc_ESM+num_Loc_ESM_C_Crack) = 
     &                  Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
              num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
            end do
          endif
          U(1:num_Loc_ESM) = c_DISP(Location_ESM(1:num_Loc_ESM))
          !***********
          ! Node Loop
          !***********
          do i_N = 1,4
              c_D     = D(Elem_Mat(i_E),:,:)     
              
              ! Elastic modulus of the Weibull distribution. 2024-06-25. NEWFTU2024062402.
              if(Flag_Weibull_E)then
                  if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
                      c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
                  endif
              endif
          
              c_v     = v(Elem_Mat(i_E),1)
              c_T_Alpha = T_Alpha(Elem_Mat(i_E))
              C_Node = c_NN(i_N)
              c_v_all(C_Node) = c_v
              c_kesi = T_Kesi(i_N)
              c_yita = T_Yita(i_N)
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !Calculate the B Matrix, Loop through each crack.
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              B(1:3,1:80) = ZR
              num_B = 0 
              !Calculate the B Matrix, Loop through each crack.
              do i_C =1,num_Crack 
                  call Cal_B_Matrix_Crack(c_kesi,c_yita,
     &                              i_C,i_E,i_N,
     &                              c_NN,c_X_NODES,c_Y_NODES,
     &                              tem_B,num_tem_B)                  
                  B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                          tem_B(1:3,1:num_tem_B)
                  num_B = num_B + num_tem_B  
              end do
              !Calculate the B Matrix, Loop through each hole.
              if(num_Hole/=0)then
                  do i_H =1,num_Hole 
                      call Cal_B_Matrix_Hl(c_kesi,c_yita,
     &                                    i_H,i_E,i_N,
     &                                    c_NN,c_X_NODES,c_Y_NODES,
     &                                    tem_B,num_tem_B)
                      B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                    tem_B(1:3,1:num_tem_B)
                      num_B = num_B + num_tem_B
                  end do
              endif
              !Calculate the B Matrix, Loop through each cross.
              if(num_Cross/=0)then
                  do i_Cross =1,num_Cross 
                      call Cal_B_Matrix_Cross(c_kesi,c_yita,
     &                                    i_Cross,i_E,i_N,
     &                                    c_NN,c_X_NODES,c_Y_NODES,
     &                                    tem_B,num_tem_B)
                      B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                    tem_B(1:3,1:num_tem_B)
                      num_B = num_B + num_tem_B
                  end do
              endif
              !Calculate the B Matrix, Loop through each circle inclusion.
              if(num_Circ_Incl/=0)then
                  do i_Incl =1,num_Circ_Incl 
                      call Cal_B_Matrix_Circ_Incl(c_kesi,c_yita,
     &                                        i_Incl,i_E,i_N,
     &                                        c_NN,c_X_NODES,c_Y_NODES,
     &                                        tem_B,num_tem_B)
                      B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                        tem_B(1:3,1:num_tem_B)
                      num_B = num_B + num_tem_B
                  end do
                  ! Determine whether the current node is within a certain mixture
                  call Tool_Yes_Point_in_Inclusions(
     &                                c_X_NODES(i_N),c_Y_NODES(i_N),
     &                                Yes_Gauss_in_Incl,c_Incl_Num)
                  if(Yes_Gauss_in_Incl)then
                      ! Obtain the current composite material matrix D
                      c_MatNum  = Circ_Inclu_Mat_Num(c_Incl_Num)
                      c_D       = D(c_MatNum,:,:)
                      c_v       = v(c_MatNum,1)
                      c_T_Alpha = T_Alpha(c_MatNum)
                  endif
              endif
              !Calculate the B Matrix, Loop through each polygon inclusion.
              if(num_Poly_Incl/=0)then
                  do i_Incl =1,num_Poly_Incl 
                      call Cal_B_Matrix_Poly_Incl(c_kesi,c_yita,
     &                                        i_Incl,i_E,i_N,
     &                                        c_NN,c_X_NODES,c_Y_NODES,
     &                                        tem_B,num_tem_B)
                      B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                        tem_B(1:3,1:num_tem_B)
                      num_B = num_B + num_tem_B
                  end do
                  ! Determine whether the current node is within a certain mixture
                  call Tool_Yes_Point_in_Inclusions(
     &                                c_X_NODES(i_N),c_Y_NODES(i_N),
     &                                Yes_Gauss_in_Incl,c_Incl_Num)
                  if(Yes_Gauss_in_Incl)then
                      ! Obtain the current composite material matrix D
                      c_MatNum = Poly_Inclu_Mat_Num(c_Incl_Num)
                      c_D     = D(c_MatNum,:,:)   
                      c_v = v(c_MatNum,1)
                      c_T_Alpha = T_Alpha(c_MatNum)
                  endif
              endif
              c_Stress=MATMUL(MATMUL(c_D,B(1:3,1:num_Loc_ESM)),
     &                        U(1:num_Loc_ESM))
              ! Thermal stress
              c_TStress(1:3) =ZR
              if(Key_Thermal_Stress==1)then
                  if(Key_Type_2D==1)then
                    c_TStress = c_T_Alpha*
     &                         Thermal_Str_Temper(Elem_Mat(i_E))*
     &                         MATMUL(c_D,[ONE,ONE,ZR])
                  elseif(Key_Type_2D==2)then
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
              endif
            c_Count(C_Node) = c_Count(C_Node) + 1
              
              ! Subtract thermal stress
              Stress_xx_Node(C_Node) = Stress_xx_Node(C_Node) + 
     &                                 c_Stress(1)-c_TStress(1)
            Stress_yy_Node(C_Node) = Stress_yy_Node(C_Node) + 
     &                                 c_Stress(2)-c_TStress(2)
            Stress_xy_Node(C_Node) = Stress_xy_Node(C_Node) + 
     &                                 c_Stress(3)-c_TStress(3)    
              
              ! Add initial stress field
              if(Key_InSitu_Strategy==2)then
                  Stress_xx_Node(C_Node) = Stress_xx_Node(C_Node) +
     &                                     Str_xx_InSitu(C_Node)
                  Stress_yy_Node(C_Node) = Stress_yy_Node(C_Node) +
     &                                     Str_yy_InSitu(C_Node)
                  Stress_xy_Node(C_Node) = Stress_xy_Node(C_Node) +
     &                                     Str_xy_InSitu(C_Node)
              endif
              ! If it is a killed element, the Gaussian stress is 0.
              if(Yes_Killed_Ele==1) then
                  Stress_xx_Node(C_Node) = ZR
                  Stress_yy_Node(C_Node) = ZR
                  Stress_xy_Node(C_Node) = ZR
              endif   
            end do
      end do  
      
c     --------------------
C     Node Stress Average
c     --------------------
      do i_N=1,num_Node   
          Stress_xx_Node(i_N) = Stress_xx_Node(i_N)/c_Count(i_N)
            Stress_yy_Node(i_N) = Stress_yy_Node(i_N)/c_Count(i_N)
            Stress_xy_Node(i_N) = Stress_xy_Node(i_N)/c_Count(i_N)
     
          call Tool_von_Mises(Stress_xx_Node(i_N),
     &                            Stress_yy_Node(i_N),
     &                            Stress_xy_Node(i_N),
     &                            c_v_all(i_N),Stress_vm_Node(i_N))
          if(Key_Thermal_Stress==1)then
            TStress_xx_Node(i_N) = TStress_xx_Node(i_N)/c_Count(i_N)
              TStress_yy_Node(i_N) = TStress_yy_Node(i_N)/c_Count(i_N)
              TStress_xy_Node(i_N) = TStress_xy_Node(i_N)/c_Count(i_N)
              TStress_vm_Node(i_N) = sqrt(TStress_xx_Node(i_N)**2+
     &                                  TStress_yy_Node(i_N)**2-
     &                                  TStress_xx_Node(i_N)*
     &                                  TStress_yy_Node(i_N)   +
     &                                  THR*TStress_xy_Node(i_N)**2) 
          call Tool_von_Mises(TStress_xx_Node(i_N),
     &                            TStress_yy_Node(i_N),
     &                            TStress_xy_Node(i_N),
     &                            c_v_all(i_N),TStress_vm_Node(i_N))
          endif
      end do
      
      
      RETURN
      END SUBROUTINE Get_Node_Stress_XFEM
