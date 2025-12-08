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
 
      SUBROUTINE Get_Gauss_Stress_XFEM(isub,c_DISP)
c     Calculate Gauss point stress
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
      integer i_E,i_G,c_NN(4),c_Num_Gauss_Point,G_Counter
      real(kind=FT) c_thick,c_D(3,3),U(80)
      real(kind=FT) c_X_NODES(4),c_Y_NODES(4),
     &                 c_kesi,c_yita,c_Stress(3),B(3,80),
     &                 tem_B(3,80)
      integer i_C,num_B,num_tem_B
      integer:: Location_ESM(MDOF_2D)
      integer Location_ESM_C_Cr_NoFEM(60)
      integer num_Loc_ESM_C_Cr_NoFEM
      integer num_Loc_ESM
      integer::Location_ESM_C_Crack(80)
      integer num_Loc_ESM_C_Crack
      real(kind=FT) kesi_Enr(Num_Gauss_Points),
     &                 yita_Enr(Num_Gauss_Points),
     &                 weight_Enr(Num_Gauss_Points)
      real(kind=FT) kesi_N_Enr(Num_Gauss_P_FEM),
     &                 yita_N_Enr(Num_Gauss_P_FEM)
      real(kind=FT)weight_N_Enr(Num_Gauss_Points)
      real(kind=FT) kesi(900),yita(900)
      integer i_H,i_Incl
      real(kind=FT) detJ,c_N(2,8)
      real(kind=FT) c_G_x,c_G_y
      logical Yes_Gauss_in_Incl
      integer c_Incl_Num,c_MatNum
      real(kind=FT) kesi_Enr_64(64),
     &              yita_Enr_64(64),
     &              weight_Enr_64(64)    
      real(kind=FT) c_T_Alpha,c_TStress(3),c_v 
      real(kind=FT) c_Ele_Gauss(3)
      integer i_Cross
      integer num_Ele_Killed,Yes_Killed_Ele
      real(kind=FT) c_S_1,c_S_3,c_theta_stress
      print *,'    Calculating stress of all Gauss points...'
      
c     -----------------
C     Initialized to 0
c     -----------------
      G_Counter  =0 
      
      ! Standard 64 Gauss points
      if (Key_Integral_Sol  == 2)then
          call Cal_Gauss_Points_QUAD(Num_Gauss_Points,kesi_Enr,yita_Enr,
     &                           weight_Enr)
          call Cal_Gauss_Points_QUAD(64,kesi_Enr_64,yita_Enr_64,
     &                         weight_Enr_64)
          call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,kesi_N_Enr,
     &                           yita_N_Enr,weight_N_Enr)
      ! Subdivision of quadrilaterals to calculate Gauss points (Gauss points are unrelated to cracks,
      ! they consist of several regular quadrilaterals)
      elseif (Key_Integral_Sol  == 3)then
          call Cal_Gauss_Points_QUAD_for_SUBQUAD(Num_Sub_Quads,
     &                           kesi_Enr,yita_Enr,
     &                           weight_Enr)
          Num_Gauss_Points = Num_Sub_Quads*4
          Num_Gauss_P_Inc  = Num_Sub_Quads*4
          call Cal_Gauss_Points_QUAD(64,kesi_Enr_64,yita_Enr_64,
     &                         weight_Enr_64)
          call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,kesi_N_Enr,
     &                           yita_N_Enr,weight_N_Enr)
      endif
      ! Life-and-Death element preparation
      if(Key_EKILL==1)then
          num_Ele_Killed = count(Ele_Killed_Each_Load_Step(1:isub,:)>0)
      endif      
     
c     -----------------
C     Inter-unit cycle
c     -----------------
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
          c_Ele_Gauss(1:3) = ZR
          
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
          if(num_Crack/=0)then
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
          endif
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
        !U(1:50) = ZR
          U(1:num_Loc_ESM) = c_DISP(Location_ESM(1:num_Loc_ESM))
          ! ------------------------------
          ! Point Scheme 1: Triangulation
          ! ------------------------------
          if(Key_Integral_Sol.eq.1)then
              !call Cal_Gauss_Points_Subtriangle(kesi,yita,weight)	
              !c_Num_Gauss_Point  = size(kesi,2)
          ! ----------------------------------------------------------------------------------
          ! Integration scheme 2 or 3: a fixed number of integration points, Num_Gauss_Points
          ! ----------------------------------------------------------------------------------
          elseif(Key_Integral_Sol.eq.2 .or.
     &           Key_Integral_Sol.eq.3)then
              !if the current element are enriched element, then 8x8 gauss points is suggested:
              if((num_Crack/= 0) .and. 
     &            maxval(Enriched_Node_Type(c_NN,1:num_Crack)).ne.0)then
                  kesi(1:Num_Gauss_Points)   = kesi_Enr
                  yita(1:Num_Gauss_Points)   = yita_Enr
                  c_Num_Gauss_Point = Num_Gauss_Points
              ! If it is a Hole enhanced node, then 8x8 Gauss points are suggested:
              elseif(num_Hole/= 0 .and.
     &            (maxval(Enriched_Node_Type_Hl(c_NN,1:num_Hole)).ne.0))
     &                                                             then
                  kesi(1:Num_Gauss_Points)   = kesi_Enr
                  yita(1:Num_Gauss_Points)   = yita_Enr
                  c_Num_Gauss_Point = Num_Gauss_Points
              ! If it is a Cross enhanced node, then 8x8 Gauss points are suggested:
              elseif(num_Cross/= 0 .and.
     &        (maxval(Enriched_Node_Type_Cross(c_NN,1:num_Cross)).ne.0))
     &                                                             then
                  kesi(1:Num_Gauss_Points)   = kesi_Enr
                  yita(1:Num_Gauss_Points)   = yita_Enr
                  c_Num_Gauss_Point = Num_Gauss_Points
              ! If it is a hybrid reinforcement node, it should be divided into two situations. For reinforcement
              ! elements that include material interfaces, use at least 400 Gauss points.
              ! For a general unit, 64 integration points are used.
              elseif(num_Inclusion/= 0 .and.
     &               (maxval(Enriched_Node_Type_Incl
     &                                   (c_NN,1:num_Inclusion)).ne.0))
     &                                                             then
                  kesi(1:Num_Gauss_Points)   = kesi_Enr
                  yita(1:Num_Gauss_Points)   = yita_Enr
                  c_Num_Gauss_Point = Num_Gauss_Points
              !if the current element are not enriched element, then 2x2 gauss points:
              else 
                  kesi(1:Num_Gauss_P_FEM)   = kesi_N_Enr
                  yita(1:Num_Gauss_P_FEM)   = yita_N_Enr
                  c_Num_Gauss_Point = Num_Gauss_P_FEM
              end if             
          end if    
          !******************
          ! Gauss Point Loop
          !******************
          do i_G = 1,c_Num_Gauss_Point
              c_D     = D(Elem_Mat(i_E),:,:) 
              
              ! Elastic modulus of the Weibull distribution. 2024-06-25. NEWFTU2024062402.
              if(Flag_Weibull_E)then
                  if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
                      c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
                  endif
              endif
          
              c_v     = v(Elem_Mat(i_E),1)
              c_T_Alpha = T_Alpha(Elem_Mat(i_E))
              G_Counter =  G_Counter +1    
              c_kesi = kesi(i_G)
              c_yita = yita(i_G)
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !Calculate the B Matrix, Loop through each crack.
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              B(1:3,1:80) = ZR
              num_B = 0 
              call Cal_N(kesi(i_G),yita(i_G),c_N)
              call Cal_detJ(kesi(i_G),yita(i_G),
     &                          c_X_NODES,c_Y_NODES,detJ)  
              !Global coordinates of the gauss point.
              c_G_x = DOT_PRODUCT(c_N(1,1:7:2),c_X_NODES(1:4))
              c_G_y = DOT_PRODUCT(c_N(1,1:7:2),c_Y_NODES(1:4)) 
              !Calculate the B Matrix, Loop through each crack.
              if(num_Crack/=0)then
                do i_C =1,num_Crack 
                  call Cal_B_Matrix_Crack(c_kesi,c_yita,i_C,i_E,i_G,
     &                              c_NN,c_X_NODES,c_Y_NODES,
     &                              tem_B,num_tem_B)                  
                  B(1:3,num_B+1:num_B+num_tem_B)=tem_B(1:3,1:num_tem_B)
                  num_B = num_B + num_tem_B                      
                end do
              endif
              !Calculate the B Matrix, Loop through each hole.
              if(num_Hole/=0)then
                  do i_H =1,num_Hole 
                      call Cal_B_Matrix_Hl(kesi(i_G),yita(i_G),
     &                                    i_H,i_E,i_G,
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
                      call Cal_B_Matrix_Cross(kesi(i_G),yita(i_G),
     &                                    i_Cross,i_E,i_G,
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
                      call Cal_B_Matrix_Circ_Incl(kesi(i_G),yita(i_G),
     &                                        i_Incl,i_E,i_G,
     &                                        c_NN,c_X_NODES,c_Y_NODES,
     &                                        tem_B,num_tem_B)
                      B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                        tem_B(1:3,1:num_tem_B)
                      num_B = num_B + num_tem_B
                  end do
                  ! Determine whether the current Gauss point is within a certain inclusion
                  call Tool_Yes_Point_in_Inclusions(c_G_x,c_G_y,
     &                                Yes_Gauss_in_Incl,c_Incl_Num)
                  if(Yes_Gauss_in_Incl)then
                      ! Obtain the current composite material matrix D
                      c_MatNum = Circ_Inclu_Mat_Num(c_Incl_Num)
                      c_D     = D(c_MatNum,:,:)   
                      c_T_Alpha = T_Alpha(c_MatNum)
                      c_v = v(c_MatNum,1)
                  endif
              endif
              !Calculate the B Matrix, Loop through each polygon inclusion.
              if(num_Poly_Incl/=0)then
                  do i_Incl =1,num_Poly_Incl 
                      call Cal_B_Matrix_Poly_Incl(kesi(i_G),yita(i_G),
     &                                        i_Incl,i_E,i_G,
     &                                        c_NN,c_X_NODES,c_Y_NODES,
     &                                        tem_B,num_tem_B)
                      B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                        tem_B(1:3,1:num_tem_B)
                      num_B = num_B + num_tem_B
                  end do
                  ! Determine whether the current Gauss point is within a certain inclusion
                  call Tool_Yes_Point_in_Inclusions(c_G_x,c_G_y,
     &                                Yes_Gauss_in_Incl,c_Incl_Num)
                  if(Yes_Gauss_in_Incl)then
                      ! Obtain the current composite material matrix D
                      c_MatNum = Poly_Inclu_Mat_Num(c_Incl_Num)
                      c_D     = D(c_MatNum,:,:)   
                      c_T_Alpha = T_Alpha(c_MatNum)
                      c_v = v(c_MatNum,1)
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

                  TStress_xx_Gauss(G_Counter) = c_TStress(1)
                TStress_yy_Gauss(G_Counter) = c_TStress(2)
                TStress_xy_Gauss(G_Counter) = c_TStress(3)  
                    TStress_vm_Gauss(G_Counter) = sqrt(c_TStress(1)**2+
     &                                     c_TStress(2)**2-
     &                                     c_TStress(1)*c_TStress(2)+
     &                                     THR*c_TStress(3)**2) 
              endif
              ! Stress
              Stress_xx_Gauss(G_Counter) = c_Stress(1)-c_TStress(1)
            Stress_yy_Gauss(G_Counter) = c_Stress(2)-c_TStress(2)
            Stress_xy_Gauss(G_Counter) = c_Stress(3)-c_TStress(3) 
              
              ! Accumulation of stresses at each Gauss point of the current element (used to calculate the average
              ! Gauss stress of each element)
              c_Ele_Gauss(1) = c_Ele_Gauss(1)+c_Stress(1)-c_TStress(1)
              c_Ele_Gauss(2) = c_Ele_Gauss(2)+c_Stress(2)-c_TStress(2)
              c_Ele_Gauss(3) = c_Ele_Gauss(3)+c_Stress(3)-c_TStress(3)
              
              ! Add the initial stress field (including the average geostress of the four Gauss points of the
              ! current FEM element)
              if(Key_InSitu_Strategy==2)then
                  Stress_xx_Gauss(G_Counter)=Stress_xx_Gauss(G_Counter)
     &              + sum(InSitu_Strs_Gaus_xx(i_E,1:Num_Gauss_P_FEM))/
     &                Num_Gauss_P_FEM
                  Stress_yy_Gauss(G_Counter)=Stress_yy_Gauss(G_Counter)
     &              + sum(InSitu_Strs_Gaus_yy(i_E,1:Num_Gauss_P_FEM))/
     &                Num_Gauss_P_FEM
                  Stress_xy_Gauss(G_Counter)=Stress_xy_Gauss(G_Counter)
     &              + sum(InSitu_Strs_Gaus_xy(i_E,1:Num_Gauss_P_FEM))/
     &                Num_Gauss_P_FEM 
                 ! Accumulation of stresses at each Gauss point of the current element (used to calculate the average
                 ! Gauss stress of each element)
                 c_Ele_Gauss(1) = c_Ele_Gauss(1)+ 
     &                  sum(InSitu_Strs_Gaus_xx(i_E,1:Num_Gauss_P_FEM))/
     &                  Num_Gauss_P_FEM
                 c_Ele_Gauss(2) = c_Ele_Gauss(2)
     &              + sum(InSitu_Strs_Gaus_yy(i_E,1:Num_Gauss_P_FEM))/
     &                Num_Gauss_P_FEM
                 c_Ele_Gauss(3) = c_Ele_Gauss(3)
     &              + sum(InSitu_Strs_Gaus_xy(i_E,1:Num_Gauss_P_FEM))/
     &                Num_Gauss_P_FEM             
              endif
              ! If the element is killed or destroyed, the Gauss stress is 0.
              if(Yes_Killed_Ele==1 .or. Elem_Break(i_E)) then
                  Stress_xx_Gauss(G_Counter) = ZR
                  Stress_yy_Gauss(G_Counter) = ZR
                  Stress_xy_Gauss(G_Counter) = ZR
                  ! Accumulation of stresses at each Gauss point of the current element (used to calculate the average
                  ! Gauss stress of each element)
                  c_Ele_Gauss(1) = ZR
                  c_Ele_Gauss(2) = ZR
                  c_Ele_Gauss(3) = ZR                  
              endif
              
              call Tool_von_Mises(Stress_xx_Gauss(G_Counter),
     &                                Stress_yy_Gauss(G_Counter),
     &                                Stress_xy_Gauss(G_Counter),
     &                                c_v,Stress_vm_Gauss(G_Counter)) 
              !if (Stress_vm_Gauss(G_Counter)>1000) then
              !endif
            end do
          ! Average Gauss stress of each element
          Elem_Ave_Gauss_Stress(i_E,1:3) = c_Ele_Gauss(1:3)/
     &                                     c_Num_Gauss_Point
          ! Calculate the average principal stress for each element
          call Tool_Principal_Stresses_2D(Elem_Ave_Gauss_Stress(i_E,1),
     &                                    Elem_Ave_Gauss_Stress(i_E,2),
     &                                    Elem_Ave_Gauss_Stress(i_E,3),
     &                                    c_S_1,c_S_3,c_theta_stress)
          Elem_Ave_Gauss_S1(i_E) =  c_S_1
      end do  
      RETURN
      END SUBROUTINE Get_Gauss_Stress_XFEM
