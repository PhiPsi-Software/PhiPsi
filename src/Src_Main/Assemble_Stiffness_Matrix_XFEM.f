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
 
      SUBROUTINE Assemble_Stiffness_Matrix_XFEM(isub,
     &                c_globalK,c_Total_Freedom,Total_Num_G_P)
c     Assemble the stiffness matrix.

      !**********************
      ! Read public variable
      !**********************
      use Global_Float_Type
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Filename
      use Global_Common
      use Global_Material
      use Global_HF
      use Global_Contact
      use Global_Inclusion
      use Global_Cross
      
      !**********************
      ! Variable Declaration
      !**********************
      implicit none
      integer,intent(in)::isub,c_Total_Freedom
      integer,intent(out)::Total_Num_G_P
      real(kind=FT) ,intent(out)::c_globalK(c_Total_Freedom,
     &                                         c_Total_Freedom)
      integer i_C,i_E,i_G,i_Incl,i_H
      real(kind=FT) c_thick,c_D(3,3)
      real(kind=FT) c_X_NODES(4),c_Y_NODES(4)
      integer c_NN(4)
      real(kind=FT) kesi_Enr(Num_Gauss_Points),
     &                 yita_Enr(Num_Gauss_Points),
     &                 weight_Enr(Num_Gauss_Points)    
      real(kind=FT) kesi_N_Enr(Num_Gauss_P_FEM),
     &                 yita_N_Enr(Num_Gauss_P_FEM),
     &                 weight_N_Enr(Num_Gauss_P_FEM)     
      real(kind=FT) kesi(900),yita(900),weight(900)
      integer:: Location_ESM(MDOF_2D)
      integer::Location_ESM_C_Crack(80),Location_ESM_C_Cr_NoFEM(60)
      integer num_Loc_ESM_C_Cr_NoFEM
      real(kind=FT) B(3,80),tem_B(3,80)
      integer num_B,num_tem_B
      integer num_Loc_ESM,num_Loc_ESM_C_Crack
      integer c_Num_Gauss_Point
      real(kind=FT) detJ,c_N(2,8)
      real(kind=FT) c_G_x,c_G_y
      logical Yes_Gauss_in_Incl
      integer c_Incl_Num,c_MatNum
      real(kind=FT) kesi_Enr_64(64),
     &              yita_Enr_64(64),
     &              weight_Enr_64(64)    
      real(kind=FT) penalty_k
      integer i_CP_node,c_CP_node,first_node
      integer first_node_DOF,c_node_DOF
      integer i_Cross
      integer i_Boux_nonzero,i_Bouy_nonzero,c_node
      integer c_i
      integer num_Ele_Killed
      
      ! Test variable
      real(kind=FT) localK(MDOF_2D,MDOF_2D)
      
      
      !*************************
      ! Variable Initialization
      !*************************
      kesi(1:900)  = ZR
      yita(1:900)  = ZR
      weight(1:900)= ZR
      
      ! Life-and-Death element preparation
      if(Key_EKILL==1)then
          num_Ele_Killed = count(Ele_Killed_Each_Load_Step(1:isub,:)>0)
      endif      
      

      !****************************************************************************************
      ! Initialize the stiffness matrix and calculate the local coordinates and weights of the
      ! element Gauss points
      !****************************************************************************************
      c_globalK(1:c_Total_Freedom,1:c_Total_Freedom) = ZR
      Total_Num_G_P = 0
      
      ! Standard 64 Gauss points
      if (Key_Integral_Sol  == 2)then
          call Cal_Gauss_Points_QUAD(Num_Gauss_Points,kesi_Enr,yita_Enr,
     &                           weight_Enr)
          call Cal_Gauss_Points_QUAD(64,kesi_Enr_64,yita_Enr_64,
     &                         weight_Enr_64)
          call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,kesi_N_Enr,
     &                           yita_N_Enr,weight_N_Enr)
      ! Subdivision of quadrilaterals to calculate Gauss points (the Gauss points are unrelated to cracks,
      ! consisting of several regular quadrilaterals)
      elseif (Key_Integral_Sol  == 3)then
          call Cal_Gauss_Points_QUAD_for_SUBQUAD(Num_Sub_Quads,
     &                           kesi_Enr,yita_Enr,
     &                           weight_Enr)
          call Cal_Gauss_Points_QUAD(64,kesi_Enr_64,yita_Enr_64,
     &                         weight_Enr_64)
          call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,kesi_N_Enr,
     &                           yita_N_Enr,weight_N_Enr)
      endif
      !**************************************************************
      ! Calculate the stiffness matrix, looping through each element
      !**************************************************************
      EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gauss_Points)= .False.

      do i_E = 1,Num_Elem
          
          localK(1:MDOF_2D,1:MDOF_2D) = ZR
          c_thick   = thick(Elem_Mat(i_E))
          c_D       = D(Elem_Mat(i_E),:,:)  
          
          ! Elastic modulus of the Weibull distribution. 2024-06-25. NEWFTU2024062402.
          if(Flag_Weibull_E)then
              if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
                  c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
              endif
          endif
          
          ! Life and Death element Handling: Weaken D Matrix (Elastic Modulus)
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
          ! ------------------------------
          ! Point Scheme 1: Triangulation
          ! ------------------------------
          if(Key_Integral_Sol.eq.1)then
          
          ! -------------------------------------
          ! Point Plan 2 or 3: Even Distribution
          ! -------------------------------------
          elseif(Key_Integral_Sol.eq.2 .or. Key_Integral_Sol.eq.3)then
              ! Starting number of Gauss points for each element
              Ele_GP_Start_Num(i_E) = Total_Num_G_P + 1
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Determine the Gauss points of the current element
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !if the current element are enriched element, then 8x8 gauss points is suggested:
              if (num_Crack/= 0 .and. 
     &            (maxval(Enriched_Node_Type(c_NN,1:num_Crack)).ne.0))
     &                                                             then
                  kesi(1:Num_Gauss_Points)    = kesi_Enr
                  yita(1:Num_Gauss_Points)    = yita_Enr
                  weight(1:Num_Gauss_Points)  = weight_Enr
                  c_Num_Gauss_Point = Num_Gauss_Points
                  Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
              ! If it is a Hole enhanced node, then 8x8 Gauss points are suggested:
              elseif(num_Hole/= 0 .and.
     &            (maxval(Enriched_Node_Type_Hl(c_NN,1:num_Hole)).ne.0))
     &                                                             then
                  kesi(1:Num_Gauss_Points)    = kesi_Enr
                  yita(1:Num_Gauss_Points)    = yita_Enr
                  weight(1:Num_Gauss_Points)  = weight_Enr
                  c_Num_Gauss_Point = Num_Gauss_Points
                  Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
              ! If it is a Cross-enhanced node, then 8x8 Gauss points are suggested:
              elseif(num_Cross/= 0 .and.
     &        (maxval(Enriched_Node_Type_Cross(c_NN,1:num_Cross)).ne.0))
     &                                                             then
                  kesi(1:Num_Gauss_Points)    = kesi_Enr
                  yita(1:Num_Gauss_Points)    = yita_Enr
                  weight(1:Num_Gauss_Points)  = weight_Enr
                  c_Num_Gauss_Point = Num_Gauss_Points
                  Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
              ! If there are embedded reinforcement nodes, then 8x8 Gauss points are suggested:
              elseif(num_Inclusion/= 0 .and.
     &            (maxval(Enriched_Node_Type_Incl
     &                        (c_NN,1:num_Inclusion)).ne.0))then
                  kesi(1:Num_Gauss_Points)    = kesi_Enr
                  yita(1:Num_Gauss_Points)    = yita_Enr
                  weight(1:Num_Gauss_Points)  = weight_Enr
                  c_Num_Gauss_Point = Num_Gauss_Points
                  Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
              !if the current element are not enriched element, then 2x2 gauss points:
              else 
                  kesi(1:Num_Gauss_P_FEM)    = kesi_N_Enr
                  yita(1:Num_Gauss_P_FEM)    = yita_N_Enr
                  weight(1:Num_Gauss_P_FEM)  = weight_N_Enr
                  c_Num_Gauss_Point = Num_Gauss_P_FEM
                  Total_Num_G_P     = Total_Num_G_P +c_Num_Gauss_Point
              end if
              ! Save the number of Gauss points for each element to a global variable
              num_GP_Elem(i_E) = c_Num_Gauss_Point
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !Decide the location of each element stiffness 
              !matrix in the global stiffness matrix.   
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              Location_ESM(1:MDOF_2D)   = 0
              num_Loc_ESM           = 0
              
              if(num_Crack/=0)then
                do i_C =1,num_Crack
                  call Location_Element_Stiff_Matrix(i_E,i_C,
     &                                c_POS(1:Num_Node,i_C),
     &                                Location_ESM_C_Crack,
     &                                num_Loc_ESM_C_Crack,
     &                                Location_ESM_C_Cr_NoFEM,
     &                                num_Loc_ESM_C_Cr_NoFEM)
                  ! Includes FEM degrees of freedom
                  Location_ESM(num_Loc_ESM+1:
     &                      num_Loc_ESM+num_Loc_ESM_C_Crack) = 
     &                      Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
                  num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
                end do
              endif
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !Decide the location of each element stiffness 
              !matrix in the global stiffness matrix for hole.   
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if(num_Hole/=0)then
                do i_H =1,num_Hole
                  call Location_Element_Stiff_Matrix_Hl(i_E,i_H,
     &                                c_POS_Hl(1:Num_Node,i_H),
     &                                Location_ESM_C_Crack,
     &                                num_Loc_ESM_C_Crack,
     &                                Location_ESM_C_Cr_NoFEM,
     &                                num_Loc_ESM_C_Cr_NoFEM)
                  ! Includes FEM degrees of freedom
                  Location_ESM(num_Loc_ESM+1:
     &                      num_Loc_ESM+num_Loc_ESM_C_Crack) = 
     &                      Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
                  num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
                end do
              endif
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !Decide the location of each element stiffness 
              !matrix in the global stiffness matrix for cross.   
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if(num_Cross/=0)then
                do i_Cross =1,num_Cross
                  call Location_Element_Stiff_Matrix_Cross(i_E,i_Cross,
     &                                c_POS_Cross(1:Num_Node,i_Cross),
     &                                Location_ESM_C_Crack,
     &                                num_Loc_ESM_C_Crack,
     &                                Location_ESM_C_Cr_NoFEM,
     &                                num_Loc_ESM_C_Cr_NoFEM)
                  ! Includes FEM degrees of freedom
                  Location_ESM(num_Loc_ESM+1:
     &                      num_Loc_ESM+num_Loc_ESM_C_Crack) = 
     &                      Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
                  num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
                end do
              endif
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !Decide the location of each element stiffness 
              !matrix in the global stiffness matrix for inclusion.   
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if(num_Inclusion/=0)then
                do i_Incl =1,num_Inclusion
                  call Location_Element_Stiff_Matrix_Incl(i_E,i_Incl,
     &                                   c_POS_Incl(1:Num_Node,i_Incl),
     &                                   Location_ESM_C_Crack,
     &                                   num_Loc_ESM_C_Crack,
     &                                   Location_ESM_C_Cr_NoFEM,
     &                                   num_Loc_ESM_C_Cr_NoFEM)
                  ! Includes FEM degrees of freedom
                  Location_ESM(num_Loc_ESM+1:
     &                      num_Loc_ESM+num_Loc_ESM_C_Crack) = 
     &                      Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
                  num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
                end do
              endif
              
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Loop through each Gauss point to assemble the stiffness matrix of the current
              ! element
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              do i_G = 1,c_Num_Gauss_Point
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
                      call Cal_B_Matrix_Crack(kesi(i_G),yita(i_G),
     &                                    i_C,i_E,i_G,
     &                                    c_NN,c_X_NODES,c_Y_NODES,
     &                                    tem_B,num_tem_B)
                      B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                          tem_B(1:3,1:num_tem_B)
                      
                      num_B = num_B + num_tem_B
                    end do
                  endif
                  !Calculate the B Matrix, Loop through each hole.
                  if(num_Hole/=0)then
                      do i_H =1,num_Hole 
                          call Cal_B_Matrix_Hl(kesi(i_G),yita(i_G),
     &                                        i_H,i_E,i_G,
     &                                        c_NN,c_X_NODES,c_Y_NODES,
     &                                        tem_B,num_tem_B)
                          B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                        tem_B(1:3,1:num_tem_B)
                          num_B = num_B + num_tem_B
                      end do
                  endif
                  !Calculate the B Matrix, Loop through each cross.
                  if(num_Cross/=0)then
                      do i_Cross =1,num_Cross
                          call Cal_B_Matrix_Cross(kesi(i_G),yita(i_G),
     &                                        i_Cross,i_E,i_G,
     &                                        c_NN,c_X_NODES,c_Y_NODES,
     &                                        tem_B,num_tem_B)
                          B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                        tem_B(1:3,1:num_tem_B)
                          num_B = num_B + num_tem_B
                      end do
                  endif
                  !Calculate the B Matrix, Loop through each circle inclusion.
                  if(num_Circ_Incl/=0)then
                      do i_Incl =1,num_Circ_Incl 
                          call Cal_B_Matrix_Circ_Incl(kesi(i_G),
     &                                       yita(i_G),i_Incl,i_E,i_G,
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
                      endif
                  endif
                  !Calculate the B Matrix, Loop through each polygon inclusion.
                  if(num_Poly_Incl/=0)then
                      do i_Incl =1,num_Poly_Incl
                          call Cal_B_Matrix_Poly_Incl(kesi(i_G),
     &                                       yita(i_G),i_Incl,i_E,i_G,
     &                                        c_NN,c_X_NODES,c_Y_NODES,
     &                                        tem_B,num_tem_B)
                          B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                        tem_B(1:3,1:num_tem_B)
                          num_B = num_B + num_tem_B
                      end do
                      ! Determine whether the current Gauss point is within an inclusion (including circular and polygonal
                      ! inclusions)
                      call Tool_Yes_Point_in_Inclusions(c_G_x,c_G_y,
     &                                Yes_Gauss_in_Incl,c_Incl_Num)
                      if(Yes_Gauss_in_Incl)then
                          ! Obtain the current composite material matrix D
                          c_MatNum = Poly_Inclu_Mat_Num(c_Incl_Num)
                          c_D     = D(c_MatNum,:,:)   
                      endif
                  endif
                  localK(1:num_B,1:num_B) = localK(1:num_B,1:num_B) +
     &                           c_thick*weight(i_G)*detJ*
     &                           MATMUL(MATMUL(transpose
     &                               (B(1:3,1:num_B)),c_D),
     &                                B(1:3,1:num_B))
              enddo
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !Assemble the global stiffness matrix. 
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              c_globalK(Location_ESM(1:num_Loc_ESM),
     &                  Location_ESM(1:num_Loc_ESM)) = 
     &            c_globalK(Location_ESM(1:num_Loc_ESM),
     &                      Location_ESM(1:num_Loc_ESM)) + 
     &            localK(1:num_B,1:num_B)
          end if      
      end do     

      ! Filter the calculation results and remove smaller values
      !option 1
C     do i=1,c_Total_Freedom
C         do j=1,c_Total_Freedom
C             if (abs(c_globalK(i,j))< = Tol_15) then
C                  c_globalK(i,j) = ZR 
C             endif
C         enddo
C     end do

      ! Option 2, speed greatly increased
      where (abs(c_globalK)< = Tol_15)
          c_globalK = ZR 
      end where
      
      
      !.................................................................
      ! Handling of non-zero displacement boundary conditions:
      ! Finite Element Methods in Engineering, Zeng Pan, Fourth Edition
      !.................................................................
      if(Num_Boux_nonzero >0)then
          !penalty_k =  1.0D4*maxval(c_globalK)
          penalty_k =  penalty_k_bou_nonzero 
          do i_Boux_nonzero =1,Num_Boux_nonzero
              c_node = int(Bou_x_nonzero(i_Boux_nonzero,1))
              c_node_DOF = 2*c_node-1
              c_globalK(c_node_DOF,c_node_DOF) =
     &                c_globalK(c_node_DOF,c_node_DOF)+penalty_k              
          enddo
      endif
      if(Num_Bouy_nonzero >0)then
          !penalty_k =  1.0D4*maxval(c_globalK)
          penalty_k =  penalty_k_bou_nonzero 
          do i_Bouy_nonzero =1,Num_Bouy_nonzero
              c_node = int(Bou_y_nonzero(i_Bouy_nonzero,1))
              c_node_DOF = 2*c_node
              c_globalK(c_node_DOF,c_node_DOF) =
     &                c_globalK(c_node_DOF,c_node_DOF)+penalty_k              
          enddo
      endif
      
      
      !...................................................................................
      ! Penalty stiffness method for node coupling, theoretical basis:
      !http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch09.d/IFEM.Ch09.pdf
      !...................................................................................
      ! Used for option 1 - directly defining coupled nodes through a keyword file (in this method, only
      ! one set of couplings can be defined for each direction)
      if (num_CP_x_nodes >0)then
          penalty_k =  1.0D4*maxval(c_globalK)
          first_node = CP_x_nodes(1)
          do i_CP_node =2,num_CP_x_nodes
              c_CP_node = CP_x_nodes(i_CP_node)
              first_node_DOF = 2*first_node-1
              c_node_DOF = 2*c_CP_node-1
              c_globalK(first_node_DOF,first_node_DOF) =
     &                c_globalK(first_node_DOF,first_node_DOF)+penalty_k
              c_globalK(c_node_DOF,c_node_DOF) =
     &                c_globalK(c_node_DOF,c_node_DOF)+penalty_k
              c_globalK(first_node_DOF,c_node_DOF) =
     &                c_globalK(first_node_DOF,c_node_DOF)-penalty_k
              c_globalK(c_node_DOF,first_node_DOF) =
     &                c_globalK(c_node_DOF,first_node_DOF)-penalty_k
          enddo
      endif
      if (num_CP_y_nodes > 0)then
          penalty_k =  1.0D4*maxval(c_globalK)
          first_node = CP_y_nodes(1)
          do i_CP_node =2,num_CP_y_nodes
              c_CP_node = CP_y_nodes(i_CP_node)
              first_node_DOF = 2*first_node
              c_node_DOF = 2*c_CP_node
              c_globalK(first_node_DOF,first_node_DOF) =
     &                c_globalK(first_node_DOF,first_node_DOF)+penalty_k
              c_globalK(c_node_DOF,c_node_DOF) =
     &                c_globalK(c_node_DOF,c_node_DOF)+penalty_k
              c_globalK(first_node_DOF,c_node_DOF) =
     &                c_globalK(first_node_DOF,c_node_DOF)-penalty_k
              c_globalK(c_node_DOF,first_node_DOF) =
     &                c_globalK(c_node_DOF,first_node_DOF)-penalty_k
          enddo
      endif
      ! Used for option 2 - defining coupled nodes through dofx, dofy, dofz files (this method can only
      ! define multiple sets of coupling for each direction)
      if (num_CP_set_x>0)then
          penalty_k =  1.0D4*maxval(c_globalK)
          do c_i = 1,num_CP_set_x
              first_node = CP_nodes_x(c_i,1)
              do i_CP_node =2, num_nodes_CP_set_x(c_i)
                  c_CP_node = CP_nodes_x(c_i,i_CP_node)
                  first_node_DOF = 2*first_node-1
                  c_node_DOF = 2*c_CP_node-1
                  c_globalK(first_node_DOF,first_node_DOF) =
     &                c_globalK(first_node_DOF,first_node_DOF)+penalty_k
                  c_globalK(c_node_DOF,c_node_DOF) =
     &                c_globalK(c_node_DOF,c_node_DOF)+penalty_k
                  c_globalK(first_node_DOF,c_node_DOF) =
     &                c_globalK(first_node_DOF,c_node_DOF)-penalty_k
                  c_globalK(c_node_DOF,first_node_DOF) =
     &                c_globalK(c_node_DOF,first_node_DOF)-penalty_k
              enddo
          enddo
      endif      
      if (num_CP_set_y>0)then
          penalty_k =  1.0D4*maxval(c_globalK)
          do c_i = 1,num_CP_set_y
              first_node = CP_nodes_y(c_i,1)
              do i_CP_node =2, num_nodes_CP_set_y(c_i)
                  c_CP_node = CP_nodes_y(c_i,i_CP_node)
                  first_node_DOF = 2*first_node
                  c_node_DOF = 2*c_CP_node
                  c_globalK(first_node_DOF,first_node_DOF) =
     &                c_globalK(first_node_DOF,first_node_DOF)+penalty_k
                  c_globalK(c_node_DOF,c_node_DOF) =
     &                c_globalK(c_node_DOF,c_node_DOF)+penalty_k
                  c_globalK(first_node_DOF,c_node_DOF) =
     &                c_globalK(first_node_DOF,c_node_DOF)-penalty_k
                  c_globalK(c_node_DOF,first_node_DOF) =
     &                c_globalK(c_node_DOF,first_node_DOF)-penalty_k
              enddo
          enddo
      endif   
      RETURN
      END SUBROUTINE Assemble_Stiffness_Matrix_XFEM
