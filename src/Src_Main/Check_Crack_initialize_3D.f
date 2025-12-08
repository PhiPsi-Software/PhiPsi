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
 
      subroutine Check_Crack_initialize_3D(ifra,iter,New_Crack_Flag)
c     Check if any new cracks have developed (3D problem).
c     2020-03-12      

c     **************
c     Public Module
c     **************
      use Global_Float_Type
      use Global_Crack_Common
      use Global_Crack_3D
      use Global_DISP
      use Global_Elem_Area_Vol
      use Global_Model
      use Global_HF
      use Global_Common
      use Global_Material
      use Global_Dynamic 
      
      implicit none
      integer,intent(in)::ifra,iter
      logical,intent(out)::New_Crack_Flag
      real(kind=FT) Stress_El_Center(num_Elem,7)  
      real(kind=FT) S1_El(num_Elem),S2_El(num_Elem),S3_El(num_Elem) 
      real(kind=FT) Shear_El(num_Elem)
      real(kind=FT) Vector_S1_El(num_Elem,3),Vector_S2_El(num_Elem,3),
     &              Vector_S3_El(num_Elem,3)
      integer i_E
      real(kind=FT) S_xx,S_yy,S_zz,S_yz,S_xz,S_xy
      real(kind=FT) S_1,S_2,S_3
      real(kind=FT) Vector_S1(3),Vector_S2(3),Vector_S3(3)
      integer max_S1_El,mat_num,c_num_Crack,max_Shear_El
      real(kind=FT) St_Critical,max_S1,max_Shear
      real(kind=FT) rotated_vector_a(3),vector_a(3),vector_b(3)
      real(kind=FT) c_Function_Concrete(Num_Elem)
      real(kind=FT) f_t(Max_Materials),f_c(Max_Materials)
      real(kind=FT) x,
     &             S1, S2,S3, a0(Max_Materials), 
     &             a1(Max_Materials), a2(Max_Materials),
     &             b0(Max_Materials), b1(Max_Materials),
     &             b2(Max_Materials)
      real(kind=FT) cos_yita,SQRT2,fcb,f1,f2,s_ha,Temp_Value
      real(kind=FT) tem_Matrix1(3,3),tem_F1(3),tem_Sol1(3),SQRT15 
      real(kind=FT) tem_Matrix2(3,3),tem_F2(3),tem_Sol2(3)
      real(kind=FT) Kesi_t,Kesi_cb,Kesi_1,Kesi_2,Kesi_0,a,b,c
      real(kind=FT) Root1,Root2,FuncF1,FuncF2,FuncF3,FuncF4     
      real(kind=FT) FuncS1,FuncS2,FuncS3,FuncS4,p1,p2    
      real(kind=FT) FuncS2_tem1,FuncS2_tem2,FuncS2_tem3,FuncS2_tem4
      integer max_F_El
      real(kind=FT) max_F,s_h,kesi,r1,r2
      real(kind=FT) FuncS1_tem1,FuncS1_tem2,FuncS1_tem3,FuncS1_tem4
      integer i_Mat,c_Mat
      real(kind=FT) Shear_Strength
      

c     *****************
c     Data Preparation
c     *****************
      New_Crack_Flag = .False.
     
c     ********************
c     Interactive display
c     ********************
      print *,'    Checking initialization of crack......'  
      
c     ***************************************************
c     Calculate the stress at the center of each element
c     ***************************************************
      call Get_Element_Center_Stress_FEM_3D(.True.,
     &        iter,DISP,Stress_El_Center)

c     *****************************************************
c     If it is the Maximum Tensile Stress Criterion (MTSC)
c     *****************************************************
      if(Key_Ini_Rule==1)then  
          !--------------------------------------------------------------
          ! Calculate the principal stress at the center of each element
          !--------------------------------------------------------------
          do i_E=1,Num_Elem   
              S_xx = Stress_El_Center(i_E,1) 
              S_yy = Stress_El_Center(i_E,3) 
              S_zz = Stress_El_Center(i_E,3) 
              S_xy = Stress_El_Center(i_E,4) 
              S_yz = Stress_El_Center(i_E,5) 
              S_xz = Stress_El_Center(i_E,6) 
                     
              call Tool_Principal_Stresses_3D(S_xx,S_yy,S_zz,
     &              S_yz,S_xz,S_xy,
     &              S1_El(i_E),S2_El(i_E),S3_El(i_E),
     &              Vector_S1_El(i_E,1:3),Vector_S2_El(i_E,1:3),
     &              Vector_S3_El(i_E,1:3))
          enddo
          !---------------------------
          !Find the element of max S1
          !---------------------------
          max_S1_El = maxloc(S1_El,1)
          max_S1    = S1_El(max_S1_El)
          mat_num   = Elem_Mat(max_S1_El)
          St_Critical  = Material_Para(mat_num,5)
          if(max_S1 >= St_Critical .and.
     &       MAT_ALLOW_CRACK_Initiation(mat_num)==1)then
              print *,'    >>>>>>>>> A Crack Emerged at Elememt ',
     &               max_S1_El, ' of material type ',
     &              mat_num, '<<<<<<<<<'
              New_Crack_Flag = .True.
          endif
          !--------------
          !Add the crack
          !--------------
          !Circle.
          if(Key_Ini_Cr_3D_Type==1 .and.
     &       MAT_ALLOW_CRACK_Initiation(mat_num)==1)then
              if(Num_Initiation_Cracks <
     &           Key_Max_Num_Initiation_Cracks)then
              !num_Crack = num_Crack + 1
              c_num_Crack = num_Crack + 1
              Num_Initiation_Cracks=  Num_Initiation_Cracks +1
              Crack3D_Cir_Coor(c_num_Crack,1:3) =
     &                         Elem_Centroid(max_S1_El,1:3)
              Crack3D_Cir_Coor(c_num_Crack,4:6) =
     &                         Vector_S1_El(max_S1_El,1:3)              
              Crack3D_Cir_Coor(c_num_Crack,7) =Size_Ini_Crack/TWO      
              endif
          !Rectangle.
          elseif(Key_Ini_Cr_3D_Type==2)then
          endif
          
      endif
      
c     ********************************************
c     If it is the maximum shear stress criterion
c     ********************************************
      if(Key_Ini_Rule==2)then  
          !--------------------------------------------------------------
          ! Calculate the principal stress at the center of each element
          !--------------------------------------------------------------
          do i_E=1,Num_Elem   
              S_xx = Stress_El_Center(i_E,1) 
              S_yy = Stress_El_Center(i_E,3) 
              S_zz = Stress_El_Center(i_E,3) 
              S_xy = Stress_El_Center(i_E,4) 
              S_yz = Stress_El_Center(i_E,5) 
              S_xz = Stress_El_Center(i_E,6) 
                     
              call Tool_Principal_Stresses_3D(S_xx,S_yy,S_zz,
     &              S_yz,S_xz,S_xy,
     &              S1_El(i_E),S2_El(i_E),S3_El(i_E),
     &              Vector_S1_El(i_E,1:3),Vector_S2_El(i_E,1:3),
     &              Vector_S3_El(i_E,1:3))
              Shear_El(i_E) = (S1_El(i_E)-S3_El(i_E))/TWO
          enddo
          !-------------------------------------
          !Find the element of max Shear stress
          !-------------------------------------
          max_Shear_El = maxloc(Shear_El,1)
          max_Shear    = Shear_El(max_Shear_El)
          mat_num   = Elem_Mat(max_Shear_El)
          ! St_Critical = Material_Para(mat_num,5)          !Tensile Strength
          Shear_Strength = Material_Para_Added(mat_num,12)
          if(max_Shear>=Shear_Strength .and.
     &       MAT_ALLOW_CRACK_Initiation(mat_num)==1)then   
              print *,'    >>>>>>>>> A Crack Emerged at Elememt ',
     &               max_Shear_El,' of material type ',
     &              mat_num, '<<<<<<<<<' 
              New_Crack_Flag = .True.
          endif
          !--------------
          !Add the crack
          !--------------
          !Circle.
          if(Key_Ini_Cr_3D_Type==1 .and.
     &       MAT_ALLOW_CRACK_Initiation(mat_num)==1)then
              if(Num_Initiation_Cracks <
     &           Key_Max_Num_Initiation_Cracks)then
              !num_Crack = num_Crack + 1
              c_num_Crack = num_Crack + 1
              Num_Initiation_Cracks=  Num_Initiation_Cracks +1
              Crack3D_Cir_Coor(c_num_Crack,1:3) =
     &                         Elem_Centroid(max_Shear_El,1:3)

              vector_a = Vector_S1_El(max_Shear_El,1:3)
              vector_b = Vector_S2_El(max_Shear_El,1:3)
              ! The normal vector of the plane with the maximum shear stress: obtained by rotating the direction
              ! of the maximum principal stress 45 degrees around the direction of the intermediate principal
              ! stress.
              call Tool_Vector_a_Rotate_around_Vector_b(
     &                   vector_a,vector_b,pi/FOU,rotated_vector_a)
              Crack3D_Cir_Coor(c_num_Crack,4:6) = rotated_vector_a
              Crack3D_Cir_Coor(c_num_Crack,7)   = Size_Ini_Crack/TWO      
              endif
          !Rectangle.
          elseif(Key_Ini_Cr_3D_Type==2)then
          
          endif
          
c     *******************************
c        Concrete rule.
c     Refe: ANSYS theory manual-4.11
c     *******************************
      elseif(Key_Ini_Rule==3)then        
          ! First calculate the relevant parameter values for various materials to avoid duplicate
          ! calculations. 2024-03-13. IMPROV2024031301.
          do i_Mat = 1,NUM_OF_MATERIAL
              ! Concrete Material Parameters
              f_t(i_Mat) = MATERIAL_PARA(i_Mat,5)
              f_c(i_Mat) = MATERIAL_PARA(i_Mat,7)
              ! Constant
              s_ha= 0.0D6
              SQRT2  = SQRT(TWO)
              SQRT15 = SQRT(15.0D0)
              fcb = 1.2D0*f_c(i_Mat)
              f1  = 1.45D0*f_c(i_Mat)
              f2  = 1.725D0*f_c(i_Mat)
              Kesi_t = f_t(i_Mat)/THR/f_c(i_Mat)
              Kesi_cb = -TWO*fcb/THR/f_c(i_Mat)
              Kesi_1 = -s_ha/f_c(i_Mat)-TWO*f1/THR/f_c(i_Mat)
              Kesi_2 = -s_ha/f_c(i_Mat)-f2/THR/f_c(i_Mat)
              c_Function_Concrete(1:Num_Elem) = -Con_Big_20
              !**************
              !Get a0,a1,a2.
              !**************
              !tem_Matrix1(3,3),tem_F1(3),tem_Sol1(3)
              ! Calculate tem_F1(ansys_thry.pdf<4-340)
              tem_F1(1) = sqrt(f_t(i_Mat)**2+f_t(i_Mat)**2)/
     &                    SQRT15/f_c(i_Mat)
              tem_F1(2) = sqrt(fcb**2+fcb**2)/SQRT15/f_c(i_Mat)
              tem_F1(3) = sqrt(f1**2+f1**2)/SQRT15/f_c(i_Mat)
              ! Calculate tem_Matrix1(ansys_thry.pdf<4-340)
              tem_Matrix1(1,1:3) =[ONE,Kesi_t, Kesi_t**2 ]
              tem_Matrix1(2,1:3) =[ONE,Kesi_cb,Kesi_cb**2]
              tem_Matrix1(3,1:3) =[ONE,Kesi_1, Kesi_1**2 ]
              ! Calculate tem_Sol1, that is a0, a1, a2
              call Matrix_Solve_LSOE(1,0,5,tem_Matrix1,
     &                               tem_F1,tem_Sol1,3)
              a0(i_Mat) = tem_Sol1(1)
              a1(i_Mat) = tem_Sol1(2)
              a2(i_Mat) = tem_Sol1(3)
              !**************
              !Get b0,b1,b2.
              !**************
              ! Calculate tem_F2(ansys_thry.pdf<4-342)
              !tem_F2(1) = sqrt(f_c**2+fc**2)/SQRT15/f_c
              tem_F2(1) = sqrt(f_c(i_Mat)**2+f_c(i_Mat)**2)
     &                  /SQRT15/f_c(i_Mat)
              tem_F2(2) = sqrt(f2**2+f2**2)/SQRT15/f_c(i_Mat)
              tem_F2(3) = ZR
              ! Calculate Kesi_0(ansys_thry.pdf<4-343)
              a = a2(i_Mat)
              b = a1(i_Mat)
              c = a0(i_Mat)
              Root1=(-b+sqrt(b*b-4.0D0*a*c))/(2.0D0*a)
              Root2=(-b-sqrt(b*b-4.0D0*a*c))/(2.0D0*a)
              if(Root1>ZR)then
                  Kesi_0 = Root1
              elseif(Root2>ZR) then
                  Kesi_0 = Root2
              else
                  print *,'    ERROR :: WRONG Kesi_0!'
                  call Warning_Message('S',Keywords_Blank)
              endif
              ! Calculate tem_Matrix2(ansys_thry.pdf<4-342)
              tem_Matrix2(1,1:3) =[ONE,-ONE/THR,ONE/NIN]
              tem_Matrix2(2,1:3) =[ONE,Kesi_2,Kesi_2**2]
              tem_Matrix2(3,1:3) =[ONE,Kesi_0,Kesi_0**2]
              ! Calculate tem_Sol2, that is b0, b1, b2
              call Matrix_Solve_LSOE(1,0,5,
     &                               tem_Matrix2,tem_F2,tem_Sol2,3)
              b0(i_Mat) = tem_Sol2(1)
              b1(i_Mat) = tem_Sol2(2)
              b2(i_Mat) = tem_Sol2(3)   
          enddo
          
          
          !--------------------------------------------------------------
          ! Calculate the principal stress at the center of each element
          !--------------------------------------------------------------
          do i_E=1,Num_Elem   
              S_xx = Stress_El_Center(i_E,1) 
              S_yy = Stress_El_Center(i_E,3) 
              S_zz = Stress_El_Center(i_E,3) 
              S_xy = Stress_El_Center(i_E,4) 
              S_yz = Stress_El_Center(i_E,5) 
              S_xz = Stress_El_Center(i_E,6) 
                     
              call Tool_Principal_Stresses_3D(S_xx,S_yy,S_zz,
     &              S_yz,S_xz,S_xy,
     &              S1_El(i_E),S2_El(i_E),S3_El(i_E),
     &              Vector_S1_El(i_E,1:3),Vector_S2_El(i_E,1:3),
     &              Vector_S3_El(i_E,1:3))
     
              S1 = S1_El(i_E)
              S2 = S2_El(i_E)
              S3 = S3_El(i_E)
              
              c_Mat   = Elem_Mat(i_E)
              
              ! Find the cells where S1 > 0, S2 < 0, and S3 < 0
              if(S1>ZR .and. S2<ZR .and. S3<ZR)then
                  Temp_Value = SQRT((S1-S2)**2+(S2-S3)**2+(S1-S3)**2)
                  x = (S2+S3)/(THR*f_c(c_Mat))
                  p1=a0(c_Mat)+a1(c_Mat)*x+a2(c_Mat)*x**2
                  p2=b0(c_Mat)+b1(c_Mat)*x+b2(c_Mat)*x**2
                  cos_yita = (TWO*S1-S2-S3)/(SQRT2*Temp_Value)
                  FuncF2 = SQRT((S2-S3)**2+S2**2+S3**2)/SQRT15 
                  FuncS2_tem1 = ONE-S1/f_t(c_Mat)
                  FuncS2_tem2 = TWO*p2*(p2**2-p1**2)*cos_yita 
                  FuncS2_tem3 = p2*(TWO*p1-p2)*
     &                          SQRT(FOU*(p2**2-p1**2)*cos_yita**2 +
     &                               FIV*p1**2-FOU*P1*P2)
                  FuncS2_tem4 = FOU*(p2**2-p1**2)*cos_yita**2+
     &                          (p2-TWO*p1)**2
                  FuncS2 = FuncS2_tem1*(FuncS2_tem2+FuncS2_tem3)
     &                     /FuncS2_tem4
                  c_Function_Concrete(i_E) =FuncF2/f_c(c_Mat)-FuncS2
              elseif(S1>ZR .and. S2>ZR .and. S3<ZR)then
              elseif(S1<ZR .and. S2<ZR .and. S3<ZR)then
                  Temp_Value = SQRT((S1-S2)**2+(S2-S3)**2+(S1-S3)**2)
                  FuncF1 = Temp_Value/SQRT15
                  s_h = (S1+S2+S3)/THR
                  kesi = s_h/f_c(c_Mat)
                  r1 = a0(c_Mat)+a1(c_Mat)*kesi+a2(c_Mat)*kesi**2
                  r2 = b0(c_Mat)+b1(c_Mat)*kesi+b2(c_Mat)*kesi**2
                  cos_yita = (TWO*S1-S2-S3)/(SQRT2*Temp_Value)
                  FuncS1_tem2 = TWO*r2*(r2**2-r1**2)*cos_yita 
                  FuncS1_tem3 = r2*(TWO*r1-r2)*
     &                          SQRT(FOU*(r2**2-r1**2)*cos_yita**2 +
     &                               FIV*r1**2-FOU*r1*r2)
                  FuncS1_tem4 = FOU*(r2**2-r1**2)*cos_yita**2+
     &                          (r2-TWO*r1)**2
                  FuncS1 = (FuncS1_tem2+FuncS1_tem3)/FuncS1_tem4
                  c_Function_Concrete(i_E) =FuncF1/f_c(c_Mat)-FuncS1
              endif
          enddo
          
          !--------------------------
          !Find the element of max F
          !--------------------------
          max_F_El = maxloc(c_Function_Concrete,1)
          max_F    = c_Function_Concrete(max_F_El)
          mat_num   = Elem_Mat(max_F_El)
          
          if(max_F >=ZR .and.
     &       MAT_ALLOW_CRACK_Initiation(mat_num)==1)then
              print *,'    >>>>>>>>> A Crack (TYPE T-C-C) Emerged at'
     &               //' Elememt ',max_F_El,' of material type ',
     &              mat_num, '<<<<<<<<<'
              New_Crack_Flag = .True.
          endif
          
          !--------------
          !Add the crack
          !--------------
          !Circle.
          if(Key_Ini_Cr_3D_Type==1 .and.
     &       MAT_ALLOW_CRACK_Initiation(mat_num)==1)then
              if(Num_Initiation_Cracks <
     &           Key_Max_Num_Initiation_Cracks)then
              !num_Crack = num_Crack + 1
              c_num_Crack = num_Crack + 1
              Num_Initiation_Cracks=  Num_Initiation_Cracks +1
              Crack3D_Cir_Coor(c_num_Crack,1:3) =
     &                         Elem_Centroid(max_F_El,1:3)
              Crack3D_Cir_Coor(c_num_Crack,4:6) =
     &                         Vector_S1_El(max_F_El,1:3)
              Crack3D_Cir_Coor(c_num_Crack,7) =Size_Ini_Crack/TWO      
              endif
          !Rectangle.
          elseif(Key_Ini_Cr_3D_Type==2)then
          endif
          
      endif  
              
      RETURN
      END SUBROUTINE Check_Crack_initialize_3D      