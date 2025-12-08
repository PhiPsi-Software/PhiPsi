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
 
subroutine Check_Crack_initialize_2D(ifra,iter,New_Crack_Flag)
! Check for the formation of new cracks (2D problem).
! 2024-06-22.      
!
!***************
! Public Module
!***************
use Global_Float_Type
use Global_Crack_Common
use Global_DISP
use Global_Elem_Area_Vol
use Global_Model
use Global_HF
use Global_Common
use Global_Material
use Global_Dynamic 
use Global_Crack
use Global_Plasticity

implicit none
integer,intent(in)::ifra,iter
logical,intent(out)::New_Crack_Flag
real(kind=FT) S1_El(num_Elem),S2_El(num_Elem),S3_El(num_Elem) 
real(kind=FT) Shear_El(num_Elem)
real(kind=FT) S_1,S_2,S_3
integer max_S1_El,mat_num,c_num_Crack,max_Shear_El
real(kind=FT) St_Critical,max_S1,max_Shear
real(kind=FT) Shear_Strength
real(kind=FT) Sxx_Center(num_Elem), Syy_Center(num_Elem), Sxy_Center(num_Elem)
integer i_E
real(kind=FT) c_D(3,3),U(8)
real(kind=FT) c_X_NODES(4),c_Y_NODES(4),c_kesi,c_yita,c_Stress(3)
integer c_NN(4)
real(kind=FT) c_v
real(kind=FT) theta_stress(num_Elem)
real(kind=FT) crack_angle,crack_angle_360,crack_center(2),crack_half_length
real(kind=FT) ele_x,ele_y
      
!******************
! Data Preparation
!******************
New_Crack_Flag = .False.

!*********************
! Interactive display
!*********************
print *,'    Checking initialization of crack......'  

!********************************************************
!                                                      *
! Calculate the stress at the center of each element   *
!                                                      *
!********************************************************
!call Get_Element_Center_Stress_FEM_3D(.True.,iter,DISP,Stress_El_Center)
do i_E=1,Num_Elem 
  !//////////////////////////////////////////////////////
  ! 2D Nonlinear Analysis PhiPsi2D_Static_NonLinear.
  ! Obtain stress data from Gauss point state variables.
  !//////////////////////////////////////////////////////
  if(Key_Analysis_Type==7) then
      c_Stress(1) = sum(STATEV(i_E,1:Num_Gauss_P_FEM,2))/Num_Gauss_P_FEM
      c_Stress(2) = sum(STATEV(i_E,1:Num_Gauss_P_FEM,3))/Num_Gauss_P_FEM
      c_Stress(3) = sum(STATEV(i_E,1:Num_Gauss_P_FEM,4))/Num_Gauss_P_FEM
  !///////////////////
  ! General analysis.
  !///////////////////
  else
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
    U =  [DISP(c_NN(1)*2-1),DISP(c_NN(1)*2),DISP(c_NN(2)*2-1),&
            DISP(c_NN(2)*2),DISP(c_NN(3)*2-1),&
            DISP(c_NN(3)*2),DISP(c_NN(4)*2-1),DISP(c_NN(4)*2)]
    ! Stress at the center of the element.
    if (Yes_XFEM .eqv. .True.) then
          call Cal_Any_Point_Stress_KesiYita(i_E,ZR,ZR,1,DISP,Sxx_Center(i_E), Syy_Center(i_E), Sxy_Center(i_E))
    elseif (Yes_XFEM .eqv. .False.) then
          call Cal_Ele_Stress_N4(i_E,1,c_X_NODES,c_Y_NODES,c_D,ZR,ZR,U,c_Stress) 
    endif
  
  endif
  
  Sxx_Center(i_E) = c_Stress(1)
  Syy_Center(i_E) = c_Stress(2)
  Sxy_Center(i_E) = c_Stress(3)

  ! If an initial crack generation region is defined.
  if(Key_Ini_Crack_Zone==1)then
      ele_x = Elem_Centroid(i_E,1)
      ele_y = Elem_Centroid(i_E,2)
      ! If the element center is outside the rupture zone, the stress is recorded as 0.
      if((ele_x < Ini_Crack_Zone_MinX) .or. (ele_x > Ini_Crack_Zone_MaxX) .or. &
         (ele_y < Ini_Crack_Zone_MinY) .or. (ele_y > Ini_Crack_Zone_MaxY))then
          Sxx_Center(i_E) = ZR; Syy_Center(i_E) = ZR; Sxy_Center(i_E) = ZR
      endif
  endif  
enddo
  
!**********************************************************
!                                                        *
! If it is the Maximum Tensile Stress Criterion (MTSC).  *
!                                                        *
!**********************************************************
if(Key_Ini_Rule==1)then  
  do i_E=1,Num_Elem   
      call Tool_Principal_Stresses_2D(Sxx_Center(i_E),Syy_Center(i_E),Sxy_Center(i_E),S1_El(i_E),S3_El(i_E),theta_stress(i_E))
      Shear_El(i_E) = (S1_El(i_E)-S3_El(i_E))/TWO
  enddo
  !-------------------------------------
  !Find the element of max Shear stress
  !-------------------------------------
  max_S1_El = maxloc(S1_El,1)
  max_S1    = S1_El(max_S1_El)
  mat_num   = Elem_Mat(max_S1_El)
  St_Critical  = Material_Para(mat_num,5)
  if(max_S1>=St_Critical .and.MAT_ALLOW_CRACK_Initiation(mat_num)==1)then   
      print *,'    >>> A Crack Emerged at Elememt ',max_S1_El,' of material type ',mat_num, '<<<' 
      New_Crack_Flag = .True.
  endif
  
  !--------------
  !Add the crack
  !--------------
  ! 2D line segment.
  if(MAT_ALLOW_CRACK_Initiation(mat_num)==1)then
      if(Num_Initiation_Cracks <Key_Max_Num_Initiation_Cracks)then
          num_Crack = num_Crack + 1
          Num_Initiation_Cracks=  Num_Initiation_Cracks +1
          crack_angle = theta_stress(max_S1_El)
          crack_angle_360 = crack_angle*180.0D0/pi
          crack_half_length = Size_Ini_Crack/TWO 
          crack_center(1:2)  = Elem_Centroid(max_S1_El,1:2) 
          ! Fine-tune the position of the crack center to prevent it from being located at the center of the
          ! element.
          crack_center(1) =  crack_center(1) - Ave_Elem_L/3.0D0
          crack_center(2) =  crack_center(2) + Ave_Elem_L/5.0D0
          ! The normal vector of the plane with the maximum shear stress: obtained by rotating the direction
          ! of the maximum principal stress 45 degrees around the direction of the intermediate principal
          ! stress.
          Crack_Coor(num_Crack,1,1) = crack_center(1) + crack_half_length*cos(crack_angle)
          Crack_Coor(num_Crack,1,2) = crack_center(2) + crack_half_length*sin(crack_angle)
          Crack_Coor(num_Crack,2,1) = crack_center(1) - crack_half_length*cos(crack_angle)
          Crack_Coor(num_Crack,2,2) = crack_center(2) - crack_half_length*sin(crack_angle)
          Each_Cr_Poi_Num(num_Crack)  = 2
          Cracks_Allow_Propa(num_Crack) = 1
      endif
  endif
endif

!*************************************************
!                                               *
! If it is the maximum shear stress criterion.  *
!                                               *
!*************************************************
if(Key_Ini_Rule==2)then  
!  !-------------------------------
! !Calculate the principal stress at the center of each element
!  !-------------------------------
  do i_E=1,Num_Elem   
      call Tool_Principal_Stresses_2D(Sxx_Center(i_E),Syy_Center(i_E),Sxy_Center(i_E),S1_El(i_E),S3_El(i_E),theta_stress(i_E))
      Shear_El(i_E) = (S1_El(i_E)-S3_El(i_E))/TWO
  enddo
  !-------------------------------------
  !Find the element of max Shear stress
  !-------------------------------------
  max_Shear_El = maxloc(Shear_El,1)
  max_Shear    = Shear_El(max_Shear_El)
  mat_num      = Elem_Mat(max_Shear_El)
  ! St_Critical = Material_Para(mat_num,5)          !Tensile Strength
  Shear_Strength = Material_Para_Added(mat_num,12)
  ! The shear strength is defined as the 12th parameter in MATERIAL_PARA_ADDED_1.
  if(max_Shear>=Shear_Strength .and. MAT_ALLOW_CRACK_Initiation(mat_num)==1)then   
      print *,'    >>> A Crack Emerged at Elememt ',max_Shear_El,' of material type ',mat_num, '<<<' 
      New_Crack_Flag = .True.
  endif
  
  !--------------
  !Add the crack
  !--------------
  ! 2D line segment.
  if(MAT_ALLOW_CRACK_Initiation(mat_num)==1)then
      if(Num_Initiation_Cracks <Key_Max_Num_Initiation_Cracks)then
          num_Crack = num_Crack + 1
          Num_Initiation_Cracks=  Num_Initiation_Cracks +1
          crack_angle = theta_stress(max_Shear_El) + pi/FOU
          crack_angle_360 = crack_angle*180.0D0/pi
          crack_half_length = Size_Ini_Crack/TWO 
          crack_center(1:2)  = Elem_Centroid(max_Shear_El,1:2) 
          ! Fine-tune the position of the crack center to prevent it from being located at the center of the
          ! element.
          crack_center(1) =  crack_center(1) - Ave_Elem_L/3.0D0
          crack_center(2) =  crack_center(2) + Ave_Elem_L/5.0D0
          !crack_angle = 45.0D0*pi/180.0D0
          ! The normal vector of the plane with the maximum shear stress: obtained by rotating the direction
          ! of the maximum principal stress 45 degrees around the direction of the intermediate principal
          ! stress.
          Crack_Coor(num_Crack,1,1) = crack_center(1) + crack_half_length*cos(crack_angle)
          Crack_Coor(num_Crack,1,2) = crack_center(2) + crack_half_length*sin(crack_angle)
          Crack_Coor(num_Crack,2,1) = crack_center(1) - crack_half_length*cos(crack_angle)
          Crack_Coor(num_Crack,2,2) = crack_center(2) - crack_half_length*sin(crack_angle)
          Each_Cr_Poi_Num(num_Crack)  = 2
          Cracks_Allow_Propa(num_Crack) = 1
      endif
  endif
  
!*********************************************
!                                           *    
!           Concrete rule.                  *
! Refe: ANSYS theory manual-4.11            *
! The following is three-dimensional code.  *
!*********************************************
elseif(Key_Ini_Rule==3)then  

  
endif  

      
RETURN
END SUBROUTINE Check_Crack_initialize_2D      