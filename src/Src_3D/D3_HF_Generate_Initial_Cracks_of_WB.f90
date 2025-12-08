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
 
SUBROUTINE D3_HF_Generate_Initial_Cracks_of_WB(i_Prop,i_WB,i_Stage)
! Determine the initial fractures of each stage of the segmented fracturing based on the wellbore
! setup parameters.
! At the same time, this subroutine also updates the status of the previous cracks in each section.
! NEWFTU2022053101.
! 2022-05-31. 

! The global variable Crack_Type_Status_3D(i_C,10) is used to indicate the type and status of
! cracks.
! Column 1 (Fracture Type): =1, HF fracture; =2, natural fracture; =3, post-fracturing hydraulic
! fracture
! Note: Natural fractures and propped hydraulic fractures may potentially turn into HF fractures.
! Column 2 (Fracture Status): =1, HF fracturing not completed; =2, HF fracturing completed
! Column 3 (Can the crack continue to propagate): =1, yes; =0, no
! Column 4 (Whether the fracture has obtained a fluid node): =1, Yes; =0, No
! Column 5 (Did the crack propagate in the previous step?): =1, Yes; =0, No
! The global variable Cracks_Stages_Wellbores(i_WB, i_Stage, i_C) is used for the crack numbers
! corresponding to each segment of each wellbore.

!#############################
! Read public variable module
!#############################
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack_Common
use Global_Crack_3D
use Global_HF
use Global_Filename
use Global_Stress,only:InSitu_Strs_Gaus_xx,&
                       InSitu_Strs_Gaus_yy,&
                       InSitu_Strs_Gaus_zz,&
                       InSitu_Strs_Gaus_yz,&
                       InSitu_Strs_Gaus_xz,&
                       InSitu_Strs_Gaus_xy

!###########################
! Variable Type Declaration
!###########################
implicit none
integer,intent(in)::i_Prop,i_WB,i_Stage
integer num_Stages,num_Cr_Stage 
real(kind=FT) Start_Point(3),End_Point(3)
real(kind=FT) c_Crack_Center(3),delta_L,L_Stage
integer i_C 
real(kind=FT) Tool_Function_2Point_Dis_3D
real(kind=FT) L_WB_Fracturing,L_Crack_Gap
real(kind=FT) Crack_L,n_Vector(3)
integer Stage_Last_WB,num_Crs_Last_Stage,c_Cracks(Max_Clust)
character(200) Filename_1
real(kind=FT) c_Point(3)
integer c_i_WB,c_i_Stage
integer Center_Elem_num
integer Ele_Num_Cache 
real(kind=FT) S_1,S_2,S_3
real(kind=FT) Vector_S1(3),Vector_S2(3),Vector_S3(3)


!#########################################################
! Generate only at the first expansion step: if it is not 
! the first expansion step, exit immediately
!#########################################################
if (i_Prop>=2)then
  return
endif

!##################
! Data Preparation
!##################
Ele_Num_Cache = 1
! Current number of wellbore sections
num_Stages = num_Stages_Wellbores(i_WB)
! Number of fractures in the current section of the current wellbore
num_Cr_Stage = num_Crs_Stages_Wellbores(i_WB,i_Stage)
! Starting coordinates of staged fracturing in the wellbore (used for automatically generating
! initial fractures)
Start_Point = Wellbores_Start_Point(i_WB,1:3)   
! Endpoint coordinates of staged fracturing in the wellbore (used for automatically generating
! initial fractures)
End_Point   = Wellbores_End_Point(i_WB,1:3)    
! Length of the wellbore fracturing section
L_WB_Fracturing=Tool_Function_2Point_Dis_3D(Start_Point,End_Point)
! The length of each stage
L_Stage = L_WB_Fracturing/num_Stages
! Current wellbore current section fracture spacing
L_Crack_Gap = L_Stage/(num_Cr_Stage+1)
! Rectangular crack size
Crack_L = Size_Ini_Crack_Wellbores
!////////////////
! Normal vector.
!////////////////
! The fracture face is perpendicular to the wellbore.
if (Key_Gen_Ini_Crack_Rule_Wellbores==1)then
  n_Vector = End_Point - Start_Point
  call Vector_Normalize(3,n_Vector)
endif
! User-defined.
if (Key_Gen_Ini_Crack_Rule_Wellbores==3)then
  n_Vector = Normal_Vector_Gen_Ini_Crack_Wellbores
  call Vector_Normalize(3,n_Vector)
endif

!#######################
! Obtain initial cracks
!#######################
do i_C = 1,num_Cr_Stage
  delta_L = (i_Stage-1)*L_Stage + i_C*L_Crack_Gap
  !**************************************************
  ! Obtain crack centroid coordinates c_Crack_Center
  !**************************************************
  call Tool_Offset_Point_A_to_Point_B_3D(Start_Point,End_Point,delta_L,c_Crack_Center)
  
  !**************************************************************
  ! Assuming the crack orientation criterion is perpendicular to
  ! the direction of the minimum principal stress.
  ! NEWFTU2024121401.
  !**************************************************************
  if (Key_Gen_Ini_Crack_Rule_Wellbores==2)then 
      !//////////////////////////////////////
      ! If initial stress is not considered.
      !//////////////////////////////////////
      if(Key_InSitu_Strategy ==0)then
          ! The cracks are still vertical to the wellbore direction.
          n_Vector = End_Point - Start_Point
          call Vector_Normalize(3,n_Vector)
      !//////////////////////////////////
      ! If initial stress is considered.
      !//////////////////////////////////
      elseif(Key_InSitu_Strategy /=0)then
          ! Obtain the element number based on the crack center coordinates.
          call Cal_Ele_Num_by_Coors_3D(c_Crack_Center(1),c_Crack_Center(2),c_Crack_Center(3),Ele_Num_Cache,Center_Elem_num)
          if (Center_Elem_num ==0) then
            print *,'    Error :: failed to get element number in D3_HF_Generate_Initial_Cracks_of_WB.f90'
            call Warning_Message('S',Keywords_Blank)
          endif
          ! Calculate the principal stress using the stress at the first Gauss point of the element.
          ! Note that the stress here is positive in compression.
          call Tool_Principal_Stresses_3D(InSitu_Strs_Gaus_xx(Center_Elem_num,1),&
                                     InSitu_Strs_Gaus_yy(Center_Elem_num,1),&
                                     InSitu_Strs_Gaus_zz(Center_Elem_num,1),&
                                     InSitu_Strs_Gaus_yz(Center_Elem_num,1),&
                                     InSitu_Strs_Gaus_xz(Center_Elem_num,1),&
                                     InSitu_Strs_Gaus_xy(Center_Elem_num,1),&
                                     S_1,S_2,S_3,&
                                     Vector_S1,Vector_S2,Vector_S3)
          n_Vector = Vector_S3
          call Vector_Normalize(3,n_Vector)
      endif
  endif
  
  !*************************
  ! Generate initial cracks
  !*************************
  num_Crack = num_Crack + 1
  
  !//////////////////////////////////////////////////////////////////////////////////////
  ! Obtain 3D rectangular HF fractures based on the normal vector and center coordinates
  !//////////////////////////////////////////////////////////////////////////////////////
  if(Key_Gen_Ini_Crack_Wellbores==1)then
      call Tool_Get_3D_Rect_by_n_Vector_and_Center(n_Vector(1:3),c_Crack_Center(1:3),Crack_L,Crack3D_Coor(num_Crack,1:4,1:3))
      Each_Cr_Poi_Num(num_Crack)   = 4
  
  !///////////////////
  ! Circular HF crack
  !///////////////////
  !NEWFTU2022061702.
  elseif(Key_Gen_Ini_Crack_Wellbores==2)then
      ! Update circular crack coordinates
      Crack3D_Cir_Coor(num_Crack,1:3) = c_Crack_Center(1:3)
      Crack3D_Cir_Coor(num_Crack,4:6) = n_Vector(1:3)
      Crack3D_Cir_Coor(num_Crack,7)   = Crack_L/TWO
  
  !//////////////////////
  ! Polygonal HF cracks.
  !//////////////////////
  !NEWFTU2022061702.    
  elseif(Key_Gen_Ini_Crack_Wellbores==3)then
      call Tool_Get_3D_Polygon_by_n_Vector_and_Center(&
             n_Vector(1:3),Num_Poly_Edges_PolyCr_WB,&
             c_Crack_Center(1:3),Crack_L,&
             Crack3D_Coor(num_Crack,1:Num_Poly_Edges_PolyCr_WB,1:3))    
      Each_Cr_Poi_Num(num_Crack)   = Num_Poly_Edges_PolyCr_WB
  endif
  !****************************************
  ! Related variables updated, 2022-06-02.
  !****************************************
  ! The fracture numbers corresponding to each section of each wellbore
  Cracks_Stages_Wellbores(i_WB,i_Stage,i_C) = num_Crack
  
  ! Update parameters such as crack status
  Crack_Type_Status_3D(num_Crack,1) = 1
  Crack_Type_Status_3D(num_Crack,2) = 1
  Crack_Type_Status_3D(num_Crack,3) = 1
enddo

!##############################################################
! Update the status of the previous crack section. 2022-06-02.
!##############################################################
! If it is the first section of the first shaft, no operation is needed.
if (i_Stage==1 .and. i_WB==1) then
  return
endif

! If it is the first section of a subsequent wellbore, then:
if (i_Stage==1 .and. i_WB>=2) then
  ! Get the number of sections of the previous shaft
  Stage_Last_WB = num_Stages_Wellbores(i_WB-1)
  ! Obtain the number of fractures in the last section of the previous wellbore
  num_Crs_Last_Stage = num_Crs_Stages_Wellbores(i_WB-1,Stage_Last_WB)
  ! Obtain the fracture number of the last section of the previous wellbore
  c_Cracks(1:num_Crs_Last_Stage)  =  Cracks_Stages_Wellbores(i_WB-1,Stage_Last_WB,1:num_Crs_Last_Stage)

  ! Update the active HF fracture status of the previous crack segment
  do i_C =1,num_Crs_Last_Stage
      Crack_Type_Status_3D(c_Cracks(i_C),1) = 3
      Crack_Type_Status_3D(c_Cracks(i_C),2) = 2
      Crack_Type_Status_3D(c_Cracks(i_C),3) = 0
  enddo
  
  ! Update the naturally occurring fractures communicated by HF to post-fracture cracks.
  ! IMPROV2022061101.
  if(Key_Random_NaCr>=1) then      
      !BUGFIX2023021902.
      do i_C = 1,num_Crack 
          ! If it is a natural fracture communicated by HF
          if (Crack_Type_Status_3D(i_C,6)/=0) then
              Crack_Type_Status_3D(i_C,1) = 3
              Crack_Type_Status_3D(i_C,2) = 2
              Crack_Type_Status_3D(i_C,3) = 0
          endif
      enddo
  endif
  

  
  return
endif

! If it is a subsequent shaft section, not the first one, then:
if (i_Stage>=2) then
  ! Get the number of cracks in the previous segment
  num_Crs_Last_Stage = num_Crs_Stages_Wellbores(i_WB,i_Stage-1)
  ! Get the number of cracks from the previous section
  c_Cracks(1:num_Crs_Last_Stage)  =  Cracks_Stages_Wellbores(i_WB,i_Stage-1,1:num_Crs_Last_Stage)
  ! Update the active HF fracture status of the previous crack segment
  do i_C =1,num_Crs_Last_Stage
      Crack_Type_Status_3D(c_Cracks(i_C),1) = 3
      Crack_Type_Status_3D(c_Cracks(i_C),2) = 2
      Crack_Type_Status_3D(c_Cracks(i_C),3) = 0
  enddo
  ! Update the naturally occurring fractures communicated by HF to post-fracture cracks.
  ! IMPROV2022061101.
  if(Key_Random_NaCr>=1) then      
      !BUGFIX2023021902.
      do i_C = 1,num_Crack  
          if (Crack_Type_Status_3D(i_C,6)/=0) then
              Crack_Type_Status_3D(i_C,1) = 3
              Crack_Type_Status_3D(i_C,2) = 2
              Crack_Type_Status_3D(i_C,3) = 0
          endif
      enddo              
  endif          
  return
endif    

RETURN
END SUBROUTINE D3_HF_Generate_Initial_Cracks_of_WB
