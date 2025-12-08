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
 
SUBROUTINE Input_Check_Display
! Check and display input information

!-----------------------------
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Filename
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack
use Global_Crack_Common
use Global_Crack_3D
use Global_HF
use Global_POST
use Global_Field_Problem
use Global_Inclusion
use Global_Dynamic
use Global_MD
use Global_Cohesive
use Module_Tool_Logger
use Global_Surface_Load
use Global_3D_HF_Experiment 
use Global_Material

!---------------------------
! Variable Type Declaration
!---------------------------
implicit none
integer i_C,i_S,i_P,i_Node,i_NC,i_Seg
character*2 tem_char
real(kind=FT) crack_p1(2),crack_p2(2)
logical Yes_Seg_P1_in,Yes_Seg_P2_in
integer EndByEnd_Outline(size(Outline,1)+1)
integer c_Ele
logical Yes_In,Yes_on_Ele_Edge
real(kind=FT) Modifi_Factor,Check_Point(2),c_Node_x,c_Node_y
logical Yes_Node_OnSeg,Yes_Inj_On_Cr
real(kind=FT) Inj_Poin_x,Inj_Poin_y,Cr_A(2),Cr_B(2)
real(kind=FT) Inj_P_x,Inj_P_y
integer i,c_C,i_H,iii
logical Yes_InjP_On_Cr,Yes_Point_in
real(kind=FT) Hole_point(4,2),x0,y0,R0,x1,y1,R1
logical Yes_Intersect
integer j_H
character(len=8) :: i_char
integer i_Incl,j_Incl,c_num_Edges
real(kind=FT) Incl_point(4,2)
integer i_Ac_node
integer i_CP_node
real(kind=FT) Tool_Function_2Point_Dis
integer c_num_Inter,c_State
real(kind=FT) c_Inter(2,2)
real(kind=FT) c_Arc_End_1(2),c_Arc_End_2(2)
real(kind=FT) c_Arc_r,c_Arc_Radian_S,c_Arc_Radian_E,c_Arc_Radian
real(kind=FT) c_Arc_Center(2)
logical c_Yes_feasible
integer num_Centroid
integer i_E
real(kind=FT) c_Centroid(2),c_Arc_Direction,a_Hole,b_Hole
character(200) temp_log
integer max_material_type
integer j
character(200) Filename_1
integer c_num_points
logical Yes_on
logical Yes_Dup,Yes_Cross
integer num_Dup
logical c_Yes_on
logical c_Yes_Cr_In_Model
real(kind=FT) c_P_x,c_P_y,c_P_z
integer i_Mat
real(kind=FT) Ran_Num
integer i_wb,i_Point_wb
character(5) s_Key_Num_Process,s_num_of_Material
character(14) s_Min_X_Coor,s_Max_X_Coor,s_Min_Y_Coor,s_Max_Y_Coor,s_Min_Z_Coor,s_Max_Z_Coor
real(kind=FT) Dis_Tol
external Warning_Message

5001 FORMAT('     !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!')

!-----------------------------
! Check the input information
!-----------------------------
print *
print *, '    ----------------------------------------------'
print *, '        <  Input    Information    Summary  >'
print *, '    ----------------------------------------------'

!-------------------------------
! Display all input information
!-------------------------------
if(len_trim(Name_of_Keywords)>=1) then
    print *, '    Keywords file:             ',trim(Name_of_Keywords)
endif
print *, '    Filename:                  ',trim(Filename)
print *, '    Working directory:         ',trim(Work_Directory)
!log file
if(len_trim(Name_of_Keywords)>=1) then
    temp_log = trim('Keywords file:             '//Name_of_Keywords)
    call log_msg(temp_log)
endif
temp_log = trim('Filename:                  '//Filename)
call log_msg(temp_log)
temp_log = trim('Working directory:         '//Work_Directory)
call log_msg(temp_log)

!************
! float type
!************
#ifndef Silverfrost
if (FT == 4)then
  print *, '    Float type:                single precision'
elseif (FT == 8)then
  print *, '    Float type:                double precision'
else
  print *,'    Wrong float type (should be 4 or 8)!'
  call Warning_Message('S',Keywords_Blank)
end if
#endif

#ifdef Silverfrost
if (FT == 1)then
  print *, '    Float type:                single precision'
elseif (FT == 2)then
  print *, '    Float type:                double precision'
else
  print *,'    Wrong float type (should be 1 or 2)!'
  call Warning_Message('S',Keywords_Blank)
end if
#endif

!**************************
! Analysis Type (2D or 3D)
!**************************
if (Key_Dimension == 2)then
    print *, '    Dimension of the problem:  2D'
    call log_msg('Dimension of the problem:  2D')
    ! 2D problem type
    if (Key_Type_2D == 1)then
        print *, '    2D problem type:           plane stress'
        call log_msg('2D problem type:           plane stress')
    elseif (Key_Type_2D == 2)then
        print *, '    2D problem type:           plane strain'
        call log_msg('2D problem type:           plane strain')
    else
        call Warning_Message('E','Key_Type_2D')
    end if
elseif (Key_Dimension == 3)then
    print *, '    Dimension of the problem:  3D'
    call log_msg('Dimension of the problem:  3D')
else
    call Warning_Message('E','Key_Dimension')
end if

!**************************************
! Is it a solid model (3D), 2022-04-14
!**************************************
!1020 FORMAT(5X,50('-'))
  if (Key_Block_Model == 1)then
      print *, '    3D Block model:            yes'
      call log_msg('3D Block model:            yes')
      print *, '    Attention :: 3D Block model is assigned!'
      ! Only works in 3D problems
      if (Key_Dimension == 2)then
          print *, '    Error :: *Key_Block_Model = 1 for 3D problem only!'
          call Warning_Message('S',Keywords_Blank)
      endif
  elseif (Key_Block_Model == 0)then
      print *, '    3D Block model:            no'
      call log_msg('3D Block model:            no')
  else
      call Warning_Message('E','Key_Block_Modeln')
  end if

!*****************************
! Big Deformation, 2021-07-24
!*****************************
  if (Key_Large_Deform == 0) then
      print *, '    Large strain formulation:  no'
      call log_msg('Large strain formulation:  no')
  elseif (Key_Large_Deform == 1)then
      print *, '    Large strain formulation:  yes'
      call log_msg('Large strain formulation:  yes')
      ! Large deformation only applies to Key_Analysis_Type=8, HYPLAS
      if (Key_Analysis_Type /= 8)then
          print *, '    Error :: Large strain formulation is available only when *Key_Analysis_Type = 8!'
          call Warning_Message('S',Keywords_Blank)
      end if
      print *, '    Attention :: Large deformation is active!'
  else
      call Warning_Message('E','Key_Dimension')
  end if

  !*************
  ! Unit system
  !*************
  ! The No. 5 iterator for hydraulic fracturing analysis and pure hydrostatic analysis supports the
  ! mm-ton-s unit system.
  if(Key_Unit_System==1)then
          print *, '    Unit system:               SI      '
  elseif(Key_Unit_System==2)then
      if(Key_Analysis_Type == 3)then
          print *, '    Unit system:               mm-ton-s'
      elseif(Key_Analysis_Type == 4)then
          print *, '    Unit system:               mm-ton-s'
      else
          print *, '    Error :: unit system mm-ton-s only for solver 5 for HF analysis and uniform pressure analysis!'
          call Warning_Message('S',Keywords_Blank)
      endif
  else
      call Warning_Message('E','Key_Unit_System')
  endif

  !**************************************************************************************************
  ! Problem types: static analysis, dynamic analysis, hydraulic fracturing analysis, test procedures
  ! applying pure water pressure, etc.
  !**************************************************************************************************
  if (Key_Analysis_Type == 1) then
      print *, '    Analysis type:             static analysis'
      call log_msg('Analysis type:             static analysis')
  elseif (Key_Analysis_Type == 2) then
      print *, '    Analysis type:             implicit dynamic analysis'
      call log_msg('Analysis type:             implicit dynamic')
  elseif (Key_Analysis_Type == 3) then
      print *, '    Analysis type:             hydraulic fracturing analisys'
      call log_msg('Analysis type:             hydraulic fracturing analisys')
  elseif (Key_Analysis_Type == 4 .and. Key_Dimension==2) then
      print *, '    Analysis type:             test program for uniform water pressure'
      call log_msg('Analysis type:             test program for uniform water pressure')
  elseif (Key_Analysis_Type == 4 .and. Key_Dimension==3) then
      print *, '    Analysis type:             HF with slipwater'
      call log_msg('Analysis type:             HF with slipwater')
      ! Local refinement with Key_Local_Mesh_Refine=1 is not supported
      if(Key_Local_Mesh_Refine==1)then
          print *, '    Error :: *Key_Local_Mesh_Refine=1 is not available for 3D *Key_Analysis_Type=4!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Partial encryption with key_front_segmentation=1 is not supported
      if(key_front_segmentation==1)then
          print *, '    Error :: *key_front_segmentation=1 is not available for 3D *Key_Analysis_Type=4!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Cylindrical coordinate system not supported
      if(Key_Post_CS_N_Strs_Cylindr==1)then
          print *, '    Error :: *Key_Post_CS_N_Strs_Cylindr=1 is not available for 3D *Key_Analysis_Type=4!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Only supports EBE-PCG solver No. 11
!      if(Key_SLOE/=11)then
!          call Warning_Message('S',Keywords_Blank)
!      endif
      !NEWFTU-2025120201.
      if(Key_SLOE/=11)then
          if(Key_K_Sparse==0)then  
              print *, '    Error :: 3D *Key_Analysis_Type =4 supports only sparse K! Please set *Key_K_Sparse=1'
              call Warning_Message('S',Keywords_Blank)
          endif
      endif
      
      ! Internal water pressure in the seam is not supported when Key_Crack_Inner_Pressure=1
      if(Key_Crack_Inner_Pressure==1)then
          print *, '    Error :: *Key_Crack_Inner_Pressure=1 is not available for 3D *Key_Analysis_Type=4!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! It must be staged fracturing.
      if(Key_HF_Multistage_3D==0)then
          print *, '    Error :: *Key_HF_Multistage_3D=0 is not available for 3D *Key_Analysis_Type=4!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Touch-based calculation is not supported because the computation is too intensive. 2022-06-29
      !if(Key_Contact/=0)then
      !    call Warning_Message('S',Keywords_Blank)
      !endif
      ! Only supports crack propagation criteria with CFCP=3. 2022-06-30.
      if(CFCP/=3)then
          print *, '    Error :: *CFCP/=3 is not available for 3D *Key_Analysis_Type=4!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Only supports crack propagation criteria with Key_CFCP_3_Type=2. 2022-07-02. NEWFTU2022070201.
      if(Key_CFCP_3_Type/=2)then
          print *, '    Error :: *Key_CFCP_3_Type/=2 is not available for 3D *Key_Analysis_Type=4!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Only supports the crack extension criterion for Key_3D_Cr_Update_Scheme=2. 2022-07-12.
      ! NEWFTU2022071201.
      if(Key_3D_Cr_Update_Scheme/=2)then
          print *, '    Error :: *Key_3D_Cr_Update_Scheme/=2 is not available for 3D *Key_Analysis_Type=4!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Key_Adj_Ini_Crack_3D=1 is not supported. July 8, 2022.
      if(Key_Adj_Ini_Crack_3D==1)then
          print *, '    Error :: *Key_Adj_Ini_Crack_3D=1 is not available for 3D *Key_Analysis_Type=4!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Time-stepping iteration method.
      if(Key_3D_HF_Time_Step_Method==1) then
          print *, '    Time-step method:          Newton-raphson method'
      elseif(Key_3D_HF_Time_Step_Method==2)then
          print *, '    Time-step method:          Bisection method'
      else
          call Warning_Message('E','Key_3D_HF_Time_Step_Method')
      endif
      ! 3D Clean Water Fracturing Crack Tip fk Function Calculation Type. 2023-08-08.
      if(Key_3D_HF_SlipWater_fk_Type==1) then
          print *, '    fk type         :          Max_KI_eq_3D'
      elseif(Key_3D_HF_SlipWater_fk_Type==2) then
          print *, '    fk type         :          Ave_KI_eq_3D'
      elseif(Key_3D_HF_SlipWater_fk_Type==3) then
          print *, '    fk type         :          Min_KI_eq_3D'
      else
          call Warning_Message('E','Key_3D_HF_SlipWater_fk_Type')
      endif
  elseif (Key_Analysis_Type == 5) then
      print *, '    Analysis type:             viscosity-dominated theoretical water pressure'
      call log_msg('Analysis type:             viscosity-dominated theoretical water pressure')
  elseif (Key_Analysis_Type == 6) then
      print *, '    Analysis type:             explicit dynamic analysis'
      call log_msg('Analysis type:             explicit dynamic analysis')
  elseif (Key_Analysis_Type == 7) then
      print *, '    Analysis type:             non-linear analysis'
      call log_msg('Analysis type:             non-linear analysis')
  elseif (Key_Analysis_Type == 8) then
!      call log_msg('Analysis type:             static HYPLAS')
      print *, '    Error :: HYPLAS is no longer supported!'
      call Warning_Message('S',Keywords_Blank)
  elseif (Key_Analysis_Type == 11) then
      print *, '    Analysis type:             contact analysis of rigid balls'
      call log_msg('Analysis type:             contact analysis of rigid balls')
      ! Since this analysis has no input data, no check is needed and it can proceed directly.
      goto 9999
  elseif (Key_Analysis_Type == 12) then
      print *, '    Analysis type:             Impact explicit dynamic analysis'
      call log_msg('Analysis type:             Impact explicit dynamic analysis')
      ! Since this analysis has no input data, no check is needed and it can proceed directly.
      goto 9999
  elseif (Key_Analysis_Type == 15) then
      print *, '    Analysis type:             static field problem'
      call log_msg('Analysis type:             static field problem')
  elseif (Key_Analysis_Type == 16) then
      print *, '    Analysis type:             transient field problem'
      call log_msg('Analysis type:             transient field problem')
  elseif (Key_Analysis_Type == 17) then
      print *, '    Analysis type:             transient field problem (gas production coupled with deformation)'
      call log_msg('Analysis type:             transient field problem (gas production coupled with deformation)')
  elseif (Key_Analysis_Type == 21) then
      print *, '    Analysis type:             molecular dynamics (OpenMP version)'
      call log_msg('Analysis type:             molecular dynamics (OpenMP version)')
      ! Since this analysis has no input data, no check is needed and it can proceed directly.
      goto 9999
  elseif (Key_Analysis_Type == 31) then
      print *, '    Analysis type:             peridynamic'
      call log_msg('Analysis type:             peridynamic')
      ! Since this analysis has no input data, no check is needed and it can proceed directly.
      goto 9999
  elseif (Key_Analysis_Type == 41) then
      print *, '    Analysis type:             CFD'
      call log_msg('Analysis type:             CFD')
      ! Since this analysis has no input data, no check is needed and it can proceed directly.
      goto 9999
  elseif (Key_Analysis_Type == 51) then
      print *, '    Analysis type:             PFEM Structure'
      call log_msg('Analysis type:             PFEM Structure')
      ! Since there is currently no input data for this analysis, no checks are needed, and you can
      ! proceed directly.
      goto 9999
  elseif (Key_Analysis_Type == 52) then
      print *, '    Analysis type:             PFEM with Umat'
      call log_msg('Analysis type:             PFEM with Umat')
      ! Since there is currently no input data for this analysis, no checks are needed, and you can
      ! proceed directly.
      goto 9999
  elseif (Key_Analysis_Type == 53) then
      print *, '    Analysis type:             PFEM FS coupling'
      call log_msg('Analysis type:             PFEM FS coupling')
      ! Since there is currently no input data for this analysis, no checks are needed, and you can
      ! proceed directly.
      goto 9999
  elseif (Key_Analysis_Type == 54) then
      print *, '    Analysis type:             PFEM FS-EBE'
      call log_msg('Analysis type:             PFEM FS-EBE')
      ! Since there is currently no input data for this analysis, no checks are needed, and you can
      ! proceed directly.
      goto 9999
  elseif (Key_Analysis_Type == 55) then
      print *, '    Analysis type:             PFEM FS-MC'
      call log_msg('Analysis type:             PFEM FS-MC')
      ! Since there is currently no input data for this analysis, no checks are needed, and you can
      ! proceed directly.
      goto 9999
  elseif (Key_Analysis_Type == 61) then
      print *, '    Analysis type:             Hydraulic fracturing experiment'
      call log_msg('Analysis type:             Hydraulic fracturing experiment')
      ! Since there is currently no input data for this analysis, no detection is needed, and you can
      ! proceed directly.
  else
      call Warning_Message('E','Key_Analysis_Type')
  end if

  !************************************************
  ! 2D Structure-Field Coupled Problem Computation
  !************************************************
  if(Key_FS_Seque_Coupling==1)then
      ! For 2D problems only
      if(Key_Dimension/=2)then
          print *, '    Error :: *Key_FS_Seque_Coupling=1 for 2D problems only!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Only for statics problems
      if(Key_Analysis_Type/=1)then
          print *, '    Error :: *Key_FS_Seque_Coupling=1 for *Key_Analysis_Type=1 only!'
          call Warning_Message('S',Keywords_Blank)
      endif
      print *,'    WARNING :: *Key_FS_Seque_Coupling:    active!'
      print *,'                solid-field sequential coupling.'
  endif
  
  !*****************************************************************
  ! 2D Structure - Thermal Coupling Calculation Related, 2021-11-03
  !*****************************************************************
  if(Key_TS_Seque_Coupling==1)then
      ! For 2D problems only
      if(Key_Dimension/=2)then
          print *, '    Error :: *Key_TS_Seque_Coupling=1 for 2D problems only!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Only for statics problems
      if(Key_Analysis_Type/=1)then
          print *, '    Error :: *Key_TS_Seque_Coupling=1 for *Key_Analysis_Type=1 only!'
          call Warning_Message('S',Keywords_Blank)
      endif
      print *,'    WARNING :: *Key_TS_Seque_Coupling:    active!'
      print *,'                solid-thermo sequential coupling.'
  endif

  !************************************
  ! Random number generation algorithm
  !************************************
1721 FORMAT(5X,'Seed is taken as:      ',T27,I13)
1722 FORMAT(5X,'Randomly generated seed:',I8)
  if (Key_Random == 0) then
      print *, '    Random number algorithm:   call random_number and never change!'
  elseif (Key_Random == 1) then
      print *, '    Random number algorithm:   call random_number and change every time!'
  elseif (Key_Random == 2) then
      print *, '    Random number algorithm:   generate with a seed'
      write (*,1721) Seed
      print *, '    Saving seed file...'
      call Save_File_Seed(Seed)
  ! Random Number Generation Scheme 3. 2022-10-26. NEWFTU2022102601.
  elseif (Key_Random == 3) then
      print *, '    Random number algorithm:   generate with a random seed'
      call Init_Random_Seed()
      call random_number (Ran_Num)
      Seed = int(Ran_Num*1000000.0D0)
      write (*,1722) Seed
      print *, '    Saving seed file...'
      call Save_File_Seed(Seed)
  else
      call Warning_Message('E','Key_Random')
  end if

  !****************************************************************************
  ! Initial inclusions: including circular inclusions and polygonal inclusions
  !****************************************************************************
  if (Key_Rand_Poly_Incl ==1 .and.Key_Rand_Poly_Incl_Irregular  == 1)then
      print *, '    Warning: *Key_Rand_Poly_Incl and *Key_Rand_Poly_Incl_Irregular cannot be'  &
               //' set to 1 at the same time!'
      call Warning_Message('E','Key_Rand_Poly_Incl')
  endif


  !********************************
  ! Ground Stress Treatment Method
  !********************************
  if(Key_InSitu_Strategy==0)then
      print *, '    Insitu stress strategy:    no'
  elseif(Key_InSitu_Strategy==1)then
      print *, '    Insitu stress strategy:    linear superposition'
  elseif(Key_InSitu_Strategy==2)then
      print *, '    Insitu stress strategy:    Zienkiewicz_7ed_P216'
      ! For clean water fracturing analysis (analysis type = 4), the geostress calculation should be 0.
      if (Key_Analysis_Type==4) then
          print *,'    *Key_InSitu_Strategy should be 0 if *Key_Analysis_Type =4'
          call Warning_Message('S',Keywords_Blank)
      endif
  elseif(Key_InSitu_Strategy==3)then
      print *, '    Insitu stress strategy:    neglected!'
  elseif(Key_InSitu_Strategy==4)then
      print *, '    Insitu stress strategy:    Zienk_with_fixed_bou'
      !2022-12-21.
      if(Key_Read_Initial_Node_Stress_File==1) then
          print *, '    Initial nodal stresses:    read from *.istn file.'
      endif
  else
      call Warning_Message('E','Key_InSitu_Strategy')
  endif
  ! Geostress processing does not support elliptical holes (2020-08-09)
  if(Key_InSitu_Strategy/=0)then
      if(num_Ellip_Hole>=1)then
          print *,'    *Key_InSitu_Strategy should be 0 if *Key_Analysis_Type =4'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif
  !2022-12-21.
  if(Key_Read_Initial_Node_Stress_File==1) then
      if(Key_InSitu_Strategy/=4) then
          print *,'    *Key_Read_Initial_Node_Stress_File =1 is valid only when *Key_InSitu_Strategy=4!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif

  !*************************************************************
  ! Applying Prestress to Specific Material Number (2019-09-21)
  !*************************************************************
1776 FORMAT(5X,'Mat contains prestress:',I3)
1777 FORMAT(5X,'Mat ',I3, ' not defined!')
1778 FORMAT(5X,'Parameters of mat ',I3, ' not defined!')
  if(Key_InStress_for_Mat==1)then
      if(sum(Mat_Number_of_InStress(1:100))==0)then
          print *,'    *Material number not defined for *Key_InStress_for_Mat!'
          call Warning_Message('S',Keywords_Blank)
      else
        do i=1,100
          if(Mat_Number_of_InStress(i)>=1)then
              print *, '    Prestress for material:    yes'
              write(*,1776) Mat_Number_of_InStress(i)
              ! Check whether the parameters for this material number have been correctly defined
              if(num_of_Material <Mat_Number_of_InStress(i))then
                  write(*,1777) Mat_Number_of_InStress(i)
                  call Warning_Message('S',Keywords_Blank)
              endif
              if(Material_Para(Mat_Number_of_InStress(i),1)==ZR)then
                  write(*,1778) Mat_Number_of_InStress(i)
                  call Warning_Message('S',Keywords_Blank)
              endif
          endif
        enddo
      endif
      if(Mat_InStress_x==ZR .and. Mat_InStress_y==ZR .and. Mat_InStress_z==ZR) then
          print *,'    *Mat_InStress_x or _y or _z not defined for *Key_InStress_for_Mat!'
          call Warning_Message('S',Keywords_Blank)
      endif
  elseif (Key_InStress_for_Mat==0)then
      print *, '    Prestress for material:    no'
  endif

  !***********************************
  ! Check num_of_Material, 2021-07-18
  !***********************************
  ! Get the maximum material number in the *.elem file
  3776 FORMAT(14X,'Maximum material type number in *.elem file is',I4,'.')
  3777 FORMAT(14X,'Defined material type number is',I4,'.')
  if (Key_Cpp_Call_Fortran_Lib/=1) then
    max_material_type = maxval(Elem_Mat(1:Num_Elem))
    if (num_of_Material < max_material_type)then
          print *, '    Error :: Number of material types is not enough!'
          write(*,3776) max_material_type
          write(*,3777) num_of_Material
          call Warning_Message('S',Keywords_Blank)
    endif
  endif

  !*********************
  ! Mianli (2023-01-21)
  !*********************
  3778 FORMAT(5X,'Number of surface loads:',I4)
  if(Key_Dimension==3) then
      if (Num_Surface_Loads>=1)then
          print *, '    Surface loads:             active'
          write(*,3778) Num_Surface_Loads
      endif
      if (Num_Surface_Loads>100)then
          print *, '    Error :: Num_Surface_Loads cannot be larger than 100!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif

  !**************************************************
  ! Regarding the testing related to proppant-supported fracture aperture calculation, 2016-08-25
  !--------------------------------------------------
  !Key_Propped_Width=1 and Key_Analysis_Type = 1
  !**************************************************
  if(Key_Propped_Width==1)then
      if(Key_Analysis_Type == 1)then
          print *, '    Warning: special static analysis, calculate width of propped fractures!'
          ! If calculating the fracture width supported by proppant, then only one step is calculated.
          if(Num_Substeps>1)then
              print *, '    Warning :: only 1 substep required,*Num_Substeps is reset to 1!'
              Num_Substeps = 1
          endif
          ! Quasi-static analysis of support crack opening requires turning on the contact analysis switch.
          if(Key_Contact==0)then
              print *, '    Error :: Key_Contact can not be 0 when Key_Propped_Width=1!'
              call Warning_Message('S',Keywords_Blank)
          endif
      else
      endif
  endif

  !********************************************
  ! Initial pore pressure, added on 2016-08-14
  !********************************************
  501 FORMAT(5X,'Pore pressure value:  ',F10.3,' MPa')
  ! Check for traffic sources within the area
  if(Key_PoreP==1)then
      if (Key_Cpp_Call_Fortran_Lib/=1) then
          print *, '    Initial pore pressure:     yes    '
          if(Initial_PoreP==-999.0D0)then
              print *, '    Error :: *Initial_PoreP not defined!'
              call Warning_Message('S',Keywords_Blank)
          else
              write(*,501) Initial_PoreP/1.0D6
          endif
      endif
  elseif(Key_PoreP==0)then
      print *, '    Initial pore pressure:     no'
  else
      call Warning_Message('E','Key_PoreP')
  endif

  !**************************
  ! Fatigue Analysis Related
  !**************************
  if (Key_Static_Fatigue==1)then
      if(Key_Analysis_Type == 1)then
          print *, '    Fatigue analysis:          yes'
          call log_msg('Fatigue analysis:          yes')
          ! Fatigue criteria
          if(Key_Fatigue_Cri==1)then
              print *, '    Fatigue criterion:         Paris'
          else
              call Warning_Message('E','Key_Fatigue_Cri')
          endif
      else
          print *, '    Error :: fatigue analysis is usable only for static analysis, check *Key_Analysis_Type!'
          call Warning_Message('S',Keywords_Blank)
      endif
  elseif(Key_Static_Fatigue==0)then
      print *, '    Fatigue analysis:          no'
  else
      call Warning_Message('E','Key_Static_Fatigue')
  end if

  !*************
  ! Data Format
  !*************
  if (Key_Data_Format==1) then
      print *, '    Data format:               formatted'
      call log_msg('Data format:               formatted')
  elseif (Key_Data_Format==2) then
      ! Note: Binary storage does not support the mm-ton-s unit system
      if(Key_Unit_System==2)then
          print *, '    Error :: binary data is not avaliable when *Key_Unit_System=2!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Inspection passed
      print *, '    Data format:               binary'
      call log_msg('Data format:               binary')
  else
      call Warning_Message('E','Key_Data_Format')
  end if

  !*******************************
  ! Number of cores set by OpenMP
  !*******************************
  1771 FORMAT(5X,'Number of threads:          ',T31,I3)
  if (Key_Num_Process==99) then
      print *, '    Number of threads:         all'
      call log_msg('Number of threads:         all')
  elseif (Key_Num_Process<=0) then
      call Warning_Message('E','Key_Num_Process')
  else
      write (*,1771) Key_Num_Process
      call Tool_chrpak_i4_to_s_left (Key_Num_Process, s_Key_Num_Process)
      temp_log = trim('Number of threads:         '//s_Key_Num_Process)
      call log_msg(temp_log)
  endif

  !*****************
  ! Program Control
  !*****************
  if (Key_Play_Sounds==0) then
      print *, '    Play sounds:               no'
  elseif (Key_Play_Sounds==1) then
      print *, '    Play sounds:               yes'
  else
      call Warning_Message('E','Key_Play_Sounds')
  end if
  if (Key_Close_Window==0) then
      print *, '    Close terminal when done:   no'
  elseif (Key_Close_Window==1) then
      print *, '    Close terminal when done:   yes'
  else
      call Warning_Message('E','Key_Close_Window')
  end if
  if (Key_Window_Log==0) then
      print *, '    Save log file:             no'
  elseif (Key_Window_Log==1) then
      print *, '    Save log file:             yes'
  else
      call Warning_Message('E','Key_Window_Log')
  end if
  if (Key_Clear_All==0) then
      print *, '    Clear files when starting: no'
  elseif (Key_Clear_All==1) then
      print *, '    Clear files when starting: yes (default)'
  else
      call Warning_Message('E','Key_Clear_All')
  end if
  if (Key_Visit_PhiPsi_top==0) then
      print *, '    Visit http://phipsi.top:   no (default)'
  elseif (Key_Visit_PhiPsi_top==1) then
      print *, '    Visit http://phipsi.top:   yes'
  else
      call Warning_Message('E','Key_Visit_PhiPsi_top')
  end if

  !***********************
  ! Check junction points
  !***********************
  if (Key_Junction_Check/=0 .and. Key_Junction_Check/=1) then
      call Warning_Message('E','Key_Play_Sounds')
  end if

  !***********************************************************************************************
  ! Stiffness matrix storage format: sparse or full matrix format (currently only supports static
  ! non-contact problems)
  !***********************************************************************************************
1101 FORMAT(5X,'Assumptive sparse ratio:  ',F5.2,'%')
1109 FORMAT(5X,'Assumptive sparse ratio:  ',F5.2,'%, given by PhiPsi')
  if (Key_K_Sparse==0) then
      print *, '    Sparse matric format:      no'
      ! For three-dimensional problems, if the number of elements exceeds 10,000, it is recommended to use
      ! sparse matrix storage.
      if (Key_SLOE/=11)then
        if(Key_Dimension==3 .and. Num_Elem > 10000)then
          print *, '      WARNING :: Key_K_Sparse=1 is suggested!  '
          print *, '                 Sparse matrix can save memory!'
        endif
      endif
  elseif (Key_K_Sparse==1) then
      print *, '    Sparse matric format:      yes'
      
      !IMPROV2024110801.
      ! For 2D problems, if the number of elements is less than 5,000, the sparsity ratio is set to 0.1.
      ! 2024-11-08.
      if(Num_Elem<=5000 .and. Key_Dimension==2) then
          if (Sparse_Ratio<0.1D0) then
              Sparse_Ratio = 0.1D0
          endif
          print *,'    Warning: *Sparse_Ratio has been set as 0.1 for small model!'
      endif
      
      if(Sparse_Ratio > ONE)then
          print *,'    *Sparse_Ratio not defined or illegal!'
          call Warning_Message('E','Sparse_Ratio')
      else
          if(Sparse_Ratio > ZR)then
              write(*,1101) Sparse_Ratio*100.0D0
          else
              write(*,1109) Sparse_Ratio*100.0D0
          endif
      endif
      ! Detect Sparse_Store_Method
      if(Key_Dimension==2)then
          if(Sparse_Store_Method==3)then
          print *,'    ERROR::Sparse_Store_Method=3 not avaliavble for 2D problems!'
          call Warning_Message('S',Keywords_Blank)
          endif
      endif
      
      ! Temporarily not supported for contact issues
      if (Key_Contact /=0) then
          print *,'    ERROR :: Sparse matric format does not support contact yet!'
          print *,'             Set *Key_K_Sparse=1 and try again.'
          call Warning_Message('S',Keywords_Blank)
      endif
      
      ! For 2D, Only static analysis is supported.
      if(Key_Dimension==2)then
          if (Key_Analysis_Type /=1) then
              print *,'    ERROR :: Sparse matric format only supports static analysis yet for 2D problems!'
              call Warning_Message('S',Keywords_Blank)
          endif
      endif
      
      ! For 3D, Only static analysis or slipwater HF is supported. 2025-12-02.
      if(Key_Dimension==3)then
          if (Key_Analysis_Type /=1 .and. Key_Analysis_Type /=4) then
              print *,'    ERROR :: For 3D problems, sparse matric format only supports static analysis or slipwater HF yet!'
              call Warning_Message('S',Keywords_Blank)
          endif
      endif
      
      ! Only supports solvers 6, 7, 8, 9, and 12 (such as UMFPACK and Lis)
      if (Key_SLOE /=6 .and. Key_SLOE /=7 .and. Key_SLOE/=8 .and. Key_SLOE/=9 .and. Key_SLOE/=12) then
          print *,'    Sparse matric format only supports solver 6 or 7 or 8 or 9 or 12!'
          call Warning_Message('S',Keywords_Blank)
      endif
  else
      call Warning_Message('E','Key_K_Sparse')
  end if

  !***************************************************
  ! Whether embedding of the crack surface is allowed
  !***************************************************
1191 FORMAT(5X,'Normal penalty paremeter:',E12.4,' GPa/m')
1192 FORMAT(5X,'Shear penalty paremeter: ',E12.4,' GPa/m')
1193 FORMAT(5X,'Friction coefficient:     ',F4.1)
  if(Key_Contact==1) then
      print *, '    Crack contact:             yes'
      call log_msg('Crack contact:             yes')
      print *, '    Contact scheme:            1, penalty method'
      ! Next, check whether the tangential and normal penalty parameters are defined
      if(kn_Cont_Penalty<0) then
          print *, '    Error :: *kn_Cont_Penalty is wrong!'
          call Warning_Message('E','kn_Cont_Penalty')
      elseif(kt_Cont_Penalty<0) then
          print *, '    Error :: *kt_Cont_Penalty is wrong!'
          call Warning_Message('E','kt_Cont_Penalty')
      else
          write(*,1191) kn_Cont_Penalty/1.0D9
          write(*,1192) kt_Cont_Penalty/1.0D9
      endif
      ! 3D Contact Analysis
      if(Key_Dimension==3)then
          ! The contact penalty stiffness in 3D contact analysis should not be too high.
          if(kn_Cont_Penalty >=1.0D13) then
            print *, '    WARNING :: *kn_Cont_Penalty is too big!'
            print *, '                1.0D12 is suggested for 3D problems!'
          endif
          ! The contact penalty stiffness in 3D contact analysis should not be too high.
          if(kt_Cont_Penalty >=1.0D13) then
            print *, '    WARNING :: *kt_Cont_Penalty is too big!'
            print *, '                1.0D12 is suggested for 3D problems!'
          endif
      endif
  ! Penalty function method (reduction).
  elseif(Key_Contact==2) then
      print *, '    Embedded crack:            not allowed'
      print *, '    Contact scheme:            2, reduced penalty method'
  ! Only keep the aperture open; if the aperture is less than 0, it will automatically be adjusted to
  ! 0. 2022-10-20. NEWFTU2022102001.
  elseif(Key_Contact==5) then
      print *, '    Contact scheme:            5, keep only positive aperture!'
  ! Two-step calculation: if the opening is negative, the corresponding element is marked as a
  ! shear-compression element, and the stiffness matrix is corrected using the penalty function
  ! method. NEWFTU2022102101.
  elseif(Key_Contact==6) then
      print *, '    Crack contact:             yes'
      print *, '    Contact scheme:            6, CS K modification!'
      ! Only supports Solver 11. 2022-10-21.
      if(Key_SLOE /=11)then
          print *,'    Contact scheme 6 only supports solver 11!'
          call Warning_Message('S',Keywords_Blank)
      endif
  elseif(Key_Contact==7) then
      print *, '    Crack contact:             yes'
      call log_msg('Crack contact:             yes')
      print *, '    Contact scheme:            7, CS penalty method'
      if(Key_SLOE /=11)then
          print *,'    Contact scheme 7 only supports solver 11!'
          call Warning_Message('S',Keywords_Blank)
      endif
      if(Key_EBE_Precondition ==1)then
          print *,'    Contact scheme 7 is not available when *Key_EBE_Precondition =1!'
          call Warning_Message('S',Keywords_Blank)
      endif
  elseif(Key_Contact==0) then
      print *, '    Embedded crack:            allowed'
  else
      call Warning_Message('E','Key_Contact')
  end if
  if(Key_Contact/=0) then
      ! Check the number of contact integration points (should be 1 or 2)
      if(Conta_Integ_Point/=1 .and. Conta_Integ_Point/=2)then
          call Warning_Message('E','Conta_Integ_Point')
      endif
      ! Check if the coefficient of friction is greater than or equal to 0
      if(fric_mu_Cont<0) then
          print *, '    Error :: *fric_mu_Cont cannot be negative!'
          call Warning_Message('E','fric_mu_Cont')
      else
          write(*,1193) fric_mu_Cont
      endif
  endif

  !*********************************************************************************
  ! Penalty function treatment for compression-shear natural fractures. 2022-10-22.
  ! NEWFTU2022102201.
  !*********************************************************************************
  if(Key_Dimension==3) then
    if(Key_CS_Natural_Crack==1) then
      print *, '    CS natural crack penalty:  active'
    end if
  endif

  !********************************
  ! Tip-splitting enhancement type
  !********************************
1901 FORMAT(5X,'Tip enriched node radius factor: ',F8.4)
  ! 3D can only use Key_TipEnrich==0 and Key_TipEnrich==1
  if(Key_Dimension==3) then
    if(Key_TipEnrich==0) then
      print *, '    Tip enrichment method:     no tip enrichment'
      call log_msg('Tip enrichment method:     no tip enrichment')
    elseif(Key_TipEnrich==1) then
      print *, '    Tip enrichment method:     common 4 F-functions'
      call log_msg('Tip enrichment method:     common 4 F-functions')
    elseif(Key_TipEnrich==2) then
      print *, '    Tip enrichment method:     only 1 F-function'
      call log_msg('Tip enrichment method:     only 1 F-function')
    else
      print *, '    Error :: tip enrichment error for 3D problems!'
      print *, '             PLease check *Key_TipEnrich!'
      call Warning_Message('S',Keywords_Blank)
    end if
  endif

  ! 2D problem
  if(Key_Dimension==2) then
      if(Key_TipEnrich==0) then
        print *, '    Tip enrichment method:     no tip enrichment'
        call log_msg('Tip enrichment method:     no tip enrichment')
      elseif(Key_TipEnrich==1) then
        print *, '    Tip enrichment method:     common 4 F-functions'
        call log_msg('Tip enrichment method:     common 4 F-functions')
      elseif(Key_TipEnrich==2) then
        print *, '    Tip enrichment method:     only 1 F-function'
        call log_msg('Tip enrichment method:     only 1 F-function')
      elseif(Key_TipEnrich==3) then
        print *, '    Tip enrichment method:     smoothed Heaviside'
        call log_msg('Tip enrichment method:     smoothed Heaviside')
      elseif(Key_TipEnrich==4) then
        print *, '    Tip enrichment method:     cohesive crack'
        call log_msg('Tip enrichment method:     cohesive crack')
        !Cohesive crack only for Key_Analysis_Type = 1
        if(Key_Analysis_Type /= 1)then
          print *, '    Error :: Cohesive crack is available only for static analysis!'
          call Warning_Message('E','Key_TipEnrich')
        endif
        ! The following is about detection related to cohesive cracks
        !The following checks are about the cohesive cracks, 2017-04-15.

        if(Coh_Constitutive_type==1)then
          if(Coh_Width_Critical1==-999.0D0)then
            print *, '    Error :: *Coh_Width_Critical1 not given!'
            call Warning_Message('S',Keywords_Blank)
          endif
          if(Coh_Width_Critical2==-999.0D0)then
            print *, '    Error :: *Coh_Width_Critical2 not given!'
            call Warning_Message('S',Keywords_Blank)
          endif
          if(Coh_Width_Critical1 > Coh_Width_Critical2) then
              print *, '    Error :: *Coh_Width_Critical2 must > *Coh_Width_Critical1!'
              call Warning_Message('S',Keywords_Blank)
          endif
          if(Coh_f_Ultimate==-999.0D0)then
              print *, '    Error :: *Coh_f_Ultimate not given!'
              call Warning_Message('S',Keywords_Blank)
          endif
          print *, '    Cohesive crack type:       bilinear'
        elseif(Coh_Constitutive_type==2)then
          if(Coh_Width_Critical2==-999.0D0)then
            print *, '    Error :: *Coh_Width_Critical2 not given!'
            call Warning_Message('S',Keywords_Blank)
          endif
          if(Coh_f_Ultimate==-999.0D0)then
              print *, '    Error :: *Coh_f_Ultimate not given!'
              call Warning_Message('S',Keywords_Blank)
          endif
          print *, '    Cohesive crack type:       linear decline'
        elseif(Coh_Constitutive_type==3)then
          if(Coh_Width_Critical2==-999.0D0)then
            print *, '    Error :: *Coh_Width_Critical2 not given!'
            call Warning_Message('S',Keywords_Blank)
          endif
          if(Coh_f_Ultimate==-999.0D0)then
              print *, '    Error :: *Coh_f_Ultimate not given!'
              call Warning_Message('S',Keywords_Blank)
          endif
          print *, '    Cohesive crack type:       linear decline'
        elseif(Coh_Constitutive_type==-999)then
          print *, '    Error :: *Coh_Constitutive_type not given!'
          call Warning_Message('S',Keywords_Blank)
        else
          call Warning_Message('E','Coh_Constitutive_type')
        endif
        !Tangential
        if (Coh_Tangential_Key ==0)then
              print *, '    Cohesive in tangential:    deactivated'
        elseif (Coh_Tangential_Key ==1)then
          if(Coh_Constitutive_type==1)then
              if(Coh_Width_Critical1_T==-999.0D0)then
                print *, '    Error :: *Coh_Width_Critical1_T not given!'
                call Warning_Message('S',Keywords_Blank)
              endif
              if(Coh_Width_Critical2_T==-999.0D0)then
                print *, '    Error :: *Coh_Width_Critical2_T not given!'
                call Warning_Message('S',Keywords_Blank)
              endif
              if(Coh_Width_Critical1_T>Coh_Width_Critical2_T)then
                print *, '    Error :: *Coh_Width_Critical2_T must > *Coh_Width_Critical1_T!'
                call Warning_Message('S',Keywords_Blank)
              endif
              if(Coh_f_Ultimate_T==-999.0D0)then
                print *, '    Error :: *Coh_f_Ultimate_T not given!'
                call Warning_Message('S',Keywords_Blank)
              endif
          elseif(Coh_Constitutive_type==2)then
              if(Coh_Width_Critical2_T==-999.0D0)then
                print *, '    Error :: *Coh_Width_Critical2_T not given!'
                call Warning_Message('S',Keywords_Blank)
              endif
              if(Coh_f_Ultimate_T==-999.0D0)then
                print *, '    Error :: *Coh_f_Ultimate_T not given!'
                call Warning_Message('S',Keywords_Blank)
              endif
          elseif(Coh_Constitutive_type==3)then
              if(Coh_Width_Critical2_T==-999.0D0)then
                print *, '    Error :: *Coh_Width_Critical2_T not given!'
                call Warning_Message('S',Keywords_Blank)
              endif
              if(Coh_f_Ultimate_T==-999.0D0)then
                print *, '    Error :: *Coh_f_Ultimate_T not given!'
                call Warning_Message('S',Keywords_Blank)
              endif
          endif
          print *, '    Cohesive in tangential:    activated'
        elseif(Coh_Tangential_Key ==-999)then
              print *, '    Error :: *Coh_Tangential_Key not given!'
              call Warning_Message('S',Keywords_Blank)
        else
              call Warning_Message('E','Coh_Tangential_Key')
        endif
      else
          call Warning_Message('E','Key_TipEnrich')
      end if
  endif

  if(Key_TipEnrich>=1 .and. Key_Multi_TipEnrNode==1) then
      print *, '    Multi-tip-enriched-nodes:  activated!'
      if(Key_Dimension==2)then
          if(Key_TipEnrich_Radius_Factor<=ZR) then
            call Warning_Message('E','Key_TipEnrich_Radius_Factor')
          else
            write(*,1901) Key_TipEnrich_Radius_Factor
          endif
      endif
  endif
  
  !*************************************************
  ! Method for Calculating Crack Width. 2023-08-15.
  !*************************************************
  if (Key_Crack_Aperture_Method==1) then
      print *, '    Method to calculate crack aperture:       1 - equation.'
  elseif (Key_Crack_Aperture_Method==2) then
      print *, '    Method to calculate crack aperture:       2 - offset.'
  else
    call Warning_Message('E','Key_Crack_Aperture_Method')
  endif
  
  
  !***************************************
  !Key_Ele_Max_Related_Cracks,2023-02-25.
  !***************************************
  if(Key_Dimension==3) then
      print *, '    Key_Ele_Max_Related_Cracks: ',Key_Ele_Max_Related_Cracks
      if(Key_Ele_Max_Related_Cracks<=1)then
          print *, '    Error :: *Key_Ele_Max_Related_Cracks cannot <=1!'
          call Warning_Message('S',Keywords_Blank)
      elseif(Key_Ele_Max_Related_Cracks>=20)then
          print *, '    Error :: *Key_Ele_Max_Related_Cracks cannot >=10!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif

  !**********************************************
  ! Enhanced Unit Local Optimization, 2021-06-27
  !**********************************************
  if(Key_Local_Mesh_Refine==0) then
    print *, '    Local refinement of enriched element:     no'
  elseif(Key_Local_Mesh_Refine==1) then
    print *, '    Local refinement of enriched element:     all enriched elements!'
    print *, '    Attention :: Local refinement of enriched element is active!'
  elseif(Key_Local_Mesh_Refine==2) then
    print *, '    Local refinement of enriched element:     only tip enriched elements!'
    print *, '    Attention :: Local refinement of enriched element is active!'
  else
    call Warning_Message('E','Key_Local_Mesh_Refine')
  endif
  if(Key_Local_Mesh_Refine >0) then
      !Now for static problems only.
      if (Key_Analysis_Type /= 1)then
          print *, '    Error :: Local mesh refinement for static problems only!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif


  !********************************
  !Key_Adj_Ini_Crack_3D,2021-08-19
  !********************************
  if (Key_Dimension ==3)then
      if(Key_Adj_Ini_Crack_3D==0) then
        print *, '    Adjust initial 3D Crack:                  no'
      elseif(Key_Adj_Ini_Crack_3D==1) then
        print *, '    Adjust initial 3D Crack:                  yes'
        print *, '    Attention :: Adjust initial crack is active!'
      else
        call Warning_Message('E','Key_Adj_Ini_Crack_3D')
      endif
  endif

  if(Key_Adj_Ini_Crack_3D==1) then
      !Now for 3D problems only.
      if (Key_Dimension ==2)then
          print *, '    Error :: Adjust initial crack for 3D problems only!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif

  !******************************************
  !Key_Check_and_Adjust_Cracks_3D,2022-08-01
  !******************************************
  if (Key_Dimension ==3)then
      if(Key_Check_and_Adjust_Cracks_3D==0) then
        print *, '    Check and adjust 3D Cracks:               no'
      elseif(Key_Check_and_Adjust_Cracks_3D>=1) then
        print *, '    Check and adjust 3D Cracks:               yes'
        
        print *, '    Attention :: Adjust initial crack is active!'
        print *, '                 Scheme is     ',Key_Check_and_Adjust_Cracks_3D
        if(Key_Check_and_Adjust_Cracks_3D>=2 .and. Key_Check_and_Adjust_Cracks_3D/=5) then
            print *, '                 Resolution is ',Adjust_Cracks_Resolution
        endif
        
        ! Resolution 6, 18, 26
        if(Adjust_Cracks_Resolution==6 .or. Adjust_Cracks_Resolution==18 .or. Adjust_Cracks_Resolution==26) then
            !Do nothing.
        else
            call Warning_Message('E','Adjust_Cracks_Resolution')
        endif
      else
        call Warning_Message('E','Key_Check_and_Adjust_Cracks_3D')
      endif
  endif


  !************************************
  !Crack front segmentation,2021-08-19
  !************************************
5101 FORMAT(5X,'Number of front segments:  ',I4)
  if (Key_Dimension ==3)then
      if(Key_Front_Segmentation==0) then
        print *, '    Crack front segmentation:                 no'
      elseif(Key_Front_Segmentation==1) then
        print *, '    Crack front segmentation:                 yes'
        
        print *, '    Attention :: Crack front segmentation is active!'
        
      else
        call Warning_Message('E','Key_Front_Segmentation')
      endif
  endif

  if(Key_Front_Segmentation ==1) then
      !Now for 3D problems only.
      if (Key_Dimension ==2)then
          print *, '    Error :: Crack front segmentation for 3D problems only!'
          call Warning_Message('S',Keywords_Blank)
      endif
      !Now for static problems only.
      if (Key_Analysis_Type /= 1)then
          print *, '    Error :: Crack front segmentation for static problems only!'
          call Warning_Message('S',Keywords_Blank)
      endif
      !Check Number_front_Segments (2021-08-20).
      if(Number_front_Segments>=2) then
          write(*,5101) Number_front_Segments
      else
          print *, '    Error :: *Number_front_Segments must >=2!'
          call Warning_Message('E','Number_front_Segments')
      endif
      ! There must be an initial crack
      if(num_Crack<=0) then
          print *, '    Error :: Crack front segmentation needs initial crack!'
          print *, '             num_Crack equals 0 now!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Only supports initial circular cracks
      if(sum(abs(Crack3D_Cir_Coor(1:num_Crack,1:7)))<=Tol_10)then
          print *, '    Error :: Crack front segmentation needs initial circular crack!'
          print *, '             Can not find Crack3D_Cir_Coor!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif

  !************************************
  !Crack front segmentation,2021-08-19
  !************************************
  if (Key_Dimension ==3)then
    if(Key_Smooth_Front==0) then
      print *,'    Crack front smooth:                       no'
    elseif(Key_Smooth_Front==1) then
      print *,'    Crack front smooth:                       yes'
      print *,'    Smooth method:                            line'
      
      print *,'    Attention :: Crack front smooth is active!'
      
    elseif(Key_Smooth_Front==2) then
      print *,'    Crack front smooth:                       yes'
      print *,'    Smooth method:                            circle'
      
      print *,'    Attention :: Crack front smooth is active!'
      
    elseif(Key_Smooth_Front==3) then
      print *,'    Crack front smooth:                       yes'
      print *,'    Smooth method:                            CSCPS'
      
      print *,'    Attention :: Crack front smooth is active!'
      print *,'                 Call Smooth_3D_Points_CSAPS.py!'
      
    elseif(Key_Smooth_Front==4) then
      print *,'    Crack front smooth:                       yes'
      print *,'    Smooth method:                            Bspline'
      
      print *,'    Attention :: Crack front smooth is active!'
      print *,'                 call Smooth_3D_Points_Scipy.py!'
      
    elseif(Key_Smooth_Front==5) then
      print *,'    Crack front smooth:                       yes'
      print *,'    Smooth method:                            Ellipse'
      
      print *,'    Attention :: Crack front smooth is active!'
      
    elseif(Key_Smooth_Front==6) then
      print *,'    Crack front smooth:                       yes'
      print *,'    Smooth method:                            Taubin'
      
      print *,'    Attention :: Crack front smooth is active!'
      
    else
      call Warning_Message('E','Key_Smooth_Front')
    endif
  endif

  if(Key_Smooth_Front >=1) then
      !Now for 3D problems only.
      if (Key_Dimension ==2)then
          print *, '    Error :: Crack front smooth for 3D problems only!'
          call Warning_Message('S',Keywords_Blank)
      endif
      !Now for static and HF problems only.
      if (Key_Analysis_Type /= 1 .and. Key_Analysis_Type /= 4 .and. Key_Analysis_Type /= 61)then
          print *, '    Error :: Crack front smooth for static problems only!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif

  if(Key_Smooth_Front_Twice ==1) then
      print *,'    Crack front smooth twice:                 yes'
      
      print *,'    Attention :: Crack front smooth-twice is active!'
      print *,'                 For closed crack front only!'
      
      !Now for 3D problems only.
      if (Key_Dimension ==2)then
          print *, '    Error :: *Key_Smooth_Front_Twice for 3D problems only!'
          call Warning_Message('S',Keywords_Blank)
      endif
      !Now for static and HF problems only.
      if (Key_Analysis_Type /= 1 .and. Key_Analysis_Type /= 4)then
          print *, '    Error :: *Key_Smooth_Front_Twice for static problems only!'
          call Warning_Message('S',Keywords_Blank)
      endif
  elseif(Key_Smooth_Front_Twice >=2)then
      call Warning_Message('E','Key_Smooth_Front_Twice')
  endif

  !************************************
  ! Junction Enhanced. 3D. 2022-08-25.
  !************************************
  if (Key_Dimension ==3)then
    if(Key_Junction_Enrich==1) then
        print *,'    Junction enrichment for 3D problems:      yes'
        
        print *,'    Attention :: Junction enrichment for 3D problems is active!'
        
    elseif(Key_Junction_Enrich==0) then
        print *,'    Junction enrichment for 3D problems:      no'
    else
        call Warning_Message('E','Key_Junction_Enrich')
    endif
  endif

  !*********************************************
  ! Crack Front Edge Data Denoising, 2022-04-25
  !*********************************************
!1127 FORMAT(5X,55('/'))
  if (Key_Dimension ==3)then
    if(Key_Denoise_Vertex_Value==0) then
      print *,'    Crack front vertex value denoising:       no'
    elseif(Key_Denoise_Vertex_Value==1) then
      print *,'    Crack front vertex value denoising:       yes'
      
      print *,'    Attention :: Crack front value denoising is active!'
      
    else
      call Warning_Message('E','Key_Denoise_Vertex_Value')
    endif
  endif

  !*******************************************
  ! Smoothing of crack front data, 2022-04-25
  !*******************************************
!1128 FORMAT(5X,55('\'))
  if (Key_Dimension ==3)then
    if(Key_Smooth_Vertex_Value==0) then
      print *,'    Crack front vertex value smoothing:       no'
    elseif(Key_Smooth_Vertex_Value==1) then
      print *,'    Crack front vertex value smoothing:       yes'
      print *,'    Smoothing method:                         moving average method'
      
      print *,'    Attention :: Crack front value smoothing is active!'
      
    else
      call Warning_Message('E','Key_Smooth_Vertex_Value')
    endif
  endif

  !***************************************
  ! Pre-crack Theta denoising, 2022-07-14
  !***************************************
  if (Key_Dimension ==3)then
    if(Key_Denoise_Theta_Value==0) then
      print *,'    Crack front Theta angle denoising:       no'
    elseif(Key_Denoise_Theta_Value==1) then
      ! Only applicable for CFCP=3
      if(CFCP/=3) then
          print *, '    Error :: *Key_Denoise_Theta_Value=1 is valid only when *CFCP=3!'
          call Warning_Message('S',Keywords_Blank)
      endif
      print *,'    Crack front Theta angle denoising:       yes'
      
      print *,'    Attention :: Crack front Theta angle denoising is active!'
      
    else
      call Warning_Message('E','Key_Denoise_Theta_Value')
    endif
  endif

  !************************************************************
  ! Smooth treatment of the crack front edge Theta, 2022-07-14
  !************************************************************
  if (Key_Dimension ==3)then
    if(Key_Smooth_Theta_Value==0) then
      print *,'    Crack front Theta angle smoothing:       no'
    elseif(Key_Smooth_Theta_Value==1) then
      ! Only applicable for CFCP=3
      if(CFCP/=3) then
          print *, '    Error :: *Key_Smooth_Theta_Value=1 is valid only when *CFCP=3!'
          call Warning_Message('S',Keywords_Blank)
      endif
      print *,'    Crack front Theta angle smoothing:       yes'
      print *,'    Smoothing method:                        moving average method'
      print *,'    Attention :: Crack front value smoothing is active!'
      
    else
      call Warning_Message('E','Key_Smooth_Theta_Value')
    endif
  endif

  !******************************************************
  ! Crack Front Edge Growth Factor Denoising, 2022-07-14
  !******************************************************
  if (Key_Dimension ==3)then
    if(Key_Denoise_GF_Value==0) then
      print *,'    Crack front Growth Factor denoising:       no'
    elseif(Key_Denoise_GF_Value==1) then
      ! Only applicable for CFCP=3
      if(CFCP/=3) then
          print *, '    Error :: *Key_Denoise_GF_Value=1 is valid only when *CFCP=3!'
          call Warning_Message('S',Keywords_Blank)
      endif
      print *,'    Crack front GF angle denoising:       yes'
      
      print *,'    Attention :: Crack front GF angle denoising is active!'
      
    else
      call Warning_Message('E','Key_Denoise_GF_Value')
    endif
  endif

  !**********************************************************************
  ! Front edge of the crack GrowthFactor smoothing treatment, 2022-07-14
  !**********************************************************************
  if (Key_Dimension ==3)then
    if(Key_Smooth_GF_Value==0) then
      print *,'    Crack front GF smoothing:       no'
    elseif(Key_Smooth_GF_Value==1) then
      ! Only applicable for CFCP=3
      if(CFCP/=3) then
          print *, '    Error :: *Key_Smooth_GF_Value=1 is valid only when *CFCP=3!'
          call Warning_Message('S',Keywords_Blank)
      endif
      print *,'    Crack front Growth Factor smoothing:       yes'
      print *,'    Smoothing method:                          moving average method'
      
      print *,'    Attention :: Crack front Growth Factor smoothing is active!'
      
    else
      call Warning_Message('E','Key_Smooth_GF_Value')
    endif
  endif

  !*************************************************
  ! Load Control Mode Key_Force_Control, 2017-04-21
  !*************************************************
  if(Key_Force_Control==1) then
      print *, '    Control of applied force:  1-constant'
  elseif(Key_Force_Control==2) then
      print *, '    Control of applied force:  2-linear increase'
  elseif(Key_Force_Control==3) then
      print *, '    Control of applied force:  3-reduce force after propagation'
  elseif(Key_Force_Control==4) then
      ! Check: Only applicable for quasi-static analysis
      if(Key_Analysis_Type/=1) then
          print *, '    Error :: *Key_Force_Control=4 is avaliable only when *Key_Analysis_Type=1!!!'
          call Warning_Message('S',Keywords_Blank)
      endif
     ! Check: Only applicable to cohesive cracks
      if(Key_TipEnrich/=4) then
          print *, '    Error :: *Key_Force_Control=4 is avaliable only when *Key_TipEnrich=4!!!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Check: Supports only one crack at most
      if(num_Crack > 1)then
          print *, '    Error :: *Key_Force_Control=4 is avaliable only for one crack at most!!!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Temporarily not available for contact
      if(Key_Contact /=0)then
          print *, '    Error :: *Key_Force_Control=4 is not avaliable when contact is activated!!!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Inspection passed
      print *, '    Control of applied force:  4-cohesive-scheme 1'
  elseif(Key_Force_Control==5) then
      ! Check: Only applicable for quasi-static analysis
      if(Key_Analysis_Type/=1) then
          print *, '    Error :: *Key_Force_Control=5 is avaliable only when *Key_Analysis_Type=1!!!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Check: Random generation of initial cracks (natural fractures) is not supported
      if(Key_Random_NaCr == 1)then
          print *, '    Error :: *Key_Force_Control=5 is not avaliable when *Key_Random_NaCr = 1!!!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Inspection: Must contain at least one crack
      if(num_Crack == 0)then
          print *, '    Error :: *Key_Force_Control=5 is avaliable only when crack exists!!!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Temporarily not available for contact
      if(Key_Contact /=0)then
          print *, '    Error :: *Key_Force_Control=5 is not avaliable when contact is enabled!!!'
          call Warning_Message('S',Keywords_Blank)
      endif
      print *, '    Control of applied force:  5- only 1 crack propagate'
  else
      call Warning_Message('E','Key_Force_Control')
  end if

  !*****************************
  ! Life-and-Death element Test
  !*****************************
3011 FORMAT(5X,'ERROR :: One of the elements to be killed in step',I3,' does not exists!')
  if(Key_EKILL==1) then
    if(Key_Dimension==3)then
      print *, '    *Key_EKILL=1 is avaliable for only 2D problems!'
      call Warning_Message('E','Key_EKILL')
    endif
    if(Key_Analysis_Type/=1)then
      print *, '    *Key_EKILL=1 is avaliable for only 2D static problems!'
      call Warning_Message('E','Key_EKILL')
    endif
    print *, '    Elements kill:             yes'
    call log_msg('Elements kill:             yes')
    ! Check whether the elements to be killed at each step exist
    do i=1,Num_Substeps
      if(maxval(Ele_Killed_Each_Load_Step(i,:)) >= Num_Elem)then
          write(*,3011) i
          call Warning_Message('E','*Ele_Killed_Each_Load_Step')
      endif
    enddo
  endif

  !**********************************************************
  ! Normal Water Pressure in Seam Crack_Pressure, 2019-07-04
  !**********************************************************
6001 FORMAT(16X,'Inner pressure of crack ',I3, ' is ', F9.3 ,' MPa')
  if(Key_Crack_Inner_Pressure==1)then
      if (sum(abs(Crack_Pressure(1:Max_Num_Cr)))>Tol_10) then
          print *,'    Warning :: crack inner pressure defined!'
          do i_C=1,num_Crack
              write (*,6001) i_C,Crack_Pressure(i_C)/Cn_M
          enddo
          ! Only applicable for Key_Analysis_Type=1 (2D and 3D)
          if(Key_Analysis_Type/=1)then
              print *,'    Error :: Crack inner pressure is avaliable only when *Key_Analysis_Type = 1!'
              call Warning_Message('S',Keywords_Blank)
          endif
      endif
  endif

  !****************************************
  ! Gravitational acceleration, 2017-12-28
  !****************************************
  if (Key_Gravity==1)then
      if(abs(sum(g_X_Y_Z))<=Tol_20)then
          print *, '    Gravitational acceleration not defined!'
          call Warning_Message('S',Keywords_Blank)
      else
          print *, '    Gravity:                   active'
          call log_msg('Gravity:                   active')
          
          print *, '    Attention :: Gravitational acceleration is active!'
          
      endif
      ! For 2D problems, the gravitational acceleration in the y-direction should be 9.8.
      if(Key_Dimension     == 2)then
          if(g_X_Y_Z(2)<=9.0D0)then
              print *, '    WARNING :: g in y direction is not 9.8, check keyword *g_X_Y_Z!'
              
          endif
      endif
  else
      print *, '    Gravity:                   inactive'
  endif

  !*********************************************
  ! Load-displacement curve related, 2017-04-23
  !*********************************************
  if(Key_Save_f_d_Curve==1) then
      ! Check if the node number is defined
      if(f_d_Curve_node_num == -999)then
          print *, '    Error :: *Coh_f_d_Curve_node_num not defined!'
          call Warning_Message('S',Keywords_Blank)
      ! Check if the maximum node number is exceeded
      else
          if(f_d_Curve_node_num > Num_Node)then
              print *, '    Error :: *f_d_Curve_node_num is larger than Max node number!'
              call Warning_Message('S',Keywords_Blank)
          elseif(f_d_Curve_node_num <= 0) then
              print *, '    Error :: illegal *f_d_Curve_node_num!'
              call Warning_Message('S',Keywords_Blank)
          endif
      endif
      print *, '    Key_Save_f_d_Curve:        active'
      print *, '    Node number to get f-d curve: ',f_d_Curve_node_num
  elseif(Key_Save_f_d_Curve==0)then
      !do nothing
  else
      call Warning_Message('E','Key_Save_f_d_Curve')
  endif

  !*******************************************************
  ! Load-COD Curve Related, 2022-10-15. NEWFTU2022101502.
  !*******************************************************
  if(Key_Save_f_COD_Curve==1) then
      ! Check if the node number is defined
      if(f_COD_Curve_Crack_Num ==0)then
          print *, '    Error :: *f_COD_Curve_Crack_Num not defined!'
          call Warning_Message('S',Keywords_Blank)
      ! Check if the maximum node number is exceeded
      else
          if(f_COD_Curve_Crack_Num > Num_Crack)then
              print *, '    Error :: *f_d_Curve_node_num is larger than Max crack number!'
              call Warning_Message('S',Keywords_Blank)
          elseif(f_COD_Curve_Crack_Num < 0) then
              print *, '    Error :: illegal *f_d_Curve_node_num!'
              call Warning_Message('S',Keywords_Blank)
          endif
      endif
      print *, '    Key_Save_f_COD_Curve:      active'
      print *, '    Crack number to get f-COD curve: ',f_COD_Curve_Crack_Num
  elseif(Key_Save_f_COD_Curve==0)then
      !do nothing
  else
      call Warning_Message('E','Key_Save_f_COD_Curve')
  endif

  !******************************************************************
  ! Should the crack inflection point be set as a calculation point?
  !******************************************************************
  if(Key_Kink_Point==0) then
      print *, '    Kink point as fluid node:  no'
  elseif(Key_Kink_Point==1) then
      print *, '    Kink point as fluid node:  yes'
  else
      call Warning_Message('E','Key_Kink_Point')
  end if

  !******************************************************
  ! Hydraulic Fracturing Analysis (Key_Analysis_Type==3)
  !******************************************************
1110 FORMAT(5X,'Injection crack number:  ',I3)
1111 FORMAT(5X,'x coor of the entry:      ',F9.3,' m')
1112 FORMAT(5X,'y coor of the entry:      ',F9.3,' m')
  if (Key_Analysis_Type==3) then
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! For 2D hydraulic fracturing analysis, there must be an initial crack; this is being checked
      ! here.
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (num_Crack ==0 .and. Key_Dimension==2) then
          WRITE(*,5001)
          print *, '    Error :: Initial crack was not found for hydraulic fracturing analysis'
          call Warning_Message('S',Keywords_Blank)
      end if
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Hydraulic fracturing analysis does not support binary data storage.
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (Key_Data_Format == 2 ) then
          WRITE(*,5001)
          print *, '    Error :: binary data is not avaliable for hydraulic fracturing analysis'
          call Warning_Message('S',Keywords_Blank)
      end if

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Types of hydraulic fracturing analysis models (whether it is a symmetrical model or a full
      ! model)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (Key_Symm_HF == 1) then
          print *, '    HF model type:             symmetric'
      elseif (Key_Symm_HF == 0) then
          print *, '    HF model type:             full'
          ! If it is not staged fracturing, check whether the water injection point is defined at this
          ! location.
          if(Key_HF_Multistage==0 .and. Key_Dimension==2)then
              if(sum(abs(Inj_Point_Loc(1:2))) >ZR) then
                  ! Check whether the water injection point corresponding to the crack is on the initial crack
                  Inj_Poin_x = Inj_Point_Loc(1)
                  Inj_Poin_y = Inj_Point_Loc(2)
                  Cr_A =  Crack_Coor(Inject_Crack_Num,1,1:2)
                  Cr_B =  Crack_Coor(Inject_Crack_Num,2,1:2)
                  call Tool_Yes_On_Line(Inj_Poin_x,Inj_Poin_y,Cr_A,Cr_B,Yes_Inj_On_Cr)
                  if(Yes_Inj_On_Cr.eqv..True.)then
                      write(*,1110) Inject_Crack_Num
                      write(*,1111) Inj_Poin_x
                      write(*,1112) Inj_Poin_y
                  elseif(Yes_Inj_On_Cr.eqv..False.)then
                      print *, '    Error :: injection point is not on the corresponding crack line!'
                      call Warning_Message('S',Keywords_Blank)
                  end if
              else
                  print *, '    Error :: no injection point found!'
                  call Warning_Message('S',Keywords_Blank)
              end if
          endif
      else
          call Warning_Message('E','Key_Symm_HF')
      end if
      !~~~~~~~~~~~~
      ! Burst Step
      !~~~~~~~~~~~~
1231  FORMAT(5X,'Number of frac steps:    ',I3)
      write(*,1231) Num_Frac
      !~~~~~~~~~~~~~~~~~~~~~~~
      ! Enable online search?
      !~~~~~~~~~~~~~~~~~~~~~~~
      if(Key_HF_LineSearch==1)then
          print *, '    Line searching:            yes'
      elseif(Key_HF_LineSearch==0)then
          print *, '    Line searching:            no (warning:: turn on line searching!!!!)'
      else
          call Warning_Message('E','Key_HF_LineSearch')
      endif
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Inheritance Check:
      ! Whether each fracture step inherits the final water pressure, opening, etc., from the
      ! previous fracture step during the initial iteration
      ! (When the number of rupture steps > 1)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (Num_Frac>1) then
          if (Key_IniPre_PassOn == 0)then
              print *, '    Succeeded initial pres:    no'
          elseif (Key_IniPre_PassOn == 1)then
              print *, '    Succeeded initial pres:    yes'
          else
              call Warning_Message('E','Key_IniPre_PassOn')
          end if
      elseif(Num_Frac <=0) then
          call Warning_Message('E','Num_Frac')
      end if
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Have you considered the proppant and proppant diameter?
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1232  FORMAT(5X,'Diameter of proppant:    ',F9.5,'mm')
      if (Key_Proppant==1) then
          print *, '    Proppant:                  yes'
          if(Dia_of_Proppant >= 0.5D-3 .and. Dia_of_Proppant <= 20.0D-3)then
              write(*,1232) Dia_of_Proppant*1000.0D0
          elseif(Dia_of_Proppant < 1.0D-3) then
              print *, '    Diameter of proppant is too small!'
              call Warning_Message('S',Keywords_Blank)
          elseif(Dia_of_Proppant > 50.0D-3) then
              print *, '    Diameter of proppant is too big!'
              call Warning_Message('S',Keywords_Blank)
          endif
      elseif (Key_Proppant==0) then
          print *, '    Proppant:                  no'
      else
          call Warning_Message('E','Key_Proppant')
      endif
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Check if it is staged fracturing
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1343  FORMAT(5X,'Number of MS fracturing: ',I3)
1344  FORMAT(5X,'Error :: initial crack of MS fracture ',I3,' has not been defined!')
1345  FORMAT(5X,'Error :: injection point of MS fracture ',I3,' has not been defined!')
1346  FORMAT(5X,'Number of natural fractures of MS fracturing: ',I3)
      if(Key_HF_Multistage==0 .and. Key_Dimension==2)then
          print *, '    Multistage(MS) fracturing: no'
      elseif(Key_HF_Multistage==1 .and. Key_Dimension==2)then
          print *, '    Multistage(MS) fracturing: yes'
          !-----------------------------------------------------------
          ! Segmented fracturing does not support a symmetrical model
          !-----------------------------------------------------------
          if( Key_Symm_HF==1)then
              print *, '    MS fracturing do not support symmetry model!'
              call Warning_Message('S',Keywords_Blank)
          endif

          !----------------------------------------------------------------------------------
          ! Read the number of fracturing stages (obtained from the MS_Crack_Order variable)
          !----------------------------------------------------------------------------------
          do i=1,Max_MS_Num_Crack
              if(MS_Crack_Order(i)>0)then
                  MS_Crack_Num = MS_Crack_Num + 1
              elseif(MS_Crack_Order(i)<0)then
                  print *, '    Error :: Illegal value for *MS_Crack_Order'
                  call Warning_Message('S',Keywords_Blank)
              endif
          enddo
          ! Check whether the number of fractures in each stage of the segmented fracturing exceeds the total
          ! number of fractures.
          if(num_Crack < MS_Crack_Num)then
              print *, '    Error :: number of MS fracturing can not be larger than *num_Crack!'
              call Warning_Message('S',Keywords_Blank)
          endif
          write(*,1343) MS_Crack_Num
          MS_NaturalCr_Num = num_Crack - MS_Crack_Num
          write(*,1346) MS_NaturalCr_Num
          !-----------------------------------------------------------------------------------------
          ! Check whether the initial fracture coordinate data for each stage of the fracturing has
          ! been defined
          !-----------------------------------------------------------------------------------------
          do i=1,MS_Crack_Num
              c_C = MS_Crack_Order(i)
              if(sum(Crack_Coor(c_C,:,:))==0)then
                  write (*,1344) i
                  call Warning_Message('S',Keywords_Blank)
              endif
          end do
          !----------------------------------------------------------------------------------------
          ! Check whether the initial fracture water injection points for each stage of fracturing
          ! have been defined
          !----------------------------------------------------------------------------------------
          do i=1,MS_Crack_Num
              c_C = MS_Crack_Order(i)
              if(sum(MS_InP_Loc(c_C,:))==0)then
                  write (*,1345) i
                  call Warning_Message('S',Keywords_Blank)
              endif
          end do
          !--------------------------------------------------------------------------------------
          ! Check whether the initial water injection point coordinates for each stage of staged
          ! fracturing are correct (whether they are on the corresponding fractures).
          !--------------------------------------------------------------------------------------
          do i_C=1,MS_Crack_Num
              c_C = MS_Crack_Order(i_C)
              Yes_InjP_On_Cr = .False.
              Inj_P_x = MS_InP_Loc(i_C,1)
              Inj_P_y = MS_InP_Loc(i_C,2)
              do i_S = 1,Each_Cr_Poi_Num(C_C)-1
                  ! Crack fragment endpoint coordinates
                  crack_p1 = [Crack_Coor(c_C,i_S,1),Crack_Coor(c_C,i_S,2)]
                  crack_p2 = [Crack_Coor(c_C,i_S+1,1),Crack_Coor(c_C,i_S+1,2)]
                  call Tool_Yes_On_Line(Inj_P_x,Inj_P_y,crack_p1,crack_p2,Yes_InjP_On_Cr)
                  if(Yes_InjP_On_Cr.eqv..True.)exit
              end do
              if(Yes_InjP_On_Cr.eqv..False.) then
                  print *, '    Error :: injection point is not on the specified crack!'
                  call Warning_Message('S',Keywords_Blank)
              endif
          end do
          !--------------------------------------------------------------------------------------
          ! Check keyword Key_HF_MS_Contc_Check: Whether the contact iteration of fracture faces
          ! containing proppant is performed during staged fracturing
          !--------------------------------------------------------------------------------------
          if(Key_HF_MS_Contc_Check==0)then
              print *, '    *Key_HF_MS_Contc_Check:    0 (default)'
          elseif(Key_HF_MS_Contc_Check==1)then
              print *, '    Warning :: *Key_HF_MS_Contc_Check is 1 and this is not a recommended choice!!!'
          else
              Key_HF_MS_Contc_Check=0
                  print *, '    Caution :: *Key_HF_MS_Contc_Check was reset to 0 (default)!'
          endif
          !----------------------------------------------------------------------------------------
          ! Check the relevant data of special cases of staged fracturing in Paper1 (special cases
          ! in Paper1 must consider the support provided by the proppant)
          !----------------------------------------------------------------------------------------
          if(Key_Paper1==1)then
              print *,'    CAUTION :: SPECIAL ROUTINE FOR PAPER 1.'
              Key_Proppant=1
          endif
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Check the definition of the water injection point flow
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(sum(Inject_Q_Time)>0)then
          do i=2,20
              ! The time points must be arranged in ascending order.
              if(Inject_Q_Time(i) /=0)then
                  if(Inject_Q_Time(i)<Inject_Q_Time(i-1))then
                      print *, '    Injection quality time is wrong, check **Inject_Q_Time!'
                      call Warning_Message('S',Keywords_Blank)
                  endif
              endif
          enddo
          ! If a flow value is defined
          if(sum(Inject_Q_Val)>0)then
              ! Check if the traffic value is negative
              if(minval(Inject_Q_Val)<ZR)then
                  print *, '    Negative Inject_Q_Val, check *Inject_Q_Val!'
                  call Warning_Message('S',Keywords_Blank)
              endif
          else
              print *, '    Injection quality value curve is not defined, check *Inject_Q_Val!'
                  call Warning_Message('S',Keywords_Blank)
          endif
      else
          print *, '    Injection quality time curve is not defined, check *Inject_c_Time!'
          call Warning_Message('S',Keywords_Blank)
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Whether to consider the migration of proppant
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1233  FORMAT(5X,'Upper limit of c:         ',F5.2)
      if (Key_Propp_Trans==1) then
          print *, '    Proppant transport:        yes'
          ! Check the maximum proppant concentration value
          if(Max_c>ZR .and. Max_c<=ZP8) then
              write (*,1233)  Max_c
          elseif(Max_c==ZR)then
              print *, '    Upper limit of proppant concentration is not defined, check *Max_c'
              call Warning_Message('S',Keywords_Blank)
          elseif(Max_c>ZP8)then
              print *, '    Upper limit of proppant concentration is larger than 0.8, check *Max_c'
              call Warning_Message('S',Keywords_Blank)
          end if
          ! If proppant transport is considered, it is necessary to provide the curve of proppant
          ! concentration at the injection point over time.
          if(sum(Inject_c_Time)>0)then
              do i=2,20
                  ! The time points must be arranged in ascending order
                  if(Inject_c_Time(i) /=0)then
                      if(Inject_c_Time(i)<Inject_c_Time(i-1))then
                        print *, '    Proppant concentration time is wrong, check **Inject_c_Time!'
                        call Warning_Message('S',Keywords_Blank)
                      endif
                  endif
              enddo
              ! If a concentration value is defined
              if(sum(Inject_c_Val)>0)then
                  ! Check if it exceeds the maximum concentration value
                  if(maxval(Inject_c_Val)>max_c)then
                      print *, '    max(Inject_c_Val) is larger than max_c, check *Inject_c_Val!'
                      call Warning_Message('S',Keywords_Blank)
                  elseif(minval(Inject_c_Val)<ZR)then
                      print *, '    Negative Inject_c_Val, check *Inject_c_Val!'
                      call Warning_Message('S',Keywords_Blank)
                  endif
              else
                  print *, '    Proppant concentration value curve is not defined, check *Inject_c_Val!'
                  call Warning_Message('S',Keywords_Blank)
              endif
          else
              print *, '    Proppant concentration time curve is not defined, check *Inject_c_Time!'
              call Warning_Message('S',Keywords_Blank)
          endif
      elseif (Key_Propp_Trans==0) then
          print *, '    Proppant transport:        no'
      else
          call Warning_Message('E','Key_Propp_Trans')
      endif
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Check whether it is static viscosity or dynamic viscosity
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(Key_Visco_Type==1)then
          print *,'    Viscosity type:            static'
      elseif(Key_Visco_Type==2)then
          ! Dynamic viscosity is only applicable in situations where the transport of the proppant is
          ! considered (because the proppant concentration needs to be calculated)
          if(Key_Propp_Trans==0)then
              print *,'    Warning :: Dynamic viscosity is only for the case that *Key_Propp_Trans=1'
              print *,'               *Key_Visco_Type is reset to 1'
              ! Correct the calculation of static viscosity
              Key_Visco_Type=1
              print *,'    Viscosity type:            static'
          elseif(Key_Propp_Trans==1)then
              print *,'    Viscosity type:            dynamic'
              ! Check whether the dynamic viscosity parameter m is within the allowable range (1, 3)
              if(Viscosity_Par_m > THR .or.Viscosity_Par_m < ONE) then
                  print *, '    Error :: value of *Viscosity_Par_m is out of range (1,3)!'
                  call Warning_Message('S',Keywords_Blank)
              endif
              ! Check if there is a maximum allowable dynamic viscosity amplification factor Visco_Zoom_Factor;
              ! this value should be greater than 1.
              if(Visco_Zoom_Factor > 1)then
                  print *,'    Viscosity zoom factor:  ',Visco_Zoom_Factor
              else
                  print *, '    Error :: value of *Viso_Zoom_Factor must be larger than 1!'
                  call Warning_Message('S',Keywords_Blank)
              endif

          endif
      elseif(Key_Visco_Type==3)then
          print *,'    Viscosity type:            time-dependent'
          print *,'    Viscosity_td_m:            ',Viscosity_td_m
          print *,'    Viscosity_td_n:            ',Viscosity_td_n
      else
          call Warning_Message('E','Key_Propp_Trans')
      endif
      !~~~~~~
      ! Leak
      !~~~~~~
1243  FORMAT(5X,'Leak-off coefficient:     'F9.5,' ms^(1/2)')
      if (Key_Leakoff==1) then
          print *, '    Leak-off:                  yes'
          ! Leakage coefficient
          write (*,1243)  Coeff_Leak
      elseif (Key_Leakoff==0) then
          print *, '    Leak-off:                  no'

      else
          call Warning_Message('E','Key_Leakoff')
      endif
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! delta_Time Calculation Methods (Rectangular Method and Trapezoidal Method)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (Key_Cal_deltaTime == 1) then
          print *, '    Method to get delta time:  rectangle method'
      elseif (Key_Cal_deltaTime == 2) then
          print *, '    Method to get delta time:  trapezoidal method'
      else
          call Warning_Message('E','Key_Cal_deltaTime')
      end if
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Water Injection Crack Number
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1133  FORMAT(5X,'Crack number of injection:',I3.1)
      ! If it is not staged fracturing, check whether the water injection point at that location
      ! corresponds to the fracture.
      if(Key_HF_Multistage==0 .and. Key_Dimension==2)then
          if ((Inject_Crack_Num <= num_Crack) .and. (Inject_Crack_Num >0))    then
              write(*,1133) Inject_Crack_Num
              ! Check whether the water injection point is on the injection crack
              if(Key_Symm_HF==0)then
                  c_C=Inject_Crack_Num
                  Yes_InjP_On_Cr = .False.
                  do i_S = 1,Each_Cr_Poi_Num(c_C)-1
                      ! Crack fragment endpoint coordinates
                      crack_p1 = [Crack_Coor(c_C,i_S,1),Crack_Coor(c_C,i_S,2)]
                      crack_p2 = [Crack_Coor(c_C,i_S+1,1),Crack_Coor(c_C,i_S+1,2)]
                      Inj_P_x=Inj_Point_Loc(1)
                      Inj_P_y=Inj_Point_Loc(2)
                      call Tool_Yes_On_Line(Inj_P_x,Inj_P_y,crack_p1,crack_p2,Yes_InjP_On_Cr)
                      if(Yes_InjP_On_Cr)exit
                  end do
                  if(Yes_InjP_On_Cr.eqv..False.) then
                      print *, '    Error :: injection point is not on the specified crack!'
                      call Warning_Message('S',Keywords_Blank)
                  endif
              endif
          else
              call Warning_Message('E','Inject_Crack_Num')
          end if
      end if
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Check whether the maximum rupture steps exceed the program's allowance
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (Num_Frac > Max_Num_Frac)then
          WRITE(*,5001)
          print *, '    Max_Num_Frac is too big！'
          call Warning_Message('E','Max_Num_Frac')
      end if
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Cracks containing water at the initial moment
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1151  FORMAT(5X,'Num of fluid-driven crack:',I3)
      do i_C =  1,num_Crack
          if(Cracks_HF_State(i_C)==1) then
              write(*,1151) i_C
          endif
      end do
      if(sum(Cracks_HF_State(1:num_Crack))>=2) then
          print *, '    Warning :: more than one crack have initial water!'
          
      endif
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Hydraulic Fracturing Contact Iteration Settings
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(Key_HF_Cont_Scheme==1)then
          print *, '    Contact scheme for HF:     only one time for each time step.'
      elseif(Key_HF_Cont_Scheme==2)then
          print *, '    Contact scheme for HF:     check contact at the first two Fluid-solid iterations (default)!'
      elseif(Key_HF_Cont_Scheme==3)then
          print *, '    Contact scheme for HF:     check contact at every Fluid-solid iteration!'
      else
          call Warning_Message('E','Key_HF_Cont_Scheme')
      endif
  end if


  !******************************************************************************
  ! Hydraulic fracturing analysis, fresh water fracturing (Key_Analysis_Type==4)
  !******************************************************************************
  if (Key_Analysis_Type==4) then
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! 2D hydraulic fracturing analysis must have an initial crack; a check is performed here.
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (num_Crack ==0 .and. Key_Dimension==2) then
          WRITE(*,5001)
          print *, '    Error :: Initial crack was not found for hydraulic fracturing analysis'
          call Warning_Message('S',Keywords_Blank)
      end if

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Types of hydraulic fracturing analysis models (whether it is a symmetrical model or a full
      ! model)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (Key_Symm_HF == 1) then
          print *, '    HF model type:             symmetric'
      elseif (Key_Symm_HF == 0) then
          print *, '    HF model type:             full'
      else
          call Warning_Message('E','Key_Symm_HF')
      end if
      !~~~~~~~~~~~~
      ! Burst Step
      !~~~~~~~~~~~~
      write(*,1231) Num_Frac
      !~~~~~~~~~~~~~~~~~~~~~~~
      ! Enable online search?
      !~~~~~~~~~~~~~~~~~~~~~~~
      if(Key_HF_LineSearch==1)then
          print *, '    Line searching:            yes'
      elseif(Key_HF_LineSearch==0)then
          print *, '    Line searching:            no (warning:: turn on line searching!!!!)'
      else
          call Warning_Message('E','Key_HF_LineSearch')
      endif
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Inheritance Check:
      ! Whether each fracture step inherits the final water pressure, opening, etc., from the
      ! previous fracture step during the initial iteration
      ! (When the number of rupture steps is greater than 1)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (Num_Frac>1) then
          if (Key_IniPre_PassOn == 0)then
              print *, '    Succeeded initial pres:    no'
          elseif (Key_IniPre_PassOn == 1)then
              print *, '    Succeeded initial pres:    yes'
          else
              call Warning_Message('E','Key_IniPre_PassOn')
          end if
      elseif(Num_Frac <=0) then
          call Warning_Message('E','Num_Frac')
      end if
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Whether to consider the proppant and proppant diameter
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1732  FORMAT(5X,'Diameter of proppant:    ',F9.5,'mm')
      if (Key_Proppant==1) then
          print *, '    Proppant:                  yes'
          if(Dia_of_Proppant >= 0.5D-3 .and. Dia_of_Proppant <= 20.0D-3)then
              write(*,1732) Dia_of_Proppant*1000.0D0
          elseif(Dia_of_Proppant < 1.0D-3) then
              print *, '    Diameter of proppant is too small!'
              call Warning_Message('S',Keywords_Blank)
          elseif(Dia_of_Proppant > 50.0D-3) then
              print *, '    Diameter of proppant is too big!'
              call Warning_Message('S',Keywords_Blank)
          endif
      elseif (Key_Proppant==0) then
          print *, '    Proppant:                  no'
      else
          call Warning_Message('E','Key_Proppant')
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Check the definition of the water injection point flow
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(sum(Inject_Q_Time)>0)then
          do i=2,20
              ! The time points must be arranged in ascending order
              if(Inject_Q_Time(i) /=0)then
                  if(Inject_Q_Time(i)<Inject_Q_Time(i-1))then
                      print *, '    Injection quality time is wrong, check *Inject_Q_Time!'
                      call Warning_Message('S',Keywords_Blank)
                  endif
              endif
          enddo
          ! If a flow value is defined
          if(sum(Inject_Q_Val)>0)then
              ! Check if the traffic value is negative
              if(minval(Inject_Q_Val)<ZR)then
                  print *, '    Negative Inject_Q_Val, check *Inject_Q_Val!'
                  call Warning_Message('S',Keywords_Blank)
              endif
          else
              print *, '    Injection quality value curve is not defined, check *Inject_Q_Val!'
                  call Warning_Message('S',Keywords_Blank)
          endif
      else
          if(Key_Dimension==2) then
            print *, '    Injection quality time curve is not defined, check *Inject_c_Time!'
            call Warning_Message('S',Keywords_Blank)
          endif
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Water Injection Crack Number
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1733  FORMAT(5X,'Crack number of injection:',I3.1)
      ! If it is not staged fracturing, check whether the water injection point at that location
      ! corresponds to the fracture.
      if(Key_HF_Multistage==0 .and. Key_Dimension==2)then
          if ((Inject_Crack_Num <= num_Crack) .and. (Inject_Crack_Num >0))    then
              write(*,1733) Inject_Crack_Num
              ! Check whether the water injection point is on the injection crack
              if(Key_Symm_HF==0)then
                  c_C=Inject_Crack_Num
                  Yes_InjP_On_Cr = .False.
                  do i_S = 1,Each_Cr_Poi_Num(c_C)-1
                      ! Crack fragment endpoint coordinates
                      crack_p1 = [Crack_Coor(c_C,i_S,1),Crack_Coor(c_C,i_S,2)]
                      crack_p2 = [Crack_Coor(c_C,i_S+1,1),Crack_Coor(c_C,i_S+1,2)]
                      Inj_P_x=Inj_Point_Loc(1)
                      Inj_P_y=Inj_Point_Loc(2)
                      call Tool_Yes_On_Line(Inj_P_x,Inj_P_y,crack_p1,crack_p2,Yes_InjP_On_Cr)
                      if(Yes_InjP_On_Cr)exit
                  end do
                  if(Yes_InjP_On_Cr.eqv..False.) then
                      print *, '    Error :: injection point is not on the specified crack!'
                      call Warning_Message('S',Keywords_Blank)
                  endif
              endif
          else
              call Warning_Message('E','Inject_Crack_Num')
          end if
      end if
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Check the definition of the water injection point pressure
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(sum(Inject_P_Time)>0)then
          do i=2,20
              ! The time points must be arranged in ascending order
              if(Inject_P_Time(i) /=0)then
                  if(Inject_P_Time(i)<Inject_P_Time(i-1))then
                      print *, '    Injection pressure time is wrong, check **Inject_P_Time!'
                      call Warning_Message('S',Keywords_Blank)
                  endif
              endif
          enddo
          ! If a traffic value is defined
          if(sum(Inject_P_Val)>0)then
              ! Check if the traffic value is negative
              if(minval(Inject_P_Val)<ZR)then
                  print *, '    Negative Inject_P_Val, check *Inject_P_Val!'
                  call Warning_Message('S',Keywords_Blank)
              endif
          else
              print *, '    Injection pressure value curve is not defined, check *Inject_P_Val!'
                  call Warning_Message('S',Keywords_Blank)
          endif
      else
          if(Key_Dimension==2) then
              print *, '    Injection pressure time curve is not defined, check *Inject_P_Time!'
              call Warning_Message('S',Keywords_Blank)
          endif
      endif
  end if

  !**********************************************************************************************
  ! Apply viscosity-dominated theoretical water pressure (analysis type 5), supporting up to one
  ! fracture
  !**********************************************************************************************
  if (Key_Analysis_Type==5) then
      ! The test procedure applying pure water pressure must have an initial crack, which is checked here.
      if (num_Crack ==0) then
          WRITE(*,5001)
          print *, '    Error :: Initial crack was not found for test program for analysis type 4 or 5!'
          call Warning_Message('S',Keywords_Blank)
      end if
      if (num_Crack >1) then
          WRITE(*,5001)
          print *, '    Error :: analysis type 5 only supports one crack!'
          call Warning_Message('S',Keywords_Blank)
      end if
  end if

  !******************************************************************
  ! Post-processing related: Do not save any data. NEWFTU2022090601.
  !******************************************************************
  if (Key_Save_Nothing ==1) then
      print *, '    Save nothing:              yes'
      print *, '    Warning :: nothing will be saved for post-processing!'
      print *, '               Keyword *Key_Save_Nothing=1!'
      
  elseif(Key_Save_Nothing ==0) then
      print *, '    Save nothing:               no'
  else
      call Warning_Message('E','Key_Save_Nothing')
  end if

  !**************************************************************************
  ! Post-processing related: whether to perform streamlined post-processing.
  !**************************************************************************
  if (Key_Simple_Post ==1) then
      print *, '    Simple post:               yes'
      
      print *, '    Warning :: only basic results will be saved for post-processing!'
      print *, '               Keyword *Key_Simple_Post=1!'
      
  elseif(Key_Simple_Post ==0) then
      print *, '    Simple post:               no'
  else
      call Warning_Message('E','Key_Simple_Post')
  end if

  !*************************************************************************
  ! Post-processing related: Calculating fracture permeability. 2022-11-26.
  !*************************************************************************
  if (key_Dimension==3) then
    if (Key_Get_Permeability ==1) then
          print *, '    Warning :: Permeability of cracks will be calculated!'
          print *, '               Keyword *Key_Get_Permeability = 1!'
    elseif(Key_Get_Permeability ==0) then
          print *, '    Warning :: Permeability of cracks will not be calculated!'
          print *, '               Keyword *Key_Get_Permeability = 0!'
    else
          call Warning_Message('E','Key_Get_Permeability')
    end if
  elseif(key_Dimension==2) then
    if (Key_Get_Permeability ==1) then
          print *, '    ERROR :: Key_Get_Permeability =1 for 3D only'
          call Warning_Message('S',Keywords_Blank)
    end if
  endif

  !************************************************************************************************
  ! Post-processing related: Whether to save the number of Gauss points for each element. Used for
  ! post-processing. 2022-07-16.
  !************************************************************************************************
  if (Key_Post_Elements_Gauss_Num ==1) then
      ! Currently only applicable to 3D.
      if(Key_Dimension/=3)then
          print *, '    Warning :: Key_Post_Elements_Gauss_Num =1 for 3D only!'
          
      endif
      ! Currently, it is only applicable to analysis types Key_Analysis_Type = 4, =1, and =61.
      if(Key_Analysis_Type/=4 .and. Key_Analysis_Type/=1 .and. Key_Analysis_Type/=61)then
          print *, '    Error :: Key_Post_Elements_Gauss_Num =1 is valid only when Key_Analysis_Type =1 or 4!'
          call Warning_Message('S',Keywords_Blank)
      endif
      print *, '    Save Gauss num of elememts: yes'
      
      print *, '    Warning :: Gauss num of ELes will be saved for post-processing!'
      print *, '               Keyword *Key_Post_Elements_Gauss_Num =1!'
      
  elseif(Key_Post_Elements_Gauss_Num ==0) then
      print *, '    Save Gauss num of elememts:no'
  else
      call Warning_Message('E','Key_Post_Elements_Gauss_Num')
  end if

  !***********************************************************
  ! Post-processing related: Whether to calculate node stress
  !***********************************************************
  if (Key_Post_CS_N_Strs ==0) then
      print *, '    Calculate stress:          no'
      
      print *, '    Warning :: PhiPsi will not computes stress of nodes!'
      print *, '               Keyword *key_post_cs_n_strs=0!'
      
  elseif(Key_Post_CS_N_Strs ==1) then
      print *, '    Calculate stress:          yes'
  else
      call Warning_Message('E','Key_Post_CS_N_Strs')
  end if

  !************************************************************************
  ! Post-processing related: Whether to calculate nodal strain, 2021-09-10
  !************************************************************************
  if (Key_Post_CS_N_Stra ==0) then
      print *, '    Calculate strain:          no'
  elseif(Key_Post_CS_N_Strs ==1) then
      print *, '    Calculate strain:          yes'
      
      print *, '    Warning :: PhiPsi will computes strain of nodes!'
      print *, '               Keyword *key_post_cs_n_stra=1!'
      
      ! For 3D use only
      if(Key_Dimension /=3) then
          print *, '    Error :: *key_post_cs_n_stra=1 is available only for 3D problems!'
          call Warning_Message('S',Keywords_Blank)
      endif
  else
      call Warning_Message('E','Key_Post_CS_N_Stra')
  end if

  !**************************************************************************************************
  ! Post-processing related: Whether to calculate node displacement (cylindrical coordinate system),
  ! 2021-09-11
  !**************************************************************************************************
  if (Key_Post_CS_N_Disp_Cylindr ==0) then
  elseif(Key_Post_CS_N_Disp_Cylindr ==1) then
      print *, '    Disppalcement (cylind):    yes'
      
      print *, '    Warning :: PhiPsi will computes displacement of nodes in cylindrical CS!'
      print *, '               Keyword *Key_Post_CS_N_Disp_Cylindr=1!'
      
      ! For 3D use only
      if(Key_Dimension /=3) then
          print *, '    Error :: *Key_Post_CS_N_Disp_Cylindr=1 is available only for 3D problems!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! To calculate strains in cylindrical coordinates (such as circumferential strain), the following
      ! parameters need to be provided
      ! Post_cylinder_Coor_Center(1:3) = [0.0,0.0,0.0] !Origin of the cylindrical material coordinate
      ! system
      ! Post_cylinder_Coor_Vector_z(1:3) = [0.0,0.0,1.0] ! Z-axis vector of the cylindrical material
      ! coordinate system
      if (abs(Post_cylinder_Coor_Center(1)) > Con_Big_15) then
          print *, '    Error :: *Post_cylinder_Coor_Center is not available!'
          call Warning_Message('S',Keywords_Blank)
      endif
      if (abs(Post_cylinder_Coor_Vector_x(1)) > Con_Big_15 .or. &
          abs(Post_cylinder_Coor_Vector_y(1)) > Con_Big_15 .or. &
          abs(Post_cylinder_Coor_Vector_z(1)) > Con_Big_15) then
          print *, '    Error :: *Post_cylinder_Coor_Vector_x is not available !'
          call Warning_Message('S',Keywords_Blank)
      endif
  else
      call Warning_Message('E','Key_Post_CS_N_Disp_Cylindr')
  end if

  !********************************************************************************************
  ! Post-processing related: Whether to calculate node stress (cylindrical coordinate system),
  ! 2021-09-11
  !********************************************************************************************
  if (Key_Post_CS_N_Strs_Cylindr ==0) then
  elseif(Key_Post_CS_N_Strs_Cylindr ==1) then
      print *, '    Stress (cylindrical CS):   yes'
      
      print *, '    Warning :: PhiPsi will computes stress of nodes in cylindrical CS!'
      print *, '               Keyword *Key_Post_CS_N_Strs_Cylindr=1!'
      
      ! For 3D use only
      if(Key_Dimension /=3) then
          print *, '    Error :: *Key_Post_CS_N_Strs_Cylindr=1 is available only for 3D problems!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! To calculate strains in cylindrical coordinates (such as circumferential strain), the following
      ! parameters need to be provided
      ! Post_cylinder_Coor_Center(1:3) = [0.0,0.0,0.0] !Origin of the cylindrical material coordinate
      ! system
      ! Post_cylinder_Coor_Vector_z(1:3) = [0.0,0.0,1.0] ! Z-axis vector of the cylindrical material
      ! coordinate system
      if (abs(Post_cylinder_Coor_Center(1)) > Con_Big_15) then
          print *, '    Error :: *Post_cylinder_Coor_Center is not available !'
          call Warning_Message('S',Keywords_Blank)
      endif
      if (abs(Post_cylinder_Coor_Vector_x(1)) > Con_Big_15 .or.   &
          abs(Post_cylinder_Coor_Vector_y(1)) > Con_Big_15 .or.   &
          abs(Post_cylinder_Coor_Vector_z(1)) > Con_Big_15) then
          print *, '    Error :: *Post_cylinder_Coor_Vector_x is not available !'
          call Warning_Message('S',Keywords_Blank)
      endif
  else
      call Warning_Message('E','Key_Post_CS_N_Strs_Cylindr')
  end if

  !********************************************************************************************
  ! Post-processing related: Whether to calculate node strain (cylindrical coordinate system),
  ! 2021-09-10
  !********************************************************************************************
  if (Key_Post_CS_N_Stra_Cylindr ==0) then
  elseif(Key_Post_CS_N_Stra_Cylindr ==1) then
      print *, '    Strain (cylindrical CS):   yes'
      
      print *, '    Warning :: PhiPsi will computes strain of nodes in cylindrical CS!'
      print *, '               Keyword *Key_Post_CS_N_Stra_Cylindr=1!'
      
      ! For 3D use only
      if(Key_Dimension /=3) then
          print *, '    Error :: *Key_Post_CS_N_Stra_Cylindr=1 is available only for 3D problems!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! To calculate strains in cylindrical coordinates (such as circumferential strain), the following
      ! parameters need to be provided
      ! Post_cylinder_Coor_Center(1:3) = [0.0,0.0,0.0] !Origin of the cylindrical material coordinate
      ! system
      ! Post_cylinder_Coor_Vector_z(1:3) = [0.0,0.0,1.0] ! Z-axis vector of the cylindrical material
      ! coordinate system
      if (abs(Post_cylinder_Coor_Center(1)) > Con_Big_15) then
          print *, '    Error :: *Post_cylinder_Coor_Center is not available !'
          call Warning_Message('S',Keywords_Blank)
      endif
      if (abs(Post_cylinder_Coor_Vector_x(1)) > Con_Big_15 .or.   &
          abs(Post_cylinder_Coor_Vector_y(1)) > Con_Big_15 .or.   &
          abs(Post_cylinder_Coor_Vector_z(1)) > Con_Big_15) then
          print *, '    Error :: *Post_cylinder_Coor_Vector_x is not available !'
          call Warning_Message('S',Keywords_Blank)
      endif
  else
      call Warning_Message('E','Key_Post_CS_N_Stra_Cylindr')
  end if

  !**************************************************
  ! Keyword Key_Initiation (Generate Initial Cracks)
  !**************************************************
  if (Key_Initiation .eq. 0)then
      print *, '    Crack initiation:          not allowed'
  elseif (Key_Initiation .eq. 1)then
      print *, '    Crack initiation:          allowed'
      ! Note, Key_Initiation is only used when num_Crack=0
      if(num_Crack>=1)then
          print *, '    Error :: *Key_Initiation=1 is available only when *num_Crack=0!'
          call Warning_Message('S',Keywords_Blank)
      endif
      !2024-06-23.
      if(Key_Dimension==2)then
            if(Key_Analysis_Type /=1 .and. Key_Analysis_Type /=7) then
            print *, '    Error :: for 2D problems, *Key_Initiation=1 is valid only when Key_Analysis_Type=1 or 7!'
            call Warning_Message('S',Keywords_Blank)
        endif
      endif
        
      if (Key_Ini_Rule==1)then
        print *, '    Crack initiation rule:     max tensile stress stress'
      elseif (Key_Ini_Rule==2)then
        print *, '    Crack initiation rule:     max shear stress stress'
      elseif (Key_Ini_Rule==3)then
        print *, '    Crack initiation rule:     concrete'
        !2024-06-23.
        if(Key_Dimension==2)then
            print *, '    Error :: *Key_Ini_Rule=3 for 3D only!'
            call Warning_Message('S',Keywords_Blank)
        endif
      else
        call Warning_Message('E','Key_Ini_Rule')
      endif
      if(Size_Ini_Crack>Tol_11)then
          print *, '    Size of initial crack:    ',Size_Ini_Crack
          !Chech if *Size_Ini_Crack is reasonable or not.
          if(Size_Ini_Crack>Max_Model_Range) then
              print *, '    Error :: *Size_Ini_Crack too large!'
              print *, '             Max_Model_Range is',Max_Model_Range
              print *, '             Min_Model_Range is',Min_Model_Range
              call Warning_Message('S',Keywords_Blank)
          endif
      else
          print *, '    Error :: *Size_Ini_Crack not defined or illegal!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! 3D problems require specifying the type of initial crack to be generated.
      if(Key_Dimension==3)then
          if(Key_Ini_Cr_3D_Type==1)then
              print *, '    Type of initial crack:     circle'
          elseif(Key_Ini_Cr_3D_Type==2)then
              print *, '    Error :: *Key_Ini_Cr_3D_Type=2 not avaliable yet!'
              call Warning_Message('S',Keywords_Blank)
          else
              print *, '    Error :: *Key_Ini_Cr_3D_Type not defined or illegal!'
              call Warning_Message('S',Keywords_Blank)
          endif

      endif
  else
      call Warning_Message('E','Key_Initiation')
  end if

  !*************************************************
  ! Keyword Key_Propagation, Crack Growth Criterion
  !*************************************************
  2821 FORMAT(5X,'Crack growth factor:   ',T30,F8.3)
  2822 FORMAT(5X,'Crack growth length:   ',T30,F8.3,' m')
  2823 FORMAT(5X,'Max theta of Schollmann:   ',T30,F8.3,' degrees')
  !1027 FORMAT(5X,70('-'))
  if (Key_Propagation .eq. 0)then
      print *, '    Warning :: Crack propagation:         not allowed'
  elseif (Key_Propagation .eq. 1)then
      print *, '    Crack propagation:         allowed'
      ! Keyword CFCP
      if (CFCP .eq. 1)then
          if(Key_Dimension==3)then
              print *, '    Error :: CFCP=1 for 2D problems only!'
              call Warning_Message('E','CFCP')
          else
              print *, '    Crack growth criterion:    maximum circumferential stress criterion'
          endif
      elseif (CFCP .eq. 2)then
          if(Key_Ave_Stress==1) then
              print *, '    Crack growth criterion:    maximum principal tensile stress criterion (SHI)'
          elseif(Key_Ave_Stress==2) then
              ! For 2D only
              if(Key_Dimension==3)then
                print *, '    Error :: *Key_Ave_Stress=2 for 2D problems only!'
                call Warning_Message('E','Key_Ave_Stress')
              endif
              print *, '    Crack growth criterion:    maximum principal tensile stress criterion (Ordinary)'
          ! Do not use weighting; calculate stress at a single point
          elseif(Key_Ave_Stress==3) then
              ! Only supports 3D
              if(Key_Dimension==2)then
                print *, '    Error :: *Key_Ave_Stress=3 for 3D problems only!'
                call Warning_Message('E','Key_Ave_Stress')
              endif
              print *, '    Crack growth criterion:    maximum principal tensile stress criterion (single point)'
          ! Do not use weighting; take the average over a finite number of points within the sphere
          ! (2021-08-24)
          elseif(Key_Ave_Stress==4) then
              ! Only supports 3D
              if(Key_Dimension==2)then
                print *, '    Error :: *Key_Ave_Stress=4 for 3D problems only!'
                call Warning_Message('E','Key_Ave_Stress')
              endif
              print *, '    Crack growth criterion:    maximum principal tensile stress criterion (ball 9 points)'
          else
              call Warning_Message('E','Key_Ave_Stress')
          endif
          if(Key_Ave_Stress==3 .or.Key_Ave_Stress==4) then
            
            print *, '    Warning :: pay attention to  *Key_Ave_Stress!'
            
          endif
      elseif (CFCP .eq. 3)then
          if(Key_Dimension==2)then
              print *, '    Error :: CFCP=3 for 3D problems only!'
              call Warning_Message('E','CFCP')
          else
              print *,'    Crack growth criterion:    Schollmann criterion'
          endif
          ! Schollm_Max_Theta.! Maximum Theta deflection angle allowed by the 3D Schollmann’s criterion (in
          ! degrees, recommended to be less than 75, default is 55). NEWFTU2022071001.
          write(*,2823) Schollm_Max_Theta
          if (Schollm_Max_Theta < -Tol_3) then
              print *, '    Error :: cannot < 0!'
              call Warning_Message('E','Schollm_Max_Theta')
          elseif (Schollm_Max_Theta > 75.0D0) then
              print *, '    Error :: cannot > 75.0!'
              call Warning_Message('E','Schollm_Max_Theta')
          endif
      ! 2024-06-24. Maximum shear stress fracture propagation criterion. Maximum shear stress fracture
      ! propagation criterion. NEWFTU2024062401.
      elseif (CFCP .eq. 4)then
          if(Key_Dimension==3)then
              print *, '    Error :: CFCP=4 for 2D problems only!'
              call Warning_Message('E','CFCP')
          else
              print *,'    Crack growth criterion:    maximum shear stress criterion'
              print *,'    Attention :: The shear strength needs to be defined using *MATERIAL_PARA_ADDED_n!'
          endif
      else
          call Warning_Message('E','CFCP')
      end if

      !Check Factor_Propagation.
      write (*,2821) Factor_Propagation
      if (Factor_Propagation<=ONE) then
        
        print *, '    Warning :: *Factor_Propagation is too small!'
        
      elseif (Factor_Propagation >=TWO) then
        
        print *, '    Warning :: *Factor_Propagation is too big!'
        
      endif

      !Check Propagation_Length (2021-07-27).
      if (Propagation_Length > Tol_10) then
        write (*,2822) Propagation_Length
        
        print *, '    Warning :: a specified propagation length  is given!'
        print *, '               Keyword *Factor_Propagation will not take effect any more!'
        
      endif
  else
      call Warning_Message('E','Key_Propagation')
  end if

  !****************************************************
  ! Key_InPlane_Growth, in-plane expansion, 2022-04-18
  !****************************************************
  if(Key_InPlane_Growth==1) then
      if(Key_Dimension==2)then
          print *, '    Warning :: Key_InPlane_Growth=1 for 3D only!'
          !call Warning_Message('E','Key_InPlane_Growth')
      elseif(Key_Dimension==3)then
          print *, '    Warning :: Key_InPlane_Growth=1 is specified!'
          print *, '               All cracks propagate in plane!'
      endif
      
    !2024-05-02.  
    if(Key_Scheme_Signed_Dis_InPlane==2) then
        print *, '    Warning :: Key_Scheme_Signed_Dis_InPlane=2 is specified!'
        print *, '               This scheme to calculate signed distance is fast but not robust!'
    elseif(Key_Scheme_Signed_Dis_InPlane==3) then
        print *, '    Warning :: Key_Scheme_Signed_Dis_InPlane=3 is specified!'
        print *, '               This scheme to calculate signed distance is fast but not robust!'
    elseif(Key_Scheme_Signed_Dis_InPlane==1)then
        !do nothing
    else
          call Warning_Message('E','Key_Scheme_Signed_Dis_InPlane')
    endif
  elseif(Key_InPlane_Growth==0)then
      !do nothing
  else
      call Warning_Message('E','Key_InPlane_Growth')
  endif

  !***********************************************************************
  ! Key_Stop_Outside_Crack, model exterior crack termination. 2022-10-02.
  !***********************************************************************
  if(Key_Stop_Outside_Crack==1) then
      if(Key_Dimension==2)then
          print *, '    Error :: Key_Stop_Outside_Crack=1 for 3D only!'
          call Warning_Message('E','Key_Stop_Outside_Crack')
      endif
      print *, '    Key_Stop_Outside_Crack:         1'
      
  elseif(Key_Stop_Outside_Crack==0)then
      print *, '    Key_Stop_Outside_Crack:         0'
  else
      call Warning_Message('E','Key_Stop_Outside_Crack')
  endif

  !**********************************************************************************
  ! Key_3D_Cr_Update_Scheme, 3D fracture surface update algorithm. NEWFTU2022071201.
  !**********************************************************************************
  if(Key_Dimension==3)then
      if(Key_3D_Cr_Update_Scheme==1) then
          print *, '    Crack surface update scheme:    schmene 1'
      elseif(Key_3D_Cr_Update_Scheme==2)then
          print *, '    Crack surface update scheme:    schmene 2 (default)'
      else
          call Warning_Message('E','Key_3D_Cr_Update_Scheme')
      endif
  endif

  !************************************************************
  ! Is the allowable crack propagation angle within the range?
  !************************************************************
  3001 FORMAT(5X,'Prop_Angle_Alowed:       ',E13.5,' degrees')
  if(Prop_Angle_Alowed<ZR .or. Prop_Angle_Alowed>180.0D0)then
      print *, '    Error :: Range of Prop_Angle_Alowe is (0,180)!'
      call Warning_Message('E','Prop_Angle_Alowed')
  else
      write(*,3001) Prop_Angle_Alowed
  endif

  !*********************************************************************
  ! Keyword Key_SIFs_Method, Stress Intensity Factor Calculation Method
  !*********************************************************************
  if (Key_SIFs_Method .eq. 1)then
      print *, '    SIFs:                      displacement interpolation method'
  elseif (Key_SIFs_Method .eq. 2)then
      ! For 2D only
      if(Key_Dimension== 3)then
          print *, '    ERROR :: Key_SIFs_Method==2 for 2D problems only!'
          call Warning_Message('S',Keywords_Blank)
      endif
      print *, '    SIFs:                      interaction integral method'
  else
      call Warning_Message('E','Key_SIFs_Method')
  end if


  if(Key_K_Sparse==0)then
      if     (Key_SLOE .eq. 1)then
          print *, '    Linear solver:             direct solver'
          call log_msg('Linear solver:             direct solver')
      elseif (Key_SLOE .eq. 2)then
          print *, '    SLOE:                      Gauss method'
          call log_msg('Linear solver:             Gauss method')
      elseif (Key_SLOE .eq. 3)then
          print *, '    SLOE:                      Pardiso'
          call log_msg('Linear solver:             Pardiso')
      elseif (Key_SLOE .eq. 4)then
          print *, '    SLOE:                      ITPACK'
          call log_msg('Linear solver:             ITPACK')
          
      elseif (Key_SLOE .eq. 5)then
          print *, '    SLOE:                      LAPACK'
          call log_msg('Linear solver:             LAPACK')
      elseif (Key_SLOE .eq. 6)then
          print *, '    SLOE:                      MUMPS'
          call log_msg('Linear solver:             MUMPS')
          !call log_msg('Error :: Solver 6 (MUMPS) is no longer supported!')
          
      elseif (Key_SLOE .eq. 7)then
          print *, '    SLOE:                      UFMPACK'
          call log_msg('Linear solver:             UFMPACK')
      elseif (Key_SLOE .eq. 8)then
          print *, '    SLOE:                      Lis (openMP)'
          call log_msg('Linear solver:             Lis (openMP)')
      elseif (Key_SLOE .eq. 9)then
          print *, '    SLOE:                      SuperLU'
          call log_msg('Linear solver:             SuperLU')
      elseif (Key_SLOE .eq. 10)then
          print *, '    SLOE:                      y12m'
          call log_msg('Linear solver:             y12m')
          !call log_msg('Error :: Solver 10 (y12m) is no longer supported!')
          
      elseif (Key_SLOE .eq. 11)then
          print *, '    SLOE:                      PCG'
          call log_msg('Linear solver:             PCG')
#ifdef sffortran
      elseif (Key_SLOE .eq. 12)then
          print *, '    SLOE:                      STRUMPACK'
          call log_msg('Linear solver:             STRUMPACK')
#endif          
      else
          call Warning_Message('E','Key_SLOE')
      end if
  ! Sparse matrix storage only supports solvers 3, 6, 7, and 8.
  elseif(Key_K_Sparse==1)then
      if (Key_SLOE .eq. 3)then
          print *, '    SLOE:                      Pardiso'
          call log_msg('Linear solver:             Pardiso')
      elseif (Key_SLOE .eq. 6)then
          print *, '    SLOE:                      MUMPS'
          call log_msg('Linear solver:             MUMPS')
          !call log_msg('Error :: Solver 6 (MUMPS) is no longer supported!')
          
      elseif (Key_SLOE .eq. 7)then
          print *, '    SLOE:                      UFMPACK'
          call log_msg('Linear solver:             UFMPACK')
      elseif (Key_SLOE .eq. 8)then
          print *, '    SLOE:                      Lis'
          call log_msg('Linear solver:             Lis')
      elseif (Key_SLOE .eq. 9)then
          print *, '    SLOE:                      SuperLU'
          call log_msg('Linear solver:             SuperLU')
      elseif (Key_SLOE .eq. 12)then
          print *, '    SLOE:                      STRUMPACK'
          call log_msg('Linear solver:             STRUMPACK')
      else
          print *,'    Error :: sparse K supports only 3,6,7,8,9,12 solvers!'
          call log_msg('Error :: sparse K supports only 3,6,7,8,9,12 solvers!')
          call Warning_Message('E','Key_SLOE')
      end if
  endif
  
  ! ifort and ifx do not support the SuperLU solver. 2024-03-08.
#ifdef ifort
  if (Key_SLOE .eq. 9)then
    !call Warning_Message('E','Key_SLOE')
    print *, '    Warning :: Intel ifort compiler does not support SuperLU solver!'
    print *, '               The solver will be switched to Pardiso!' 
    call log_msg('Warning :: Intel ifort compiler does not support SuperLU solver!')
    call log_msg('The solver will be switched to Pardiso!')
    Key_SLOE = 3
  endif
#endif
#ifdef ifx
    if (Key_SLOE .eq. 9)then
        !call Warning_Message('E','Key_SLOE')
        print *, '    Warning :: Intel ifx compiler does not support SuperLU solver!'
        print *, '               The solver will be switched to Pardiso!'
        call log_msg('Warning :: Intel ifx compiler does not support SuperLU solver!')
        call log_msg('The solver will be switched to Pardiso!')
        Key_SLOE = 3
    endif
#endif
#ifdef macos
  if (Key_SLOE .eq. 9)then
    print *, '    Warning :: MacOS does not support SuperLU solver!'
    print *, '               The solver will be switched to LAPACK!' 
    call log_msg('Warning :: MacOS does not support SuperLU solver!')
    call log_msg('The solver will be switched to LAPACK!')
    Key_SLOE = 5 
  endif
#endif

! ifort and ifx do not support the UFMPACK solver. 2024-03-08.
#ifdef ifort
  if (Key_SLOE .eq. 7)then
    !call Warning_Message('E','Key_SLOE')
    print *, '    Warning :: Intel ifort compiler does not support UFMPACK solver!'
    print *, '               The solver will be switched to Pardiso!'
    call log_msg('Warning :: Intel ifort compiler does not support UFMPACK solver!')
    call log_msg('The solver will be switched to Pardiso!')
    Key_SLOE = 3
  endif
#endif
#ifdef ifx
    if (Key_SLOE .eq. 7)then
        !call Warning_Message('E','Key_SLOE')
        print *, '    Warning :: Intel ifx compiler does not support UFMPACK solver!'
        print *, '               The solver will be switched to Pardiso!'
        call log_msg('Warning :: Intel ifx compiler does not support UFMPACK solver!')
        call log_msg('The solver will be switched to Pardiso!')
        Key_SLOE = 3
    endif
#endif  
#ifdef macos
  if (Key_SLOE .eq. 7)then
    print *, '    Warning :: MacOS does not support UFMPACK solver!'
    print *, '               The solver will be switched to LAPACK!' 
    call log_msg('Warning :: MacOS does not support UFMPACK solver!')
    call log_msg('The solver will be switched to LAPACK!')
    Key_SLOE = 5 
  endif
#endif
  ! ifort and ifx do not support the Lis solver. 2024-03-08.
#ifdef ifort
  if (Key_SLOE .eq. 8)then
    !call Warning_Message('E','Key_SLOE')
    print *, '    Warning :: Intel ifort compiler does not support Lis solver!'
    print *, '               The solver will be switched to Pardiso!'
    call log_msg('Warning :: Intel ifort compiler does not support Lis solver!')
    call log_msg('The solver will be switched to Pardiso!')
    Key_SLOE = 3
  endif
#endif
#ifdef ifx
    if (Key_SLOE .eq. 8)then
        !call Warning_Message('E','Key_SLOE')
        print *, '    Warning :: Intel ifx compiler does not support Lis solver!'
        print *, '               The solver will be switched to Pardiso!'
    call log_msg('Warning :: Intel ifx compiler does not support Lis solver!')
    call log_msg('The solver will be switched to Pardiso!')
        Key_SLOE = 3
    endif
#endif 
#ifdef ifort
    if (Key_SLOE .eq. 12)then
        print *, '    Warning :: Intel ifx compiler does not support STRUMPACK solver!'
        print *, '               The solver will be switched to Pardiso!'
    call log_msg('Warning :: Intel ifx compiler does not support STRUMPACK solver!')
    call log_msg('The solver will be switched to Pardiso!')
        Key_SLOE = 3
    endif
#endif  
#ifdef ifx
    if (Key_SLOE .eq. 12)then
        print *, '    Warning :: Intel ifx compiler does not support STRUMPACK solver!'
        print *, '               The solver will be switched to Pardiso!'
    call log_msg('Warning :: Intel ifx compiler does not support STRUMPACK solver!')
    call log_msg('The solver will be switched to Pardiso!')
        Key_SLOE = 3
    endif
#endif 
#ifdef macos
  if (Key_SLOE .eq. 12)then
    print *, '    Warning :: MacOS does not support STRUMPACK solver!'
    print *, '               The solver will be switched to LAPACK!' 
    call log_msg('Warning :: MacOS does not support STRUMPACK solver!')
    call log_msg('The solver will be switched to LAPACK!')
    Key_SLOE = 5 
  endif
#endif
  ! gfortran does not support the Pardiso solver. 2024-03-08.
#ifdef gfortran
    if (Key_SLOE .eq. 3)then
        print *,'    Error :: gfortran compiler does not support Pardiso solver!'
        call log_msg('Error :: gfortran compiler does not support Pardiso solver!')
        call Warning_Message('E','Key_SLOE')
    endif
#endif
  
  ! The Element by Element solver for item 11 only supports static analysis, 2021-06-01
  if (Key_Dimension==2) then
    if (Key_SLOE==11 .and. key_analysis_type /= 6) then
      if (Key_Analysis_Type /= 1 .and. Key_Analysis_Type /= 4)then
          print *, '    Error :: Key_SLOE=11 is available only when Key_Analysis_Type = 1 or 4!'
          call log_msg('Error :: Key_SLOE=11 is available only when Key_Analysis_Type = 1 or 4!')
          call Warning_Message('E','Key_SLOE')
      endif
    endif
  endif
  if (Key_Dimension==3) then
    if (Key_SLOE==11 .and. key_analysis_type /= 6) then
      if (Key_Analysis_Type /= 1 .and. Key_Analysis_Type /= 3 .and. Key_Analysis_Type /= 4  .and. Key_Analysis_Type /= 61)then
          print *, '    Error :: Key_SLOE=11 is available only  when Key_Analysis_Type = 1 or 4 or 61!'
          call log_msg('Error :: Key_SLOE=11 is available only  when Key_Analysis_Type = 1 or 4 or 61!')
          call Warning_Message('E','Key_SLOE')
      endif
    endif
  endif

  !Single precision float supports only 1,2,5,6 solvers.
#ifndef Silverfrost
  if(FT==4)then
     if (Key_SLOE /= 1 .and.Key_SLOE /= 2 .and. Key_SLOE /= 5 .and.Key_SLOE /= 6 .and. Key_SLOE /= 11)then
          print *,'    Error :: single float supports only 1,2,5,6,11 solvers!'
          call log_msg('Error :: single float supports only 1,2,5,6,11 solvers!')
          call Warning_Message('E','Key_SLOE')
      end if
  endif
#endif

#ifdef Silverfrost
  !2024-09-11.
  if (Key_Analysis_Type /= 6) then
    if (Key_SLOE /= 1 .and. Key_SLOE /= 2 .and. Key_SLOE /= 4 .and. Key_SLOE /= 10 .and. Key_SLOE /= 11)then
          print *,'    Error :: Silverfrost compiler supports only 1,2,4,10,11 solvers!'
          call log_msg('Error :: Silverfrost compiler supports only 1,2,4,10,11 solvers!')
          call Warning_Message('E','Key_SLOE')
    endif
  endif
#endif

  !****************************************
  !About Key_EBE_Precondition. 2023-01-17.
  !****************************************
  if(Key_SLOE==11)then
    if(Key_Dimension==3)then
        if(Key_EBE_Precondition==0) then
            print *, '    Precontioner for EBE:      Inactive'
        elseif(Key_EBE_Precondition==1)then
            print *, '    Precontioner for EBE:      Diagonal (default)'
        elseif(Key_EBE_Precondition==2)then
            print *, '    Precontioner for EBE:      HW (Huges-Winget)'
            ! Only supports Key_EBE_Sym_Storage_K==0
            if(Key_EBE_Sym_Storage_K/=0) then
               print *, '    Error :: Key_EBE_Precondition=2 is available only when Key_EBE_Sym_Storage_K=0!'
               call Warning_Message('S',Keywords_Blank)
            endif
            ! Only FEM is supported, XFEM is not supported.
            if(num_Crack>0)then
               print *, '    Error :: Key_EBE_Precondition=2 for FEM only, not for XFEM!'
               call Warning_Message('S',Keywords_Blank)
            endif
        else
            call Warning_Message('E','Key_EBE_Precondition')
        endif
    endif
  endif

  !*********************************************************************************
  ! Symmetrically store the EBE-PCG element Stiffness Matrix storK_FEM. 2022-11-10.
  !*********************************************************************************
  if(Key_SLOE ==11)then
    if (Key_EBE_Sym_Storage_K ==1)then
        print *, '    Symmetric storage of storK_FEM:  active!'
    elseif(Key_EBE_Sym_Storage_K ==0) then
        print *, '    Symmetric storage of storK_FEM:  negative (default)!'
    else
        call Warning_Message('E','Key_EBE_Sym_Storage_K')
    endif
  endif

  !************************************************
  ! Calculate matrix condition number, 2020-02-10.
  !************************************************
  if (Key_Cond_Number ==1)then
      print *, '    Calculate condition number of K: active!'
      ! Only LAPACK No. 5 and Lis solver No. 8 support condition number calculation.
      if((Key_SLOE /= 5) .and. (Key_SLOE /= 8)) then
          print *,'    Error :: condition number supports only 5,8 solvers!'
          call Warning_Message('S',Keywords_Blank)
      endif
  elseif(abs(Key_Cond_Number)>=2) then
      call Warning_Message('E','Key_Cond_Number')
  endif

  !****************************************************
  ! Calculate the determinant of a matrix, 2022-08-19.
  !****************************************************
  if (Key_Determinant ==1)then
      print *, '    Calculate determinant of K: active!'
      ! The EBE-PCG solver for No. 11 is not supported
      if(Key_SLOE == 11) then
          print *,'    Error :: determinant does not support 11 solvers!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Sparse solver not supported
      if(Key_K_Sparse== 1) then
          print *,'    Error :: determinant does not support sparse solvers!'
          call Warning_Message('S',Keywords_Blank)
      endif
  elseif(abs(Key_Cond_Number)>=2) then
      call Warning_Message('E','Key_Cond_Number')
  endif

  !******************************************
  ! Calculate matrix eigenvalues, 2022-08-19
  !******************************************
  if (Key_Eigenvalue ==1)then
      print *, '    Calculate eigenvalue of K: active!'
      ! The EBE-PCG solver for No. 11 is not supported
      if(Key_SLOE == 11) then
          print *,'    Error :: eigenvalue does not support solver 11!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Sparse solver not supported
      if(Key_K_Sparse== 1) then
          print *,'    Error :: eigenvalue does not support sparse solvers!'
          call Warning_Message('S',Keywords_Blank)
      endif
  elseif(abs(Key_Cond_Number)>=2) then
      call Warning_Message('E','Key_Cond_Number')
  endif

  !**********************************
  ! BLAS library related, 2022-08-22
  !**********************************
  if (Key_BLAS ==1)then
      print *, '    Using BLAS library (Parallel for Intel Fortran): active!'
      ! Only supports EBE-PCG solver No. 11
      if(Key_SLOE /= 11) then
          print *,'    Error :: Key_BLAS=1 supports only solver 11!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Only analysis type 4 is supported.
      if(Key_Analysis_Type/= 4) then
          print *,'    Error :: Key_BLAS=1 supports only analysis type 4!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! 3D only.
      if(Key_Dimension/= 3) then
          print *,'    Error :: Key_BLAS=1 supports only 3D problems!'
          call Warning_Message('S',Keywords_Blank)
      endif
  elseif(abs(Key_BLAS)>=2) then
      call Warning_Message('E','Key_BLAS')
  endif

  !****************************************************
  ! Number of elements and nodes, average element area
  !****************************************************
  1011 FORMAT(5X,'Number of elements:        ',A8)
  1012 FORMAT(5X,'Number of nodes:           ',A8)
  1013 FORMAT(5X,'Average area of elements:',E13.5,' m^2')
  1014 FORMAT(5X,'Average vol of elements :',E13.5,' m^3')
  write (i_char, '(i8)') Num_Elem
  write(*,1011) adjustl(i_char)
  call log_msg('Number of elements:        '//adjustl(i_char))
  write (i_char, '(i8)') Num_Node
  write(*,1012) adjustl(i_char)
  call log_msg('Number of nodes:           '//adjustl(i_char))
  if (Key_Dimension==2) then
      WRITE(*,1013) Ave_Elem_Area
  elseif (Key_Dimension==3) then
      WRITE(*,1014) Ave_Elem_Vol
  endif
  
  !*********************************************************
  ! Log output of the number of material types. 2024-03-17.
  !*********************************************************
  call Tool_chrpak_i4_to_s_left (num_of_Material, s_num_of_Material)
  temp_log = trim('Number of material types:  '//s_num_of_Material(1:len_trim(s_num_of_Material)))
  call log_msg(temp_log)
  
  !************************************************************************
  ! Material Parameter Weibull Distribution. 2024-06-24. NEWFTU2024062402.
  !************************************************************************
  do i_Mat = 1,num_of_Material
    if (Key_Weibull_E(i_Mat)==1)then
          if(Key_Dimension==3)then
              print *,'    Error :: *Key_Weibull_E_n =1 supports only 2D problems!'
              call Warning_Message('S',Keywords_Blank)
          endif
          print *, '    Weibull distribution of E: active!'
          print *, '    Shape parameter of Weibull distribution of E: ',Weibull_Parameters_E(i_Mat,1)
          ! The shape parameter must be greater than 0.
          if(Weibull_Parameters_E(i_Mat,1)<ZR) then
              print *,'    Error :: Weibull_Parameters_E(i_Mat,1) must > 0!'
              call Warning_Message('S',Keywords_Blank)
          endif
          ! Only supports elastic materials.
          if (Material_Type(i_mat) /= 1 )then   
              print *,'    Error :: *Key_Weibull_E_n =1 supports only elastic material!'
              call Warning_Message('S',Keywords_Blank)
          endif
    endif
  enddo
  
  !*************************
  ! Gauss Quadrature Scheme
  !*************************
  ! Triangulation, temporarily unavailable
  if (Key_Integral_Sol == 1)then
      print *, '    Integration scheme:        triangular partitioning method'
      print *, '    NOT AVALIABLE YET!!!'
      call Warning_Message('S',Keywords_Blank)
  ! Standard Multi-Point Gauss Integration Scheme
  elseif (Key_Integral_Sol == 2)then
      print *, '    Integration scheme:        standard Gauss integration'
      !*******************************
      ! Related to Gaussian integrals
      !*******************************
1301     FORMAT(5X,'Num of GP of enriched ele: ',A8)
1302     FORMAT(5X,'Num of GP of enriched inclusion ele: ',A8)
      if (Key_Dimension == 2)then
          if ((Num_Gauss_Points == 4*4)  .or.(Num_Gauss_Points == 6*6)  .or.(Num_Gauss_Points == 8*8) .or.   &
              (Num_Gauss_Points == 10*10).or.(Num_Gauss_Points == 12*12).or.(Num_Gauss_Points == 14*14).or.  &
              (Num_Gauss_Points == 20*20).or.(Num_Gauss_Points == 26*26).or.(Num_Gauss_Points == 30*30)) then
              write (i_char, '(i8)') Num_Gauss_Points
              write(*,1301) adjustl(i_char)
           else
               call Warning_Message('E','Num_Gauss_Points')
           endif
           ! The number of Gauss points for the mixed enhancement element must be at least 64.
          if ((Num_Gauss_P_Inc == 8*8).or.(Num_Gauss_P_Inc == 10*10).or. (Num_Gauss_P_Inc == 12*12).or.   &
              (Num_Gauss_P_Inc == 14*14).or.(Num_Gauss_P_Inc == 20*20).or. (Num_Gauss_P_Inc == 26*26).or. &
              (Num_Gauss_P_Inc == 30*30)) then
              write (i_char, '(i8)') Num_Gauss_P_Inc
              write(*,1302) adjustl(i_char)
           else
               call Warning_Message('E','Num_Gauss_P_Inc')
           endif
      end if
1401     FORMAT(5X,'Num of GP of common 3D element: ',A8)
      if (Key_Dimension == 3)then
          if ((Num_Gau_Points_3D == 2*2*2) .or. (Num_Gau_Points_3D == 3*3*3) .or. (Num_Gau_Points_3D == 4*4*4) .or. &
             (Num_Gau_Points_3D == 5*5*5)  .or. (Num_Gau_Points_3D == 6*6*6).or.(Num_Gau_Points_3D == 7*7*7) .or.     &
             (Num_Gau_Points_3D == 8*8*8)  .or. (Num_Gau_Points_3D == 9*9*9) .or.(Num_Gau_Points_3D == 10*10*10).or. &
             (Num_Gau_Points_3D == 15*15*15).or. (Num_Gau_Points_3D == 18*18*18).or.(Num_Gau_Points_3D == 20*20*20)) then
              write (i_char, '(i8)') Num_Gau_Points_3D
              write(*,1401) adjustl(i_char)
           else
               call Warning_Message('E','Num_Gauss_Points')
           endif
           ! Num_Gau_Points_3D must be less than or equal to Num_Gau_Points_3D_MC. 2022-07-26.
           if(Num_Gau_Points_3D > Num_Gau_Points_3D_MC)then
             print *, '    ERROR :: Num_Gau_Points_3D must <= Num_Gau_Points_3D_MC!'
             call Warning_Message('S',Keywords_Blank)
           endif
      end if
  ! Regular Quadrilateral Subdivision Plan
  elseif (Key_Integral_Sol == 3)then
1309     FORMAT(5X,'Number of sub-quads:                  ',A8)
1310     FORMAT(5X,'Number of sub-cubes:                  ',A8)
1325     FORMAT(5X,'Number of Gauss points of sub-cubes:  ',A8)
      if (Key_Dimension == 2)then
          print *, '    Integration scheme:        rectangular sub-quads method'
          if ((Num_Sub_Quads == 2*2).or.(Num_Sub_Quads == 3*3).or.(Num_Sub_Quads == 4*4).or.       &
              (Num_Sub_Quads == 5*5).or.(Num_Sub_Quads == 6*6).or.(Num_Sub_Quads == 7*7).or.       &
              (Num_Sub_Quads == 8*8).or.(Num_Sub_Quads == 9*9).or.(Num_Sub_Quads == 10*10).or.     &
              (Num_Sub_Quads == 11*11).or.(Num_Sub_Quads == 12*12).or.(Num_Sub_Quads == 13*13).or. &
              (Num_Sub_Quads == 14*14).or.(Num_Sub_Quads == 15*15).or.(Num_Sub_Quads == 16*16).or. &
              (Num_Sub_Quads == 17*17)) then
              write (i_char, '(i8)') Num_Sub_Quads
              write(*,1309) adjustl(i_char)
           else
               call Warning_Message('E','Num_Sub_Quads')
           endif
       elseif (Key_Dimension == 3)then
          print *, '    Integration scheme:        sub-cubes method'
          if ((Num_Sub_3D_Cubes == 1*1*1).or.(Num_Sub_3D_Cubes == 2*2*2).or.(Num_Sub_3D_Cubes == 3*3*3).or. &
              (Num_Sub_3D_Cubes == 4*4*4).or.(Num_Sub_3D_Cubes == 5*5*5).or.(Num_Sub_3D_Cubes == 6*6*6).or. &
              (Num_Sub_3D_Cubes == 7*7*7).or.(Num_Sub_3D_Cubes == 8*8*8).or.(Num_Sub_3D_Cubes == 9*9*9).or. &
              (Num_Sub_3D_Cubes == 10*10*10)) then
              write (i_char, '(i8)') Num_Sub_Quads
              write(*,1310) adjustl(i_char)
           else
               call Warning_Message('E','Num_Sub_3D_Cubes')
           endif

           write (i_char, '(i8)') Num_Gau_Points_3D_Cube
           write(*,1325) adjustl(i_char)

       endif

  else
      call Warning_Message('E','Key_Integral_Sol')
  end if


  !**************************************
  ! Dynamic analysis-related examination
  !**************************************
  if (Key_Analysis_Type ==2)then
      if(IDy_Num_Iteras ==0)then
          print *,'    Error :: *IDy_Num_Iteras not defined!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif
  if (Key_Analysis_Type ==6)then
      if(EDy_Num_Iteras ==0)then
          print *,'    Error :: *EDy_Num_Iteras not defined!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! The display analysis, due to the large number of calculation steps, by default only calculates
      ! node stress and displacement; stress, displacement, and coordinates at Gauss points are not
      ! calculated.
      if(Key_Post_CS_G_Coor ==1)then
          print *,'    Warning :: *Key_Post_CS_G_Coor=1 is not recommended for explicit dynamic analysis!'
      endif
      if(Key_Post_CS_G_Disp ==1)then
          print *,'    Warning :: *Key_Post_CS_G_Disp=1 is not recommended for explicit dynamic analysis!'
      endif
      if(Key_Post_CS_G_Strs ==1)then
          print *,'    Warning :: *Key_Post_CS_G_Strs=1 is not recommended for explicit dynamic analysis!'
      endif
  endif

  !*************************
  ! HYPLAS Analysis Related
  !2021-07-06.
  !*************************
  if (Key_Analysis_Type ==8)then
      ! XFEM not supported
      if(num_Crack >0 .or. num_Circ_Hole >0 .or. num_Ellip_Hole >0 .or. num_Circ_Incl >0 .or. num_Poly_Incl >0 .or. &
              num_Ellip_Incl >0)then
        print *,'    Error :: Analysis Type 8 do not support XFEM!'
        call Warning_Message('S',Keywords_Blank)
      endif

      ! Contact not supported
      if(Key_Contact==1)then
        print *,'    Error :: Analysis Type 8 do not support contact!'
        call Warning_Message('S',Keywords_Blank)
      endif
  endif

  !************************
  ! Model coordinate range
  !************************
  1021 FORMAT(5X,'x coordinates range of model: ',E13.5,' m to ',E13.5,' m')
  1022 FORMAT(5X,'y coordinates range of model: ',E13.5,' m to ',E13.5,' m')
  1023 FORMAT(5X,'z coordinates range of model: ',E13.5,' m to ',E13.5,' m')
  WRITE(*,1021) Min_X_Coor,Max_X_Coor
  WRITE(*,1022) Min_Y_Coor,Max_Y_Coor
  call Tool_chrpak_r8_to_s_left (Min_X_Coor, s_Min_X_Coor)
  temp_log = trim('Minimum x coordinate of model: '//s_Min_X_Coor(1:len_trim(s_Min_X_Coor))//' m')
  call log_msg(temp_log)
  call Tool_chrpak_r8_to_s_left (Max_X_Coor, s_Max_X_Coor)
  temp_log = trim('Maximum x coordinate of model: '//s_Max_X_Coor(1:len_trim(s_Max_X_Coor))//' m')
  call log_msg(temp_log)
  call Tool_chrpak_r8_to_s_left (Min_Y_Coor, s_Min_Y_Coor)
  temp_log = trim('Minimum y coordinate of model: '//s_Min_Y_Coor(1:len_trim(s_Min_Y_Coor))//' m')
  call log_msg(temp_log)
  call Tool_chrpak_r8_to_s_left (Max_Y_Coor, s_Max_Y_Coor)
  temp_log = trim('Maximum y coordinate of model: '//s_Max_Y_Coor(1:len_trim(s_Max_Y_Coor))//' m')
  call log_msg(temp_log)
  if (Key_Dimension==3) then
      WRITE(*,1023) Min_Z_Coor,Max_Z_Coor
      call Tool_chrpak_r8_to_s_left (Min_Z_Coor, s_Min_Z_Coor)
      temp_log = trim('Minimum z coordinate of model: '//s_Min_Z_Coor(1:len_trim(s_Min_Z_Coor))//' m')
      call log_msg(temp_log)
      call Tool_chrpak_r8_to_s_left (Max_Z_Coor, s_Max_Z_Coor)
      temp_log = trim('Maximum z coordinate of model: '//s_Max_Z_Coor(1:len_trim(s_Max_Z_Coor))//' m')
      call log_msg(temp_log)
  end if

  !*************************************
  ! Definition related to ruptured area
  !*************************************
  if(Key_Fracture_Zone==1)then
      print *, '    Fracture zone:             yes'
      ! Check if the rupture area exceeds the limit
      if(Frac_Zone_MinX < Min_X_Coor)then
          print *,'    Error :: illegal *Frac_Zone_MinX!'
          call Warning_Message('S',Keywords_Blank)
      elseif(Frac_Zone_MaxX > Max_X_Coor)then
          print *,'    Error :: illegal *Frac_Zone_MaxX!'
          call Warning_Message('S',Keywords_Blank)
      elseif(Frac_Zone_MinY < Min_Y_Coor)then
          print *,'    Error :: illegal *Frac_Zone_MinY!'
          call Warning_Message('S',Keywords_Blank)
      elseif(Frac_Zone_MaxY > Max_Y_Coor)then
          print *,'    Error :: illegal *Frac_Zone_MaxY!'
          call Warning_Message('S',Keywords_Blank)
      endif
      if(Key_Dimension ==3) then
          if(Frac_Zone_MinZ < Min_Z_Coor)then
              print *,'    Error :: illegal *Frac_Zone_MinZ!'
              call Warning_Message('S',Keywords_Blank)
          elseif(Frac_Zone_MaxZ > Max_Z_Coor)then
              print *,'    Error :: illegal *Frac_Zone_MaxZ!'
              call Warning_Message('S',Keywords_Blank)
          endif
      endif
      ! If the inspection is correct, save the crack area coordinate file for post-processing and
      ! plotting.
      call Save_Coors_of_Fracture_Zone
  elseif(Key_Fracture_Zone==0)then
      !do nothing
  else
      call Warning_Message('E','Key_Fracture_Zone')
  endif

  !************************************************************
  ! Definition of initial crack generation region. 2024-06-23.
  !************************************************************
  if(Key_Ini_Crack_Zone==1)then
      ! Currently only used for 2D problems.
      if(Key_Dimension ==3) then
          print *,'    Error :: *Key_Ini_Crack_Zone=1 if valid only for 2D problems!'
          call Warning_Message('S',Keywords_Blank)
      endif
      
      print *, '    Crack Initiation zone:     yes'
      ! Is it out of range?
      if(Ini_Crack_Zone_MinX < Min_X_Coor)then
          print *,'    Error :: illegal *Ini_Crack_Zone_MinX!'
          call Warning_Message('S',Keywords_Blank)
      elseif(Ini_Crack_Zone_MaxX > Max_X_Coor)then
          print *,'    Error :: illegal *Ini_Crack_Zone_MaxX!'
          call Warning_Message('S',Keywords_Blank)
      elseif(Ini_Crack_Zone_MinY < Min_Y_Coor)then
          print *,'    Error :: illegal *Ini_Crack_Zone_MinY!'
          call Warning_Message('S',Keywords_Blank)
      elseif(Ini_Crack_Zone_MaxY > Max_Y_Coor)then
          print *,'    Error :: illegal *Ini_Crack_Zone_MaxY!'
          call Warning_Message('S',Keywords_Blank)
      endif
  elseif(Key_Ini_Crack_Zone==0)then
      !do nothing
  else
      call Warning_Message('E','Key_Ini_Crack_Zone')
  endif
  
  !**********************************************
  ! Randomly generated natural crack related. 2D
  !**********************************************
1511 FORMAT(5X,'      -----      number: ',I5)
1521 FORMAT(5X,'      -----  orintation: ',F10.3)
1531 FORMAT(5X,'      -----      length: ',F10.3)
1541 FORMAT(5X,'      ----- △orintation: ',F10.3)
1551 FORMAT(5X,'      -----     △length: ',F10.3)
  if(Key_Dimension==2)then
    if(Key_Random_NaCr==1)then
      print *, '    Random natural fracture:   yes'
      WRITE(*,1511) num_Rand_Na_Crack
      WRITE(*,1521) NaCr_Orientation
      WRITE(*,1531) NaCr_Length
      WRITE(*,1541) NaCr_Ori_Delta
      WRITE(*,1551) NaCr_Len_Delta
      ! Check if the length is defined
      if(NaCr_Length <ZR)then
          print *,'    Error :: illegal *NaCr_Length!'
          call Warning_Message('E','NaCr_Length')
      endif
      ! Check if a length increment is defined (must be a positive number)
      if(NaCr_Len_Delta <ZR)then
          print *,'    Error :: illegal *NaCr_Len_Delta!'
          call Warning_Message('E','NaCr_Len_Delta')
      endif
      ! Check if the azimuth increment is defined (must be a positive number)
      if(NaCr_Ori_Delta <ZR)then
          print *,'    Error :: illegal *NaCr_Ori_Delta!'
          call Warning_Message('E','NaCr_Ori_Delta')
      endif
      if(num_Na_Crack /= 0)then
          print *,'    Num of initial NaCr *num_Na_Crack must be 0 when *Key_Random_NaCr=1!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! The number of randomly generated natural fractures cannot exceed the maximum number of fractures
      ! allowed by the program, 2021-06-01.
      if (num_Rand_Na_Crack > Max_Num_Cr)then
          print *,'    num_Rand_Na_Crack cannot be larger than Max_Num_Cr!'
          call Warning_Message('S',Keywords_Blank)
      endif
    elseif(Key_Random_NaCr==0)then
      !do nothing
    else
      call Warning_Message('E','Key_Random_NaCr')
    endif
  endif

  !************************************************************
  ! Randomly generated natural cracks related. 3D. 2022-06-11.
  !************************************************************
1611 FORMAT(5X,'      -----       number: ',I5)
1621 FORMAT(5X,'      -----         Type:  Rectangle')
1622 FORMAT(5X,'      -----         Type:  Circle (21-gon)')
1623 FORMAT(5X,'      -----         Type:  Polygon')
1624 FORMAT(5X,'      ----- num of edges:  ',I5)
1631 FORMAT(5X,'      -----   orintation: ',3F10.3)
1632 FORMAT(5X,'      -----         size: ',F10.3)
1641 FORMAT(5X,'      -----  △orintation: ',F10.3, ' degrees')
1642 FORMAT(5X,'      -----        △size: ',F10.3)
1651 FORMAT(5X,'      -----         Type:  Long and narrow rectangle')
1652 FORMAT(5X,'      -----  Rect L size: ',F10.3)
1653 FORMAT(5X,'      -----  Rect W size: ',F10.3)
1654 FORMAT(5X,'      ----- L orintation: ',3F10.3)
  if(Key_Dimension==3)then
    if(Key_Random_NaCr==1)then
      print *, '    Random natural fracture:   yes'
      WRITE(*,1611) num_Rand_Na_Crack
      WRITE(*,1631) NaCr_3D_n_Vector(1:3)
      WRITE(*,1632) NaCr_3D_Size
      WRITE(*,1641) NaCr_3D_n_Vector_Delta
      WRITE(*,1642) NaCr_3D_Sz_Delta
      ! Types of Natural Fracture Shapes
      if(Key_NaCr_Type_3D==1)then
          write(*,1621)
      elseif(Key_NaCr_Type_3D==2)then
          write(*,1622)
      elseif(Key_NaCr_Type_3D==3)then
          write(*,1623)
          write(*,1624) Num_Poly_Edges_NaCr
          if(Num_Poly_Edges_NaCr<=4) then
            print *,'    Error :: illegal *Num_Poly_Edges_NaCr must > 4!'
            call Warning_Message('E','Num_Poly_Edges_NaCr')
          endif
      elseif(Key_NaCr_Type_3D==4)then
          write(*,1651)
          if(NaCr_3D_Rect_L <ZR)then
              print *,'    Error :: illegal *NaCr_3D_Rect_L!'
              call Warning_Message('E','NaCr_3D_Rect_L')
          endif
          if(sum(abs(NaCr_3D_Rect_Longside_Vector)) <=Tol_10)then
              print *,'    Error :: illegal *NaCr_3D_Rect_Longside_Vector!'
              call Warning_Message('E','NaCr_3D_Rect_Longside_Vector')
          endif
          if(NaCr_3D_Rect_L <ZR)then
              print *,'    Error :: illegal *NaCr_3D_Rect_L!'
              call Warning_Message('E','NaCr_3D_Rect_L')
          endif
          WRITE(*,1652) NaCr_3D_Rect_L
          WRITE(*,1653) NaCr_3D_Rect_W
          WRITE(*,1654) NaCr_3D_Rect_Longside_Vector
      else
          call Warning_Message('E','Key_NaCr_Type_3D')
      endif
      !2024-02-26.
      if(Num_Poly_Edges_NaCr>200) then
          print *,'    Error :: illegal *Num_Poly_Edges_NaCr must < 200!'
          call Warning_Message('E','Num_Poly_Edges_NaCr')
      endif
      ! Check if the natural fracture size is defined
      if(NaCr_3D_Size <ZR .and. Key_NaCr_Type_3D/=4)then
          print *,'    Error :: illegal *NaCr_3D_Size!'
          call Warning_Message('E','NaCr_3D_Size')
      endif
      if(num_Na_Crack /= 0)then
          print *,'    Num of initial NaCr *num_Na_Crack must be 0 when *Key_Random_NaCr=1!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! The number of randomly generated natural fractures cannot exceed the maximum number of fractures
      ! allowed by the program, 2021-06-01.
      if (num_Rand_Na_Crack > Max_Num_Cr_3D)then
          print *,'    num_Rand_Na_Crack cannot be larger than Max_Num_Cr_3D!'
          call Warning_Message('S',Keywords_Blank)
      endif
    elseif(Key_Random_NaCr==0)then
      !do nothing
    elseif(Key_Random_NaCr==2)then
      print *, '    Random natural fracture:   yes (FracMan *.fab)'
    elseif(Key_Random_NaCr==3)then
      print *, '    Random natural fracture:   given in kpp file'
    else
      call Warning_Message('E','Key_Random_NaCr')
    endif
  endif

  !*************************************************************************************************
  ! Key_NaCr_Active_Scheme_3D, related to the 3D natural fracture activation algorithm. 2023-01-07.
  !*************************************************************************************************
  if(Key_Dimension==3)then
    if(Key_Random_NaCr >= 1)then
        if(Key_NaCr_Active_Scheme_3D==1) then
            print *, '    Key_NaCr_Active_Scheme_3D: 1 (default)'
        elseif(Key_NaCr_Active_Scheme_3D==2) then
            print *, '    Key_NaCr_Active_Scheme_3D: 2'
        elseif(Key_NaCr_Active_Scheme_3D==3) then
            print *, '    Key_NaCr_Active_Scheme_3D: 3'
            
            ! Key_NaCr_Active_Scheme_3D=3 does not support circular initial cracks. 2024-02-26.
            ! NEWFTU2024022601.
            if(Key_NaCr_Type_3D ==2)then
              print *,'    Error :: *Key_NaCr_Type_3D=2 is illegal when Key_NaCr_Active_Scheme_3D=3!'
              print *,'             Set *Key_NaCr_Type_3D=3 (polygon NF) instead!'
              call Warning_Message('S',Keywords_Blank)
            endif
        else
            call Warning_Message('E','Key_NaCr_Active_Scheme_3D')
        endif
    endif
    ! Logic Check. 2023-01-09.
    if (Key_NaCr_Active_Scheme_3D>=2 .and. Key_Cpp_Call_Fortran_Lib/=1) then
        if(Key_Random_NaCr == 0)then
          print *,'    Key_NaCr_Active_Scheme_3D >=2 is illegal when Key_Random_NaCr=0!'
          call Warning_Message('S',Keywords_Blank)
        endif
    endif
  endif

  !********************************
  ! Initial Crack Basic Inspection
  !********************************
1031 FORMAT(5X,'Number of initial cracks:  ',A8)
1051 FORMAT(5X,'Error :: Coors of crack ',I3,' was not defined!')
1171 FORMAT(5X,'Error :: Points of rectangular crack ',I3,' are not on the same plane!')
1172 FORMAT(5X,'Error :: Reduplicative points of rectangular crack ',I3, ' are found!')
1173 FORMAT(5X,'Error :: Corss lines of rectangular crack ',I3, ' are found!')
1175 FORMAT(5X,'Error :: Number of points of 3D crack ',I3,' < 4!')
1176 FORMAT(5X,'         Please check array Each_Cr_Poi_Num!')
1177 FORMAT(5X,'Error :: Number of points of polygon crack ',I3,' > 10!')
  write (i_char, '(i8)') num_Crack
  write(*,1031) adjustl(i_char)
  call log_msg('Number of initial cracks:        '//adjustl(i_char))

  if (num_Crack >0) then
      ! If it is a 2-dimensional problem
      if (Key_Dimension==2)then
        do i_C = 1,num_Crack
          if(sum(abs(Crack_Coor(i_C,:,1:2)))==0)then
              write (*,1051) i_C
              call Warning_Message('S',Keywords_Blank)
          end if
        end do
      ! If it is a three-dimensional problem
      elseif(Key_Dimension==3)then
        ! If it is a multi-stage fracturing analysis and the initial hydraulic fractures are generated
        ! automatically, this stage of detection will be skipped, 2022-05-31.
        if(Key_HF_Multistage_3D==1 .and. Key_Gen_Ini_Crack_Wellbores==1)then
            goto 301
        endif
        do i_C = 1,num_Crack
          if(sum(abs(Crack3D_Coor(i_C,:,1:3)))<=Tol_20 .and.   &
             sum(abs(Crack3D_Cir_Coor(i_C,1:7)))<=Tol_20 .and. &
             sum(abs(Crack3D_Ellip_Coor(i_C,1:8)))<=Tol_20) then
              write (*,1051) i_C
              call Warning_Message('S',Keywords_Blank)
          end if
        end do
        ! Check rectangular cracks or polygonal cracks. Added on 2022-04-23. NEWFTU2022042301.
        do i_C = 1,num_Crack
          if(sum(abs(Crack3D_Coor(i_C,:,1:3)))>Tol_20) then
            c_num_points = Each_Cr_Poi_Num(i_C)

            ! If c_num_points is less than 4, there should be at least 4 points. 2022-06-10.
            if(c_num_points<4)then
              write (*,1175) i_C
              write (*,1176)
              call Warning_Message('S',Keywords_Blank)
            endif

            ! If c_num_points is greater than 10, a maximum of 10 points is supported for defining a polygonal
            ! crack. 2022-06-10.
            if(c_num_points>10)then
              write (*,1177) i_C
              write (*,1176)
              call Warning_Message('S',Keywords_Blank)
            endif

            ! Check for duplicate points.
            call Matrix_Yes_Duplicated_Row_Dou(c_num_points,3,Crack3D_Coor(i_C,1:c_num_points,1:3),Yes_Dup,num_Dup)
            if(Yes_Dup)then
              write (*,1172) i_C
              call Warning_Message('S',Keywords_Blank)
            endif

            ! Check whether all points lie on the same plane. This is required for four-point polygons, but not
            ! for polygons with five or more points.
            c_num_points = Each_Cr_Poi_Num(i_C)
            if (c_num_points==4) then
                call Tool_3D_Yes_Points_on_Same_Plane(c_num_points,Crack3D_Coor(i_C,1:c_num_points,1:3),Yes_on)
                if(Yes_on .eqv. .false.)then
                    write (*,1171) i_C
                    call Warning_Message('S',Keywords_Blank)
                endif
            endif

            ! Check for any intersections.
            call Tool_3D_Yes_Points_Lines_Cross(c_num_points, Crack3D_Coor(i_C,1:c_num_points,1:3),Yes_Cross)
            if(Yes_Cross .eqv. .true.)then
              write (*,1173) i_C
              call Warning_Message('S',Keywords_Blank)
            endif
          end if
        enddo
      endif
  end if
  301 continue

  !**************************************************************************************************
  ! Check whether the number of initial crack coordinate points (Each_Cr_Poi_Num) is the same as the
  ! actual given number (Crack_Coor)
  !**************************************************************************************************
  ! First, ensure that the number of coordinate points for each crack is >= 2.
  2352 FORMAT(5X,'Error :: Each_Cr_Poi_Num(',I3,') wrong or not defined!')
  if (Key_Dimension==2)then
    do i_C=1,num_Crack
      if(Each_Cr_Poi_Num(i_C)<2)then
          write (*,2352) i_C
          call Warning_Message('S',Keywords_Blank)
      endif
    end do
  endif
  1052 FORMAT(5X,'Error :: Point ',I3,' of crack ',I3,' was not defined!')
  if (Key_Dimension==2)then
    if (num_Crack.ne. 0) then
      ! Individual crack cycles
      do i_C=1,num_Crack
          ! Current cracks cycle between each crack segment
          do i_P = 1,Each_Cr_Poi_Num(i_C)
              if(Crack_Coor(i_C,i_P,1)==ZR.and.Crack_Coor(i_C,i_P,2)==ZR )then
                  write (*,1052) i_P,i_C
                  call Warning_Message('S',Keywords_Blank)
              end if
          end do
      end do
    end if
  endif
  ! If it is a three-dimensional problem, the initial crack must contain at least four points in
  ! space.
  if (Key_Dimension==3)then
    do i_C=1,num_Crack
      if(sum(abs(Crack3D_Coor(i_C,1:4,1:3)))>=Tol_20)then
          if(Each_Cr_Poi_Num(i_C)<4)then
              write (*,2352) i_C
              call Warning_Message('S',Keywords_Blank)
          endif
      endif
    end do
  endif
  if (Key_Dimension==3)then
    if (num_Crack.ne. 0) then
      ! Individual crack cycles
      do i_C=1,num_Crack
          if(sum(abs(Crack3D_Coor(i_C,:,1:3)))>=Tol_20)then
              ! Current cracks cycle between each crack segment
              do i_P = 1,Each_Cr_Poi_Num(i_C)
                  if(Crack3D_Coor(i_C,i_P,1)==ZR .and.Crack3D_Coor(i_C,i_P,2)==ZR .and.Crack3D_Coor(i_C,i_P,3)==ZR)then
                      write (*,1052) i_P,i_C
                      call Warning_Message('S',Keywords_Blank)
                  end if
              end do
          endif
      end do
    end if
  endif

  !**************************************************************************************************
  ! Check whether the initial cracks are within the model range (only check rectangular or polygonal
  ! cracks). 2022-08-27. NEWFTU2022082701.
  !**************************************************************************************************
  3551 FORMAT(5X,'Error :: crack ',I3,' is outside the model!')
  if (Key_Dimension==3)then
    if (num_Crack.ne. 0) then
      ! Individual crack cycles
      do i_C=1,num_Crack
          c_Yes_Cr_In_Model = .False.
          ! Rectangular or polygonal cracks
          if(sum(abs(Crack3D_Coor(i_C,:,1:3)))>=Tol_20)then
              ! Current cracks cycle between each crack segment
              do i_P = 1,Each_Cr_Poi_Num(i_C)
                  c_P_x = Crack3D_Coor(i_C,i_P,1)
                  c_P_y = Crack3D_Coor(i_C,i_P,2)
                  c_P_z = Crack3D_Coor(i_C,i_P,3)
                  ! As long as any point of a crack is within the model, the crack is considered to be inside the
                  ! model.
                  if(c_P_x < Max_X_Coor .and. c_P_x > Min_X_Coor .and.&
                     c_P_y < Max_Y_Coor .and. c_P_y > Min_Y_Coor .and.&
                     c_P_z < Max_Z_Coor .and. c_P_z > Min_Z_Coor) then
                      c_Yes_Cr_In_Model = .True.
                      exit
                  endif
              end do
              ! If it is not in the model, an error will pop up.
              !if (c_Yes_Cr_In_Model .eqv. .False.) then
              if ((c_Yes_Cr_In_Model .eqv. .False.) .and. Key_Allow_3D_Outside_Crack==0) then
                  write (*,3551) i_C
                  call Warning_Message('S',Keywords_Blank)
              endif
          endif
      end do
    end if
  endif

  !***********************************************************************
  ! Basic Inspection of Initial Natural Fractures in Hydraulic Fracturing
  !***********************************************************************
  1331 FORMAT(5X,'Caution :: number of natural cracks:',I3)
  1351 FORMAT(5X,'Error :: Coors of natural crack ',I3,' not defined!')
  ! Natural fractures are limited to hydraulic fracturing failure analysis
  if(Key_Analysis_Type /= 3 .and. Key_Analysis_Type /= 4.and. num_Na_Crack/=0)then
      print *, '    Error :: natural fractures only used for HF analysis!'
      call Warning_Message('S',Keywords_Blank)
  elseif(Key_Analysis_Type == 3 .and. num_Na_Crack>0)then
      print *,'    *******************************************'
      WRITE(*,1331) num_Crack
      print *,'    *******************************************'
      ! Types of natural fractures
      if(Key_Na_Crack_Type==1)then
          print *, '    Natural crack type:        cemented, bilateral'
      elseif(Key_Na_Crack_Type==2)then
          print *, '    Natural crack type:        cemented, unilateral'
      elseif(Key_Na_Crack_Type==3)then
          print *, '    Natural crack type:        frictional'
      else
          call Warning_Message('E','Key_Na_Crack_Type')
      endif
  endif
  ! Check whether each natural fracture has been defined
  if (num_Na_Crack >0) then
      do i_NC = 1,num_Na_Crack
          if(sum(abs(Na_Crack_Coor(i_NC,:,1:2)))==0)then
              write (*,1351) i_NC
              call Warning_Message('S',Keywords_Blank)
          end if
      end do
  end if
  ! Check whether the number of natural crack coordinate points (Each_Na_Cr_Poi_Num) is the same as
  ! the actual given number (Na_Crack_Coor)
  1352 FORMAT(5X,'Error :: Point ',I3,' of natural crack ',I3,' was not defined!')
  if (num_Na_Crack.ne. 0) then
      ! Individual crack cycles
      do i_NC=1,num_Na_Crack
          ! Current cracks are cycling between each crack segment
          do i_P = 1,Each_Na_Cr_Poi_Num(i_NC)
              if(Na_Crack_Coor(i_NC,i_P,1)==ZR.and.Na_Crack_Coor(i_NC,i_P,2)==ZR )then
                  write (*,1352) i_P,i_NC
                  call Warning_Message('S',Keywords_Blank)
              end if
          end do
      end do
  end if

  !****************************************************************************
  ! Inspection and Related Correction of Initial Crack Position in 2D Problems
  !****************************************************************************
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (1) Check whether the initial crack is within the model (allowing one point to be outside the
  ! model for defining boundary cracks)
  ! (2) Check whether the crack point falls on the element boundary. If it falls exactly on the
  ! element boundary, it may cause
  ! Disorder, needs to be corrected (uniformly shift toward the lower-left corner by
  ! Modifi_Factor*Ave_Elem_L)
  ! (3) Check whether any nodes are exactly on the crack segment; if so, they need to be corrected.
  ! Method same as (2)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Modifi_Factor = 1.0D-8
  if (Key_Dimension == 2)then
    if (num_Crack.ne. 0) then
     ! Obtain a model boundary that is connected end to end
      EndByEnd_Outline(1:size(Outline,1)) = Outline(:,1)
      EndByEnd_Outline(size(Outline,1)+1) = Outline(1,1)
      ! Individual crack cycles
      do i_C=1,num_Crack
          ! Current cracks are cycling between each crack segment
          do i_S = 1,Each_Cr_Poi_Num(i_C)-1
              ! Crack fragment endpoint coordinates
              crack_p1=[Crack_Coor(i_C,i_S,1),Crack_Coor(i_C,i_S,2)]
              crack_p2=[Crack_Coor(i_C,i_S+1,1),Crack_Coor(i_C,i_S+1,2)]
              !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              ! (1) First, check whether the initial crack is within the model (allowing one point
              ! to be outside the model for defining boundary cracks)
              !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              ! Check whether crack_p1 is inside the model
              Call Tool_Yes_In_Poly(crack_p1(1),crack_p1(2),  &
                 Coor(EndByEnd_Outline,1),Coor(EndByEnd_Outline,2),size(EndByEnd_Outline),Yes_Seg_P1_in)
              ! Check whether crack_p2 is in the model
              Call Tool_Yes_In_Poly(crack_p2(1),crack_p2(2),  &
                  Coor(EndByEnd_Outline,1),Coor(EndByEnd_Outline,2),size(EndByEnd_Outline),Yes_Seg_P2_in)
              ! If both crack_p1 and crack_p2 are outside the model, an error is prompted.
              if ((Yes_Seg_P1_in .eqv. .False.) .and. (Yes_Seg_P2_in .eqv. .False.)) then
                  print *, '    Error :: wrong location of initial cracks!'
                  call Warning_Message('S',Keywords_Blank)
              end if
              !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              ! (2) Check whether it is on the element boundary (including the node positions, of
              ! course)
              !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              Yes_on_Ele_Edge = .False.
              Check_Point = crack_p1
              call Cal_Ele_Num_by_Coors_YesInOn(Check_Point(1),Check_Point(2), c_Ele,Yes_In,Yes_on_Ele_Edge)
              if (Yes_on_Ele_Edge .eqv. .True.) then
                  Crack_Coor(i_C,i_S,1) = Crack_Coor(i_C,i_S,1) -Modifi_Factor*Ave_Elem_L
                  Crack_Coor(i_C,i_S,2) = Crack_Coor(i_C,i_S,2) -Modifi_Factor*Ave_Elem_L
                  print *, '    Warning :: location of initial crack has been modified!'
                  !write (*,1041) i_C
              end if
              if(i_S ==Each_Cr_Poi_Num(i_C)-1) then
                  Check_Point = crack_p2
                  call Cal_Ele_Num_by_Coors_YesInOn(Check_Point(1),Check_Point(2),c_Ele,Yes_In,Yes_on_Ele_Edge)
                  if (Yes_on_Ele_Edge .eqv. .True.) then
                      Crack_Coor(i_C,i_S+1,1) = Crack_Coor(i_C,i_S+1,1) -Modifi_Factor*Ave_Elem_L
                      Crack_Coor(i_C,i_S+1,2) =Crack_Coor(i_C,i_S+1,2) -Modifi_Factor*Ave_Elem_L
                      print *, '    Warning :: location of initial crack ha been modified!'
                      !write (*,1041) i_C
                  end if
              end if
              !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              ! (3) Check whether there are nodes exactly on the crack fragment
              !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              Yes_Node_OnSeg = .False.
              Modifi_Factor = 1.0D-8
              ! Loop through all nodes (Note: this scheme has room for improvement, only nodes near the cracks
              ! need to be checked in a loop, TDL)
              do i_Node = 1,Num_Node
                  ! Coordinates of the current node
                  c_Node_x = Coor(i_Node,1)
                  c_Node_y = Coor(i_Node,2)
                  ! Check if the current node is on the current crack segment
                  call Tool_Yes_On_Line(c_Node_x,c_Node_y,crack_p1,crack_p2,Yes_Node_OnSeg)
                  ! If on the current crack segment, then correct the left endpoint.
                  if (Yes_Node_OnSeg.eqv..True.) then
                      Crack_Coor(i_C,i_S,1)=Crack_Coor(i_C,i_S,1) -Modifi_Factor*Ave_Elem_L
                      Crack_Coor(i_C,i_S,2)=Crack_Coor(i_C,i_S,2) -TWO*Modifi_Factor*Ave_Elem_L
                      print *, '    Warning :: location of initial crack had been modified!'
                      exit
                  end if
              end do
          end do
      end do
    end if
  endif

  !*************************************************************
  ! Initial arc-shaped crack basic inspection, added on 2017-07-15
  !  -----------------------------------------
  ! If a certain crack segment is detected to be an arc, then calculate Arc_Crack_Coor
  ! Variables 3, 4, 5, and 6 are r, Radian_Start, Radian_End, and Radian, respectively.
  ! Variables 1-2 (x, y) have already been provided.
  !*************************************************************
1054 FORMAT(5X,'Warning :: crack ',I3,' has arc segment (seg',I3,')!')
1055 FORMAT(5X,'Error :: illegal arc segment ',I3,' of crack ',I3,'!')
1041 FORMAT(16X,'Radian_Start of arc segment is ',I3,' degrees')
1042 FORMAT(16X,'Radian_End of arc segment is ',I3,' degrees')
1056 FORMAT(16X,'Radian of arc segment is ',I3,' degrees')
1057 FORMAT(16X,'Radius of arc segment is',E13.5,' m')
1058 FORMAT(16X,'Direction arc segment is anti-clockwise')
1059 FORMAT(16X,'Direction of arc segment is clockwise')
  if (num_Crack >0) then
      do i_C = 1,num_Crack
          ! Current cracks cycle between each crack segment
          do i_Seg = 1,Each_Cr_Poi_Num(i_C)-1
              if(sum(abs(Arc_Crack_Coor(i_C,i_Seg,1:10))) > ZR)then
                  write (*,1054) i_C,i_Seg
                  Yes_Arc_Crack = .True.
                  ! Calculate the radius and the starting angle of an arc (3 o'clock is 0 degrees, 0-360 degrees
                  ! counterclockwise) based on the coordinates of the two endpoints of the arc and the center of the
                  ! arc.
                  c_Arc_End_1(1:2) =[Crack_Coor(i_C,i_Seg,1),Crack_Coor(i_C,i_Seg,2)]
                  c_Arc_End_2(1:2) =[Crack_Coor(i_C,i_Seg+1,1),Crack_Coor(i_C,i_Seg+1,2)]
                  c_Arc_Center(1:2)= Arc_Crack_Coor(i_C,i_Seg,1:2)
                  c_Arc_Direction  = Arc_Crack_Coor(i_C,i_Seg,3)
                  call Tool_Arc_r_and_Radian_Given_Coors(c_Arc_Direction,c_Arc_End_1,c_Arc_End_2,c_Arc_Center,  &
                            c_Yes_feasible,c_Arc_r,c_Arc_Radian_S,c_Arc_Radian_E,c_Arc_Radian)
                  if (c_Yes_feasible .EQV. .False.) then
                      write (*,1055) i_Seg,i_C
                      call Warning_Message('S',Keywords_Blank)
                  else
                      Arc_Crack_Coor(i_C,i_Seg,4)=c_Arc_r
                      Arc_Crack_Coor(i_C,i_Seg,5)=c_Arc_Radian_S
                      Arc_Crack_Coor(i_C,i_Seg,6)=c_Arc_Radian_E
                      Arc_Crack_Coor(i_C,i_Seg,7)=c_Arc_Radian
                      Arc_Crack_Coor(i_C,i_Seg,8)=c_Arc_End_1(1)
                      Arc_Crack_Coor(i_C,i_Seg,9)=c_Arc_End_1(2)
                      Arc_Crack_Coor(i_C,i_Seg,10)=c_Arc_End_2(1)
                      Arc_Crack_Coor(i_C,i_Seg,11)=c_Arc_End_2(2)
                      write (*,1056) int(c_Arc_Radian)
                      write (*,1057) c_Arc_r
                      write (*,1041) int(c_Arc_Radian_S)
                      write (*,1042) int(c_Arc_Radian_E)
                      if (c_Arc_Direction > HLF)then
                          write (*,1058)
                      else
                          write (*,1059)
                      endif
                  endif
              end if
          end do
      end do
  end if

  !****************************************************************************
  ! Inspection and Related Correction of Initial Natural Crack Positions in 2D
  ! The same as the inspection of the initial cracks
  !****************************************************************************
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (1) Check whether the initial crack is within the model (allowing one point to be outside the
  ! model for defining boundary cracks)
  ! (2) Check whether the crack point falls on the element boundary. If it happens to fall on the
  ! element boundary, it will cause
  ! Disorder, needs to be corrected (uniformly shift toward the lower-left corner by
  ! Modifi_Factor*Ave_Elem_L)
  ! (3) Check whether any nodes are exactly on the crack segment; if so, they need to be corrected.
  ! Method same as (2)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Modifi_Factor = 1.0D-5
  if (Key_Dimension == 2)then
    if (num_Na_Crack.ne. 0) then
     ! Obtain a model boundary that is connected end to end
      EndByEnd_Outline(1:size(Outline,1)) = Outline(:,1)
      EndByEnd_Outline(size(Outline,1)+1) = Outline(1,1)
      ! Individual crack cycles
      do i_C=1,num_Na_Crack
          ! Current cracks are cycling between each crack segment
          do i_S = 1,Each_Na_Cr_Poi_Num(i_C)-1
              ! Crack fragment endpoint coordinates
              crack_p1=[Na_Crack_Coor(i_C,i_S,1),Na_Crack_Coor(i_C,i_S,2)]
              crack_p2=[Na_Crack_Coor(i_C,i_S+1,1),Na_Crack_Coor(i_C,i_S+1,2)]
              !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              ! (1) First, check whether the initial crack is within the model (allowing one point
              ! to be outside the model for defining boundary cracks)
              !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              ! Check whether crack_p1 is inside the model
              Call Tool_Yes_In_Poly (crack_p1(1),crack_p1(2),  &
                 Coor(EndByEnd_Outline,1),Coor(EndByEnd_Outline,2),size(EndByEnd_Outline),Yes_Seg_P1_in)
              ! Check whether crack_p2 is in the model
              Call Tool_Yes_In_Poly(crack_p2(1),crack_p2(2),   &
                 Coor(EndByEnd_Outline,1),Coor(EndByEnd_Outline,2),size(EndByEnd_Outline),Yes_Seg_P2_in)
              ! If both crack_p1 and crack_p2 are outside the model, an error is prompted.
              if ((Yes_Seg_P1_in .eqv. .False.) .and.(Yes_Seg_P2_in .eqv. .False.)) then
                  print *, '    Error :: wrong location of initial natural cracks!'
                  call Warning_Message('S',Keywords_Blank)
              end if
              !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              ! (2) Check whether it is on the element boundary (including the node position, of
              ! course)
              !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              Yes_on_Ele_Edge = .False.
              Check_Point = crack_p1
              call Cal_Ele_Num_by_Coors_YesInOn(Check_Point(1),Check_Point(2),c_Ele,Yes_In,Yes_on_Ele_Edge)
              if (Yes_on_Ele_Edge .eqv. .True.) then
                  Na_Crack_Coor(i_C,i_S,1) = Na_Crack_Coor(i_C,i_S,1) -Modifi_Factor*Ave_Elem_L
                  Na_Crack_Coor(i_C,i_S,2) = Na_Crack_Coor(i_C,i_S,2) -Modifi_Factor*Ave_Elem_L
                  print *, '      +++ Warning :: location of initial natural crack had been modified!'
                  !write (*,1041) i_C
              end if
              if(i_S ==Each_Na_Cr_Poi_Num(i_C)-1) then
                  Check_Point = crack_p2
                  call Cal_Ele_Num_by_Coors_YesInOn(Check_Point(1),Check_Point(2), c_Ele,Yes_In,Yes_on_Ele_Edge)
                  if (Yes_on_Ele_Edge .eqv. .True.) then
                      Na_Crack_Coor(i_C,i_S+1,1) =  Na_Crack_Coor(i_C,i_S+1,1) -Modifi_Factor*Ave_Elem_L
                      Na_Crack_Coor(i_C,i_S+1,2) =  Na_Crack_Coor(i_C,i_S+1,2) -Modifi_Factor*Ave_Elem_L
                      print *, '      +++ Warning :: location of initial natural crack had been modified!'
                      !write (*,1041) i_C
                  end if
              end if
              !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              ! (3) Check whether there are nodes exactly on the crack fragment
              !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              Yes_Node_OnSeg = .False.
              Modifi_Factor = 1.0D-4
              ! Loop through all nodes (Note: this scheme has room for improvement, only nodes near the cracks
              ! need to be checked in a loop, TDL)
              do i_Node = 1,Num_Node
                  ! Coordinates of the current node
                  c_Node_x = Coor(i_Node,1)
                  c_Node_y = Coor(i_Node,2)
                  ! Check if the current node is on the current crack segment
                  call Tool_Yes_On_Line(c_Node_x,c_Node_y,crack_p1,crack_p2,Yes_Node_OnSeg)
                  ! If on the current crack segment, then correct the left endpoint.
                  if (Yes_Node_OnSeg.eqv..True.) then
                      Na_Crack_Coor(i_C,i_S,1)=Na_Crack_Coor(i_C,i_S,1) -Modifi_Factor*Ave_Elem_L
                      Na_Crack_Coor(i_C,i_S,2)=Na_Crack_Coor(i_C,i_S,2) -TWO*Modifi_Factor*Ave_Elem_L
                      print *, '      +++ Warning :: location of initial natural crack had been modified!'
                      exit
                  end if
              end do
          end do
      end do
    end if
  endif

  !*****************************************
  ! Inspection of the initial circular hole
  !*****************************************
2031 FORMAT(5X,'Number of initial circle holes:   ',A8)
2051 FORMAT(5X,'Error :: Location of circle hole ',I3,' not defined!')
2061 FORMAT(5X,'Warning :: circle hole ',I3, ' outside the model')
  write (i_char, '(i8)') num_Circ_Hole
  write(*,2031) adjustl(i_char)
  call log_msg('Number of initial circle holes:  '//adjustl(i_char))

  if (num_Circ_Hole >0) then
      do i_H = 1,num_Circ_Hole
          if(sum(abs(Hole_Coor(i_H,1:3)))==0)then
              write (*,2051) i_H
              call Warning_Message('S',Keywords_Blank)
          end if
      end do
      ! Holes are only applicable to static analysis or dynamic analysis (including explicit and
      ! implicit).
      ! or static field analysis
      if(Key_Analysis_Type /= 1 .and. Key_Analysis_Type /= 2 .and. Key_Analysis_Type /= 6.and. Key_Analysis_Type /= 15)then
          print *, '    Error :: holes are avaliable only when Key_Analysis_Type =1 or 2 or 6!'
          call Warning_Message('S',Keywords_Blank)
      endif
  end if
  ! Check whether the holes exceed the bounds of the model (warning only, does not terminate the
  ! program)
  if (num_Circ_Hole .ne. 0) then
     ! Obtain a model boundary that is continuous from start to end
      EndByEnd_Outline(1:size(Outline,1)) = Outline(:,1)
      EndByEnd_Outline(size(Outline,1)+1) = Outline(1,1)
      ! Circulation of each hole
      do i_H=1,num_Circ_Hole
          ! Left
          Hole_point(1,1) = Hole_Coor(i_H,1) - Hole_Coor(i_H,3)
          Hole_point(1,2) = Hole_Coor(i_H,2)
          ! Right
          Hole_point(2,1) = Hole_Coor(i_H,1) + Hole_Coor(i_H,3)
          Hole_point(2,2) = Hole_Coor(i_H,2)
          ! Up
          Hole_point(3,1) = Hole_Coor(i_H,1)
          Hole_point(3,2) = Hole_Coor(i_H,2) + Hole_Coor(i_H,3)
          ! Down
          Hole_point(4,1) = Hole_Coor(i_H,1)
          Hole_point(4,2) = Hole_Coor(i_H,2) - Hole_Coor(i_H,3)
          do iii=1,4
              ! Check whether Hole_p1 is inside the model
              Call Tool_Yes_In_Poly(Hole_point(iii,1),Hole_point(iii,2),Coor(EndByEnd_Outline,1),Coor(EndByEnd_Outline,2), &
                  size(EndByEnd_Outline),Yes_Point_in)
              if(.not. Yes_Point_in)then
                  write (*,2061) i_H
                  
              endif
          enddo
      end do
  end if
  ! Check if the holes overlap (warning only, does not terminate the program)
4011 FORMAT(5X,'Warning :: hole ',I3,' and hole ',I3,' intersects!')
  do i_H=1,num_Circ_Hole
      x0 =  Hole_Coor(i_H,1)
      y0 =  Hole_Coor(i_H,2)
      R0 =  Hole_Coor(i_H,3)
      do j_H = i_H,num_Circ_Hole
          if(j_H /= i_H)then
              x1 =  Hole_Coor(j_H,1)
              y1 =  Hole_Coor(j_H,2)
              R1 =  Hole_Coor(j_H,3)
              call Tool_Yes_Two_Circles_Intersects(x0,y0,R0,x1,y1,R1,Yes_Intersect)
              if(Yes_Intersect)then
                  write(*,4011) i_H,j_H
              endif
          endif
      enddo
  enddo
  !Check if any crack intersects holes, added on 2017-05-14.
4012 FORMAT(5X,'ERROR :: hole',I3,' contains segment',I3,' of crack',I3)
4013 FORMAT(5X,'ERROR :: hole',I3,' contains wrong segment of crack',I3)
4014 FORMAT(5X,'Warning :: crack',I3,' intersects hole',I3,'!')
4015 FORMAT(5X,'Warning :: tip',I2,' of crack',I3,' modified to hole egde!')
4016 FORMAT(5X,'ERROR :: segment',I3,' of crack',I3,' goes through hole',I3,'!')
!      goto 777
  do i_H=1,num_Circ_Hole
      x0 =  Hole_Coor(i_H,1)
      y0 =  Hole_Coor(i_H,2)
      R0 =  Hole_Coor(i_H,3)
      do i_C=1,num_Crack
          do i_S = 1,Each_Cr_Poi_Num(i_C)-1
              ! Crack fragment endpoint coordinates
              crack_p1=[Crack_Coor(i_C,i_S,1),Crack_Coor(i_C,i_S,2)]
              crack_p2=[Crack_Coor(i_C,i_S+1,1),Crack_Coor(i_C,i_S+1,2)]
              ! Calculate the intersection between the current segment and the current hole
              call Tool_Intersection_Line_and_Circle(x0,y0,R0,crack_p1,crack_p2,c_num_Inter,c_State,c_Inter)
              ! If there are two intersection points, it means it is completely penetrated, and the program
              ! terminates.
              if(c_num_Inter==2)then
                  write(*,4016) i_S,i_C,i_H
                  call Warning_Message('S',Keywords_Blank)
              endif
              !Check if crack_p1 is inside the current hole.
              Yes_Seg_P1_in = .False.
              if(Tool_Function_2Point_Dis([x0,y0],crack_p1)<R0)then
                  Yes_Seg_P1_in  = .True.
              endif
              !Check if crack_p2 is inside the current hole.
              Yes_Seg_P2_in = .False.
              if(Tool_Function_2Point_Dis([x0,y0],crack_p2)<R0)then
                  Yes_Seg_P2_in  = .True.
              endif
              !if both points are inside the current hole, then stop PhiPsi.
              if(Yes_Seg_P1_in .and. Yes_Seg_P2_in)then
                  write(*,4012) i_H,i_S,i_C
                  call Warning_Message('S',Keywords_Blank)
              endif
              ! Only allow two crack tips within the current hole; otherwise, stop PhiPsi.
              if(Yes_Seg_P1_in .and. (i_S /= 1))then
                  write(*,4013) i_H,i_C
                  call Warning_Message('S',Keywords_Blank)
              endif
              if(Yes_Seg_P2_in .and.  (i_S /= (Each_Cr_Poi_Num(i_C)-1)))then
                  write(*,4013) i_H,i_C
                  call Warning_Message('S',Keywords_Blank)
              endif
              !All check passed
              if(Yes_Seg_P1_in .and. c_num_Inter==1)then
                  write(*,4014) i_C,i_H
                  write(*,4015) 1,i_C
                  !Modify the tip of the crack to the edge of the hole.
                  Crack_Coor(i_C,1,1:2) = c_Inter(1,1:2)
              endif
              if(Yes_Seg_P2_in .and. c_num_Inter==1)then
                  write(*,4014) i_C,i_H
                  write(*,4015) 2,i_C
                  !Modify the tip of the crack to the edge of the hole.
                  Crack_Coor(i_C,Each_Cr_Poi_Num(i_C),1:2)=c_Inter(1,1:2)
              endif
          enddo
      enddo
  enddo

  !*************************************
  ! Elliptical hole related, 2020-08-09
  !*************************************
2032 FORMAT(5X,'Number of initial ellipse holes:   ',A8)
2052 FORMAT(5X,'Error :: ellipse hole ',I3,' not defined!')
2152 FORMAT(5X,'Error :: a must > b for ellipse hole ',I3,'!')
2062 FORMAT(5X,'Warning :: ellipse hole ',I3, ' outside the model!')
  write (i_char, '(i8)') num_Ellip_Hole
  write(*,2032) adjustl(i_char)
  call log_msg('Number of initial ellipse holes: '//adjustl(i_char))
  if(num_Ellip_Hole>=1)then
      do i_H = 1,num_Ellip_Hole
          if(sum(abs(Ellip_Hole_Coor(i_H,1:5)))==0)then
              write (*,2052) i_H
              call Warning_Message('S',Keywords_Blank)
          end if
          x0 =  Ellip_Hole_Coor(i_H,1)
          y0 =  Ellip_Hole_Coor(i_H,2)
          a_Hole =  Ellip_Hole_Coor(i_H,3)
          b_Hole =  Ellip_Hole_Coor(i_H,4)
          !theta_Hole = Ellip_Hole_Coor(i_H,5)
          !a must > b.
          if(a_Hole<b_Hole)then
              write(*,2152) i_H
              call Warning_Message('S',Keywords_Blank)
          endif
          !Check if the center of the ellipse outside the model.
          if (x0>Max_X_Coor .or. x0<Min_X_Coor .or.y0>Max_Y_Coor .or. y0<Min_Y_Coor) then
              write(*,2062) i_H
              call Warning_Message('S',Keywords_Blank)
          endif
      end do
      ! Holes are only applicable to static analysis or dynamic analysis (including explicit and
      ! implicit).
      ! or static field analysis
      if(Key_Analysis_Type /= 1 .and. Key_Analysis_Type /= 2 .and. Key_Analysis_Type /= 6.and. Key_Analysis_Type /= 15)then
          print *, '    Error :: holes are avaliable only when Key_Analysis_Type =1 or 2 or 6!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Circular holes and elliptical holes cannot exist at the same time
      if(num_Circ_Hole>=1)then
          print *, '    Error :: circle holes and ellipse holes can not coexist!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif


  !******************************
  ! Initial Inclusion Inspection_2016-09-28
  !  --------------
  ! Spherical inclusion
  !******************************
  2131 FORMAT(5X,'Number of initial circle inclusions:   ',A8)
  2151 FORMAT(5X,'Error :: Location of circle inclusion ',I3,'!')
  2161 FORMAT(5X,'Error :: wrong location of circle inclusion ',I3)
  write (i_char, '(i8)') num_Circ_Incl
  write(*,2131) adjustl(i_char)
  call log_msg('Number of circle inclusions:     '//adjustl(i_char))

  if (num_Circ_Incl >0) then
      do i_Incl = 1,num_Circ_Incl
          if(sum(abs(Circ_Inclu_Coor(i_Incl,1:3)))==0)then
              write (*,2151) i_Incl
              call Warning_Message('S',Keywords_Blank)
          end if
      end do
      ! Holes are only applicable to static analysis or dynamic analysis.
      if(Key_Analysis_Type /= 1 .and. Key_Analysis_Type /= 2 .and. Key_Analysis_Type /= 6)then
          print *, '    Error :: inclusions are avaliable only when Key_Analysis_Type =1 or 2 or 6!'
          call Warning_Message('S',Keywords_Blank)
      endif
  end if
  ! Check whether each circular inclusion exceeds the range of the model
  if (num_Circ_Incl .ne. 0) then
     ! Obtain a model boundary that is connected end to end
      EndByEnd_Outline(1:size(Outline,1)) = Outline(:,1)
      EndByEnd_Outline(size(Outline,1)+1) = Outline(1,1)
      ! Circular Inclusion Loop
      do i_Incl=1,num_Circ_Incl
          ! Left
          Incl_point(1,1) = Circ_Inclu_Coor(i_Incl,1)-Circ_Inclu_Coor(i_Incl,3)
          Incl_point(1,2) = Circ_Inclu_Coor(i_Incl,2)
          ! Right
          Incl_point(2,1) = Circ_Inclu_Coor(i_Incl,1)+Circ_Inclu_Coor(i_Incl,3)
          Incl_point(2,2) = Circ_Inclu_Coor(i_Incl,2)
          ! Up
          Incl_point(3,1) = Circ_Inclu_Coor(i_Incl,1)
          Incl_point(3,2) = Circ_Inclu_Coor(i_Incl,2)+Circ_Inclu_Coor(i_Incl,3)
          ! Down
          Incl_point(4,1) = Circ_Inclu_Coor(i_Incl,1)
          Incl_point(4,2) = Circ_Inclu_Coor(i_Incl,2)-Circ_Inclu_Coor(i_Incl,3)
          do iii=1,4
              ! Check if it is within the model
              Call Tool_Yes_In_Poly(Incl_point(iii,1),Incl_point(iii,2),  &
                 Coor(EndByEnd_Outline,1),Coor(EndByEnd_Outline,2),size(EndByEnd_Outline),Yes_Point_in)
              if(.not. Yes_Point_in)then
                  write (*,2161) i_Incl
                  call Warning_Message('S',Keywords_Blank)
              endif
          enddo
      end do
  end if
  ! Check if circular inclusions overlap (warning only, does not terminate the program)
4111 FORMAT(5X,'Warning :: circle inclusion ',I3,' and circle inclusion ',I3,' intersects!')
  do i_Incl=1,num_Circ_Incl
      x0 =  Circ_Inclu_Coor(i_Incl,1)
      y0 =  Circ_Inclu_Coor(i_Incl,2)
      R0 =  Circ_Inclu_Coor(i_Incl,3)
      do j_Incl = i_Incl,num_Circ_Incl
          if(j_Incl /= i_Incl)then
              x1 =  Circ_Inclu_Coor(j_Incl,1)
              y1 =  Circ_Inclu_Coor(j_Incl,2)
              R1 =  Circ_Inclu_Coor(j_Incl,3)
              call Tool_Yes_Two_Circles_Intersects(x0,y0,R0,x1,y1,R1,Yes_Intersect)
              if(Yes_Intersect.eqv..True.)then
                  write(*,4111) i_Incl,j_Incl
              endif
          endif
      enddo
  enddo
  ! Check whether a material type number has been defined for each circular inclusion
4112 FORMAT(5X,'Warning :: circle inclusion ',I3,' has no mat number!')
4201 FORMAT(5X,'ERROR :: mat of circle inclusion ',I3,' not defined!')
  do i_Incl=1,num_Circ_Incl
      if (Circ_Inclu_Mat_Num(i_Incl) ==0) then
          write(*,4112) i_Incl
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Check whether the inclusion materials have been defined, 2021-06-30.
      max_material_type=maxval(Circ_Inclu_Mat_Num(1:num_Circ_Incl))
      if (num_of_Material > num_of_Material) then
          write(*,4201) i_Incl
          call Warning_Message('S',Keywords_Blank)
      endif

  enddo

  ! Check whether each circular inclusion is too small (in the presence of cracks), 2017-07-19
  ! Check the number of unit centroids contained in each circular inclusion, with at least 9 unit
  ! centroids.
4113 FORMAT(5X,'Warning :: circle inclusion ',I3,' is too small!')
4114 FORMAT(5X,'           it contains only ',I3,' element centroid!')
4115 FORMAT(5X,'           should contain 9 elements at least!')
  if (num_Crack>=1)then
      do i_Incl=1,num_Circ_Incl
          x0 =  Circ_Inclu_Coor(i_Incl,1)
          y0 =  Circ_Inclu_Coor(i_Incl,2)
          R0 =  Circ_Inclu_Coor(i_Incl,3)
          num_Centroid = 0
          do i_E = 1,Num_Elem
              c_Centroid(1:2) = Elem_Centroid(i_E,1:2)
              if(sqrt((c_Centroid(1)-x0)**2+ (c_Centroid(2)-y0)**2)<=R0)then
                  num_Centroid = num_Centroid + 1
              endif
          enddo
          if (num_Centroid<=9)then
              write(*,4113) i_Incl
              write(*,4114) num_Centroid
              write(*,4115)
              call Warning_Message('S',Keywords_Blank)
          endif
      enddo
  endif

  !******************************
  ! Initial Inclusion Inspection_2016-10-4
  !  --------------
  ! Polygonal inclusion
  !******************************
2231 FORMAT(5X,'Number of initial polygon inclusions:   ',A8)
2251 FORMAT(5X,'Error :: X coor of point',I3,' of poly inclusion',I3)
2252 FORMAT(5X,'Error :: Y coor of point',I3,' of poly inclusion',I3)
2261 FORMAT(5X,'Error :: wrong location of polygon inclusion ',I3)
2262 FORMAT(5X,'         point ',I3,' is beyond the range of model!')
  write (i_char, '(i8)') num_Poly_Incl
  write (*,2231) adjustl(i_char)
  call log_msg('Number of polygon inclusions:    '//adjustl(i_char))

  if (num_poly_Incl >0) then
      ! Interleaving is only applicable to static analysis or dynamic analysis.
      if(Key_Analysis_Type /= 1 .and. Key_Analysis_Type /= 2 .and. Key_Analysis_Type /= 6 )then
          print *, '    Error :: inclusions are avaliable only when Key_Analysis_Type =1 or 2 or 6!'
          call Warning_Message('S',Keywords_Blank)
      endif
  end if
  ! Check whether the coordinates of each polygon inclusion are provided
  if (num_poly_Incl .ne. 0) then
      ! Various polygonal interspersed loops
      do i_Incl=1,num_poly_Incl
          c_num_Edges = Poly_Inclu_Edges_Num(i_Incl)
          do iii=1,c_num_Edges
              ! Check if the x-coordinate of the polygon point has been defined
              if(Poly_Incl_Coor_x(i_Incl,iii)<-1.0D7)then
                  write (*,2251) iii, i_Incl
                  call Warning_Message('S',Keywords_Blank)
              endif
              ! Check if the y-coordinate of the polygon points has been defined
              if(Poly_Incl_Coor_y(i_Incl,iii)<-1.0D7)then
                  write (*,2252) iii, i_Incl
                  call Warning_Message('S',Keywords_Blank)
              endif
          enddo
      end do
  end if
  ! Check whether any polygons extend beyond the model's boundaries
  if (num_poly_Incl .ne. 0) then
     ! Obtain a model boundary that is connected end to end
      EndByEnd_Outline(1:size(Outline,1)) = Outline(:,1)
      EndByEnd_Outline(size(Outline,1)+1) = Outline(1,1)
      ! Various polygonal inclusions and loops
      do i_Incl=1,num_poly_Incl
          c_num_Edges = Poly_Inclu_Edges_Num(i_Incl)
          do iii=1,c_num_Edges
              Call Tool_Yes_In_Poly(Poly_Incl_Coor_x(i_Incl,iii),Poly_Incl_Coor_y(i_Incl,iii), &
                 Coor(EndByEnd_Outline,1),Coor(EndByEnd_Outline,2),size(EndByEnd_Outline),Yes_Point_in)
              if(.not. Yes_Point_in)then
                  write (*,2261) i_Incl
                  write (*,2262) iii
                  call Warning_Message('S',Keywords_Blank)
              endif
          enddo
      end do
  end if
  ! Check if polygon inclusions overlap (warning only, does not terminate the program)
  ! work to be done!
  ! Check whether each polygon inclusion has been assigned a material type number
  4212 FORMAT(5X,'Warning :: polygon inclusion ',I3,' has no mat num!')
  do i_Incl=1,num_poly_Incl
      if (Poly_Inclu_Mat_Num(i_Incl) ==0) then
          write(*,4212) i_Incl
          call Warning_Message('S',Keywords_Blank)
      endif
  enddo
  !***************************************************
  ! Randomly generated polygons with related elements
  !  Material_Interface,added on 2017-07-15
  !***************************************************
  ! The elongation of polygonal inclusions must be between 1 and 4.
  if (Rand_Poly_Incl_Irregular_Extension_Factor .LT. ONE)then
      Rand_Poly_Incl_Irregular_Extension_Factor =ONE
  elseif(Rand_Poly_Incl_Irregular_Extension_Factor .GT.FOU)then
      Rand_Poly_Incl_Irregular_Extension_Factor =FOU
  endif

  !**********************************************
  ! Randomly generated holes related, 2018-08-27
  !**********************************************
  if(Key_Random_Hole==1)then
      ! Only applicable to quasi-static analysis
      if(Key_Analysis_Type /=1)then
          print *, '    Error :: *Key_Random_Hole=1 is avaliable only when *Key_Analysis_Type=1!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! num_Rand_Hole must be greater than zero
      if(num_Rand_Hole <=0 )then
          print *, '    Error :: Illegal *num_Rand_Hole!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Rand_Hole_R must be greater than zero
      if(Rand_Hole_R <=0 )then
          print *, '    Error :: Illegal *Rand_Hole_R!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Rand_Hole_R_Delta must be thrown a lot
      if(Rand_Hole_R_Delta <0 )then
          print *, '    Error :: Illegal *Rand_Hole_R_Delta!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif

  !**********************************
  ! Holes causing cracks, 2018-08-27
  !**********************************
  if(Key_Random_Hole==1 .or. num_Circ_Hole>=1)then
      if(Key_Hole_Crack_Generate==1)then
          print *, '    Cracks emerged from hole:  allowed'
          if(num_Ellip_Hole>=1)then
              print *, '    Key_Hole_Crack_Generate =1 only for circle hole!'
              call Warning_Message('S',Keywords_Blank)
          endif
      endif
  endif

  !********************************************************************************************
  ! If both circular inclusions and cracks are present and crack propagation is allowed, it is
  ! necessary to define the interface material.
  !  Material_Interface,added on 2017-07-15
  !********************************************************************************************
1411 FORMAT(5X,'S_t of interface:        ',E13.5,' MPa')
1412 FORMAT(5X,'K_Ic of interface:       ',E13.5,' MPa.m^1/2')
  if ((num_Crack >0) .and. (num_Circ_Incl>0)) then
      ! At least one crack is allowed to extend
      if(sum(Cracks_Allow_Propa(1:num_Crack)) >= 1)then
          if(sum(abs(Material_Interface(1:2)))==ZR)then
              print *, '      +++ Warning :: Keyword *Material_Interface is not available!'
              print *, '                     *Material_Interface has been set to default values: 0.1e6,0.1e6!'
              
              Material_Interface(1) = 0.1e6
              Material_Interface(2) = 0.1e6
          else
              write(*,1411) Material_Interface(1)
              write(*,1412) Material_Interface(2)
          endif
      endif
  endif

  !******************************************
  !  Cracks_Allow_Propa, added on 2020-03-30
  !******************************************
3731 FORMAT(5X,'Warning :: crack ',I5,' is not allowed to propagate!')
  if (num_Crack >0) then
      do i_C=1,num_Crack
          if(Cracks_Allow_Propa(i_C) == 0)then
              write(*,3731) i_C
          endif
      enddo
  endif

  !*******************************
  ! Material Types and Parameters
  !*******************************
  
  1001 FORMAT(20X,'Elastic modulus    -- ',F16.8,' GPa')
  1002 FORMAT(20X,'Poisson ratio      -- ',F16.8)
  1003 FORMAT(20X,'Density            -- ',F16.8,' Kg/m^3')
  1004 FORMAT(20X,'Thickness          -- ',F16.8,' m')
  1005 FORMAT(20X,'Tensile strength   -- ',F16.8,' MPa')
  1006 FORMAT(20X,'Fracture toughness -- ',F16.8,' MPa.m^1/2')
  1007 FORMAT(20X,'Yield stress       -- ',F16.8,' MPa')
  1008 FORMAT(20X,'Tangent modulus    -- ',F16.8,' GPa')
  1201 FORMAT(20X,'E1 of fiber        -- ',F16.8,' GPa')
  1202 FORMAT(20X,'E2 of fiber        -- ',F16.8,' GPa')
  1203 FORMAT(20X,'E3 of fiber        -- ',F16.8,' GPa')
  1204 FORMAT(20X,'v12 of fiber       -- ',F16.8)
  1205 FORMAT(20X,'v23 of fiber       -- ',F16.8)
  1206 FORMAT(20X,'v13 of fiber       -- ',F16.8)
  1207 FORMAT(20X,'G12 of fiber       -- ',F16.8,' GPa')
  1208 FORMAT(20X,'G23 of fiber       -- ',F16.8,' GPa')
  1209 FORMAT(20X,'G13 of fiber       -- ',F16.8,' GPa')
  1210 FORMAT(20X,'Volume ratio       -- ',F16.8,'%')
  1211 FORMAT(20X,'Mat coodinates type-- Cartesian coordinate system')
  1212 FORMAT(20X,'Mat coodinates type-- Cylinder coordinate system')
  1311 FORMAT(20X,'Mat coor vector x  -- ',F16.8,F16.8,F16.8)
  1312 FORMAT(20X,'Mat coor vector y  -- ',F16.8,F16.8,F16.8)
  1313 FORMAT(20X,'Mat coor center    -- ',F16.8,F16.8,F16.8)
  1314 FORMAT(20X,'Mat coor vector z  -- ',F16.8,F16.8,F16.8)
  1315 FORMAT(5X,'Parameters of material type ',I3,':')
  1316 FORMAT(18X,43('*'))
  do i=1,num_of_Material
      !Check E, if zreo, then remind and stop.
      if(Material_Para(i,1)<=Tol_10)then
          print *, '    Error :: Illegal material parameter E!'
          call Warning_Message('S',Keywords_Blank)
      endif
      if (Material_Type(i) .eq. 1) then
          WRITE (tem_char,"(i2)") i
          write(*,1315) i
          write(*,1316)
          print *, '                     Material ' //tem_char// ' is isotropic material'
          write(*,1316)
          if(Key_Unit_System ==1)then
              WRITE(*,1001) Material_Para(i,1)/Cn_G
              WRITE(*,1002) Material_Para(i,2)
              WRITE(*,1003) Material_Para(i,3)
              WRITE(*,1004) Material_Para(i,4)
              WRITE(*,1005) Material_Para(i,5)/Cn_M
              WRITE(*,1006) Material_Para(i,6)/Cn_M
          elseif(Key_Unit_System ==2)then
              WRITE(*,1001) Material_Para(i,1)/1.0D3
              WRITE(*,1002) Material_Para(i,2)
              WRITE(*,1003) Material_Para(i,3)*1.0D12
              WRITE(*,1004) Material_Para(i,4)
              WRITE(*,1005) Material_Para(i,5)
              WRITE(*,1006) Material_Para(i,6)
          endif
      elseif (Material_Type(i) .eq. 2 .or.   &
              Material_Type(i) .eq. 3 .or.   &
              Material_Type(i) .eq. 4 .or.   &
              Material_Type(i) .eq. 6  ) then
          WRITE (tem_char,"(i2)") i
          write(*,1315) i
          write(*,1316)
          print *, '                     Material ' //tem_char//  ' is plastic material'
          write(*,1316)
          if(Key_Unit_System ==1)then
              WRITE(*,1001) Material_Para(i,1)/Cn_G
              WRITE(*,1002) Material_Para(i,2)
              WRITE(*,1003) Material_Para(i,3)
              WRITE(*,1004) Material_Para(i,4)
              WRITE(*,1005) Material_Para(i,5)/Cn_M
              WRITE(*,1006) Material_Para(i,6)/Cn_M
              WRITE(*,1007) Material_Para(i,12)/Cn_M
              WRITE(*,1008) Material_Para(i,13)/Cn_G
          elseif(Key_Unit_System ==2)then
              WRITE(*,1001) Material_Para(i,1)/1.0D3
              WRITE(*,1002) Material_Para(i,2)
              WRITE(*,1003) Material_Para(i,3)*1.0D12
              WRITE(*,1004) Material_Para(i,4)
              WRITE(*,1005) Material_Para(i,5)
              WRITE(*,1006) Material_Para(i,6)
              WRITE(*,1007) Material_Para(i,12)
              WRITE(*,1008) Material_Para(i,13)/1.0D3
          endif
          ! Activate plastic analysis flag
          Key_Plasticity = 1
          if(Material_Type(i) .eq. 3)then
              Key_Damage =1
          endif
          ! Plastic analysis currently only supports analysis types 7 and 8.
          if(Key_Analysis_Type /=7 .and. Key_Analysis_Type /=8)then
            if (Key_Damage ==1) then
                !no nothing
                ! Damage analysis supports quasi-static analysis and also supports cracks.

                ! Damage quasi-static analysis, only supports plane strain
                if(Key_Analysis_Type ==1 .and.Key_Type_2D       /=2) then
                    print *,'    Error ::  Damage material is avaliable only for plane strain model ' &
                     //'when *Key_Analysis_Type=1!'
                    call Warning_Message('S',Keywords_Blank)
                endif
            else
                print *,'    Error ::  Plastic material is avaliable only when *Key_Analysis_Type=7 or 8!'
                call Warning_Message('S',Keywords_Blank)
                ! Plastic analysis currently does not support cracks.
                if(num_Crack >0)then
                    print *, '    Error ::  Plastic material is avaliable only when *num_Crack=0!'
                    call Warning_Message('S',Keywords_Blank)
                endif
            endif
          endif
          ! Plastic analysis does not support cohesive cracks
          if(Key_TipEnrich==4) then
            print *, '    Error ::  Plastic material is not avaliable for cohesive crack (*Key_TipEnrich=4)!'
            call Warning_Message('S',Keywords_Blank)
          endif
          ! Plastic analysis does not currently support plane stress. For theory, see Chapter 9 of
          ! Computational Methods for Plasticity.
          if(Key_Type_2D==1) then
            print *, '    Error ::  Plastic material is not avaliable for plane stress model!'
            call Warning_Message('S',Keywords_Blank)
          endif
      elseif (Material_Type(i) .eq. 5) then
          WRITE (tem_char,"(i2)") i
          write(*,1315) i
          write(*,1316)
          print *, '                     Material ' //tem_char// ' is composite material'
          write(*,1316)
          WRITE(*,1001) Material_Para(i,1)/Cn_G
          WRITE(*,1002) Material_Para(i,2)
          WRITE(*,1003) Material_Para(i,3)
          WRITE(*,1004) Material_Para(i,4)
          WRITE(*,1005) Material_Para(i,5)/Cn_M
          WRITE(*,1006) Material_Para(i,6)/Cn_M
          WRITE(*,1007) Material_Para(i,12)/Cn_M
          WRITE(*,1008) Material_Para(i,13)/Cn_G
          WRITE(*,1201) Material_Para_Added(i,1)/Cn_G
          WRITE(*,1202) Material_Para_Added(i,2)/Cn_G
          WRITE(*,1203) Material_Para_Added(i,3)/Cn_G
          WRITE(*,1204) Material_Para_Added(i,4)
          WRITE(*,1205) Material_Para_Added(i,5)
          WRITE(*,1206) Material_Para_Added(i,6)
          WRITE(*,1207) Material_Para_Added(i,7)/Cn_G
          WRITE(*,1208) Material_Para_Added(i,8)/Cn_G
          WRITE(*,1209) Material_Para_Added(i,9)/Cn_G
          WRITE(*,1210) Material_Para_Added(i,10)*100.0D0
          if (Material_Para_Added(i,11)<=2.0) then
              WRITE(*,1211)
              WRITE(*,1311) Mat_Cartesian_Coor_Vector_x(i,1:3)
              WRITE(*,1312) Mat_Cartesian_Coor_Vector_y(i,1:3)
              if(sum(abs(Mat_Cartesian_Coor_Vector_x))<Tol_10)then
                  print *, '    Error ::  *Mat_Cartesi an_Coor_Vector_x not defined or illegal!'
                  call Warning_Message('S',Keywords_Blank)
              endif
              if(sum(abs(Mat_Cartesian_Coor_Vector_y))<Tol_10)then
                  print *, '    Error ::  *Mat_Cartesi an_Coor_Vector_y not defined or illegal!'
                  call Warning_Message('S',Keywords_Blank)
              endif
          elseif(Material_Para_Added(i,11)>=2.0) then
              WRITE(*,1212)
              WRITE(*,1313) Mat_cylinder_Coor_Center(i,1:3)
              WRITE(*,1314) Mat_cylinder_Coor_Vector_z(i,1:3)
              if(sum(abs(Mat_cylinder_Coor_Vector_z))<Tol_10)then
                  print *, '    Error ::  *Mat_cylinder_Coor_Vector_z not defined or illegal!'
                  call Warning_Message('S',Keywords_Blank)
              endif
          endif
          if(sum(abs(Material_Para_Added(i,:)))<Tol_10)then
              print *, '    Error ::  *Material_Para_Added not defined or illegal!'
              call Warning_Message('S',Keywords_Blank)
          endif
      end if
  end do


  !*****************************************************************************************
  ! Analysis Type 15 (Static Field Problem) Special Detection Area or Dynamic Field Problem
  !*****************************************************************************************
! 1601 FORMAT(5X,'kxx for field prob:    ',F10.3)
! 1603 FORMAT(5X,'kyy for field prob:    ',F10.3)
1604 FORMAT(5X,'Uniform body flux:     ',F10.3)
  if (Key_Analysis_Type == 15 .or. Key_Analysis_Type == 16) then
      ! Check for traffic sources within the area
      if(Key_Fd_Body_Source==1)then
          write(*,1604) Fd_Body_Source
      elseif(Key_Fd_Body_Source==0)then
          print *, '    Uniform body flux:     no!'
      else
          call Warning_Message('E','Key_Fd_Body_Source')
      endif
  endif

  !************************
  ! Thermal Stress Related
  !************************
  3162 FORMAT(5X,'Applied temperature of material type ',I3,' is ',F8.3,' K/℃')
  if(Key_Thermal_Stress==1)then
      ! Thermal stress calculation only supports quasi-static analysis
      if(Key_Analysis_Type/=1 .and. Key_Analysis_Type/=4)then
          print *, '    Error :: *Key_Thermal_Stress=1 only for *Key_Analysis_Type=1 or 4'
          call Warning_Message('S',Keywords_Blank)
      else
          print *,'    Thermal stress:             yes!'
          !Added on 2022-10-15.
          do i_Mat = 1,100
              if(abs(Thermal_Str_Temper(i_Mat))>=Tol_10) then
                  write(*,3162) i_Mat,Thermal_Str_Temper(i_Mat)
              endif
          enddo
      endif
      ! Temperature Stress Calculation Method. 2023-03-13.
      if(Key_Scheme_Thermal_Stress==1) then
          print *, '    Thermal stress method: according to current temperature'
      elseif(Key_Scheme_Thermal_Stress==2)then
          print *, '    Thermal stress method: according to temperature increment'
      else
          call Warning_Message('E','Key_Scheme_Thermal_Stress')
      endif

      ! Initial Temperature Definition Method. 2023-03-13.
      if(Key_Initial_Temperature==1) then
          print *, '    Initial temperature:   defined by material ID'
      elseif(Key_Initial_Temperature==2)then
          print *, '    Initial temperature:   read from external file'
          print *, '    Error :: *Key_Initial_Temperature=2 is not available yet!'
          call Warning_Message('S',Keywords_Blank)
      else
          call Warning_Message('E','Key_Initial_Temperature')
      endif

      ! Key_Initial_Temperature >= 1 is temporarily not applicable to 2D. 2023-03-13.
      if(Key_Dimension==2 .and. Key_Initial_Temperature>=1) then
          print *, '    Error :: *Key_Initial_Temperature>=1 is not available for 2D yet!'
          call Warning_Message('S',Keywords_Blank)
      endif

      ! Key_Scheme_Thermal_Stress >=1 is not yet applicable to 2D. 2023-03-13.
      if(Key_Dimension==2 .and. Key_Scheme_Thermal_Stress>=1) then
          print *, '    Error :: *Key_Scheme_Thermal_Stress>=1 is not available for 2D yet!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif

  !*****************************
  ! Earthquake Analysis Related
  !*****************************
  if(Key_EQ==1)then
      ! Seismic analysis is only available under dynamic analysis.
      if(Key_Analysis_Type/=2 .and. Key_Analysis_Type/=6)then
          print *, '    Error :: *Key_EQ=1 only for *Key_Analysis_Type=2 or 6'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Seismic acceleration data time interval
      if(EQ_Ac_Time_Gap<=0)then
          print *, '    Error :: *EQ_Ac_Time_Gap must > 0'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Whether the number of nodes for applying acceleration is defined
      if(num_EQ_Ac_nodes<=0)then
          print *,'    Error :: *num_EQ_Ac_nodes not defined!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Have all the nodes for applying acceleration been defined?
      if(num_EQ_Ac_nodes>0)then
          do i_Ac_node =1,num_EQ_Ac_nodes
              ! If undefined or out of node range, an error will be reported
              if (EQ_Ac_nodes(i_Ac_node)==0 .or.EQ_Ac_nodes(i_Ac_node)> Num_Node) then
                print *,'    Error :: *EQ_Ac_nodes not defined or illegal!'
                call Warning_Message('S',Keywords_Blank)
              endif
          enddo
      endif
      ! Earthquake acceleration and sine acceleration cannot exist at the same time
      if (Key_Sin_Accel==1)then
          print *,'    Error :: *Key_Sin_Accel cannot be 1 when *Key_EQ=1!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Output
      print *,'    Earthquake analysis:    Yes!'
  endif

  !******************************************
  ! Sine acceleration excitation correlation
  !******************************************
  if(Key_Sin_Accel==1)then
      ! Seismic analysis is only available under dynamic analysis.
      if(Key_Analysis_Type/=2 .and. Key_Analysis_Type/=6)then
          print *, '    Error :: *Key_Sin_Accel=1 only for *Key_Analysis_Type=2 or 6'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Whether the number of nodes for applying acceleration is defined
      if(Sin_Accel_num_Nodes==0)then
          print *,'    Error :: *Sin_Accel_num_Nodes not defined!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Have all the nodes for applying acceleration been defined?
      if(Sin_Accel_num_Nodes>0)then
          do i_Ac_node =1,Sin_Accel_num_Nodes
              ! If undefined or out of node range, an error will be reported
              if (Sin_Accel_Nodes(i_Ac_node)==0 .or.Sin_Accel_Nodes(i_Ac_node)> Num_Node) then
                print *,'    Error :: *Sin_Accel_Nodes not defined or illegal!'
                call Warning_Message('S',Keywords_Blank)
              endif
          enddo
      endif
      ! Check incentive direction
      if (Sin_Accel_Dire/=1 .and. Sin_Accel_Dire/=2)then
          print *,'    Error :: *Sin_Accel_Dire not defined or illegal!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Check excitation amplitude
      if (Sin_Accel_A<= ZR)then
        print *,'    Error :: *Sin_Accel_A not defined or illegle!'
        call Warning_Message('S',Keywords_Blank)
      endif
      ! Check incentive cycle
      if (Sin_Accel_T<= ZR)then
        print *,'    Error :: *Sin_Accel_T not defined or illegle!'
        call Warning_Message('S',Keywords_Blank)
      endif
      ! Output
      print *,'    Earthquake analysis:    Yes!'
  endif

  !*********************************************************************************************
  ! Analysis of Production Estimate on the 17th (Indirect Fluid-Structure Interaction Analysis)
  !*********************************************************************************************
  ! Only analyses on the 16th and 17th involve yield assessment.
  if(Key_Analysis_Type /=16 .and. Key_Analysis_Type /=17)then
      Key_Gas_Production =0
  endif
  ! Production evaluation analysis on the 17th: When Key_Proppant_Creep == 1, Key_Proppant_Active must
  ! be ensured to be 1.
  if(Key_Proppant_Creep  == 1 .and. Key_Proppant_Active /= 1)then
      print *,'    Error :: *Key_Proppant_Active must be 1 when *Key_Proppant_Creep = 1!'
      call Warning_Message('S',Keywords_Blank)
  endif
  ! The production evaluation analysis on the 17th only supports one fracture.
  if(Key_Analysis_Type ==17)then
      if(num_Crack >1)then
          print *,'    Error :: only support 1 crack when *Key_Analysis_Type = 17!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif
  ! Key_Changing_Kf and Key_Proppant_Active cannot both be set to 1
  if(Key_Analysis_Type ==17)then
      if(Key_Changing_Kf==1 .and. Key_Proppant_Active==1)then
          print *,'    Error :: *Key_Proppant_Active and *Key_Changing_Kf cannot be 1 at the same time!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif

  !*****************************************************
  ! Node degree of freedom coupling related, 2016-10-20
  !*****************************************************
  ! Check if the list of degrees of freedom coupling nodes is defined
  if(num_CP_x_nodes > 0)then
      do i_CP_node =1,num_CP_x_nodes
          ! If undefined or out of node range, an error will be reported
          if (CP_x_nodes(i_CP_node)==0 .or.CP_x_nodes(i_CP_node)> Num_Node) then
            print *,'    Error :: *CP_x_nodes not defined or illegal!'
            call Warning_Message('S',Keywords_Blank)
          endif
      enddo
      print *,'    COUPLE nodes (x):    Yes!'
  endif
  if(num_CP_y_nodes > 0)then
      do i_CP_node =1,num_CP_y_nodes
          ! If undefined or out of node range, an error will be reported
          if (CP_y_nodes(i_CP_node)==0 .or.CP_y_nodes(i_CP_node)> Num_Node) then
            print *,'    Error :: *CP_y_nodes not defined or illegal!'
            call Warning_Message('S',Keywords_Blank)
          endif
      enddo
      print *,'    COUPLE nodes (y):    Yes!'
  endif


  !************************************************************************************************
  ! 3D Analysis of Crack Tip Local Coordinate System Transformation Matrix Calculation, 2020-01-15
  !************************************************************************************************
  ! Num_Check_T_Matrix range check: acceptable range is 50 to 200
  if(Key_Dimension==3) then
      if (Num_Check_T_Matrix<50)then
          print *,'    Error :: *Num_Check_T_Matrix can not < 50!'
          call Warning_Message('S',Keywords_Blank)
      elseif(Num_Check_T_Matrix>200) then
          print *,'    Error :: *Num_Check_T_Matrix can not > 200!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif
  
  !********************************************************************
  ! Jagged crack pattern. 2024-06-25. For demonstration purposes only.
  !********************************************************************
  if (Key_Sawtooth_Crack_Path==1)then
      if(Key_Dimension==3) then
          print *,'    Error :: *Key_Sawtooth_Crack_Path=1 for 2D problems only!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! For use only with the maximum shear stress criterion.
      if(CFCP/=4) then
          print *,'    Error :: *Key_Sawtooth_Crack_Path=1 is valid only when *CFCP=4!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif
  
  !*********************************************************
  ! User-defined 2D crack extension path. NEWFTU2024111701.
  !*********************************************************
  if (Key_User_Defined_2D_Crack_Path==1)then
      if(Key_Dimension==3) then
          print *,'    Error :: *Key_User_Defined_2D_Crack_Path=1 for 2D problems only!'
          call Warning_Message('S',Keywords_Blank)
      endif
      if(Key_Analysis_Type /= 1) then
          print *,'    Error :: *Key_User_Defined_2D_Crack_Path=1 for static problem only!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif

  !************************************************
  ! Wellbore related, 2022-04-19, NEWFTU2022041901
  !************************************************
  1032 FORMAT(5X,'Number of wellbores(WB):  ',I3)
  1033 FORMAT(5X,'Point', I3,' of wellbore ',I3,' : ',F11.3,F11.3,F11.3)
  1034 FORMAT(5X,'Size of initial cracks along WB:',F11.3,' * ',F11.3, 'm')
  1035 FORMAT(5X,'Number of stages of WB', I3,' is ',I3)
  1036 FORMAT(5X,'Number of cracks of stage',I3,' of WB', I3,' is ',I3)
  1037 FORMAT(5X,'Coors of start point of fracturing of WB ',I3,': ',F12.3,F12.3,F12.3)
  1038 FORMAT(5X,'Coors of end point of fracturing of WB   ',I3,': ',F12.3,F12.3,F12.3)
  1039 FORMAT(5X,'The defined start point of fracturing of WB ',I3,' is not on WB Line!')
  1040 FORMAT(5X,'The defined end point of fracturing of WB ',I3,' is not on WB Line!')
  1141 FORMAT(5X,'Injection Q of stage',I3,' of WB', I3,' is ',F9.3, ' m^3/s')
  1142 FORMAT(5X,'Injection time of stage',I3,' of WB', I3,' is ',F11.3,' mins')
  1161 FORMAT(5X,'Injection Q of stage',I3,' of WB',I3,' is illegal!')
  1162 FORMAT(5X,'Injection time of stage',I3,' of WB',I3,' is illegal!')
  ! Stage fracturing, 2022-05-31.
  if(Key_HF_Multistage_3D==1) then
      print *,'    Staged fracturing:      Yes!'
  elseif(Key_HF_Multistage_3D==0) then
      print *,'    Staged fracturing:      No!'
  else
      call Warning_Message('E','Key_HF_Multistage_3D')
  endif

  if (num_Wellbore>=1)then
  !if (num_Wellbore>=1 .and. Key_XA/=1)then
      if(Key_Dimension==2)then
          print *,'    Error :: Wellbore for 3D problems only!'
          call Warning_Message('E','num_Wellbore')
      endif
      write(*,1032) num_Wellbore
      do i=1,num_Wellbore
          do j=1,num_Points_WB(i)
              write(*,1033) j,i,Wellbore_Coors(i,j,1:3)
          enddo
      enddo

      ! Check whether the wellbore coordinate points are within the model range. 2022-11-01.
      do i_wb = 1,num_Wellbore
          do i_Point_wb = 1,num_Points_WB(i_wb)
              c_P_x = Wellbore_Coors(i_wb,i_Point_wb,1)
              c_P_y = Wellbore_Coors(i_wb,i_Point_wb,2)
              c_P_z = Wellbore_Coors(i_wb,i_Point_wb,3)
              ! As long as any point of a crack is within the model, the crack is considered to be inside the
              ! model.
              if(c_P_x < Max_X_Coor .and. c_P_x > Min_X_Coor .and.&
                 c_P_y < Max_Y_Coor .and. c_P_y > Min_Y_Coor .and.&
                 c_P_z < Max_Z_Coor .and. c_P_z > Min_Z_Coor) then
                  !do nothing
              ! If it is not within the model.
              else
                  print *,'    ERROR-2022110109 :: illegal Wellbore_Coors!'
                  print *,'                        Maybe outside the model!'
                  print *,'                        i_wb, i_Point_wb:',i_wb, i_Point_wb
                  call Warning_Message('S',Keywords_Blank)
              endif
          enddo
      enddo

      ! Save to a wbif file for Matlab post-processing
      print *, "    Saving *.wbif file of wellbore..."
      Filename_1=trim(Full_Pathname)//'.wbif'
      open(101,file=Filename_1,status='unknown')
      write(101, '(6I3)') num_Wellbore,num_Points_WB(1),num_Points_WB(2),num_Points_WB(3),num_Points_WB(4),num_Points_WB(5)
      close(101)
      ! Save to a wbco file for Matlab post-processing
      print *, "    Saving *.wbco file of wellbore..."
      if (num_Wellbore>=1)then
          Filename_1=trim(Full_Pathname)//'.wbco_1'
          open(101,file=Filename_1,status='unknown')
          do j=1,10
              write(101, '(3E20.12)') Wellbore_Coors(1,j,1:3)
          enddo
          close(101)
      endif
      if (num_Wellbore>=2)then
          Filename_1=trim(Full_Pathname)//'.wbco_2'
          open(101,file=Filename_1,status='unknown')
          do j=1,10
              write(101, '(3E20.12)') Wellbore_Coors(2,j,1:3)
          enddo
          close(101)
      endif
      if (num_Wellbore>=3)then
          Filename_1=trim(Full_Pathname)//'.wbco_3'
          open(101,file=Filename_1,status='unknown')
          do j=1,10
              write(101, '(3E20.12)') Wellbore_Coors(3,j,1:3)
          enddo
          close(101)
      endif
      if (num_Wellbore>=4)then
          Filename_1=trim(Full_Pathname)//'.wbco_4'
          open(101,file=Filename_1,status='unknown')
          do j=1,10
              write(101, '(3E20.12)') Wellbore_Coors(4,j,1:3)
          enddo
          close(101)
      endif
      if (num_Wellbore>=5)then
          Filename_1=trim(Full_Pathname)//'.wbco_5'
          open(101,file=Filename_1,status='unknown')
          do j=1,10
              write(101, '(3E20.12)') Wellbore_Coors(5,j,1:3)
          enddo
          close(101)
      endif
      ! Number of sections for each wellbore, 2022-05-31.
      if(Key_HF_Multistage_3D==1) then
          do i=1,num_Wellbore
              write(*,1035) i,num_Stages_Wellbores(i)
          enddo
      endif
      ! Number of clusters per segment in each wellbore (initial number of fractures), 2022-05-31.
      if(Key_HF_Multistage_3D==1) then
          do i=1,num_Wellbore
              do j=1,num_Stages_Wellbores(i)
                  write(*,1036) j,i,num_Crs_Stages_Wellbores(i,j)
              enddo
          enddo
      endif
      ! Flow rate and injection time for each segment of each wellbore, May 31, 2022.
      if(Key_HF_Multistage_3D==1) then
          do i=1,num_Wellbore
            do j=1,num_Stages_Wellbores(i)
              ! If the traffic is less than 0, an error message will be displayed. 2022-07-17.
              if(Injection_Q_Stages_Wellbores(i,j)<ZR)then
                   write(*,1161) j,i
              endif
              ! If the time is less than 0, an error message will pop up. 2022-07-17.
              if(Injection_T_Stages_Wellbores(i,j)<ZR)then
                   write(*,1162) j,i
              endif
              write(*,1141) j,i,Injection_Q_Stages_Wellbores(i,j)
              write(*,1142) j,i,Injection_T_Stages_Wellbores(i,j)/60.0D0
            enddo
          enddo
      endif
      ! Wellbore fracturing section start and end coordinates, 2022-05-31.
      if(Key_HF_Multistage_3D==1) then
          do i=1,num_Wellbore
              write(*,1037) i,Wellbores_Start_Point(i,1:3)
              write(*,1038) i,Wellbores_End_Point(i,1:3)
          enddo
          ! Check whether the starting point is on WB, 2022-05-31.
          do i=1,num_Wellbore
               Dis_Tol = Max_Model_Range/TEN_3
!              call Tool_Yes_Point_on_Line_Segments_3D(Wellbores_Start_Point(i,1:3),  &
!                          Wellbore_Coors(i,1:num_Points_WB(i),1:3),num_Points_WB(i),c_Yes_on)
                          
              call Tool_Yes_Point_on_Line_Segments_3D_Tol(Wellbores_Start_Point(i,1:3),  &
                          Wellbore_Coors(i,1:num_Points_WB(i),1:3),Dis_Tol,num_Points_WB(i),c_Yes_on)
              if(c_Yes_on .eqv. .False.) then
                  write(*,1039) i
                  call Warning_Message('S',Keywords_Blank)
              endif
              call Tool_Yes_Point_on_Line_Segments_3D_Tol(Wellbores_End_Point(i,1:3),    &
                         Wellbore_Coors(i,1:num_Points_WB(i),1:3),Dis_Tol,num_Points_WB(i),c_Yes_on)
              if(c_Yes_on .eqv. .False.) then
                  write(*,1040) i
                  call Warning_Message('S',Keywords_Blank)
              endif
          enddo
      endif
      ! Automatically generated initial natural fractures, 2022-05-31.
      if(Key_Gen_Ini_Crack_Wellbores>=1)then
          print *,'    Generate initial HF cracks:      Yes!'
          if(Key_Gen_Ini_Crack_Wellbores==1)then
           print *,'    Initial HF crack type:           Rectangle'
          elseif(Key_Gen_Ini_Crack_Wellbores==2)then
           print *,'    Initial HF crack type:           Circle'
          elseif(Key_Gen_Ini_Crack_Wellbores==3)then
           print *,'    Initial HF crack type:           Polygon'
           print *,'    Number of polygon edges:        ',Num_Poly_Edges_PolyCr_WB
           if(Num_Poly_Edges_PolyCr_WB<=4)then
            print *,'    Error :: illegal *Num_Poly_Edges_PolyCr_WB!'
            print *,'             Num_Poly_Edges_PolyCr_WB must > 4!'
            call Warning_Message('E','Num_Poly_Edges_PolyCr_WB')
           endif
          endif
          if(Key_HF_Multistage_3D==0)then
              print *,'    Error :: *Key_Gen_Ini_Crack_Wellbores>=1 is available only when *Key_HF_Multistage_3D=1!'
              call Warning_Message('S',Keywords_Blank)
          endif
          ! Automatically generated initial crack size (square)
          write(*,1034) Size_Ini_Crack_Wellbores,Size_Ini_Crack_Wellbores
          ! If the initial cracks are generated automatically, there can be no cracks.
          if(num_Crack>=1)then
              print *,'    Error :: *num_Crack>=1 is illegal when *Key_Gen_Ini_Crack_Wellbores=1!'
              call Warning_Message('S',Keywords_Blank)
          endif
      endif
  endif
  ! For 3D segmented fracturing, the wellbore must be defined. 2022-05-31.
  if (Key_HF_Multistage_3D==1) then
      if(num_Wellbore<=0)then
          print *,'    Error :: *num_Wellbore<=0! Illegal for *Key_HF_Multistage_3D=1!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif

  !************************************
  ! Detection of some logical problems
  !************************************
  print *,'    Logic checking......'
  ! If it is the weighted average maximum principal stress criterion, it is necessary to calculate the
  ! Gauss point stress.
  if(CFCP==2)then
      Key_Post_CS_G_Strs=1
  endif
  ! If you need to calculate the displacement at Gauss points, you also need to calculate the
  ! coordinates of the Gauss points.
  if (Key_Post_CS_G_Disp==1) then
      Key_Post_CS_G_Coor = 1
  end if
  ! If you need to calculate the stress at Gauss points, you also need to calculate the coordinates
  ! and displacements at the Gauss points.
  if (Key_Post_CS_G_Strs==1) then
      Key_Post_CS_G_Disp = 1
      Key_Post_CS_G_Coor = 1
  end if
  ! For hydraulic fracturing analysis, if proppants are considered, it is obviously also necessary to
  ! take into account the crack surface contact iteration (using the penalty method).
  if (Key_Analysis_Type==3) then
      if (Key_Proppant==1) then
          Key_Contact = 1
      !elseif(Key_Proppant==0)then
      !    Key_Contact = 0
      endif
  endif
  ! If the transport of proppant is to be supported, it is obviously necessary to support the
  ! detection of proppant fracture closure.
  if(Key_Propp_Trans==1)then
      !Key_Contact  = 1
      !Key_Proppant = 1
  end if
  ! If it is staged fracturing, then Key_Post_Cracked_Ele = 1
  if(Key_HF_Multistage==1)then
      Key_Post_Cracked_Ele=1
  endif
  ! If it is not staged fracturing, then MS_Crack_Num = 1
  if(Key_HF_Multistage==0)then
      MS_Crack_Num=1
  endif

  ! The step size control for fatigue analysis and dynamic analysis (later removed) was set to the
  ! actual step size, rather than a dynamic step size.
  if(Key_Analysis_Type ==1 .and. Key_Static_Fatigue==1)then
      Key_Propa_Type= 2
      print *,'    --Set *Key_Propa_Type to 2 for fatigue analysis!'
  elseif(Key_Analysis_Type ==2)then
      ! Key_Propa_Type= 2 !Extension Type:=1, fixed step extension;=2, actual step extension (used for
      ! fatigue analysis and dynamic analysis, with automatic logic correction)
  endif
  ! Fatigue analysis obviously allows for crack propagation.
  if(Key_Static_Fatigue ==1)then
      Key_Propagation =1
  endif

  ! If it is the Gauss point integration scheme 3, then set the number of Gauss integration points.
  if(Key_Integral_Sol ==3)then
      Num_Gauss_Points  = Num_Sub_Quads*4
      Num_Gauss_P_Inc   = Num_Sub_Quads*4
  endif

  !Num_Gau_Points_3D must >= num_Gauss_S1 (2021-08-19)
  if (Num_Gau_Points_3D < num_Gauss_S1)then
      print *,'    Error :: *Num_Gau_Points_3D must > *num_Gauss_S1!'
      call Warning_Message('S',Keywords_Blank)
  endif

  ! The maximum number of load steps is 5000
  1674 FORMAT(5X,'*Num_Substeps should be less than',I5,'!')
  if (Num_Substeps>Max_Step) then
      write(*,1674) Max_Step
      call Warning_Message('E','Num_Substeps')
  endif
  
  ! Key_3D_HF_SlipWater_fk_Type is only applicable to 3D clean water fracturing analysis. 2023-08-08.
  if(Key_3D_HF_SlipWater_fk_Type >1) then
      if(Key_Dimension /=3)then
          print *,'    Error :: *Key_3D_HF_SlipWater_fk_Type for 3D only!'
          call Warning_Message('S',Keywords_Blank)
      endif
      if(Key_Analysis_Type /= 4) then
          print *,'    Error :: *Key_3D_HF_SlipWater_fk_Type is valid only when Key_Analysis_Type==4!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif

  ! Key_3D_HF_Time_Step_Method is only applicable to 3D hydraulic fracturing analysis. 2023-08-08.
  if(Key_3D_HF_Time_Step_Method >1) then
      if(Key_Dimension /=3)then
          print *,'    Error :: *Key_3D_HF_Time_Step_Method for 3D only!'
          call Warning_Message('S',Keywords_Blank)
      endif
      if(Key_Analysis_Type /= 4) then
          print *,'    Error :: *Key_3D_HF_Time_Step_Method is valid only when Key_Analysis_Type==4!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif
  
  ! Key_SIFs_DIM_Points==3 applies only to 2D. 2023-08-12.
  if(Key_SIFs_DIM_Points==3) then
      if(Key_Dimension ==3)then
          print *,'    Error :: *Key_SIFs_DIM_Points=3 for 2D only!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif
  
  9999 continue

  !*********************************************
  ! Molecular Dynamics Analysis Related Testing
  !*********************************************
  if(Key_Analysis_Type ==21)then
      if(MD_num_molecule <=0)then
          print *,'    Error :: *MD_num_molecule is illegal or not defined!'
          call Warning_Message('S',Keywords_Blank)
      endif
      if(MD_mss_molecule <=ZR)then
          print *,'    Error :: *MD_mss_molecule is illegal or not defined!'
          call Warning_Message('S',Keywords_Blank)
      endif
      if(MD_num_time_step <=0)then
          print *,'    Error :: *MD_num_time_step is illegal or not defined!'
          call Warning_Message('S',Keywords_Blank)
      endif
      if(MD_Delt_Time <=ZR)then
          print *,'    Error :: *MD_Delt_Time is illegal or not defined!'
          call Warning_Message('S',Keywords_Blank)
      endif
      if(MD_Dimension_x <=ZR)then
          print *,'    Error :: *MD_Dimension_x is illegal or not defined!'
          call Warning_Message('S',Keywords_Blank)
      endif
      if(MD_Dimension_y <=ZR)then
          print *,'    Error :: *MD_Dimension_y is illegal or not defined!'
          call Warning_Message('S',Keywords_Blank)
      endif
      if(Key_Dimension==3)then
        if(MD_Dimension_z <=ZR)then
          print *,'    Error :: *MD_Dimension_z is illegal or not defined!'
          call Warning_Message('S',Keywords_Blank)
        endif
      endif
  endif

  !*********************************************
  ! Near-field dynamic analysis related testing
  !*********************************************
  if(Key_Analysis_Type ==31)then
      !TO BE DONE HERE
  endif

  !******************************
  ! CFD Analysis Related Testing
  !******************************
  if(Key_Analysis_Type ==41)then
      !TO BE DONE HERE
  endif

  !**************************************************************
  ! 3D Hydraulic Fracturing Analysis Related Testing, 2021-12-03
  !**************************************************************
  if (Key_Analysis_Type ==3 .and. Key_Dimension==3)then
      ! Only supports Solver 11 Element by Element (Solver 11)
      if(Key_SLOE    /= 11)then
          print *,'    Error :: 3D HF analysis supports only solver 11!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Touch calculation is not supported
      if(Key_Contact   == 1)then
          print *,'    Error :: 3D HF analysis does not support contact!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Local refinement of tip mesh is not supported
      if(Key_Local_Mesh_Refine   == 1)then
          print *,'    Error :: 3D HF analysis does not support *Key_Local_Mesh_Refine=1!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! There must be an initial crack
      if(num_Crack   == 0)then
          print *,'    Error :: 3D HF analysis needs initial crack! Now num_Crack  = 0!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Composite materials are not supported
      if (any(Material_Type == 5))then
          print *,'    Error :: 3D HF analysis does not support composite material!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Filtration is temporarily not supported
      if (Key_Leakoff == 1) then
          print *,'    Error :: 3D HF analysis does not support fluid leakoff yet!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Proppants are not supported at the moment.
      if (Key_Proppant == 1) then
          print *,'    Error :: 3D HF analysis does not support proppant yet!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif

  !*********************************************
  ! Key_XA Analysis Related Testing, 2022-09-10
  !*********************************************
  if(Key_XA >= 1)then
      ! 3D only.
      if (Key_Dimension /= 3) then
          print *,'    Error :: Key_XA =1 supports only 3D problems!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Only Solver 11 is supported.
      if (Key_SLOE /= 11) then
          print *,'    Error :: Key_XA =1 supports only solver 11!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Only analysis type 1 is supported.
      if (Key_Analysis_Type /= 1) then
          print *,'    Error :: Key_XA =1 supports only analysis type 1!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! If Key_XA = 2. 2023-03-13.
      if(Key_XA == 2)then
          ! The thermal stress switch must be turned on.
          if (Key_Thermal_Stress /= 1) then
              print *,'    Error :: Key_Thermal_Stress must equals 1 when Key_XA = 2!'
              call Warning_Message('S',Keywords_Blank)
          endif
          ! Thermal stress must be calculated according to the temperature difference.
          if (Key_Scheme_Thermal_Stress /= 2) then
              print *,'    Error :: Key_Scheme_Thermal_Stress must equals 2 when Key_XA = 2!'
              call Warning_Message('S',Keywords_Blank)
          endif
      endif
  endif

  !*****************************************************************
  ! About 3D Fracturing Experiment Simulation Analysis. 2023-01-23.
  !*****************************************************************
  if(Key_Analysis_Type ==61)then
      ! 3D only.
      if (Key_Dimension /= 3) then
          print *,'    Error :: Key_Analysis_Type =61 supports only 3D problems!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Only Solver 11 is supported.
      if (Key_SLOE /= 11) then
          print *,'    Error :: Key_Analysis_Type =61 supports only solver 11!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Initial stress field is not supported.
      if (Key_InSitu_Strategy /= 0 ) then
          print *,'    Error :: Key_Analysis_Type =61 does not support insitu stress!'
          print *,'             Check keywrods *Key_InSitu_Strategy!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Segmented hydraulic fracturing analysis is not supported.
      if (Key_HF_Multistage_3D /= 0 ) then
          print *,'    Error :: Key_Analysis_Type =61 is not avaliable when Key_HF_Multistage_3D/=0!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! key_front_segmentation is not supported.
      if (key_front_segmentation /= 0 ) then
          print *,'    Error :: Key_Analysis_Type =61 is not avaliable when key_front_segmentation/=0!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! key_front_segmentation is not supported.
      if (key_front_segmentation /= 0 ) then
          print *,'    Error :: Key_Analysis_Type =61 is not avaliable when key_front_segmentation/=0!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Key_Local_Mesh_Refine is not supported.
      if (Key_Local_Mesh_Refine /= 0 ) then
          print *,'    Error :: Key_Analysis_Type =61 is not avaliable when Key_Local_Mesh_Refine/=0!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Key_Determinant is not supported.
      if (Key_Determinant /= 0 ) then
          print *,'    Error :: Key_Analysis_Type =61 is not avaliable when Key_Determinant/=0!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Key_Integral_Sol==4 is not supported.
      if (Key_Integral_Sol==4) then
          print *,'    Error :: Key_Analysis_Type =61 is not avaliable when Key_Integral_Sol=4!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Key_EBE_Sym_Storage_K==1 is not supported.
      if (Key_EBE_Sym_Storage_K==1) then
          print *,'    Error :: Key_Analysis_Type =61 is not avaliable when Key_EBE_Sym_Storage_K=1!'
          call Warning_Message('S',Keywords_Blank)
      endif
      ! Key_EBE_Sym_Storage_K==1 is not supported.
      if (Key_Post_Elements_Gauss_Num==1) then
          print *,'    Error :: Key_Analysis_Type =61 is not avaliable when Key_Post_Elements_Gauss_Num=1!'
          call Warning_Message('S',Keywords_Blank)
      endif

      if(HFE_Surface_Load_Num<=0 .or. HFE_Surface_Load_Num > Num_Surface_Loads) then
          print *,'    Error :: illegal *HFE_Surface_Load_Num!'
          call Warning_Message('S',Keywords_Blank)
      endif
      if(HFE_Hole_Mat_Number<=0 .or. HFE_Hole_Mat_Number > num_of_Material) then
          print *,'    Error :: illegal *HFE_Surface_Load_Num!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif
  

  !********************************
  ! Save crack radius. 2024-02-15.
  !********************************
  if(Key_Save_Crack_Radius==1)then
      if(Key_Dimension == 2)then
          print *,'    Error :: *Key_Save_Crack_Radius=1 for 3D only!'
          call Warning_Message('S',Keywords_Blank)
      endif
  endif
print *, "    Pre-Process done!"

RETURN
END SUBROUTINE Input_Check_Display
