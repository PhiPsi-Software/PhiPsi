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
 
module Global_Float_Type
  implicit none
  save
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !The following is single precision:
  !---------------------------
  ! Using single precision for calculating hydraulic fracturing problems yields extremely inaccurate
  ! results.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#ifdef SP
  !SP: Single 
  integer,PARAMETER::FT = 4
  REAL(kind=FT)::ZR   = 0.0e0
  REAL(kind=FT)::pi   = 3.1415926e0
  REAL(kind=FT)::Cn_H = 1.0e2
  REAL(kind=FT)::Cn_K = 1.0e3
  REAL(kind=FT)::Cn_M = 1.0e6
  REAL(kind=FT)::Cn_G = 1.0e9
  REAL(kind=FT)::HLF  = 0.5e0
  REAL(kind=FT)::ONE  = 1.0e0
  REAL(kind=FT)::TWO  = 2.0e0
  REAL(kind=FT)::THR  = 3.0e0
  REAL(kind=FT)::FOU  = 4.0e0
  REAL(kind=FT)::FIV  = 5.0e0
  REAL(kind=FT)::SIX  = 6.0e0
  REAL(kind=FT)::SEV  = 7.0e0
  REAL(kind=FT)::EIG  = 8.0e0
  REAL(kind=FT)::NIN  = 9.0e0
  REAL(kind=FT)::TEN  = 1.0e1
  REAL(kind=FT)::ZP1  = 0.1e0
  REAL(kind=FT)::ZP2  = 0.2e0
  REAL(kind=FT)::ZP3  = 0.3e0
  REAL(kind=FT)::ZP4  = 0.4e0
  REAL(kind=FT)::ZP5  = 0.5e0
  REAL(kind=FT)::ZP6  = 0.6e0
  REAL(kind=FT)::ZP7  = 0.7e0
  REAL(kind=FT)::ZP8  = 0.8e0
  REAL(kind=FT)::ZP9  = 0.9e0
  REAL(kind=FT)::ZPZ1 = 0.01e0
  REAL(kind=FT)::ZPZZ1= 0.001e0
  REAL(kind=FT)::ZPZZZ1= 0.0001e0
  REAL(kind=FT)::ZP25 = 0.25e0
  REAL(kind=FT)::ONEP5= 1.5e0
  REAL(kind=FT)::Time_year  = 365.0e0*30.416e0*24.0e0*3600.0e0
  REAL(kind=FT)::Time_month = 30.416e0*24.0e0*3600.0e0
  REAL(kind=FT)::Time_week  = 7.0e0*24.0e0*3600.0e0
  REAL(kind=FT)::Time_day   = 24.0e0*3600.0e0
  REAL(kind=FT)::Time_hour  = 3600.0e0
  REAL(kind=FT)::Time_min   = 60.0e0
  REAL(kind=FT)::Tol_4      = 1.0e-4
  REAL(kind=FT)::Tol_5      = 1.0e-5
  REAL(kind=FT)::Tol_6      = 1.0e-6
  REAL(kind=FT)::Tol_7      = 1.0e-7
  REAL(kind=FT)::Tol_8      = 1.0e-8
  REAL(kind=FT)::Tol_9      = 1.0e-9
  REAL(kind=FT)::Tol_10     = 1.0e-10
  REAL(kind=FT)::Tol_12     = 1.0e-12
  REAL(kind=FT)::Tol_13     = 1.0e-13
  REAL(kind=FT)::Tol_15     = 1.0e-15
  REAL(kind=FT)::Tol_20     = 1.0e-20
  REAL(kind=FT)::Tol_25     = 1.0e-25
  REAL(kind=FT)::Tol_30     = 1.0e-30
  REAL(kind=FT)::Tol_35     = 1.0e-35
  REAL(kind=FT)::Tol_40     = 1.0e-40
  REAL(kind=FT)::Con_30     = 30.0e0
  REAL(kind=FT)::Con_90     = 90.0e0
  REAL(kind=FT)::Con_180    = 180.0e0
  REAL(kind=FT)::Con_240    = 240.0e0
  REAL(kind=FT)::Con_270    = 270.0e0
  REAL(kind=FT)::Con_360    = 360.0e0
  REAL(kind=FT)::Con_Big_15 = 1.0e15
  REAL(kind=FT)::Con_Big_20 = 1.0e20
  REAL(kind=FT)::TEN_2      = 1.0e2
  REAL(kind=FT)::TEN_3      = 1.0e3
  REAL(kind=FT)::TEN_4      = 1.0e4
  REAL(kind=FT)::TEN_5      = 1.0e5
  REAL(kind=FT)::TEN_6      = 1.0e6
  REAL(kind=FT)::TEN_7      = 1.0e7
  REAL(kind=FT)::TEN_8      = 1.0e8
  REAL(kind=FT)::TEN_9      = 1.0e9
  REAL(kind=FT)::TEN_10     = 1.0e10
  REAL(kind=FT)::TEN_11     = 1.0e11
  REAL(kind=FT)::TEN_12     = 1.0e12
  REAL(kind=FT)::TEN_13     = 1.0e13
  REAL(kind=FT)::TEN_14     = 1.0e14
  REAL(kind=FT)::TEN_15     = 1.0e15
#endif
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !The following is double precision:
  !------------------------------
  ! Note: Don't forget to modify the DGAMMA function in Tool_Function_GAMMA.f!!!!!!
  ! It was later proven that no modification was needed, 2019-04-27
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifndef SP

#ifndef Silverfrost
  integer,PARAMETER::FT  = 8
  integer,PARAMETER::LIT = 8
#endif 
  
#ifdef Silverfrost
  ! Silverfrost Fortran Compiler
  integer,PARAMETER::FT  = 2
  integer,PARAMETER::LIT = 4
#endif 
 
  REAL(kind=FT)::ZR   = 0.0D0
  REAL(kind=FT)::pi   = 3.14159265358979323846264D0
  REAL(kind=FT)::Cn_H = 1.0D2
  REAL(kind=FT)::Cn_K = 1.0D3
  REAL(kind=FT)::Cn_M = 1.0D6
  REAL(kind=FT)::Cn_G = 1.0D9
  REAL(kind=FT)::HLF  = 0.5D0
  REAL(kind=FT)::ONE  = 1.0D0
  REAL(kind=FT)::TWO  = 2.0D0
  REAL(kind=FT)::THR  = 3.0D0
  REAL(kind=FT)::FOU  = 4.0D0
  REAL(kind=FT)::FIV  = 5.0D0
  REAL(kind=FT)::SIX  = 6.0D0
  REAL(kind=FT)::SEV  = 7.0D0
  REAL(kind=FT)::EIG  = 8.0D0
  REAL(kind=FT)::NIN  = 9.0D0
  REAL(kind=FT)::TEN  = 1.0D1
  REAL(kind=FT)::ZP1  = 0.1D0
  REAL(kind=FT)::ZP2  = 0.2D0
  REAL(kind=FT)::ZP3  = 0.3D0
  REAL(kind=FT)::ZP4  = 0.4D0
  REAL(kind=FT)::ZP5  = 0.5D0
  REAL(kind=FT)::ZP6  = 0.6D0
  REAL(kind=FT)::ZP7  = 0.7D0
  REAL(kind=FT)::ZP8  = 0.8D0
  REAL(kind=FT)::ZP9  = 0.9D0
  REAL(kind=FT)::ZPZ1 = 0.01D0
  REAL(kind=FT)::ZPZZ1= 0.001D0
  REAL(kind=FT)::ZPZZZ1= 0.0001D0
  REAL(kind=FT)::ZP25 = 0.25D0
  REAL(kind=FT)::ONEP5= 1.5D0
  ! Others
  REAL(kind=FT)::Time_year  = 365.0D0*30.416D0*24.0D0*3600.0D0
  REAL(kind=FT)::Time_month = 30.416D0*24.0D0*3600.0D0
  REAL(kind=FT)::Time_week  = 7.0D0*24.0D0*3600.0D0
  REAL(kind=FT)::Time_day   = 24.0D0*3600.0D0
  REAL(kind=FT)::Time_hour  = 3600.0D0
  REAL(kind=FT)::Time_min   = 60.0D0
  REAL(kind=FT)::Tol_2      = 1.0D-2
  REAL(kind=FT)::Tol_3      = 1.0D-3
  REAL(kind=FT)::Tol_4      = 1.0D-4
  REAL(kind=FT)::Tol_5      = 1.0D-5
  REAL(kind=FT)::Tol_6      = 1.0D-6
  REAL(kind=FT)::Tol_7      = 1.0D-7
  REAL(kind=FT)::Tol_8      = 1.0D-8
  REAL(kind=FT)::Tol_9      = 1.0D-9
  REAL(kind=FT)::Tol_10     = 1.0D-10
  REAL(kind=FT)::Tol_11     = 1.0D-11
  REAL(kind=FT)::Tol_12     = 1.0D-12
  REAL(kind=FT)::Tol_13     = 1.0D-13
  REAL(kind=FT)::Tol_15     = 1.0D-15
  REAL(kind=FT)::Tol_20     = 1.0D-20
  REAL(kind=FT)::Tol_25     = 1.0D-25
  REAL(kind=FT)::Tol_30     = 1.0D-30
  REAL(kind=FT)::Tol_35     = 1.0D-35
  REAL(kind=FT)::Tol_40     = 1.0D-40
  REAL(kind=FT)::Con_30     = 30.0D0
  REAL(kind=FT)::Con_60     = 60.0D0
  REAL(kind=FT)::Con_90     = 90.0D0
  REAL(kind=FT)::Con_180    = 180.0D0
  REAL(kind=FT)::Con_240    = 240.0D0
  REAL(kind=FT)::Con_270    = 270.0D0
  REAL(kind=FT)::Con_360    = 360.0D0
  REAL(kind=FT)::Con_Big_15 = 1.0D15
  REAL(kind=FT)::Con_Big_20 = 1.0D20
  REAL(kind=FT)::TEN_2      = 1.0D2
  REAL(kind=FT)::TEN_3      = 1.0D3
  REAL(kind=FT)::TEN_4      = 1.0D4
  REAL(kind=FT)::TEN_5      = 1.0D5
  REAL(kind=FT)::TEN_6      = 1.0D6
  REAL(kind=FT)::TEN_7      = 1.0D7
  REAL(kind=FT)::TEN_8      = 1.0D8
  REAL(kind=FT)::TEN_9      = 1.0D9
  REAL(kind=FT)::TEN_10     = 1.0D10
  REAL(kind=FT)::TEN_11     = 1.0D11
  REAL(kind=FT)::TEN_12     = 1.0D12
  REAL(kind=FT)::TEN_13     = 1.0D13
  REAL(kind=FT)::TEN_14     = 1.0D14
  REAL(kind=FT)::TEN_15     = 1.0D15
  REAL(kind=FT)::TEN_20     = 1.0D20
#endif
end module Global_Float_Type

!--------------------------------------------
! 01. Jagged Array Class (Real). 2022-09-02.
!--------------------------------------------
!Ref: https://en.wikipedia.org/wiki/Jagged_array
!Ref: https://stackoverflow.com/questions/18316592/multidimensional-array-with-different-lengths
!Ref:
!https://stackoverflow.com/questions/14857366/are-there-any-problems-with-using-jagged-arrays-in-fortran-with-multiple-levels
module Global_Ragged_Array_Real_Classs
  use Global_Float_Type
  implicit none
  save
  !private
  !public::Ragged_Array_1D,Ragged_Array_2D,Ragged_Array_3D,Ragged_Array_4D

  !///////////
  !  1D Array
  !///////////
  type :: Ragged_Array_1D
      real(kind=FT), allocatable :: row(:)
  end type Ragged_Array_1D

  !Usage:
  !type(Ragged_Array_1D) :: Ragged(3)
  !Ragged(1)%row = [1]
  !Ragged(2)%row = [1,2]
  !Ragged(3)%row = [1,2,3]

  !or

  !type(ragged_array),Ragged(:)
  !allocate(Ragged_Array_2D_1(3))
  !allocate(Ragged_Array_2(1)%row(1))
  !allocate(Ragged_Array_2(2)%row(2))
  !allocate(Ragged_Array_2(3)%row(3))
  !Ragged_Array_2(1)%row(:) = 1.0D0
  !Ragged_Array_2(2)%row(:) = 2.0D0
  !Ragged_Array_2(3)%row(:) = 3.0D0

  !///////////
  !  2D Array
  !///////////
  type :: Ragged_Array_2D
      real(kind=FT), allocatable :: row(:,:)
  end type Ragged_Array_2D

  !Usage:
  !allocate(Ragged_Array_2D_1(3))
  !allocate(Ragged_Array_2D_1(1)%row(1,5))
  !allocate(Ragged_Array_2D_1(2)%row(2,5))
  !allocate(Ragged_Array_2D_1(3)%row(3,5))
  !Ragged_Array_2D_1(1)%row(:,:) = 1.0D0
  !Ragged_Array_2D_1(2)%row(:,:) = 2.0D0
  !Ragged_Array_2D_1(3)%row(:,:) = 3.0D0

  !///////////
  !  3D Array
  !///////////
  type :: Ragged_Array_3D
      real(kind=FT), allocatable :: row(:,:,:)
  end type Ragged_Array_3D

  !///////////
  !  4D Array
  !///////////
  type :: Ragged_Array_4D
      real(kind=FT), allocatable :: row(:,:,:,:)
  end type Ragged_Array_4D


end module Global_Ragged_Array_Real_Classs

!-------------------------------------------
! 02. Jagged Array Class (Int). 2022-09-03.
!-------------------------------------------
!Ref: https://en.wikipedia.org/wiki/Jagged_array
!Ref: https://stackoverflow.com/questions/18316592/multidimensional-array-with-different-lengths
!Ref:
!https://stackoverflow.com/questions/14857366/are-there-any-problems-with-using-jagged-arrays-in-fortran-with-multiple-levels
module Global_Ragged_Array_Int_Classs
  implicit none
  save
  !private
  !public::Ragged_Int_Array_1D,Ragged_Int_Array_2D,Ragged_Int_Array_3D,Ragged_Int_Array_4D

  !///////////
  !  1D Array
  !///////////
  type :: Ragged_Int_Array_1D
      integer, allocatable :: row(:)
  end type Ragged_Int_Array_1D

  !///////////
  !  2D Array
  !///////////
  type :: Ragged_Int_Array_2D
      integer, allocatable :: row(:,:)
  end type Ragged_Int_Array_2D


  !///////////
  !  3D Array
  !///////////
  type :: Ragged_Int_Array_3D
      integer, allocatable :: row(:,:,:)
  end type Ragged_Int_Array_3D

  !///////////
  !  4D Array
  !///////////
  type :: Ragged_Int_Array_4D
      integer, allocatable :: row(:,:,:,:)
  end type Ragged_Int_Array_4D


end module Global_Ragged_Array_Int_Classs

!-----------------------------
! 1. General Global Variables
!-----------------------------
module Global_Common
  use Global_Float_Type
  implicit none
  save
  character(120) :: PhiPsi_Version,PhiPsi_Release_Date,Compiler
  integer Operation_System_Type
  integer Compiler_Type
  integer(kind = LIT) S_time
  integer(kind = LIT) Label_number
  integer system_bit
  real(kind=FT) Delta_Factor_Aper,Delta_Factor_Junc,Delta_Factor_Edge
  real(kind=FT) Factor_Propagation
  real(kind=FT) Propagation_Length
  integer Key_Propa_Type
  integer Max_Step,Max_Itera
  integer max_nthreads
  character(80) :: mac_Address
  character(80) :: CPU_ID,CPU_Type,DISK_ID,MAC_ID
  character(1) :: String_Connector
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! cccccccccccccccccc         Internal Program Control Parameters        ccccccccccccccccccccccc
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  parameter (Delta_Factor_Edge    = 0.001D0 )
  parameter (Delta_Factor_Aper    = 0.001D0 )
                                                    ! 0.001 * average unit length
  parameter (Delta_Factor_Junc    = 0.001D0 )
                                                    ! 0.001 * average unit length
  parameter (Max_Step        = 1000  )
  parameter (Max_Itera       = 10000 )
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! cccccccccccccccccc         Internal Program Control Parameters        ccccccccccccccccccccccc
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer Key_Heaviside_Value
  integer Key_Cond_Number
  integer Key_Determinant
  integer Key_Eigenvalue
  integer Key_BLAS
                                                    !EBE_XFEM_PCG_3D_with_K.f90
  integer Key_Hole_Value
  real(kind=FT) Delta_Factor_SIFs
  integer Key_Tip_Pres_Zero
  integer Num_Gauss_Points,Num_Substeps,Key_Integral_Sol
  integer Num_Sub_Quads
  integer Num_Sub_3D_Cubes
  integer Num_Gau_Points_3D_Cube
  integer Num_Gau_Points_3D,Num_Gauss_P_FEM_3D
  integer Num_Gau_Points_3D_MC
  integer Num_Gauss_P_Inc
  integer Num_Gauss_P_FEM
  integer Num_Gau_P_SubInteg_6
  integer Num_Gau_P_SubInteg_4
  integer Key_Initiation,Key_Propagation,Key_Type_2D,CFCP
  integer Key_CFCP_3_Type
  integer Key_3D_Cr_Update_Scheme
                                                    ! 2: Non-propagating crack tip does not extend; the crack tip extends in very small increments
                                                    ! (default).
  integer Key_Large_Deform
  integer Key_Ini_Rule
  integer Key_Ini_Cr_3D_Type
  real(kind=FT) Size_Ini_Crack
  integer Key_Schollmann_Twist
  real(kind=FT) Schollm_Max_Theta
  integer Key_Adj_Prp_Step_3D
  real(kind=FT) Adj_Prp_Step_3D_Max
  real(kind=FT) Prp_3D_Factor_m
  real(kind=FT) Prp_Bisec_Factor_3D
  integer Key_Smooth_Front
  integer Key_Smooth_Front_Twice
  integer Key_Num_Process

  integer Key_Ave_Stress
  real(kind=FT) S1_Point_Factor
  integer Key_Ave_Half_HF
  real(kind=FT) a_Ave_Shi
  real(kind=FT) Factor_Ave_R
  integer num_Gauss_S1
  real(kind=FT) Prop_Angle_Alowed
  integer Key_HF_Del_Ne_Pres
  real(kind=FT) Coff_a_AMPT
  integer,ALLOCATABLE::Material_Type(:)
  real(kind=FT),ALLOCATABLE:: Material_Para(:,:)
  real(kind=FT),ALLOCATABLE:: Material_Para_Added(:,:)
  real(kind=FT),ALLOCATABLE:: Mat_Cartesian_Coor_Vector_x(:,:)
  real(kind=FT),ALLOCATABLE:: Mat_Cartesian_Coor_Vector_y(:,:)
  real(kind=FT),ALLOCATABLE:: Mat_cylinder_Coor_Center(:,:)
  real(kind=FT),ALLOCATABLE:: Mat_cylinder_Coor_Vector_z(:,:)
  real(kind=FT) Desired_KIc

  real(kind=FT) Material_Interface(2)
  integer Key_SIFs_Method,Key_Force_Control,Key_Analysis_Type
  integer Key_SIFs_DIM_Points
  integer Key_SIFs_DIM_Method
  integer Key_FS_Seque_Coupling
  integer Key_TS_Seque_Coupling
  integer Num_Force_Divs
  integer Key_Save_f_d_Curve
  integer f_d_Curve_node_num
  integer Key_Save_f_COD_Curve
  integer f_COD_Curve_Crack_Num
  integer num_Crack_Log(Max_Itera)
  integer Key_Gravity
  real(kind=FT) g_X_Y_Z(3)
  integer Key_Geo_Stress
  integer Key_SLOE
  integer Key_EBE_Precondition
  integer Key_EBE_Condition_Number
  integer Key_EBE_Sym_Storage_K
  real(kind=FT),ALLOCATABLE::EBE_Condition_Number(:)
  integer Lis_Number_Iteration
  integer MDOF_2D
  integer MDOF_3D
  integer Old_MDOF_3D
  integer Key_Kink_Point
  integer Key_Data_Format
  integer Key_XA
  integer Key_K_Sparse
  real(kind=FT) Sparse_Ratio
  integer Sparse_Store_Method
  integer Key_Unit_System
  integer Key_Dimension
  integer Key_Contact
  real(kind=FT) Contact_Aper_Tol
  integer Key_TipEnrich
                                                    ! 2: Keep only the first item (default for dynamic analysis)
                                                    ! 3: Heaviside smooth transition scheme (see my PhD thesis for details)
                                                    !              4: Cohesive
  integer Key_CS_Natural_Crack
  real(kind=FT) Penalty_CS_Natural_Crack
  integer Key_Penalty_CS_Method
                                                    ! Then only the penalty function controls the displacement in the n direction, that is, the normal
                                                    ! direction displacement. NEWFTU2023082701.
  integer Num_F_Functions
  integer Key_Junction_Enrich
  integer Key_Tip_Fluid_Element_3D
  integer Key_Multi_TipEnrNode
  integer Key_Pre_Conditioner
  real(kind=FT) Key_TipEnrich_Radius_Factor
  integer Key_Local_Mesh_Refine
  logical Flag_Local_Refined
  integer Key_Front_Segmentation
  integer Number_front_Segments

  integer Key_Adj_Ini_Crack_3D
  integer Key_Check_and_Adjust_Cracks_3D
  real(kind=FT) Adjust_Cracks_Delta_Angle
  integer Adjust_Cracks_Resolution
  real(kind=FT) k_tt
  real(kind=FT) k_nn
  integer Key_Min_Aperture
  real(kind=FT) Min_Aperture
  ! Contact analysis related
  real(kind=FT) fric_mu_Cont
  integer Max_Contact_Iter
  integer Key_Conta_ConCrit
  real(kind=FT) Aper_Tol_Factor
  real(kind=FT) kn_Cont_Penalty,kt_Cont_Penalty
  real(kind=FT) Conve_Tol_Penalty
  ! Cohesive crack tips and cohesive cracks related
  integer Max_Cohesive_Iter
  integer Coh_Integ_Point
  integer Key_Coh_ConCrit
  integer Key_Play_Sounds
  integer Key_Memo_Moni
  integer Key_Clear_All
  integer Key_Junction_Check
  real(kind=FT) Factor_Check_Dis
  real(kind=FT) Factor_L_Jun_NewCr
  !----------------------------------------
  real(kind=FT) HF_Ini_Pressure
  !----------------------------------------
  integer Key_Post_CS_N_Strs
  integer Key_Post_CS_N_Stra
  integer Key_Save_avol_file
  integer Key_Post_CS_N_Stra_Cylindr
  integer Key_Post_CS_N_Strs_Cylindr
  integer Key_Post_CS_N_Disp_Cylindr
  ! To calculate strains in cylindrical coordinates (such as circumferential strain), the following
  ! parameters need to be provided
  real(kind=FT) Post_cylinder_Coor_Center(3)
  real(kind=FT) Post_cylinder_Coor_Vector_x(3)
  real(kind=FT) Post_cylinder_Coor_Vector_y(3)
  real(kind=FT) Post_cylinder_Coor_Vector_z(3)
  real(kind=FT),ALLOCATABLE::Theta_Cartesian_to_Cylinder(:,:)
  real(kind=FT),ALLOCATABLE::Theta_Cartesian_to_Cylinder_Node(:)
  real(kind=FT) Water_Pressure
  integer Type_Water_Pressure
  real(kind=FT) Num_Div_Elment
  integer Key_InSitu_Method
  integer Key_InSitu_Strategy
  integer Key_Read_Initial_Node_Stress_File
  real(kind=FT) InSitu_x
  real(kind=FT) InSitu_y
  real(kind=FT) InSitu_z
  real(kind=FT) InSitu_xy,InSitu_yz,InSitu_xz
  integer Key_InStress_for_Mat
  integer Mat_Number_of_InStress(100)
  real(kind=FT) Mat_InStress_x
  real(kind=FT) Mat_InStress_y
  real(kind=FT) Mat_InStress_z
  integer Key_HF_Secant_TS
  integer Key_HF_LineSearch
  integer Key_HF_Multistage
  integer Key_HF_Multistage_3D
  integer Key_HF_MS_Contc_Check
  integer Key_HF_Cont_Scheme
                                                    ! =2, perform contact detection only during the first two fluid-structure interaction iterations;
                                                    ! =3, contact iterations are performed in each fluid-structure coupling iteration.
  character*1  Keywords_Blank
  parameter (Keywords_Blank  = ' ')
  character*4 Space_4
  parameter  (Space_4  = '    ')
  character*8 Space_8
  parameter  (Space_8  = '        ')
  ! Optimization related to specific papers
  integer Key_Paper1
  integer Type_Cal_Propped
  ! integer Key_No_Tip_Enrich                        !Whether to perform tip fracture enhancement
  ! Related to ground stress
  integer State_InSitu
  ! ---------Fatigue Analysis Related--------------
  integer Key_Static_Fatigue
  integer Key_Fatigue_Cri
  integer Num_Fatigue_N
  real(kind=FT) Fatigue_Paris_C,Fatigue_Paris_m
  ! ---------Initial Pore Pressure Related----------
  integer Key_PoreP
  real(kind=FT) Initial_PoreP
  real(kind=FT) Initial_Biot_Coeff
  ! -------The following is used to support crack width calculation,
  ! Key_Propped_Width=1------------------
  integer Key_Propped_Width
                                                    ! W_P_0 (wpnp file) at each crack calculation point, used only for analyzing type 1
  real(kind=FT)Prop_Size
  real(kind=FT)Prop_Yita
  real(kind=FT)Prop_Elas
  real(kind=FT)Prop_Poss
  real(kind=FT)Prop_D
  ! real(kind=FT) Prop_C1, Prop_C2                      ! C1 and C2 in the constitutive formula, see Paper03 for details
  logical Yes_XFEM
  ! Related to thermal stress calculation
  integer Key_Thermal_Stress
  integer Key_Initial_Temperature
  integer Key_Scheme_Thermal_Stress
  real(kind=FT) Thermal_Str_Temper(100)
  integer Key_Random
                          ! 1: The system's built-in auto-generated function random_number generates different numbers each
                          ! time.
                          !             2:Handbook of Simulation: Principles, Methodology, Advances, Applications
                          ! The seed scheme in the book, each generation is controlled by the seed
  real(kind=FT) Inject_Prop_c
  integer Key_Plasticity
  ! Others
  integer Key_Close_Window
  integer Key_Visit_PhiPsi_top
  integer Key_Window_Log
  character(256) Window_log_filename
  character(256) PhiPsi_Current_Directory
  character(256) Name_of_Keywords
  integer Seed
  integer Key_Damage
  real(kind=FT),ALLOCATABLE:: Ele_Damage(:,:)
  real(kind=FT),ALLOCATABLE:: Damage_Gauss(:)
  real(kind=FT) Material_Dam2Frac_Value
  integer Crack_Gen_from_Damage
  ! real(kind=FT), ALLOCATABLE:: Ele_Damage_D(:,:,:,:)! D matrix after considering damage at each
  ! Gaussian point of each element
  logical file_Sparse_K_Location_COO_bin
  logical file_Sparse_K_Mask_COO_bin

  integer Key_InPlane_Growth
  integer Key_Stop_Outside_Crack
  integer Key_3D_HF_Time_Step_Method
  !********************************
  ! Life-and-Death element related
  !********************************
  integer Key_EKILL
  real(kind=FT) EKILL_Weaken_Factor
  integer Ele_Killed_Each_Load_Step(Max_Step,1000)
  ! Feedback related to ALLOCATA and DEALLOCATA errors
  integer allo_STAT
  CHARACTER(len=80)::err_msg
  !**********************
  !About elements break.
  !**********************
  integer Key_Element_Break
  integer Key_Element_Break_Rule
  integer Break_Mat_Num
  !**********************************
  !About EBE-PCG solver, 2021-07-31.
  !**********************************
  real(kind=FT) cg_tol
  integer       max_cg
  !parameter (cg_tol = 1.0D-7 )   !Default: 1.0D-7
  !parameter (cg_tol = 1.0D-6 )   !Default: 1.0D-7
  !parameter (cg_tol = 1.0D-5 )   !Default: 1.0D-7
  parameter (max_cg = 5000)
  !*************************************************
  ! 3D crack leading edge S1 or K smooth treatment.
  !*************************************************
  integer Key_Denoise_Vertex_Value
  integer Key_Smooth_Vertex_Value
  integer Smooth_Vertex_n
  integer Key_Smooth_Vertex_Value2
  integer Smooth_Vertex_n2
  !*********************************************
  ! 3D crack front Theta smoothing. 2022-07-14.
  !*********************************************
  integer Key_Denoise_Theta_Value
  integer Key_Smooth_Theta_Value
  integer Smooth_Theta_n
  integer Key_Smooth_Theta_Value2
  integer Smooth_Theta_n2
  !************************************************************
  ! 3D crack front GrowthFactor smoothing process. 2022-07-14.
  !************************************************************
  integer Key_Denoise_GF_Value
  integer Key_Smooth_GF_Value
  integer Smooth_GF_n
  integer Key_Smooth_GF_Value2
  integer Smooth_GF_n2
  
  !********
  ! Others
  !********
  integer Key_Crack_Inner_Pressure
  integer Key_Block_Model
  integer Flag_HF_3D
  integer Key_Cpp_Call_Fortran_Lib
  integer XA_Step_Count
  integer Key_Save_Crack_Radius
  integer Circle_3D_Eqv_Polygon_Resolution
  integer Key_Non_Negtive_Aperture_3D
  integer Key_Print_EBEPCG_Solution_Time
  integer Key_Scheme_Signed_Dis_InPlane
  integer Key_Sawtooth_Crack_Path
  integer Key_User_Defined_2D_Crack_Path
  integer Num_User_Defined_2D_Crack_Path
  
end module Global_Common

!-----------------------------------
! 2. Model-Related Global Variables
!-----------------------------------
module Global_Model
  use Global_Float_Type
  use Global_Ragged_Array_Real_Classs
  use Global_Ragged_Array_Int_Classs
  implicit none
  save
  integer Num_Node, Num_Elem, num_of_Material
  integer Num_Bou_x,Num_Bou_y,Num_Bou_z,Num_Foc_x,Num_Foc_y,Num_Foc_z
  integer Num_Boux_nonzero,Num_Bouy_nonzero,Num_Bouz_nonzero
  real(kind=FT) Min_X_Coor,Max_X_Coor,Min_Y_Coor,Max_Y_Coor,Min_Z_Coor,Max_Z_Coor
  !Model coordinates range.
  real(kind=FT) Model_X_Range,Model_Y_Range,Model_Z_Range
  real(kind=FT) Max_Model_Range,Min_Model_Range
  real(kind=FT),ALLOCATABLE::Coor(:,:)
  integer Max_Materials
  integer Model_Nx,Model_Ny,Model_Nz
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! cccccccccccccccccc         Program Size Control Parameters        ccccccccccccccccccccccc
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  parameter (Max_Materials     = 100)
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! cccccccccccccccccc         Program Size Control Parameters        ccccccccccccccccccccccc
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer,ALLOCATABLE::Bou_x(:)
  integer,ALLOCATABLE::Bou_y(:)
  integer,ALLOCATABLE::Bou_z(:)
  real(kind=FT),ALLOCATABLE::Foc_x(:,:)
  real(kind=FT),ALLOCATABLE::Foc_y(:,:)
  real(kind=FT),ALLOCATABLE::Foc_z(:,:)
  real(kind=FT),ALLOCATABLE::Bou_x_nonzero(:,:)
  real(kind=FT),ALLOCATABLE::Bou_y_nonzero(:,:)
  real(kind=FT),ALLOCATABLE::Bou_z_nonzero(:,:)
  integer,ALLOCATABLE::Elem_Node(:,:)
  integer,ALLOCATABLE::Elem_Node_Num(:,:)
  integer,ALLOCATABLE::Node_Elements(:,:)
  ! integer, ALLOCATABLE :: Node_Elements(:,:)     ! The elements surrounding each node.
  type(Ragged_Int_Array_1D),allocatable::Node_Elements_3D(:)
  integer,ALLOCATABLE::num_Node_Elements(:)
  integer,ALLOCATABLE::Ele_Elements(:,:)
  integer,ALLOCATABLE::num_Ele_Eles(:)
  integer,ALLOCATABLE::Element_Edges(:,:,:)
  integer Max_Diff_Elem_Num
  integer Max_Half_Band_Width
  integer,ALLOCATABLE::Elem_Mat(:)
  integer,ALLOCATABLE::Outline(:,:)
  integer,ALLOCATABLE::OutArea(:,:)
  integer Num_Surface_Nodes
  integer,ALLOCATABLE::Surface_Nodes(:)
  integer,ALLOCATABLE:: G_NN(:,:)
  real(kind=FT),ALLOCATABLE:: G_X_NODES(:,:),G_Y_NODES(:,:),G_Z_NODES(:,:)
  real(kind=FT),ALLOCATABLE:: Elem_Area(:),Elem_Vol(:),Elem_Max_L(:),Elem_Min_L(:),Elem_Ave_L(:)
  real(kind=FT),ALLOCATABLE:: Node_Max_L(:)
  logical,ALLOCATABLE:: Elem_Break(:)
  real(kind=FT),ALLOCATABLE:: Elem_Ave_Gauss_Stress(:,:)
  real(kind=FT),ALLOCATABLE:: Elem_Ave_Gauss_S1(:)
  real(kind=FT),ALLOCATABLE:: Elem_Centroid(:,:)
  logical*1,ALLOCATABLE:: EleGaus_yes_FEM_asemd(:,:)
  integer Total_FD,Usual_Freedom,Enrich_Freedom
  integer Total_FD_P,Total_FD_U
  integer,ALLOCATABLE:: Flag_FreeDOF(:)
  integer,ALLOCATABLE:: Location_FreeDOF(:)
  integer,ALLOCATABLE:: FreeDOF_for_Det(:)
  integer,ALLOCATABLE:: Real_FreeDOF_Index(:)
  real(kind=FT) k_PenaltyStiffness
  ! -------------Gauss Point Numbering Related---------------
  ! The following two global variables are obtained during the assembly of the stiffness matrix.
  integer,ALLOCATABLE::num_GP_Elem(:)
  integer,ALLOCATABLE::Ele_GP_Start_Num(:)
  ! -------------Gauss Point Coordinates---------------
  real(kind=FT),ALLOCATABLE:: Gauss_CoorX(:)
  real(kind=FT),ALLOCATABLE:: Gauss_CoorY(:)
  ! -------------Rupture Zone---------
  integer Key_Fracture_Zone
  real(kind=FT) Frac_Zone_MinX
  real(kind=FT) Frac_Zone_MaxX
  real(kind=FT) Frac_Zone_MinY
  real(kind=FT) Frac_Zone_MaxY
  real(kind=FT) Frac_Zone_MinZ
  real(kind=FT) Frac_Zone_MaxZ
  ! -------------Initial Crack Generation Area.2024-06-23---------
  integer Key_Ini_Crack_Zone              
  real(kind=FT) Ini_Crack_Zone_MinX,Ini_Crack_Zone_MaxX          
  real(kind=FT) Ini_Crack_Zone_MinY,Ini_Crack_Zone_MaxY       
  real(kind=FT) Ini_Crack_Zone_MinZ,Ini_Crack_Zone_MaxZ         
  ! -------------Random Natural Cracks---------
  integer Key_Random_NaCr
  integer num_Rand_Na_Crack
  real(kind=FT) NaCr_Orientation
  real(kind=FT) NaCr_Ori_Delta
  real(kind=FT) NaCr_Length
  real(kind=FT) NaCr_Len_Delta
  real(kind=FT) Random_NaCr_Rad_Factor
  ! -------------Initial ground stress of the model
  ! -------------Node Coupling Related
  ! option1 - Define coupling nodes directly through a keyword file (in this method, only one set of
  ! couplings can be defined per direction)
  integer num_CP_x_nodes,num_CP_y_nodes
  integer CP_x_nodes(1:5000)
  integer CP_y_nodes(1:5000)
  ! option2 - Define coupled nodes through dofx, dofy, dofz files (in this method, only multiple sets
  ! of couplings can be defined for each direction)
  integer num_CP_set_x,num_CP_set_y
  integer num_nodes_CP_set_x(10),num_nodes_CP_set_y(10)
  integer CP_nodes_x(10,5000),CP_nodes_y(10,5000)

  ! -------------Range of x, y, z coordinates for each element--------
  real(kind=FT),ALLOCATABLE:: x_max_Elements(:)
  real(kind=FT),ALLOCATABLE:: x_min_Elements(:)
  real(kind=FT),ALLOCATABLE:: y_max_Elements(:)
  real(kind=FT),ALLOCATABLE:: y_min_Elements(:)
  real(kind=FT),ALLOCATABLE:: z_max_Elements(:)
  real(kind=FT),ALLOCATABLE:: z_min_Elements(:)
  real(kind=FT) penalty_k_bou_nonzero

  integer D3_nband_FEM

  integer Num_Elem_Block_Bou
  integer,ALLOCATABLE::Elems_Block_Bou(:)
  integer Ele_3D_Edges_Node(12,2)
  ! -------------Related to Gaussian Integrals--------
  integer,ALLOCATABLE::Elements_Gauss_Num(:)
  integer,ALLOCATABLE::Elems_Integration_Type(:,:)
  integer,ALLOCATABLE::Elems_Num_SubEles(:,:)
  integer,ALLOCATABLE::Elems_Type_SubEles(:,:)
  integer,ALLOCATABLE::Elems_SubEles_Index(:,:)
  integer num_SubEles
  integer,ALLOCATABLE::SubEles_Integ_Num(:)
  real(kind=FT),ALLOCATABLE::SubEles_Integ_Coors(:,:,:)
  integer Num_Max_3D_gauss
  !-------------------------------------------------------------------------------------------------
  ! integer Ele_Num_Cache_by_Coors_3D ! Cache, stores the most recently obtained element numbers by
  ! coordinates. 2022-09-24.
  !-------------------------------------------------------------------------------------------------
  !2022-11-21.
  ! Divide the model into 8 regions based on coordinates. Save the element numbers of the 8 regions.
  ! This makes it easier to find element numbers based on coordinates. IMPROV2022112101.
  real(kind=FT) Model_Center_x,Model_Center_y,Model_Center_z
  integer Domain_Elements_Num(8)
  integer,ALLOCATABLE::Ele_Domain_ID(:)
  integer,ALLOCATABLE::Domain_Elements(:,:)
  integer First_XFEM_Step
  ! element list corresponding to each material. 2023-01-24. NEWFTU2023012401.
  type(Ragged_Int_Array_1D),allocatable::List_Elements_Mat(:)
  integer Elements_Num_Mat(Max_Materials)
  ! The initial and current temperature of the element. 2023-03-13.
  real(kind=FT),ALLOCATABLE:: Elem_Initial_T(:),Elem_Current_T(:),Elem_T_for_Stress(:)
  ! Core global variables. Size_Local_3D, All_Local_3D. 2023-03-15. IMPROV2023031501.
  integer,ALLOCATABLE::Size_Local_3D(:)
  integer,ALLOCATABLE::All_Local_3D(:,:)
  !
  ! Initial pore pressure of the element, current pore pressure, Biot coefficient. 2023-03-19.
  ! IMPROV2023031901.
  real(kind=FT),ALLOCATABLE::Elem_Initial_PoreP(:),Elem_Current_PoreP(:),Elem_Biots(:)
  ! element's elastic modulus, Poisson's ratio, coefficient of thermal expansion, fracture toughness,
  ! tensile strength. 2023-03-19. IMPROV2023031903.
  real(kind=FT),ALLOCATABLE::Elem_E_XA(:),Elem_Mu_XA(:),Elem_TEC_XA(:),Elem_KIc_XA(:),Elem_St_XA(:)
  ! D matrix of the element. 2023-03-19. IMPROV2023031904.
  real(kind=FT),ALLOCATABLE::Elem_D_XA(:,:,:)
  ! Whether the running crack is located outside the model (boundary crack). Default is 0.
  ! NEWFTU2023050701.
  integer Key_Allow_3D_Outside_Crack
  ! Number of independent nodes in each element. 2023-06-14. IMPROV2023061402.
  integer,ALLOCATABLE::Elem_Uniqued_Nodes(:)
  logical,ALLOCATABLE::Yes_Degenarated_Elem(:)
  integer Num_Degenarated_Elems
  integer MAT_ALLOW_CRACK_Initiation(Max_Materials)
  integer Key_Max_Num_Initiation_Cracks
  integer Num_Initiation_Cracks 
end module Global_Model

!----------------------------------------------
! 2.1 enhancement element Related. 2022-06-24.
!----------------------------------------------
module Global_XFEM_Elements
  ! ---------------FEM Element List and XFEM Enhanced Element List, 2022-04-16------------
  integer num_FEM_Elem,num_XFEM_Elem
  integer,ALLOCATABLE::FEM_Elem_List(:)
  integer,ALLOCATABLE::XFEM_Elem_List(:)
  integer,ALLOCATABLE::Elem_XFEM_Flag(:)
  integer,ALLOCATABLE::Elem_Location(:,:)
  integer,ALLOCATABLE::Elem_New_XFEM_Flag(:)
  integer,ALLOCATABLE::Elem_Update_XFEM_Flag(:)
  integer,ALLOCATABLE::Elem_XFEM_Flag_Old(:)
  integer,ALLOCATABLE::Elem_Location_Old(:,:)
  integer,ALLOCATABLE::size_local_Old(:)

  !BUGFIX2022092602.
  integer Must_Gauss_Number_3D
  parameter (Must_Gauss_Number_3D = 50)
  !2023-02-15.
  integer,ALLOCATABLE::Rollbacked_FEM_Elements(:)
  integer Num_Rollbacked_FEM_Elements
end module Global_XFEM_Elements


!-------------------------------------------
! 3. Global variables related to file names
!-------------------------------------------
module Global_Filename
  implicit none
  save
  character(256) Filename,Work_Directory
  character(256) Full_Pathname
  character(256) Python_Directory
  character(256) PhiPsi_Directory
end module Global_Filename

!---------------------------------------------------------------
! 4. Dynamic analysis related (including implicit and explicit)
!---------------------------------------------------------------
module Global_Dynamic
  use Global_Float_Type
  implicit none
  save
  integer Num_Ivex
  integer Num_Ivey
  integer Num_Ivez
  real(kind=FT),ALLOCATABLE::Ive_x(:,:)
  real(kind=FT),ALLOCATABLE::Ive_y(:,:)
  real(kind=FT),ALLOCATABLE::Ive_z(:,:)
  integer Num_Iacx
  integer Num_Iacy
  integer Num_Iacz
  real(kind=FT),ALLOCATABLE::Iac_x(:,:)
  real(kind=FT),ALLOCATABLE::Iac_y(:,:)
  real(kind=FT),ALLOCATABLE::Iac_z(:,:)
  integer IDy_Num_Iteras
  integer IDy_Num_force_Itr
  real(kind=FT) delt_time_NewMark
  integer Key_EQ
  real(kind=FT)  EQ_Ac_Time_Gap
  real(kind=FT),ALLOCATABLE::EQ_Accel_data(:)
  integer num_EQ_Accel
  integer num_EQ_Ac_nodes
  integer EQ_Ac_nodes(5000)
  ! Sine acceleration excitation related
  integer Key_Sin_Accel
  integer  Sin_Accel_Dire
  real(kind=FT) Sin_Accel_A
  real(kind=FT) Sin_Accel_T
  integer  Sin_Accel_num_Nodes
  integer Sin_Accel_Nodes(5000)

  real(kind=FT) Factor_Prop_Dy
  integer EDy_Num_Iteras
  integer EDy_Num_force_Itr
  real(kind=FT) Delt_Time_Explicit
  integer Key_Mass_Lumped
  real(kind=FT) Explicit_time_inc


  real(kind=FT),ALLOCATABLE::EDy_DISP(:),EDy_VELC(:),EDy_ACCL(:)
  real(kind=FT),ALLOCATABLE::IDy_DISP(:),IDy_VELC(:),IDy_ACCL(:)
end module Global_Dynamic

!----------------------------------------------------------------------------------------------------
! 5.2D fracture-related, including enhanced nodes, calculation points, and crack connectivity
! The largest global variables in the program; in addition, it also includes hole-related variables.
!----------------------------------------------------------------------------------------------------
module Global_Crack
  use Global_Float_Type
  implicit none
  save
  integer  Max_Num_Cr,Max_Num_Arc_Cr,Max_Num_Cr_P
  integer  Max_Num_Cr_CalP,Max_Num_Seg_CalP,Max_Num_Cone_Cr
  integer  Max_Num_Ele_QuadHF,Max_Num_Ele_CalP,Max_Num_Hl
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! cccccccccccccccccc         Program Size Control Parameters        ccccccccccccccccccccccc
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! To have gfortran support 150 cracks, you need: -fopenmp -fno-align-commons -fno-range-check
  ! -fmax-stack-var-size=100000 ! Set 100000 higher to support more cracks (2021-08-19)
  ! parameter (Max_Num_Cr = 150) !Up to 150 cracks, 100 is fine for gfortran, but causes issues with
  ! Intel Fortran; 50 is fine for all
  ! parameter (Max_Num_Cr = 100) ! Up to 150 cracks, 100 works fine in gfortran, but causes issues in
  ! Intel Fortran; 50 works fine
  parameter (Max_Num_Cr        = 50   )
  parameter (Max_Num_Cr_P      = 200   )
  parameter (Max_Num_Cr_CalP   = 300 )
  parameter (Max_Num_Seg_CalP  = 200  )
  parameter (Max_Num_Cone_Cr   = 10   )
  !---------
  parameter (Max_Num_Arc_Cr    = 100   )
  !---------
  parameter (Max_Num_Ele_QuadHF= 100  )
  !---------
  parameter (Max_Num_Ele_CalP  = 5   )
  parameter (Max_Num_Hl        = 100   )
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! cccccccccccccccccc         Program Size Control Parameters        ccccccccccccccccccccccc
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !****************************************************************************************
  ! Related to enriched nodes, the meaning of the variables can be found in the subroutine
  ! Determine_Enriched_Nodes
  !****************************************************************************************
  real(kind=FT) Crack_Coor(Max_Num_Cr,Max_Num_Cr_P,2)
  real(kind=FT) Arc_Crack_Coor(Max_Num_Cr,Max_Num_Cr_P-1,11)
                                                                ! x, y, Direction (1 for counterclockwise, -1 for clockwise), r, Radian_Start, Radian_End, Radian,
                                                                !Point_Start_x,Point_Start_y,Point_End_x,Point_End_y;
                                                                ! When entering, you only need to input the coordinates of the center (x, y) and the direction.
                                                                ! Variables 4, 5, 6, 7, 8, 9, 10, and 11 will be calculated through the
                                                                ! Tool_Cal_Arc_r_and_Radian_Given_Coors subroutine.
                                                                ! Note: The starting arc angles Radian_Start and Radian_End use 3 o'clock as 0 degrees, measured
                                                                ! counterclockwise from 0 to 360 degrees
  ! integer Arc_Crack_Passed_Ele(Max_Num_Cr, Max_Num_Cr_P-1, 1000) ! Element number passed through by
  ! the arc crack


  real(kind=FT) Hole_Coor(Max_Num_Hl,3)
  real(kind=FT) Ellip_Hole_Coor(Max_Num_Hl,5)
  real(kind=FT) Na_Crack_Coor(Max_Num_Cr,Max_Num_Cr_P,2)
  real(kind=FT) Cr_First_Tip(Max_Num_Cr,2),Cr_Second_Tip(Max_Num_Cr,2)
  real(kind=FT) Cr_First_Tip_Ori(Max_Num_Cr),Cr_Second_Tip_Ori(Max_Num_Cr)


  integer,ALLOCATABLE:: Elem_Type(:,:)
  integer ,ALLOCATABLE:: c_POS(:,:)
  integer,ALLOCATABLE:: Enriched_Node_Type(:,:)
  real(kind=FT),ALLOCATABLE::Enriched_Node_Crack_n_Vector(:,:,:)
  integer,ALLOCATABLE:: Node_Jun_elem(:,:)
  integer,ALLOCATABLE:: Jun_Ele_Negative_Cr_Num(:,:)
  integer,ALLOCATABLE:: Elem_Type_Hl(:,:)
  integer,ALLOCATABLE:: Enriched_Node_Type_Hl(:,:)
  integer ,ALLOCATABLE:: c_POS_Hl(:,:)

  real(kind=FT),ALLOCATABLE::  Coors_Element_Crack(:,:,:)
  real(kind=FT),ALLOCATABLE::  Coors_Tip(:,:)
  real(kind=FT),ALLOCATABLE::  Coors_Vertex(:,:)
  real(kind=FT),ALLOCATABLE::  Coors_Junction(:,:,:)
  real(kind=FT),ALLOCATABLE::  x_cr_tip_nodes(:,:)
  real(kind=FT),ALLOCATABLE::  y_cr_tip_nodes(:,:)
  integer,ALLOCATABLE::  Ele_Num_Tip_Enriched_Node(:,:)
  integer Crack_Tip_Type(Max_Num_Cr,2)
  real(kind=FT) Crack_Arc_Tip_A_B_C_x(Max_Num_Cr,2,3)
  real(kind=FT) Crack_Arc_Tip_A_B_C_y(Max_Num_Cr,2,3)
  integer Crack_Jun_CrNum(Max_Num_Cr,2)
  integer Crack_Jun_HoleNum(Max_Num_Cr,2)
  integer,ALLOCATABLE:: Node_Jun_Hole(:,:)
  integer,ALLOCATABLE:: Ele_Jun_Hole(:,:)
  integer,ALLOCATABLE:: TipEle_Adjacent_Ele(:,:)
  integer Crack_Jun_Elem(Max_Num_Cr,2)
  real(kind=FT) Crack_Tip_Coor(Max_Num_Cr,2,2)
  real(kind=FT) Edge_Disposed_Crack(Max_Num_Cr,Max_Num_Cr_P,2)
  logical Flag_Crack_Tip_Out_Mol(Max_Num_Cr,2)
  logical Yes_Arc_Crack
  integer,ALLOCATABLE::  Cracks_CalP_Num(:)
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_Coors(:,:,:)
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_Orient(:,:)
  integer,ALLOCATABLE::   Cracks_CalP_Seg(:,:)
  integer,ALLOCATABLE::   Cracks_CalP_Elem(:,:)
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_Aper(:,:)
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_Pres(:,:)
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_Tractions(:,:,:)
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_Pgra(:,:)
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_Velo(:,:)
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_Quan(:,:)
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_Conc(:,:)
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_Remo_Strs(:,:)
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_Contact_Strs(:,:,:)
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_Contact_Force_x(:,:)
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_Contact_Force_y(:,:)
  ! -----Added on 2016-08-25--Related to the calculation of support crack width and flow guidance
  ! capacity-------
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_wpnp(:,:)
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_wpor(:,:)
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_wdeform(:,:)
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_Conductivity(:,:)
  real(kind=FT),ALLOCATABLE::  Cracks_CalP_kf(:,:)
  ! -----------The following variables are related to staged fracturing----------
  integer,ALLOCATABLE::  MS_Cracks_CalP_Num(:,:)
  real(kind=FT),ALLOCATABLE::  MS_Cracks_CalP_Aper(:,:,:)
  real(kind=FT),ALLOCATABLE::  MS_Cracks_CalP_Conc(:,:,:)
  real(kind=FT),ALLOCATABLE::  MS_CalP_Propped_Aper(:,:)

  !**************************************************************************************************
  ! This corresponds to the crack calculation point data for the previous rupture step
  ! Do not allocate memory space initially. If the number of break steps is greater than or equal to
  ! 2, the program will allocate memory space.
  !**************************************************************************************************
  integer,ALLOCATABLE::       L_Cracks_CalP_Num(:)
  real(kind=FT),ALLOCATABLE:: L_Cracks_CalP_Pres(:,:)
  real(kind=FT),ALLOCATABLE:: L_Cracks_CalP_Aper(:,:)
  real(kind=FT),ALLOCATABLE:: L_Cracks_CalP_Coors(:,:,:)
  real(kind=FT),ALLOCATABLE:: L_Cracks_CalP_Orient(:,:)
  integer,ALLOCATABLE::       L_Cracks_CalP_Seg(:,:)
  integer,ALLOCATABLE::       L_Cracks_CalP_Elem(:,:)
  real(kind=FT),ALLOCATABLE:: L_Cracks_CalP_Pgra(:,:)
  real(kind=FT),ALLOCATABLE:: L_Cracks_CalP_Velo(:,:)
  real(kind=FT),ALLOCATABLE:: L_Cracks_CalP_Quan(:,:)
  real(kind=FT),ALLOCATABLE:: L_Cracks_CalP_Conc(:,:)
  real(kind=FT),ALLOCATABLE:: Map_L_Cracks_CalP_Conc(:,:)
  !----------
  integer Cracks_CalP_Type(Max_Num_Cr,Max_Num_Cr_CalP,2)
                                                         ! (1) Calculate point type, Cracks_CalP_Type(i_C, i_CalP, 1)
                                                         ! 0 --- Normal
                                                         ! 1 --- Crack (i_C) crack tip 1 shares a computation point with the junction of another crack
                                                         ! (other_C);
                                                         ! 2 --- Crack (i_C) crack tip 2 shares calculation points with the junction of another crack
                                                         ! (other_C);
                                                         ! 3 --- A calculation point in the middle of crack (i_C) shares a junction with another crack
                                                         ! (other_C);
                                                         ! 4 --- A calculation point shared at the intersection of a certain point in the middle of crack
                                                         ! (i_C) and another crack (other_C);
                                                         ! 5 --- Crack (i_C) ordinary calculation point when the crack tip is not connected to other cracks,
                                                         ! in this case
                                                         ! Constrain the HF pressure boundary to 0, see the subroutine Boundary_Cond_HF for details.
                                                         ! (2) Corresponding crack number other_C, Cracks_CalP_Type(i_C, i_CalP, 2)
                                                         ! Note: When Cracks_CalP_Type(i_C, i_CalP, 1) = 5, it is not required
  integer Cracks_TipJunCalpNum(Max_Num_Cr,2)
  integer Cracks_MidJunCalpNum(Max_Num_Cr,Max_Num_Cr_CalP)
  ! ----------Calculate the global numbering of points and the corresponding local
  ! numbering---------------------
  ! It includes both water-containing cracks and dry cracks that do not participate in fluid-solid
  ! coupling.
  integer Cracks_GloNumCalP(Max_Num_Cr,Max_Num_Cr_CalP)
  integer Cracks_LocalNumCalP(Max_Num_Cr*Max_Num_Cr_CalP,2)
                                                               ! 1: Crack Number
                                                               ! 2: Local Calculation Point Number
  ! It includes both water-containing cracks and dry cracks that do not participate in fluid-solid
  ! coupling.
  integer Cracks_GloNumCalP_W(Max_Num_Cr,Max_Num_Cr_CalP)
  integer Cracks_LocalNumCalP_W(Max_Num_Cr*Max_Num_Cr_CalP,2)
                                                               ! 1: Crack Number
                                                               ! 2: Local Calculation Point Number
  !----------
  integer Num_JunPair
  integer Cracks_JunPair(Max_Num_Cr,2)
  !***********************************************************************
  ! Unique to advanced hydraulic fracturing elements (secondary elements)
  !***********************************************************************
  integer Ele_Nodes_QuadHF(Max_Num_Cr,Max_Num_Ele_QuadHF,3)
  !**************************************
  ! Hydraulic fracturing element Related
  !**************************************
  real(kind=FT) Cracks_HF_Ele_L(Max_Num_Cr,Max_Num_Cr_CalP-1)
  !***************************************************************
  ! Fracture connectivity (used only during hydraulic fracturing)
  !***************************************************************
  integer  Cracks_Cone_Num(Max_Num_Cr)
  integer  Cracks_Cone_Cr(Max_Num_Cr,Max_Num_Cone_Cr)
  ! --------------Fissure Tip Communication------------
  integer Cracks_Cone_NumTipCr(Max_Num_Cr)
  integer Cracks_Cone_TipCrNum(Max_Num_Cr,2)
                                                              ! Cracks_Cone_TipType(i_C,1) -- Connected crack number at the tip of crack i_C, tip 1
                                                              ! Cracks_Cone_TipType(i_C,2) -- Connected crack number at crack tip 2 of i_C
  integer Cracks_Cone_TipJuEle(Max_Num_Cr,2)
  real(kind=FT) Cracks_Cone_TipJuCor(Max_Num_Cr,2,2)
  ! --------------Central Unicom------------
  integer Cracks_Cone_NumMidCr(Max_Num_Cr)
  integer Cracks_Cone_MidCrNum(Max_Num_Cr,Max_Num_Cone_Cr)
  integer Cracks_Cone_MidCrTip(Max_Num_Cr,Max_Num_Cone_Cr)
  integer Cracks_Cone_MidJuEle(Max_Num_Cr,Max_Num_Cone_Cr)
  real(kind=FT) Cracks_Cone_MidJuCor(Max_Num_Cr, Max_Num_Cone_Cr,2)
  !*********************************
  ! Stress intensity factor related
  !*********************************
  real(kind=FT) KI(Max_Num_Cr,2),KII(Max_Num_Cr,2)

  !*******************************************************************************
  ! Hydraulic fracturing related natural fractures (see Check_Crack_Grows_MCSC.f)
  !*******************************************************************************
  integer  Cracks_NF_JiS_Cr_Num(Max_Num_Cr)
  integer  Cracks_NF_JiS_Stat(Max_Num_Cr,2)
  integer  Cracks_NF_T_Stat(Max_Num_Cr,2)
  integer  Cracks_NF_Cement_Tip_Num(Max_Num_Cr,2)
  integer  Cracks_fric_NF_num(Max_Num_Cr)
  integer  Cracks_QinS_Stat(Max_Num_Cr)
  !*************************
  ! Cohesive cracks related
  !*************************
  integer  Cracks_Coh_Ele_Type(Max_Num_Cr,Max_Num_Cr_CalP-1)
  integer  Cracks_Tip_Num_of_Coh_Ele(Max_Num_Cr,2)

  ! Crack and polygon inclusion interaction
  integer Crack_Tip_Ploy_Inc_Info(Max_Num_Cr,2,5)
                                                          ! 1. Record the polygon containing the inclusion; 2. Record the edge number of the polygon where the
                                                          ! crack tip is located
  ! REAL(kind=FT) Cracks_CalP_Coh_FN_FT(Max_Num_Cr, Max_Num_Cr_CalP, 2) !Calculate normal and
  ! tangential cohesion for each crack point
  integer Key_Crack_Aperture_Method
  
  real(kind=FT) Crack_Max_Min_Aperture(Max_Num_Cr,3)
  !*************************************************
  ! User-defined 2D crack path related. 2024-11-17.
  !*************************************************
  real(kind=FT),ALLOCATABLE:: User_Defined_2D_Crack_Path(:,:)
  integer,ALLOCATABLE:: User_Defined_2D_Crack_Path_Num_Points(:)
end module Global_Crack

!------------------------------------------
! 5.2 2D and 3D crack sharing. 2022-09-05.
!------------------------------------------
module Global_Crack_Common
  use Global_Float_Type
  implicit none
  save

  integer num_Crack
  integer num_Arc_Crack
  integer num_Hole,num_Circ_Hole,num_Ellip_Hole
  integer Each_Cr_Poi_Num(1000)
  integer num_Na_Crack
  integer Key_Na_Crack_Type
  integer Key_NaCr_Friction
  integer n_h_Node,n_t_Node,n_j_Node,n_hl_Node,n_c_Node
  !**********************************
  ! Randomly generated holes related
  !**********************************
  integer Key_Random_Hole
  integer num_Rand_Hole
  real(kind=FT) Rand_Hole_R
  real(kind=FT) Rand_Hole_R_Delta
  !*****************************
  ! Related to crack initiation
  !*****************************
  integer Key_Hole_Crack_Generate
  integer Num_Crack_Hole_Generated
  integer num_Hole_Crack_Generated(1000)
  integer Hole_Crack_Generated_num(1000,10)

  ! ------In-seam Pressure-------
  real(kind=FT) Crack_Pressure(1000)
  integer Crack_Pressure_Type

  ! integer Key_Cr_Pressure !Whether to apply fracture pressure: =1, constant pressure; =2, linear
  ! pressure (0 at the crack tip); =3, quadratic pressure
  ! REAL(kind=FT) Cr_Pressure_Value           !Magnitude of fracture pressure for dynamic analysis

  !*********************************************************************************
  ! Whether the fracture is driven by water (used only during hydraulic fracturing)
  !*********************************************************************************
  integer Cracks_HF_State(1000)
  integer Cracks_HF_Propp(1000)
  !******************************************************
  !Location of the injection point for the full HF model
  !******************************************************
  real(kind=FT) Inj_Point_Loc(2)
  integer Current_Inj_Crack
  integer CalP_num_InjP_Local
                                                             ! For symmetrical HF fracturing, this number is obviously equal to 1, so the symmetrical model does
                                                             ! not require this parameter.
  real(kind=FT) Cracks_HF_ConPressure(1000)
  !**************************************
  ! Whether cracks are allowed to extend
  !**************************************
  ! integer Cracks_Allow_Propa(Max_Num_Cr)                     !Whether each crack is allowed to propagate
  integer,ALLOCATABLE::Cracks_Allow_Propa(:)
  ! integer Cracks_Tips_Allow_Propa(Max_Num_Cr,2) !Whether each crack tip of each crack is allowed to
  ! propagate

  integer Each_Na_Cr_Poi_Num(1000)
  integer Key_CS_Crack(1000)

  !***************************************************
  ! Related to fracture calculation points (the so-called calculation points are the nodes of the
  ! hydraulic fracture elements)
  ! For non-hydraulic fracturing problems, it is about calculating the coordinate points of the crack
  ! opening.
  !-----------
  ! This section corresponds to the data of the current rupture step.
  !***************************************************
  integer num_Tol_CalP_Water
  integer num_Tol_CalP_Water_ConP
  integer num_Tol_CalP_All
  real(kind=FT) Total_Conductivity
  real(kind=FT) Ave_Conductivity
  ! ----------Used for crack surface contact iteration
  integer,ALLOCATABLE::Ele_NumCalP(:)
  integer,ALLOCATABLE::Ele_CalPNum(:,:)
  integer Conta_Integ_Point
  real(kind=FT) Norm2_Contact_R_PSI_0

  !2022-11-14.
  real(kind=FT),ALLOCATABLE:: Crack_Coor_Range(:,:,:)
  !2023-01-07
  integer Key_NaCr_Active_Scheme_3D
                                 ! = 1, at the initial moment all natural fractures are activated, and all natural fractures are
                                 ! equal to the actual natural fracture size. After HF communication, the fractures are filled with
                                 ! fracturing fluid. Default.
                                 ! = 2, activated only after communication with HF, the fractures are filled with fracturing fluid
                                 ! after HF communication.
                                 ! = 3, only activated after communication with HF; after HF communication, only some cracks open,
                                 ! extending along the planes of natural fractures. The fracture toughness of elements penetrated by
                                 ! natural fractures is low.
  real(kind=FT) Size_Factor_of_Active_NaCr
  real(kind=FT) KIc_NaCr(50000)
  real(kind=FT) St_NaCr(50000)
  integer Key_Ele_Max_Related_Cracks
  integer Max_Related_Cracks
  integer Key_Print_SIFs_to_Screen
end module Global_Crack_Common

!----------------------
! 5.3 3D Crack Related
!----------------------
module Global_Crack_3D
  use Global_Float_Type
  use Global_Ragged_Array_Real_Classs
  use Global_Ragged_Array_Int_Classs
  implicit none
  save
  ! 3D crack-related
  integer Max_Num_Cr_3D
  integer Max_Num_El_3D
  ! ----Original Default-----
  ! parameter (Max_Num_Cr_3D = 5)                 !Maximum of 10 cracks
  ! parameter (Max_N_Node_3D = 5000) !Up to 5000 discrete nodes per fracture surface (previously used
  ! 1000)
  ! parameter (Max_N_CalP_3D        = 20000 )                 !Each crack can contain up to 20,000 calculation points
  ! -----Reduce----
  ! parameter (Max_Num_Cr_3D = 5)                 !Maximum of 10 cracks
  ! parameter (Max_N_Node_3D = 1000) !Up to 5000 discrete nodes per fracture surface (previously used
  ! 1000)
  ! parameter (Max_N_CalP_3D        = 3000 )                 !Each crack contains a maximum of 20,000 calculation points
  ! -----Increase----
  ! parameter (Max_Num_Cr_3D = 50)                 !Maximum of 50 cracks
  ! parameter (Max_N_Node_3D = 5000) !Up to 5000 discrete nodes per fracture surface (previously used
  ! 1000)
  ! parameter (Max_N_CalP_3D        = 20000 )                 !Each crack contains a maximum of 20,000 calculation points
  ! -----Increase----
  parameter (Max_Num_Cr_3D        = 10000)
  !
  ! integer :: Max_N_Node_3D = 200 ! Initially, each fracture surface has a maximum of 200 discrete
  ! nodes. IMPROV2022110501. If insufficient, each fracture will automatically expand.
  integer :: Max_N_Node_3D(1:Max_Num_Cr_3D)  = 200
  integer Max_Max_N_Node_3D
  !
  ! integer :: Max_N_FluEl_3D = 200 !Initially, each fracture face has at most 200 fluid elements.
  ! IMPROV20221105011. If insufficient, each fracture will automatically expand.
  !
  integer :: Max_N_FluEl_3D(1:Max_Num_Cr_3D)  = 200
  integer Max_Max_N_FluEl_3D
  !
  ! integer :: Max_N_CalP_3D = 200 ! Initially, each fracture surface has up to 200 fluid nodes.
  ! IMPROV20221105021. If not enough, each fracture will automatically expand.
  integer :: Max_N_CalP_3D(1:Max_Num_Cr_3D)    = 200
  integer Max_Max_N_CalP_3D
  !
  integer :: Max_ele_num_CalP     = 100
  !---------

  !---------
  ! real(kind=FT) Crack3D_Coor(Max_Num_Cr_3D,4,3) ! Crack number, crack surface coordinate point
  ! number, (x, y, z)
  real(kind=FT),allocatable::Crack3D_Coor(:,:,:)
  real(kind=FT),allocatable::Na_Crack3D_Coor(:,:,:)
  real(kind=FT),allocatable::Na_Crack3D_St(:)
  real(kind=FT),allocatable::Na_Crack3D_KIc(:)
  real(kind=FT),allocatable::Na_Crack3D_Friction(:)
  integer,ALLOCATABLE:: Each_NaCr3D_Poi_Num(:)
  integer,ALLOCATABLE:: NaCr3D_Status(:,:)
                                                            ! Column 1 stores the active status;
                                                            ! Column 2 records the number of solid cells containing the natural fracture.
  type(Ragged_Int_Array_1D),allocatable::Na_Crack3D_Ele_List(:)
  ! type(Ragged_Array_2D)::Crack3D_Coor(Max_Num_Cr_3D) !Assign a Ragged_Array_2D class object to each
  ! crack (memory allocated after allocation). 2022-09-02.

  ! real(kind=FT), ALLOCATABLE:: Dis_Node_to_FS(:,:) ! Symbolic distance from nodes to each fracture
  ! surface, FS stands for Fracture Surface
  type(Ragged_Array_1D),allocatable::Dis_Node_to_FS(:)

  real(kind=FT),ALLOCATABLE::Vector_o_Orient(:,:)
  integer,ALLOCATABLE::Sign_o_Orient(:)

  ! real(kind=FT) KI_3D(Max_Num_Cr_3D,Max_N_Node_3D) ! Stress intensity factor at each boundary point
  ! of each crack surface
  !NEWFTU2022090201.
  type(Ragged_Array_1D),allocatable::KI_3D(:)
  type(Ragged_Array_1D),allocatable::KII_3D(:)
  type(Ragged_Array_1D),allocatable::KIII_3D(:)
  type(Ragged_Array_1D),allocatable::KI_eq_3D(:)

  integer,allocatable::Crack3D_Meshed_Node_num(:)
  integer,allocatable::Crack3D_Meshed_Ele_num(:)
  integer num_Suspended_Point
  real(kind=FT),allocatable::Suspended_Points(:,:)
  integer,allocatable::Crack3D_Meshed_Outline_num(:)
  ! ------Circular Initial 3D Crack------
  real(kind=FT),allocatable::Crack3D_Cir_Coor(:,:)
  real(kind=FT),allocatable::Na_Crack3D_Cir_Coor(:,:)
  ! ------Initial 3D Elliptical Crack------
  real(kind=FT),allocatable::Crack3D_Ellip_Coor(:,:)
  type(Ragged_Array_2D),allocatable::Crack3D_Meshed_Node(:)
  type(Ragged_Int_Array_2D),allocatable::Crack3D_Meshed_Ele(:)
  type(Ragged_Array_2D),allocatable::Crack3D_Meshed_Node_Value(:)
  type(Ragged_Int_Array_1D),allocatable::Cr3D_Meshed_Node_in_Ele_Num(:)
  type(Ragged_Array_2D),allocatable::Cr3D_Meshed_Node_in_Ele_Local(:)
  type(Ragged_Array_2D),allocatable::Crack3D_Meshed_Ele_Attri(:)
  type(Ragged_Array_2D),allocatable::Crack3D_Meshed_Ele_Nor_Vector(:)
  type(Ragged_Array_2D),allocatable::Crack3D_Meshed_Node_Nor_Vector(:)
  type(Ragged_Int_Array_2D),allocatable::Crack3D_Meshed_Outline(:)
  type(Ragged_Array_2D),allocatable::Crack3D_Meshed_Vertex_x_Vector(:)
  type(Ragged_Array_2D),allocatable::Crack3D_Meshed_Vertex_y_Vector(:)
  type(Ragged_Array_2D),allocatable::Crack3D_Meshed_Vertex_z_Vector(:)
  type(Ragged_Array_3D),allocatable::Crack3D_Meshed_Vertex_T_Matrx(:)
  type(Ragged_Array_2D),allocatable::Crack3D_Meshed_Vertex_T_Theta(:)
  type(Ragged_Array_2D),allocatable::Crack3D_Vector_S1(:)

  ! ------Fluid Elements and Fluid Element Calculation Points------
  integer,allocatable::Cracks_FluidEle_num_3D(:)
  integer,allocatable::Cracks_CalP_Num_3D(:)
  integer,allocatable::Cracks_Real_CalP_Num_3D(:)
  real(kind=FT),allocatable::Cracks_Volume(:)
  real(kind=FT),allocatable::Cracks_Volume_Old(:)
  integer Key_3D_FluEle_Triang
  integer,allocatable::Cracks_FluidEle_CalP_Glo_Info(:,:)
  real(kind=FT),allocatable::Cracks_FluidEle_CalP_Glo_Insitu(:)
  real(kind=FT),allocatable::Crack3D_Centroid(:,:)
  type(Ragged_Int_Array_2D),allocatable::Cracks_FluidEle_CalP_3D(:)
  type(Ragged_Int_Array_2D),allocatable::Cracks_FluidEle_Glo_CalP_3D(:)
  type(Ragged_Int_Array_1D),allocatable::Cracks_FluidEle_num_CalP_3D(:)
  type(Ragged_Int_Array_1D),allocatable::Cracks_FluidEle_EleNum_3D(:)
  type(Ragged_Array_1D),allocatable::Cracks_FluidEle_Area_3D(:)
  type(Ragged_Array_2D),allocatable::Cracks_FluidEle_Centroid_3D(:)
  type(Ragged_Array_2D),allocatable::Cracks_FluidEle_LCS_x_3D(:)
  type(Ragged_Array_2D),allocatable::Cracks_FluidEle_LCS_y_3D(:)
  type(Ragged_Array_2D),allocatable::Cracks_FluidEle_LCS_z_3D(:)
  type(Ragged_Array_3D),allocatable::Cracks_FluidEle_LCS_T_3D(:)
  type(Ragged_Array_2D),allocatable::Cracks_FluidEle_Vector_3D(:)
  type(Ragged_Array_1D),allocatable::Cracks_FluidEle_Aper_3D(:)
  type(Ragged_Array_2D),allocatable::Cracks_CalP_Coors_3D(:)
  type(Ragged_Array_2D),allocatable::Cracks_CalP_Orient_3D(:)
  type(Ragged_Int_Array_1D),allocatable::Cracks_CalP_MeshedEl_3D(:)
  type(Ragged_Int_Array_2D),allocatable::Cracks_CalP_Elem_3D(:)
  type(Ragged_Array_1D),allocatable::Cracks_CalP_Aper_3D(:)
  type(Ragged_Array_2D),allocatable::Cracks_CalP_UpDis_3D(:)
  type(Ragged_Array_2D),allocatable::Cracks_CalP_LowDis_3D(:)
  type(Ragged_Array_1D),allocatable::Cracks_CalP_Pres_3D(:)
  type(Ragged_Array_2D),allocatable::Cracks_CalP_Tractions_3D(:)
  type(Ragged_Array_1D),allocatable::Cracks_CalP_Pgra_3D(:)
  type(Ragged_Array_1D),allocatable::Cracks_CalP_Velo_3D(:)
  type(Ragged_Array_1D),allocatable::Cracks_CalP_Quan_3D(:)
  type(Ragged_Array_1D),allocatable::Cracks_CalP_Conc_3D(:)
  type(Ragged_Array_1D),allocatable::Cracks_CalP_Remo_Strs_3D(:)
  ! ----Tip reinforcement unit and node related----
  integer Solid_El_Max_num_Crs
  parameter (Solid_El_Max_num_Crs   = 8)
  ! integer(kind=1), ALLOCATABLE :: Solid_El_num_Crs(:) ! Number of fractures related to each solid
  ! element, 2022-08-15. IMPROV2022091803.
  integer,ALLOCATABLE:: Solid_El_num_Crs(:)
  integer,ALLOCATABLE:: Solid_El_Crs(:,:)
  type(Ragged_Int_Array_1D),allocatable::Solid_El_Vertex_Num(:)
  type(Ragged_Array_3D),allocatable::Solid_El_Vertex_Coor(:)
  type(Ragged_Array_3D),allocatable::Solid_El_Vertex_Nor_Vec(:)
  type(Ragged_Array_3D),allocatable::Solid_El_Vertex_x_Vec(:)
  type(Ragged_Array_3D),allocatable::Solid_El_Vertex_y_Vec(:)
  type(Ragged_Array_3D),allocatable::Solid_El_Vertex_z_Vec(:)
  type(Ragged_Array_3D),allocatable::Solid_El_Pre_Vertex_Coor(:)
  type(Ragged_Array_3D),allocatable::Solid_El_Pre_Vertex_Nor_Vec(:)
  type(Ragged_Array_3D),allocatable::Solid_El_Pre_Vertex_x_Vec(:)
  type(Ragged_Array_3D),allocatable::Solid_El_Pre_Vertex_y_Vec(:)
  type(Ragged_Array_3D),allocatable::Solid_El_Pre_Vertex_z_Vec(:)
  type(Ragged_Array_3D),allocatable::Solid_El_Nex_Vertex_Coor(:)
  type(Ragged_Array_3D),allocatable::Solid_El_Nex_Vertex_Nor_Vec(:)
  type(Ragged_Array_3D),allocatable::Solid_El_Nex_Vertex_x_Vec(:)
  type(Ragged_Array_3D),allocatable::Solid_El_Nex_Vertex_y_Vec(:)
  type(Ragged_Array_3D),allocatable::Solid_El_Nex_Vertex_z_Vec(:)
  type(Ragged_Array_3D),allocatable::Solid_El_Tip_BaseLine(:)
  type(Ragged_Array_2D),allocatable::Solid_El_Tip_BaseLine_Nor_Vec(:)
  type(Ragged_Array_2D),allocatable::Solid_El_Tip_BaseLine_x_Vec(:)
  type(Ragged_Array_2D),allocatable::Solid_El_Tip_BaseLine_y_Vec(:)
  type(Ragged_Array_2D),allocatable::Solid_El_Tip_BaseLine_z_Vec(:)
  type(Ragged_Array_2D),allocatable::Solid_El_Tip_BaseLine_T_theta(:)
  type(Ragged_Array_3D),allocatable::Solid_El_Tip_BaseLine_T_Matrix(:)
  !--------------------
  integer Num_Check_T_Matrix
  integer,allocatable::Cracks_Initial_Meshed(:)
  integer,allocatable::Cracks_Initial_Adjusted(:)
  integer,allocatable::Cracks_Checked_and_Adjusted(:)
  !--------------------
  integer,allocatable::Crack_Type_Status_3D(:,:)
                                                        ! Column 1 (Fracture Type): =1, HF Fracture; =2, Natural Fracture; =3, Post-Fracturing Hydraulic
                                                        ! Fracture
                                                        ! Note: Natural fractures and propped hydraulic fractures may potentially turn into HF fractures.
                                                        ! Column 2 (Fracture Status): =1, HF fracturing not completed; =2, HF fracturing completed
                                                        ! Column 3 (Can the crack continue to propagate): =1, yes; =0, no
                                                        ! Column 4 (Whether the fracture has obtained a fluid node): =1, Yes; =0, No
                                                        ! Column 5 (Did the crack propagate in the previous step?): =1, Yes; =0, No
                                                        ! Column 6 (if activated from natural fractures secondarily, it corresponds to the number of the
                                                        ! respective natural fracture)
  ! 3D natural fracture related.
  integer Key_NaCr_Type_3D
  integer Num_Poly_Edges_NaCr
  integer Key_NaCr_Cross
  integer Key_NaCr_Growth
  real(kind=FT) NaCr_3D_n_Vector(3)
  real(kind=FT) NaCr_3D_n_Vector_Delta
  real(kind=FT) NaCr_3D_Size
  real(kind=FT) NaCr_3D_Sz_Delta
  real(kind=FT) NaCr_3D_Check_R
  real(kind=FT) NaCr_3D_Rect_Longside_Vector(3)
  real(kind=FT) NaCr_3D_Rect_L
  real(kind=FT) NaCr_3D_Rect_W
  real(kind=FT) NaCr_3D_Rect_L_Delta
  real(kind=FT) NaCr_3D_Rect_W_Delta
  real(kind=FT) NaCr_3D_Rect_Longside_Vector_Delta
  ! 3D Crack Intersection State.
  integer,ALLOCATABLE::Cracks_3D_Inter_Status(:,:)

  ! 3D enhancement element Related. 2022-09-05. IMPROV2022090502.
  integer,ALLOCATABLE:: Elem_Type_3D(:,:)
  integer,ALLOCATABLE:: c_POS_3D(:,:)
  integer,ALLOCATABLE:: Enriched_Node_Type_3D(:,:)
  type(Ragged_Array_2D),allocatable::Enriched_Node_Crack_n_Vector_3D(:)
  ! Related to the 3D Junction enhancement element.
  type(Ragged_Int_Array_1D),allocatable::Node_Jun_elem_3D(:)
  type(Ragged_Int_Array_1D),allocatable::Jun_Ele_Negative_Cr_Num_3D(:)
  type(Ragged_Array_2D),allocatable::Coors_Junction_3D(:)
  type(Ragged_Int_Array_1D),allocatable::Ele_Num_Tip_Enriched_Node_3D(:)

  integer n_h_Node_3D,n_t_Node_3D,n_j_Node_3D,n_hl_Node_3D,n_c_Node_3D
  ! Fluid-solid coupling matrix Q: using staggered array. 2022-09-16. IMPROV2022091601.
  type(Ragged_Array_1D),allocatable::Coupled_Q_3D(:)
  type(Ragged_Int_Array_1D),allocatable::Coupled_Q_3D_Index(:)
  ! EBE - Related to element Stiffness Matrix. 2022-09-19.
  type(Ragged_Array_2D),allocatable::storK_XFEM(:)
  type(Ragged_Array_2D),allocatable::storK_XFEM_Old(:)
  real(kind=FT),ALLOCATABLE::storK_FEM(:,:,:)
  real(kind=FT),ALLOCATABLE::storK_FEM_Sym(:,:)
  type(Ragged_Array_2D),allocatable::storK_XFEM_Updated(:)
  ! Others
  integer,ALLOCATABLE:: Elem_num_Related_Cracks(:)
  integer,ALLOCATABLE:: Elem_Related_Cracks(:,:)
  ! Permeability of solid element. 2022-11-28. NEWFTU2022112801.
  real(kind=FT),ALLOCATABLE::Ele_Permeability_3D(:,:)
  ! Daily crack boundary line area. 2023-02-22.
  real(kind=FT),allocatable::Cracks_Outline_Area(:)
  !
  ! Related to Sinopec.
  !
  integer num_XA_Input_Cracks
  real(kind=FT) XA_Min_Frac_Radius
  real(kind=FT),ALLOCATABLE:: XA_Ini_Cracks(:,:,:)
  real(kind=FT),ALLOCATABLE:: XA_Ini_Crack_St(:)
  real(kind=FT),ALLOCATABLE:: XA_Ini_Crack_FracFriction(:)
  ! solid element fracture volume ratio. 2023-03-26. NEWFTU2023032601.
  real(kind=FT),ALLOCATABLE::Ele_VolumeRatio_3D(:)
  !2023-08-10. 3D DIM SIF related.
  real(kind=FT) SIFs_DIM_3D_Offset_Delta_Factor,SIFs_DIM_3D_r_1_Factor,SIFs_DIM_3D_r_2_Factor,SIFs_DIM_3D_r_k_Factor
  !2024-05-03.
  integer,allocatable::simplified_crack_outline_num(:)
  real(kind=FT),allocatable::simplified_crack_outline_points(:,:,:)
end module Global_Crack_3D

!--------------------------
! 5.2 Cross Cracks Related
!--------------------------
module Global_Cross
  use Global_Float_Type
  implicit none
  save
  integer Max_Num_Cross
  integer n_Cross_Node
  parameter (Max_Num_Cross  = 100   )
  integer num_Cross
  integer,ALLOCATABLE:: Elem_Type_Cross(:,:)
  integer,ALLOCATABLE:: Enriched_Node_Type_Cross(:,:)
  integer,ALLOCATABLE:: c_POS_Cross(:,:)
  integer,ALLOCATABLE:: Node_Cross_elem(:,:)
  integer,ALLOCATABLE:: Cross_Point_Cr_num(:,:)
  integer,ALLOCATABLE:: Cross_Point_Ele_num(:)
  real(kind=FT),ALLOCATABLE:: Cross_Point_RABCD(:,:,:)
end module Global_Cross

!-----------------------
! 5.3 Inclusion-Related
!-----------------------
module Global_Inclusion
  use Global_Float_Type
  implicit none
  save
  integer Max_Num_Circ_Incl,Max_Num_Tri_Incl
  integer Max_Num_Quad_Incl,Max_Num_Penta_Incl
  integer Max_Num_Ellip_Incl,Max_Num_Incl,Max_Num_Poly_Incl
  integer Max_Num_Edges_Poly
  integer n_Incl_Node
  parameter (Max_Num_Incl  = 500   )
  parameter (Max_Num_Circ_Incl  = 100   )
  parameter (Max_Num_Poly_Incl = 100   )
  parameter (Max_Num_Ellip_Incl = 100   )
  parameter (Max_Num_Edges_Poly = 100   )
  integer num_Inclusion
  integer num_Circ_Incl
  integer num_Poly_Incl
  integer num_Ellip_Incl
  ! Definition rules for the position, size, and material number of circular inclusions
  real(kind=FT) Circ_Inclu_Coor(Max_Num_Circ_Incl,3)
  integer Circ_Inclu_Mat_Num(Max_Num_Circ_Incl)
  ! Definition rules for polygon inclusion position, size, and material number
  real(kind=FT) Poly_Incl_Coor_x(Max_Num_Poly_Incl,Max_Num_Edges_Poly)
  real(kind=FT) Poly_Incl_Coor_y(Max_Num_Poly_Incl,Max_Num_Edges_Poly)
  integer Poly_Inclu_Edges_Num(Max_Num_Poly_Incl)
  integer Poly_Inclu_Mat_Num(Max_Num_Poly_Incl)
  ! Closed-form polygon inclusion
  real(kind=FT) Poly_Incl_Coor_x_Cl(Max_Num_Poly_Incl,Max_Num_Edges_Poly)
  real(kind=FT) Poly_Incl_Coor_y_Cl(Max_Num_Poly_Incl,Max_Num_Edges_Poly)
  ! Related to enhancement elements
  integer,ALLOCATABLE:: Elem_Type_Incl(:,:)
  integer,ALLOCATABLE:: Enriched_Node_Type_Incl(:,:)
  integer,ALLOCATABLE:: c_POS_Incl(:,:)
  ! Randomly generated rules mixed with related elements
  integer Key_Rand_Circ_Incl
  integer num_Rand_Circ_Incl
  real(kind=FT) Rand_Circ_Incl_R
  real(kind=FT) Rand_Circ_Inc_R_Delta
  integer Key_Rand_Poly_Incl
  integer num_Rand_Poly_Incl
  integer num_Vert_Poly_Incl
  real(kind=FT)Rand_Poly_Incl_R
  real(kind=FT)Rand_Poly_Inc_R_Delta
  ! Randomly generated irregular mix related
  integer Key_Rand_Poly_Incl_Irregular
  integer num_Rand_Poly_Incl_for_Each_Type(10)
  real(kind=FT) Rand_Poly_Incl_R_Min_and_Max(10,2)
  real(kind=FT) Rand_Poly_Incl_Irregular_Extension_Factor
  real(kind=FT) Rand_Poly_Incl_Irregular_Inclination
  real(kind=FT) Rand_Poly_Incl_Irregular_R_Delta_Factor
  real(kind=FT) Rand_Poly_Incl_Irregular_Angle_Factor
end module Global_Inclusion

!--------------------------------------------------------------------------------
! 6. element side length, area, and volume (commonly used, so listed separately)
!--------------------------------------------------------------------------------
module Global_Elem_Area_Vol
  use Global_Float_Type
  implicit none
  save
  real(kind=FT) :: Max_Elem_Area,Min_Elem_Area,Ave_Elem_Area,Ave_Elem_L,Ave_Elem_Vol,Max_Elem_Vol,Min_Elem_Vol
  real(kind=FT) :: Max_Elem_Area_Enrich,Min_Elem_Area_Enrich
  real(kind=FT) :: Ave_Elem_Area_Enrich,Ave_Elem_L_Enrich
  real(kind=FT) :: Ave_Elem_Vol_Enrich
  real(kind=FT) :: Min_Ele_Edge_Length,Max_Ele_Edge_Length
  real(kind=FT) :: Ave_Elem_L_Enrich_Unlocalrefined
end module Global_Elem_Area_Vol

!------------------------
! 7. Material Parameters
!------------------------
module Global_Material
  use Global_Float_Type
  implicit none
  save
  real(kind=FT),ALLOCATABLE:: D(:,:,:)
  real(kind=FT),ALLOCATABLE:: D4(:,:,:)
  real(kind=FT),ALLOCATABLE:: D_Comp(:,:,:)
  real(kind=FT),ALLOCATABLE:: S(:,:,:)
  real(kind=FT),ALLOCATABLE:: St(:,:)
  real(kind=FT),ALLOCATABLE:: Sc(:,:)
  real(kind=FT),ALLOCATABLE:: T_Alpha(:)
  real(kind=FT),ALLOCATABLE:: KIc(:,:)
  real(kind=FT),ALLOCATABLE:: E(:,:)
  real(kind=FT),ALLOCATABLE:: v(:,:)
  real(kind=FT),ALLOCATABLE:: density(:)
  real(kind=FT),ALLOCATABLE:: thick(:)
  real(kind=FT) Global_K_m
  real(kind=FT),ALLOCATABLE:: Fd_k_MT(:,:,:)
  real(kind=FT),ALLOCATABLE:: Fd_c(:)
  real(kind=FT),ALLOCATABLE:: Lame_lambda(:)
  real(kind=FT),ALLOCATABLE:: Lame_mu(:)
  real(kind=FT),ALLOCATABLE:: Ele_ComMat_RotMatrix(:,:,:)
  real(kind=FT),ALLOCATABLE:: MC_dt(:)
  real(kind=FT),ALLOCATABLE:: MC_phi_deg(:)
  real(kind=FT),ALLOCATABLE:: MC_phi_rad(:)
  real(kind=FT),ALLOCATABLE:: MC_psi_deg(:)
  real(kind=FT),ALLOCATABLE:: MC_psi_rad(:)
  real(kind=FT),ALLOCATABLE:: MC_c(:)
  real(kind=FT),ALLOCATABLE:: D_for_cylindrical(:,:,:)
  real(kind=FT) St_KIc_Conversion
  !*****************************************************
  ! Material parameters Weibull processing. 2024-06-24.
  !*****************************************************
  logical Flag_Weibull_E 
  integer,ALLOCATABLE::Key_Weibull_E(:)
  real(kind=FT),ALLOCATABLE::Weibull_Parameters_E(:,:)
  real(kind=FT),ALLOCATABLE::Weibull_Elements_E(:)
  real(kind=FT),ALLOCATABLE::Weibull_Elements_D_Matrix(:,:,:)
end module Global_Material

!-----------------
! 8. Displacement
!-----------------
module Global_DISP
  use Global_Float_Type
  implicit none
  save
  real(kind=FT),ALLOCATABLE:: DISP(:),DISP_InSitu(:)
  real(kind=FT),ALLOCATABLE:: DISP_Cylinder(:)
  ! 2D Gauss Point Displacement
  real(kind=FT),ALLOCATABLE:: DISP_x_Gauss(:)
  real(kind=FT),ALLOCATABLE:: DISP_y_Gauss(:)
  ! 3D Gauss Point Displacement
  real(kind=FT),ALLOCATABLE:: DISP_z_Gauss(:)
end module Global_DISP

!-------------------------
! 9. Hydraulic Fracturing
!-------------------------
module Global_HF
  use Global_Float_Type
  implicit none
  save
  integer Max_Num_Frac
  integer Key_HF_num_Contact
  integer Max_Num_Inject_Crack
  integer Max_MS_Num_Crack
  integer Key_AAP
  integer Key_Proppant
  integer Key_IniPre_PassOn
                                    ! 0: Do not inherit, set total_time to zero, iterate the crack opening from 0
                                    ! 1: Inheritance, for detailed ideas see my notes, V3_P78
  real(kind=FT) SOR_factor
  integer Key_HF_Conv_Crite
                                    ! 2: Through crack width and water pressure
  integer Key_Cal_deltaTime
                                    ! The calculated delta_V is too small, so delta_Time is also too small.
                                    ! Picard iteration must use the rectangle method; otherwise, the calculated water pressure decreases
                                    ! linearly to zero.
                                    ! 2: Trapezoidal method, delta_L * 0.5 * (w1 w2)
                                    ! The calculated delta_V is too small, so delta_Time is too large.
  parameter (Max_Num_Frac  = 200 )
  parameter (Max_Num_Inject_Crack = 10)
  parameter (Max_MS_Num_Crack = 20)
  integer Num_Frac
  integer ifra_Data_iter_num(Max_Num_Frac)
                                    ! Only save the current fracture step general data at this step, such as enhanced node information.
                                    ! Crack calculation point information, etc.
  ! integer Max_Picard_Iter, Max_NR_Iter, Max_Secant_Iter ! Maximum number of Picard iterations,
  ! maximum number of NR iterations, maximum number of Secant method iterations
  ! integer Max_Picard_Iter, Max_NR_Iter! Maximum number of Picard iterations, maximum number of NR
  ! iterations
  integer Max_Picard_Iter
  real(kind=FT) Viscosity
  real(kind=FT) Viscosity_Par_m
  integer Key_Visco_Type
  real(kind=FT) Visco_Zoom_Factor
  integer Key_Leakoff
  real(kind=FT) Coeff_Leak
  real(kind=FT) Dia_of_Proppant
  real(kind=FT) Alpha_Picard
  real(kind=FT) NR_Tol
  real(kind=FT) Secant_Tol
  integer Max_PNR_Iter
  integer Max_NR_AAP_Iter
  integer Max_NR_MS_Iter
  integer Max_NR_Red_AAP_Iter
  integer Num_Pic_PNR_Iter
  integer Max_MNR_Iter
  integer Max_MNR_Red_Iter
  integer Max_Num_Lnsrch
  real(kind=FT) MNR_Tol
  real(kind=FT) PNR_Tol
  real(kind=FT) NR_AAP_Tol
  real(kind=FT) NR_Red_AAP_Tol
  real(kind=FT) NR_MS_Tol
  real(kind=FT) MNR_Red_Tol
  integer Key_Symm_HF
  integer Type_of_HF
  real(kind=FT) Viscosity_td_m,Viscosity_td_n
  ! -------Water Injection Related-------
  integer Inject_Crack_Num
  real(kind=FT) Inject_Q_Time(200)
  real(kind=FT) Inject_Q_Val(200)
  integer Key_Propp_Trans
  real(kind=FT) Max_c
  real(kind=FT) Inject_c_Time(200)
  real(kind=FT) Inject_c_Val(200)
  real(kind=FT) Inject_P_Time(200)
  real(kind=FT) Inject_P_Val(200)
  logical Propp_Trans_Start
  integer Counter_Num_iFrac(Max_Num_Frac)
  real(kind=FT),ALLOCATABLE:: F_ConP(:)
  ! --------Stage Fracturing Related--------
  integer i_MS,MS_Crack_Num
  integer MS_NaturalCr_Num
  integer MS_Crack_Order(Max_MS_Num_Crack)
  ! integer MS_Each_Cr_Poi_Num(Max_MS_Num_Crack)          !Number of initial left points for staged fracturing
  ! real(kind=FT) MS_Crack_Coor(Max_MS_Num_Crack,200,2)   ! Coordinates of cracks in each segment
  real(kind=FT) MS_InP_Loc(Max_MS_Num_Crack,1:2)
  integer MS_Finish_Counter_Iter(Max_MS_Num_Crack)
  real(kind=FT) HF_Theor_Time
  real(kind=FT) HF_Theor_Aper
  integer Key_Paper1_Alter
  integer Key_Paper1_finish
  real(kind=FT) Last_Inj_Pres
  ! Slippery Water Analysis
  real(kind=FT)  Current_SlipWater_P
  real(kind=FT)  Picard_Alpha
  real(kind=FT)  Picard_Tol
  ! Well shaft related, 2022-04-19
  integer Max_WB,Max_Stages,Max_Clust
  parameter (Max_WB       = 20)
  parameter (Max_Stages   = 20)
  parameter (Max_Clust = 20)
  integer num_Wellbore
  integer num_Points_WB(Max_WB)
  real(kind=FT) Wellbore_Coors(Max_WB,20,3)
  integer num_Stages_Wellbores(Max_WB)
  integer num_Crs_Stages_Wellbores(Max_WB,Max_Stages)
  real(kind=FT) Injection_Q_Stages_Wellbores(Max_WB,Max_Stages)
  real(kind=FT) Injection_T_Stages_Wellbores(Max_WB,Max_Stages)
  integer Key_Gen_Ini_Crack_Wellbores
                                                             ! =1, Generate an initial rectangular crack
                                                             ! =2, Generate an initial circular crack
                                                             ! =3, Generate initial cracks in the polygon
  integer Key_Gen_Ini_Crack_Rule_Wellbores
                                                             ! =1, the generated initial fractures are perpendicular to the wellbore (default)
                                                             ! =2, generate the initial cracks perpendicular to the direction of minimum principal stress.
                                                             ! =3, user-defined.
  real(kind=FT) Normal_Vector_Gen_Ini_Crack_Wellbores(3)
  integer Num_Poly_Edges_PolyCr_WB
  real(kind=FT) Size_Ini_Crack_Wellbores
  real(kind=FT) Wellbores_Start_Point(Max_WB,3)
  real(kind=FT) Wellbores_End_Point(Max_WB,3)
  integer Cracks_Stages_Wellbores(Max_WB,Max_Stages,Max_Clust)
  
  integer Key_3D_HF_SlipWater_fk_Type
                                                               ! =1, maximum value reaches KIc (default); =2, average value reaches KIc; =3, minimum value reaches
                                                               ! KIc
                                                               ! Note: Only works in PhiPsi3D_Static_HF_SlipWater.
  integer SlipWater_Max_Time_Steps_3D,SlipWater_Max_Pres_Steps_3D
  integer SlipWater_Time_Step_Conv_Check
  integer SlipWater_Pres_Step_Conv_Check
end module Global_HF

!---------------------------------------------------
! 10. Fracture Surface Contact and Proppant Related
!---------------------------------------------------
module Global_Contact
  use Global_Float_Type
  implicit none
  save
  integer,ALLOCATABLE:: Elem_Conta_Sta(:,:)
                                              ! 1: Status, =0, Not Contacted; =1, Contacted
  integer,ALLOCATABLE:: Elem_Conta_Sta_Last(:,:)
                                              ! 1: Status, =0, Not Contacted; =1, Contacted
  real(kind=FT),ALLOCATABLE:: Elem_Proppant_Coor(:,:)
end module Global_Contact

!--------------------
! 11. Stress Related
!--------------------
module Global_Stress
  use Global_Float_Type
  implicit none
  save
  ! 2D Node Stress
  real(kind=FT),ALLOCATABLE:: Stress_xx_Node(:)
  real(kind=FT),ALLOCATABLE:: Stress_yy_Node(:)
  real(kind=FT),ALLOCATABLE:: Stress_xy_Node(:)
  real(kind=FT),ALLOCATABLE:: Stress_vm_Node(:)
  ! 2D nodal stress thermal stress
  real(kind=FT),ALLOCATABLE:: TStress_xx_Node(:)
  real(kind=FT),ALLOCATABLE:: TStress_yy_Node(:)
  real(kind=FT),ALLOCATABLE:: TStress_xy_Node(:)
  real(kind=FT),ALLOCATABLE:: TStress_vm_Node(:)
  ! Supplementary 3D Node Stress and Thermal Stress
  real(kind=FT),ALLOCATABLE:: TStress_zz_Node(:)
  real(kind=FT),ALLOCATABLE:: TStress_yz_Node(:)
  real(kind=FT),ALLOCATABLE:: TStress_xz_Node(:)
  ! 3D Node Stress Supplement
  real(kind=FT),ALLOCATABLE:: Stress_zz_Node(:)
  real(kind=FT),ALLOCATABLE:: Stress_yz_Node(:)
  real(kind=FT),ALLOCATABLE:: Stress_xz_Node(:)
  ! 2D Gauss Point Stress
  real(kind=FT),ALLOCATABLE:: Stress_xx_Gauss(:)
  real(kind=FT),ALLOCATABLE:: Stress_yy_Gauss(:)
  real(kind=FT),ALLOCATABLE:: Stress_xy_Gauss(:)
  real(kind=FT),ALLOCATABLE:: Stress_vm_Gauss(:)
  ! 2D Gauss Point Thermal Stress
  real(kind=FT),ALLOCATABLE:: TStress_xx_Gauss(:)
  real(kind=FT),ALLOCATABLE:: TStress_yy_Gauss(:)
  real(kind=FT),ALLOCATABLE:: TStress_xy_Gauss(:)
  real(kind=FT),ALLOCATABLE:: TStress_vm_Gauss(:)
  ! 3D Gauss Point Stress Supplement
  real(kind=FT),ALLOCATABLE:: Stress_zz_Gauss(:)
  real(kind=FT),ALLOCATABLE:: Stress_yz_Gauss(:)
  real(kind=FT),ALLOCATABLE:: Stress_xz_Gauss(:)
  ! Unit stress state, whether σ1-σ3 > Tol is satisfied
  integer,ALLOCATABLE:: Ele_State_Stress_1_3(:)
  ! 2D Gauss Point Stress_Gaussian Stress Field Under Ground Stress Without Pore Water Pressure
  real(kind=FT),ALLOCATABLE:: Stress_xx_Gas_InSitu(:)
  real(kind=FT),ALLOCATABLE:: Stress_yy_Gas_InSitu(:)
  real(kind=FT),ALLOCATABLE:: Stress_xy_Gas_InSitu(:)
  real(kind=FT),ALLOCATABLE:: Stress_vm_Gas_InSitu(:)

  ! 2D Gauss Point Stress_Gaussian Stress Field under True Water Pressure (not Net Water Pressure)
  real(kind=FT),ALLOCATABLE:: Stress_xx_Gas_InSitu2(:)
  real(kind=FT),ALLOCATABLE:: Stress_yy_Gas_InSitu2(:)
  real(kind=FT),ALLOCATABLE:: Stress_xy_Gas_InSitu2(:)
  real(kind=FT),ALLOCATABLE:: Stress_vm_Gas_InSitu2(:)
  ! Initial stress field at Gauss points
  real(kind=FT),ALLOCATABLE:: InSitu_Strs_Gaus_xx(:,:)
  real(kind=FT),ALLOCATABLE:: InSitu_Strs_Gaus_yy(:,:)
  real(kind=FT),ALLOCATABLE:: InSitu_Strs_Gaus_xy(:,:)
  ! 3D Supplement
  real(kind=FT),ALLOCATABLE:: InSitu_Strs_Gaus_zz(:,:)
  real(kind=FT),ALLOCATABLE:: InSitu_Strs_Gaus_yz(:,:)
  real(kind=FT),ALLOCATABLE:: InSitu_Strs_Gaus_xz(:,:)
  ! Initial stress field at the node
  real(kind=FT),ALLOCATABLE::Str_xx_InSitu(:),Str_yy_InSitu(:),Str_xy_InSitu(:),Str_vm_InSitu(:)
  ! 3D Supplement
  real(kind=FT),ALLOCATABLE::Str_zz_InSitu(:),Str_yz_InSitu(:),Str_xz_InSitu(:)
  ! Cylindrical coordinate system
  real(kind=FT),ALLOCATABLE:: Stress_Crr_Node(:)
  real(kind=FT),ALLOCATABLE:: Stress_Ctt_Node(:)
  real(kind=FT),ALLOCATABLE:: Stress_Czz_Node(:)
  real(kind=FT),ALLOCATABLE:: Stress_Crt_Node(:)
  real(kind=FT),ALLOCATABLE:: Stress_Ctz_Node(:)
  real(kind=FT),ALLOCATABLE:: Stress_Crz_Node(:)
  real(kind=FT),ALLOCATABLE:: Stress_Cvm_Node(:)
  ! The following variables are used to quickly define the initial stress field (without
  ! distinguishing regions).
  real(kind=FT) InSitu_S1_3D
  real(kind=FT) InSitu_S2_3D
  real(kind=FT) InSitu_S3_3D
  real(kind=FT) InSitu_S1_nv_3D(3)
  real(kind=FT) InSitu_S2_nv_3D(3)
  real(kind=FT) InSitu_S3_nv_3D(3)
  ! Non-uniform initial stress field (in x, y, z directions). 2022-07-06. NEWFTU2022070601.
  integer Key_Nonuniform_InSitu_X_with_Z
  real(kind=FT) InSitu_Sx_3D_Seg_Strs_X_with_Z(100)
  real(kind=FT) InSitu_Sx_3D_Seg_Loca_X_with_Z(100)
  integer Key_Nonuniform_InSitu_X_with_Y
  real(kind=FT) InSitu_Sx_3D_Seg_Strs_X_with_Y(100)
  real(kind=FT) InSitu_Sx_3D_Seg_Loca_X_with_Y(100)
  integer Key_Nonuniform_InSitu_Y_with_Z
  real(kind=FT) InSitu_Sy_3D_Seg_Strs_Y_with_Z(100)
  real(kind=FT) InSitu_Sy_3D_Seg_Loca_Y_with_Z(100)

  integer Key_Nonuniform_InSitu_Y_with_X
  real(kind=FT) InSitu_Sy_3D_Seg_Strs_Y_with_X(100)
  real(kind=FT) InSitu_Sy_3D_Seg_Loca_Y_with_X(100)
  integer Key_Nonuniform_InSitu_Z_with_X
  real(kind=FT) InSitu_Sz_3D_Seg_Strs_Z_with_X(100)
  real(kind=FT) InSitu_Sz_3D_Seg_Loca_Z_with_X(100)
  integer Key_Nonuniform_InSitu_Z_with_Y
  real(kind=FT) InSitu_Sz_3D_Seg_Strs_Z_with_Y(100)
  real(kind=FT) InSitu_Sz_3D_Seg_Loca_Z_with_Y(100)
  ! Initial strain field at Gauss points
  real(kind=FT),ALLOCATABLE:: InSitu_Strain_Gaus_xx(:,:)
  real(kind=FT),ALLOCATABLE:: InSitu_Strain_Gaus_yy(:,:)
  real(kind=FT),ALLOCATABLE:: InSitu_Strain_Gaus_xy(:,:)
  ! 3D Supplement
  real(kind=FT),ALLOCATABLE:: InSitu_Strain_Gaus_zz(:,:)
  real(kind=FT),ALLOCATABLE:: InSitu_Strain_Gaus_yz(:,:)
  real(kind=FT),ALLOCATABLE:: InSitu_Strain_Gaus_xz(:,:)
  ! Initial stress of the element. XA. 2023-03-24.
  real(kind=FT),ALLOCATABLE:: XA_Ele_InSitu_S1_Vector(:,:)
  real(kind=FT),ALLOCATABLE:: XA_Ele_InSitu_S2_Vector(:,:)
  real(kind=FT),ALLOCATABLE:: XA_Ele_InSitu_S3_Vector(:,:)
  real(kind=FT),ALLOCATABLE:: XA_Ele_InSitu_S1_S2_S3(:,:)
  ! Unit Stress. XA. 2023-03-26.
  real(kind=FT),ALLOCATABLE:: XA_Ele_Stress(:,:)
end module Global_Stress

!----------------------------------
! 11.2 Strain Related (2021-09-10)
!----------------------------------
module Global_Strain
  use Global_Float_Type
  implicit none
  save
  ! 2D Node Stress
  real(kind=FT),ALLOCATABLE:: Strain_xx_Node(:)
  real(kind=FT),ALLOCATABLE:: Strain_yy_Node(:)
  real(kind=FT),ALLOCATABLE:: Strain_xy_Node(:)
  real(kind=FT),ALLOCATABLE:: Strain_vm_Node(:)
  ! Cylindrical coordinate system
  real(kind=FT),ALLOCATABLE:: Strain_Crr_Node(:)
  real(kind=FT),ALLOCATABLE:: Strain_Ctt_Node(:)
  real(kind=FT),ALLOCATABLE:: Strain_Czz_Node(:)
  real(kind=FT),ALLOCATABLE:: Strain_Crt_Node(:)
  real(kind=FT),ALLOCATABLE:: Strain_Ctz_Node(:)
  real(kind=FT),ALLOCATABLE:: Strain_Crz_Node(:)
  real(kind=FT),ALLOCATABLE:: Strain_Cvm_Node(:)

  ! 2D nodal stress thermal stress
  real(kind=FT),ALLOCATABLE:: TStrain_xx_Node(:)
  real(kind=FT),ALLOCATABLE:: TStrain_yy_Node(:)
  real(kind=FT),ALLOCATABLE:: TStrain_xy_Node(:)
  real(kind=FT),ALLOCATABLE:: TStrain_vm_Node(:)
  ! 3D Node Stress Supplement
  real(kind=FT),ALLOCATABLE:: Strain_zz_Node(:)
  real(kind=FT),ALLOCATABLE:: Strain_yz_Node(:)
  real(kind=FT),ALLOCATABLE:: Strain_xz_Node(:)
  ! 2D Gauss Point Stress
  real(kind=FT),ALLOCATABLE:: Strain_xx_Gauss(:)
  real(kind=FT),ALLOCATABLE:: Strain_yy_Gauss(:)
  real(kind=FT),ALLOCATABLE:: Strain_xy_Gauss(:)
  real(kind=FT),ALLOCATABLE:: Strain_vm_Gauss(:)
  ! 2D Gauss Point Thermal Stress
  real(kind=FT),ALLOCATABLE:: TStrain_xx_Gauss(:)
  real(kind=FT),ALLOCATABLE:: TStrain_yy_Gauss(:)
  real(kind=FT),ALLOCATABLE:: TStrain_xy_Gauss(:)
  real(kind=FT),ALLOCATABLE:: TStrain_vm_Gauss(:)
  ! 3D Gauss Point Stress Supplement
  real(kind=FT),ALLOCATABLE:: Strain_zz_Gauss(:)
  real(kind=FT),ALLOCATABLE:: Strain_yz_Gauss(:)
  real(kind=FT),ALLOCATABLE:: Strain_xz_Gauss(:)
  ! Unit stress state, whether σ1-σ3 > Tol is satisfied
  integer,ALLOCATABLE:: Ele_State_Strain_1_3(:)
  ! 2D Gauss Point Stress_Gaussian Stress Field Under Ground Stress Without Pore Water Pressure
  real(kind=FT),ALLOCATABLE:: Strain_xx_Gas_InSitu(:)
  real(kind=FT),ALLOCATABLE:: Strain_yy_Gas_InSitu(:)
  real(kind=FT),ALLOCATABLE:: Strain_xy_Gas_InSitu(:)
  real(kind=FT),ALLOCATABLE:: Strain_vm_Gas_InSitu(:)

  ! 2D Gauss Point Stress_Gaussian Stress Field under True Water Pressure (not Net Water Pressure)
  real(kind=FT),ALLOCATABLE:: Strain_xx_Gas_InSitu2(:)
  real(kind=FT),ALLOCATABLE:: Strain_yy_Gas_InSitu2(:)
  real(kind=FT),ALLOCATABLE:: Strain_xy_Gas_InSitu2(:)
  real(kind=FT),ALLOCATABLE:: Strain_vm_Gas_InSitu2(:)
  ! Initial strain field at the Gauss point. BUGFIX2022070906. Global variable redefinition.
  ! 3D Supplement
  ! Initial strain field at the node
  real(kind=FT),ALLOCATABLE::Strain_xx_InSitu(:),Strain_yy_InSitu(:),Strain_xy_InSitu(:),Strain_vm_InSitu(:)
  ! 3D Supplement
  real(kind=FT),ALLOCATABLE::Strain_zz_InSitu(:),Strain_yz_InSitu(:),Strain_xz_InSitu(:)
end module Global_Strain

!-----------------------------
! 12. Post-processing Related
!-----------------------------
module Global_POST
  use Global_Float_Type
  implicit none
  save
  integer Key_Post_CS_G_Coor
  integer Key_Post_CS_G_Disp
  integer Key_Post_CS_G_Strs
  integer Key_Post_S_Dof_F
  integer Key_Post_Cracked_Ele
  integer Key_Post_S_TanDisp
  real(kind=FT) Tol_Stress_1_3
  integer Key_Node_Value
  integer Key_Save_vtk
  integer Key_Simple_Post
  integer Key_Save_Nothing
  integer Key_Post_Elements_Gauss_Num
  integer Key_Get_Permeability
end module Global_POST

!------------------------------------------------------------------------------
! 13. Random generation of multibody dynamics simulation variables for spheres
!------------------------------------------------------------------------------
module Global_Rigid_Balls
  use Global_Float_Type
  implicit none
  save
  real(kind=FT) W_Ball_Model
  real(kind=FT) H_Ball_Model
  real(kind=FT) Ave_R_Ball
  real(kind=FT) Delta_R_Ball
  real(kind=FT) Max_R_Ball
  real(kind=FT) Min_R_Ball
  integer num_Balls
end module Global_Rigid_Balls

!--------------------------------------------------------------------------------------
! 14. Issues related to the field (including Biot's consolidation pore water pressure)
!--------------------------------------------------------------------------------------
module Global_Field_Problem
  use Global_Float_Type
  implicit none
  save
  integer Num_Fd_Bou_fixed
  integer,ALLOCATABLE::Fd_Bou_fixed(:)
  integer Num_Fd_Bou_vl
  real(kind=FT),ALLOCATABLE::Fd_Bou_vl(:,:)
  integer Num_Fd_Bou_vl_nz
  real(kind=FT),ALLOCATABLE::Fd_Bou_vl_nz(:,:)
  integer Num_Fd_Ini_vl
  real(kind=FT),ALLOCATABLE::Fd_Ini_vl(:,:)

  integer Num_Fd_Bou_qn
  real(kind=FT),ALLOCATABLE::Fd_Bou_qn(:,:)
  real(kind=FT),ALLOCATABLE::Fd_Value(:)
  real(kind=FT),ALLOCATABLE::Fd_ele_k_MT(:,:,:)
  real(kind=FT),ALLOCATABLE::Fd_Flux_x(:)
  real(kind=FT),ALLOCATABLE::Fd_Flux_y(:)
  integer Key_Fd_Body_Source
  real(kind=FT) Fd_Body_Source
  ! The following is related to Biot's consolidation pore water pressure
  real(kind=FT),ALLOCATABLE::Porous_P(:)
  real(kind=FT),ALLOCATABLE::Biot_c_MAT(:,:)
  ! Time step and time integration related
  integer Fd_IDy_Num_Iteras
  real(kind=FT) Delt_Time_Trapez
  ! The following is related to shale gas production assessment (Key_Analysis_Type=17)
  integer Key_Gas_Production
  integer GasP_num_Fractures
  real(kind=FT) GasP_Thic_Reservoir
  integer GasP_Well_Nodes(100)
  real(kind=FT) P_Langmuir
  real(kind=FT) V_Langmuir
  real(kind=FT) Density_gst
  real(kind=FT) Width_of_crack
  integer Key_Langmuir_Source
  real(kind=FT) porosity_Shale
  integer Key_Changing_BHP
  integer Num_BHP_curve_point
  real(kind=FT),ALLOCATABLE:: BHP_Curve(:,:)
  real(kind=FT) Gas_P_stress_x,Gas_P_stress_y
  integer Key_Changing_Kf
  integer Key_Proppant_Active
  integer Key_Proppant_Creep
  integer Key_Proppant_Crush
  real(kind=FT) Proppant_Strength
  real(kind=FT) Proppant_visco_factor
  real(kind=FT) Rock_visco_factor
  ! -----------Field Issues Related to XFEM-------------
  integer Fd_n_h_Node,Fd_n_t_Node,Fd_n_j_Node,Fd_n_hl_Node,Fd_n_c_Node,Fd_n_cross_Node,Fd_n_Incl_Node
  integer Fd_Usual_Freedom,Fd_Enrich_Freedom
  integer ,ALLOCATABLE:: Fd_c_POS(:,:)
  integer ,ALLOCATABLE:: Fd_c_POS_Hl(:,:)
  integer ,ALLOCATABLE:: Fd_c_POS_Cross(:,:)
  integer ,ALLOCATABLE:: Fd_c_POS_Incl(:,:)
  integer Key_Fd_TipEnrich
                            ! 1: Strong discontinuity, one item, sqrt(r) * sin(theta / 2)
                            ! 2: Weak intermittence, one term, sqrt(r)*cos(theta/2)
  integer,ALLOCATABLE:: Fd_EleGaus_yes_FEM_asemd(:,:)
  real(kind=FT),ALLOCATABLE:: Fd_Gauss_CoorX(:)
  real(kind=FT),ALLOCATABLE:: Fd_Gauss_CoorY(:)
  real(kind=FT),ALLOCATABLE:: Field_Value_Gauss(:)
  logical Fd_Yes_XFEM
end module Global_Field_Problem

!-------------------------------------------
! 15. Molecular Dynamics Simulation Related
!-------------------------------------------
module Global_MD
  use Global_Float_Type
  implicit none
  save
  integer MD_num_molecule
  real(kind=FT) MD_mss_molecule
  integer MD_num_time_step
  real(kind=FT) MD_Delt_Time
  real(kind=FT) MD_Dimension_x
  real(kind=FT) MD_Dimension_y
  real(kind=FT) MD_Dimension_z
  integer MD_step_print_num
  integer MD_step_save_num
end module Global_MD

!---------------------------------
! 16. Related to Plastic Analysis
!---------------------------------
module Global_Plasticity
  use Global_Float_Type
  implicit none
  save
  real(kind=FT),ALLOCATABLE:: dsdeEl(:,:,:,:)
  real(kind=FT),ALLOCATABLE:: dsdePl(:,:,:,:)
  real(kind=FT),ALLOCATABLE:: STATEV(:,:,:)
  real(kind=FT),ALLOCATABLE:: Last_U_of_Ele(:,:)
end module Global_Plasticity

!-----------------------------
! 17. Cohesive Cracks Related
!-----------------------------
module Global_Cohesive
  use Global_Float_Type
  implicit none
  save
  integer Coh_Constitutive_type
  real(kind=FT) Coh_Width_Critical1
  real(kind=FT) Coh_Width_Critical2
  real(kind=FT) Coh_f_Ultimate
  integer Coh_Tangential_Key
  real(kind=FT) Coh_Width_Critical1_T
  real(kind=FT) Coh_Width_Critical2_T
  real(kind=FT) Coh_f_Ultimate_T
  integer,ALLOCATABLE:: Elem_Coh_Sta(:,:)
                                            ! Status, =0, non-cohesive crack; =1, cohesive crack
end module Global_Cohesive


!------------------------------------
! 18. Near-field dynamics simulation
!------------------------------------
module Global_PD
  use Global_Float_Type
  implicit none
  save
  integer PD_num_points
  integer PD_step_print_num
  integer PD_step_save_num
  integer PD_num_time_step
  integer PD_ndivx
  integer PD_ndivy
  integer PD_nbnd
  integer PD_maxfam
  integer PD_Delt_Time
  real(kind=FT),ALLOCATABLE:: PD_coor(:,:)
  real(kind=FT),ALLOCATABLE:: PD_pforce(:,:)
  real(kind=FT),ALLOCATABLE:: PD_pforceold(:,:)
  real(kind=FT),ALLOCATABLE:: PD_bforce(:,:)
  real(kind=FT),ALLOCATABLE:: PD_stendens(:,:)
  real(kind=FT),ALLOCATABLE:: PD_fncst(:,:)
  real(kind=FT),ALLOCATABLE:: PD_disp(:,:)
  real(kind=FT),ALLOCATABLE:: PD_vel(:,:)
  real(kind=FT),ALLOCATABLE:: PD_velhalfold(:,:)
  real(kind=FT),ALLOCATABLE:: PD_velhalf(:,:)
  real(kind=FT),ALLOCATABLE:: PD_acc(:,:)
  real(kind=FT),ALLOCATABLE:: PD_massvec(:,:)
  real(kind=FT),ALLOCATABLE:: PD_enddisp(:)
  real(kind=FT),ALLOCATABLE:: PD_endtime(:)
  real(kind=FT),ALLOCATABLE:: PD_dmg(:)
  integer,ALLOCATABLE:: PD_numfam(:), PD_pointfam(:)
  integer,ALLOCATABLE:: PD_nodefam(:)
  integer,ALLOCATABLE:: PD_fail(:,:)
end module Global_PD


!-------------------------------------
! 19. Nonlinear (NL) Analysis Related
!-------------------------------------
module Global_NonLinear
      use Global_Float_Type
      implicit none
      save
      real(kind=FT) NL_TIMS(1000,5)
      INTEGER       NL_ITRA
      real(kind=FT) NL_ATOL
      INTEGER       NL_NTOL
      real(kind=FT) NL_TOL
      integer       NL_NLOAD
      real(kind=FT),ALLOCATABLE:: NL_Delta_U(:)
end module Global_NonLinear

!---------------------------------------
! 21. Surface load related. 2023-01-21.
!---------------------------------------
module Global_Surface_Load
  use Global_Float_Type
  use Global_Ragged_Array_Real_Classs
  use Global_Ragged_Array_Int_Classs
  implicit none
  save
  character(256) File_Surface_Load(100)
  integer Num_Surface_Loads
  type(Ragged_Int_Array_2D),allocatable::Surface_Load_Elements_Nodes(:)
  type(Ragged_Array_1D),allocatable::Surface_Load_Elements_Area(:)
  type(Ragged_Array_2D),allocatable::Surface_Load_Elements_Normal(:)
  real(kind=FT) Surface_Pressure(100)
end module Global_Surface_Load

!-----------------------------------------------------------------------------------
! 22. 3D Hydraulic Fracturing Experiment Simulation Related Parameters. 2023-01-23.
!-----------------------------------------------------------------------------------
module Global_3D_HF_Experiment
  use Global_Float_Type
  use Global_Ragged_Array_Real_Classs
  use Global_Ragged_Array_Int_Classs
  implicit none
  save
  integer HFE_Surface_Load_Num
  real(kind=FT) HFE_Initial_Injection_Rate
  integer HFE_Hole_Mat_Number
  real(kind=FT) HFE_Initial_Try_Pressure
  real(kind=FT) HFE_Initial_Pressure_Step_Size
end module Global_3D_HF_Experiment

!--------------------------------------
! 23. Related to Read_kpp. 2023-08-24.
!--------------------------------------
module Global_Read_kpp
  use Global_Float_Type
  implicit none
  save
  real(kind=FT) Read_kpp_Na_CRACK3D_COOR(100,10,3)
  real(kind=FT) Read_kpp_Na_Crack3D_Cir_Coor(100,7)
  integer Read_kpp_Each_NaCr3D_Poi_Num(100)
end module Global_Read_kpp

!--------------------------------------------------------------
! 24. Because the Silverfrost Compiler does not support OpenMP
! Therefore, define a virtual omp_lib.
!     2024-09-10.
!--------------------------------------------------------------
#ifdef Silverfrost
module omp_lib
  use Global_Float_Type
  implicit none
  save
  real(kind=FT) dummy_omp_lib
  
  contains
  
  integer function omp_get_thread_num()
      omp_get_thread_num = 0
  end function omp_get_thread_num
  
  integer function omp_get_max_threads()
      omp_get_max_threads = 1
  end function omp_get_max_threads
  
  integer function OMP_GET_NUM_PROCS()
      OMP_GET_NUM_PROCS = 1
  end function OMP_GET_NUM_PROCS
  
  integer function OMP_GET_NUM_THREADS()
      OMP_GET_NUM_THREADS = 1
  end function OMP_GET_NUM_THREADS
  
  
end module omp_lib
#endif

!-------------------------------------------------
! 25. ITPACK Solver Global Variables. 2024-09-13.
!-------------------------------------------------
module Global_ITPACK
  implicit none
  save
  integer ITPACK_NNZ
end module Global_ITPACK