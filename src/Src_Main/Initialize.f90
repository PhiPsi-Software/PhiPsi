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
 
SUBROUTINE Initialize
! Initialize the program and related variables.

!-----------------------------
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Filename
use Global_Crack
use Global_Crack_Common
use Global_Crack_3D
use Global_HF
use Global_Model
use Global_Field_Problem
use Global_Inclusion
use Global_Dynamic
use Global_MD
use Global_Cohesive
use Global_POST
use Global_NonLinear
use Global_Field_Problem
use Global_Elem_Area_Vol
use Global_Stress
use Global_Surface_Load
use Global_Material
use Global_3D_HF_Experiment
use Global_Cross

implicit none

!------------------------------------------
! Execute parameter initialization process
!------------------------------------------
print *, "    Initializing parameters of PhiPsi...."
! The maximum values of several global variables. 2023-08-13.
Max_Max_N_FluEl_3D = maxval(Max_N_FluEl_3D)  
Max_Max_N_CalP_3D  = maxval(Max_N_CalP_3D) 
Max_Max_N_Node_3D  = maxval(Max_N_Node_3D) 

! Crack Coordinate Initialization
Crack_Coor(1:Max_Num_Cr,1:Max_Num_Cr_P,1:2)=ZR
! Arc Crack Coordinate Initialization
Arc_Crack_Coor(1:Max_Num_Cr,1:Max_Num_Cr_P-1,1:11)  = ZR
! Initialize circular hole
Hole_Coor(1:Max_Num_Hl,1:3)=ZR
! Initialize elliptical holes
Ellip_Hole_Coor(1:Max_Num_Hl,1:5)=ZR
! Initialize circular inclusion
Circ_Inclu_Coor(1:Max_Num_Circ_Incl,1:3) =ZR
Circ_Inclu_Mat_Num(1:Max_Num_Circ_Incl)  =0
! Initialize polygon inclusions
Poly_Incl_Coor_x(1:Max_Num_Poly_Incl,1:Max_Num_Edges_Poly)=-1.0D8
Poly_Incl_Coor_y(1:Max_Num_Poly_Incl,1:Max_Num_Edges_Poly)=-1.0D8
Poly_Inclu_Mat_Num(1:Max_Num_Poly_Incl) = 0
Poly_Inclu_Edges_Num(1:Max_Num_Poly_Incl) =0
Poly_Incl_Coor_x_Cl(1:Max_Num_Poly_Incl,1:Max_Num_Edges_Poly)=-1.0D8
Poly_Incl_Coor_y_Cl(1:Max_Num_Poly_Incl,1:Max_Num_Edges_Poly)=-1.0D8
! Randomly generated irregular polygons with related inclusions
Key_Rand_Poly_Incl = 0
Key_Rand_Poly_Incl_Irregular =0
num_Rand_Poly_Incl_for_Each_Type(1:10)  =0
Rand_Poly_Incl_R_Min_and_Max(1:10,1:2)  =ZR
Rand_Poly_Incl_Irregular_Extension_Factor  = ONE
Rand_Poly_Incl_Irregular_Inclination  = ZR
Rand_Poly_Incl_Irregular_R_Delta_Factor = ZR
Rand_Poly_Incl_Irregular_Angle_Factor = ZR
! Randomly generated holes related
Key_Random_Hole     = 0
! Initialize number of cracks
num_Crack = 0
! Initialize the number of arc cracks
num_Arc_Crack = 0
! Initialize number of holes
num_Hole  = 0
num_Circ_Hole  = 0
num_Ellip_Hole  = 0
! Cracks Related to Hole Formation
Key_Hole_Crack_Generate  =0
Num_Crack_Hole_Generated    = 2
num_Hole_Crack_Generated(1:1000) = 0
Hole_Crack_Generated_num(1:1000,1:10) =0
! Initialize number of inclusions
num_Circ_Incl   = 0
num_Poly_Incl =0
num_Ellip_Incl  = 0
num_Inclusion = 0
! Initialize File Name
Filename  = '*blank*'
! Initialize file path
Work_Directory = '*blank*'
! Initialize the water content state of each crack
Cracks_HF_State(1:Max_Num_Cr)  = 0
! Initialize the state of each fracture containing proppant
Cracks_HF_Propp(1:Max_Num_Cr)  = 0
! Initialize water injection point coordinates (full HF model)
Inj_Point_Loc(1:2) =ZR
Yes_Arc_Crack      = .False.
! Post-processing related
Key_Post_CS_G_Coor = 1
Key_Post_CS_G_Disp = 1
Key_Post_CS_N_Strs = 1
Key_Post_CS_G_Strs = 1
Key_Post_S_Dof_F   = 0
Key_Post_Cracked_Ele=0
Key_Post_S_TanDisp = 0
Key_Node_Value     = 1
Key_Save_vtk       = 1
Key_Simple_Post    = 0
Key_Save_Nothing   = 0
Key_Post_Elements_Gauss_Num= 0
Key_Get_Permeability       = 0
Key_Post_CS_N_Stra  = 0
Key_Save_avol_file         = 0
Key_Post_CS_N_Stra_Cylindr= 0
Key_Post_CS_N_Strs_Cylindr= 0
Key_Post_CS_N_Disp_Cylindr =0
Post_cylinder_Coor_Center(1:3)    = -Con_Big_20
Post_cylinder_Coor_Vector_x(1:3)  = -Con_Big_20
Post_cylinder_Coor_Vector_y(1:3)  = -Con_Big_20
Post_cylinder_Coor_Vector_z(1:3)  = -Con_Big_20

! Initialize parameters related to water injection and proppant
Inject_Q_Time(1:200) = ZR
Inject_Q_Val(1:200)  = ZR
Inject_c_Time(1:200) = ZR
Inject_c_Val(1:200)  = ZR
Inject_P_Time(1:200) = ZR
Inject_P_Val(1:200)  = ZR
! Default Settings for Hydraulic Fracturing Analysis Control Parameters
MNR_Tol              = 0.05
Key_HF_Conv_Crite    = 2
                           ! 2: Through crack width and water pressure
Key_HF_LineSearch    = 1

! Maximum proppant concentration
Max_c=ZR

Key_Propped_Width = 0

! Has the proppant migrated?
Propp_Trans_Start = .False.
! The total number of iterations corresponding to the end of each rupture step
Counter_Num_iFrac(1:Max_Num_Frac) =0

! All cracks initialized as: all allowed to expand: 2022-08-21. BUGFIX2022082101.
allocate(Cracks_Allow_Propa(max(Max_Num_Cr,Max_Num_Cr_3D)))
Cracks_Allow_Propa(1:size(Cracks_Allow_Propa,1)) = 1

Visco_Zoom_Factor=20.0D0
! ----Time-Dependent Viscosity Parameters----
Viscosity_td_m    = 1.86D-2
Viscosity_td_n    = 2.5
! Initialize data related to staged fracturing
MS_Crack_Num  =0
MS_Crack_Order(1:Max_MS_Num_Crack)         = 0
! MS_Crack_Coor(1:Max_MS_Num_Crack,1:200,1:2) = ZR   ! Coordinates of cracks for each segment
MS_InP_Loc(1:Max_MS_Num_Crack,1:2)         = ZR
i_MS=0
Key_HF_MS_Contc_Check = -999
! Total iteration number corresponding to the completion of fracturing in each stage
MS_Finish_Counter_Iter(1:Max_MS_Num_Crack)= 0
! Related to the ground stress processing algorithm.
Key_InSitu_Strategy = 0
                          ! 0: Do not consider the initial stress problem
                          ! 1: Handle the displacement field caused by ground stress through the method of linear
                          ! superposition, that is, delete it during the HF iteration process
                          ! In-situ stress: When calculating the direction of crack propagation, in-situ stress is considered.
                          ! In this scheme, the water pressure at the crack tip is set to 0.
                          ! 2: Through Zienkiewicz_7ed_P216 (default)
                          ! 3: Apply the in-situ stress shielding (only for hydraulic fracturing analysis)
                          ! 4: Fixed constraints, specify initial stress or read in an initial stress file. Theory: Methods in
                          ! Xsite theory and validation examples.
Key_Read_Initial_Node_Stress_File = 0
Key_InStress_for_Mat   = 0
Mat_Number_of_InStress(1:100) =0
Mat_InStress_x         = ZR
Mat_InStress_y         = ZR
Mat_InStress_z         = ZR

Label_number = 0
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Initialization of cross-correlation variables for hydraulic fracturing HF and NF
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! The natural fracture number corresponding to the newly calculated fracture generated by the
! intersection of HF and NF
Cracks_NF_JiS_Cr_Num(1:Max_Num_Cr)  = 0
! Parasitic state at the intersection of HF and NF
Cracks_NF_JiS_Stat(1:Max_Num_Cr,1:2) = 0
! Related to the intersection of HF and NF, determine the T-shaped intersection state of newly
! generated fractures (unrelated to the crack tip)
Cracks_NF_T_Stat(1:Max_Num_Cr,1:2)   = 0
Cracks_NF_Cement_Tip_Num(1:Max_Num_Cr,1:2) =0
! Used for friction-type natural fractures, preserving the friction natural fracture number
! corresponding to each fracture
Cracks_fric_NF_num(1:Max_Num_Cr)  =0
! Used for friction-type natural fractures, to mark whether each fracture has been eroded by
! hydraulic fracturing.
Cracks_QinS_Stat(1:Max_Num_Cr)  =0
! In-situ stress
InSitu_x =ZR
InSitu_y =ZR
! element Gauss point Numbering Related
!num_GP_Elem(1:Num_Elem) = 0
!Ele_GP_Start_Num(1:Num_Elem) = 0
! Related to rupture zone
Key_Fracture_Zone = 0
Frac_Zone_MinX    = ZR
Frac_Zone_MaxX    = ZR
Frac_Zone_MinY    = ZR
Frac_Zone_MaxY    = ZR
Frac_Zone_MinZ    = ZR
Frac_Zone_MaxZ    = ZR
! Random generation of natural cracks
Key_Random_NaCr   = 0
num_Rand_Na_Crack = 0
NaCr_Orientation  = ZR
NaCr_Ori_Delta    = ZR
NaCr_Length       = ZR
NaCr_Len_Delta    = ZR
! Maximum sparsity ratio of sparse matrix K
Sparse_Ratio      = 0.03
! Data storage methods for sparse matrix intermediate variables (Location_COO, MASK_COO): 1:
! Entirely on disk (default); 2: Entirely in memory; 3: Partially on disk as binary files
Sparse_Store_Method  = 2
! Initialize the number of coordinate points for each crack to 0
Each_Cr_Poi_Num(1:Max_Num_Cr) = 0

! Crack surface contact
Key_Contact =0
Contact_Aper_Tol  = ZR
fric_mu_Cont      = ZR
! Maximum number of iterations for fracture surface contact (including proppant issues) >= 50
Max_Contact_Iter  = 50
! Number of contact integration points (1 or 2, default is 1)
Conta_Integ_Point =1
! Contact Analysis Convergence Criteria (1: Determined by residuals; 2: Determined by displacements;
! 3: Determined by crack openings)
Key_Conta_ConCrit = 2
! Related to penalty function contact algorithm
kn_Cont_Penalty   = 1.0D12
kt_Cont_Penalty   = 1.0D12
Conve_Tol_Penalty = 1.0D-3
Key_HF_Cont_Scheme= 0
Key_CS_Natural_Crack =0
Penalty_CS_Natural_Crack  = 1.0D11
Key_Penalty_CS_Method  = 1

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Field-related initialization (including initial pore pressure, etc.)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Fd_kxx            = -999.0D0      !Coefficient kxx
! Fd_kxy            = -999.0D0      !Coefficient kxy
! Fd_kyy            = -999.0D0      !Coefficient kyy
! Initial pore pressure
! Key_PoreP         = 1             !Whether there is initial pore water pressure
Key_PoreP         = 0
Initial_PoreP     = -999.0D0
Initial_Biot_Coeff= 0.75D0

! Material Parameter Initialization
ALLOCATE(Material_Type(Max_Materials))
Material_Type(1:Max_Materials) = 1
ALLOCATE(Material_Para(Max_Materials,30))
Material_Para(1:Max_Materials,1:30) = ZR
ALLOCATE(Material_Para_Added(Max_Materials,30))
Material_Para_Added(1:Max_Materials,1:30) = ZR
ALLOCATE(Mat_Cartesian_Coor_Vector_x(Max_Materials,3))
ALLOCATE(Mat_Cartesian_Coor_Vector_y(Max_Materials,3))
ALLOCATE(Mat_cylinder_Coor_Center(Max_Materials,3))
ALLOCATE(Mat_cylinder_Coor_Vector_z(Max_Materials,3))
Mat_Cartesian_Coor_Vector_x(1:Max_Materials,1:3) = ZR
Mat_Cartesian_Coor_Vector_y(1:Max_Materials,1:3) = ZR
Mat_cylinder_Coor_Center(1:Max_Materials,1:3) =ZR
Mat_cylinder_Coor_Vector_z(1:Max_Materials,1:3)  =ZR
GasP_Well_Nodes(1:100) = -999
! Initialization of Gauss points
Num_Gauss_Points  = 64
                          ! Supported number of Gauss points: 16, 36, 64, 100, 144 (12*12), 196 (14*14), 400 (20*20), 676
                          ! (26*26), 900 (30*30)
! Num_Gauss_Points = 900 !Number of Gauss integration points for enhanced elements, default is 64,
! cannot be an odd number, e.g., 9x9=81
                         
Num_Gauss_P_Inc   = 400
Num_Gau_Points_3D = 512
                          ! Supported number of Gauss points: 16, 36, 64, 100, 144 (12*12), 196 (14*14), 400 (20*20), 676
                          ! (26*26), 900 (30*30)
Num_Gau_Points_3D_MC= 1000
Num_Gau_P_SubInteg_6= 64
Num_Gau_P_SubInteg_4= 4
Num_Gauss_P_FEM   = 4
Num_Gauss_P_FEM_3D= 8
! Related to dynamic analysis
Num_Ivex =0
Num_Ivey =0
Num_Iacx =0
Num_Iacy =0
! Gravitational acceleration in all directions
! g_X_Y_Z(1:3) = ZR     ! Gravity acceleration in each direction
g_X_Y_Z           = [0.0D0,9.8D0,0.0D0]
! Earthquake-related
Key_EQ   = 0
EQ_Ac_Time_Gap  =-999.0D0
num_EQ_Ac_nodes = 0
EQ_Ac_nodes(1:5000) = 0
! Sine acceleration excitation related
Key_Sin_Accel       = 0
Sin_Accel_Dire      = 1
Sin_Accel_A         = -999.0D0
Sin_Accel_T         = -999.0D0
Sin_Accel_num_Nodes = 0
Sin_Accel_Nodes(1:5000)= 0
! Node Coupling Related
num_CP_x_nodes = 0
num_CP_y_nodes = 0
CP_x_nodes(1:5000)  = 0
CP_y_nodes(1:5000)  = 0
num_nodes_CP_set_x(1:10)=0
num_nodes_CP_set_y(1:10)=0
CP_nodes_x(1:10,1:5000)=0
CP_nodes_y(1:10,1:5000)=0
! 3D crack-related
allocate(Crack3D_Coor(Max_Num_Cr_3D,200,3))
Crack3D_Coor(1:Max_Num_Cr_3D,1:200,1:3) = ZR

allocate(Crack3D_Cir_Coor(Max_Num_Cr_3D,7))
Crack3D_Cir_Coor(1:Max_Num_Cr_3D,1:7) = ZR
allocate(Crack3D_Ellip_Coor(Max_Num_Cr_3D,8))
Crack3D_Ellip_Coor(1:Max_Num_Cr_3D,1:8) = ZR
allocate(Crack3D_Meshed_Node_num(Max_Num_Cr_3D))
Crack3D_Meshed_Node_num(1:Max_Num_Cr_3D)                       =0
allocate(Crack3D_Meshed_Ele_num(Max_Num_Cr_3D))
Crack3D_Meshed_Ele_num(1:Max_Num_Cr_3D)                        =0
Crack_Pressure(1:Max_Num_Cr) = ZR
Crack_Pressure_Type = 1
allocate(Cracks_CalP_Num_3D(Max_Num_Cr_3D))
Cracks_CalP_Num_3D(1:Max_Num_Cr_3D)                        = 0
! Number of contact iterations performed in each fracture step (by default, contact iterations are
! only performed in the first fluid-structure interaction iteration step, i.e., =1)
Key_HF_num_Contact= 1
! Related to dynamic analysis
IDy_Num_Iteras    = 0
IDy_Num_force_Itr = 0
Factor_Prop_Dy    = 1.63D0
Key_Mass_Lumped   = 0
! Molecular dynamics analysis related, analysis type number: 21
MD_num_molecule   = -999
MD_mss_molecule   = -999.0D0
MD_num_time_step  = -999
MD_Delt_Time      = -999.0D0
MD_Dimension_x    = -999.0D0
MD_Dimension_y    = -999.0D0
MD_Dimension_z    = -999.0D0
MD_step_print_num = 10
MD_step_save_num  = 10
! Related to plastic analysis
Key_Plasticity    = 0
! Cohesive cracks related
Coh_Constitutive_type   = -999
Coh_Width_Critical1     = -999.0D0
Coh_Width_Critical2     = -999.0D0
Coh_f_Ultimate          = -999.0D0
Coh_Tangential_Key      = -999
Coh_Width_Critical1_T   = -999.0D0
Coh_Width_Critical2_T   = -999.0D0
Coh_f_Ultimate_T        = -999.0D0
Key_Save_f_d_Curve      = -999
f_d_Curve_node_num      = -999
Key_Save_f_COD_Curve    = 0
f_COD_Curve_Crack_Num   = 1
Cracks_Coh_Ele_Type(1:Max_Num_Cr,1:Max_Num_Cr_CalP-1) =0
Cracks_Tip_Num_of_Coh_Ele(1:Max_Num_Cr,1:2)  =0
Max_Cohesive_Iter       = 50
Coh_Integ_Point         = 2
Key_Coh_ConCrit         = 1
Key_Save_f_d_Curve      = 0
! Others
Key_Close_Window        = 0
Key_Play_Sounds         = 0
Key_Memo_Moni           = 0
Key_Window_Log          = 0
Key_Clear_All           = 0

!Key_OpenMP              = 0          !Close openmp
Key_Num_Process         = 1
Key_Data_Format         = 1
Key_XA                  = 0
Key_Integral_Sol        = 2
Num_Sub_Quads           = 16
Num_Sub_3D_Cubes        = 125
                                   !(8=2^3,27=3^3,64=4^3,125=5^3,216=6^3,343=7^3,512=8^3,729=9^3,1000=10^3)
Num_Gau_Points_3D_Cube  = 8
Seed                    = 123456789
! Internal key control parameters (cannot be modified)
Key_Heaviside_Value    = -1
Key_Hole_Value         =  0
Key_Visit_PhiPsi_top   =  0
Key_Cond_Number        =  0
Key_Determinant        =  0
Key_Eigenvalue         =  0
Key_BLAS               =  0
                                   !EBE_XFEM_PCG_3D_with_K.f90
! Material parameters related
Material_Interface(1:2) = ZR
! Nonlinear Analysis Related (2018-01-19)
NL_TIMS(1,1:5) = [ZR,ONE,ZP1,ZR,ONE]
NL_TIMS(2:1000,1:5) = ZR
NL_ITRA        = 30
NL_ATOL        = 1.0D8
NL_NTOL        = 6
NL_TOL         = 1.0D-6

! Related to field issues
Key_Fd_Body_Source = 0
! Damaged Material Label
Key_Damage =0
Material_Dam2Frac_Value   = TWO
Crack_Gen_from_Damage =0

! Related to non-zero displacement boundary conditions
Num_Boux_nonzero = 0
Num_Bouy_nonzero = 0
Num_Bouz_nonzero = 0
penalty_k_bou_nonzero = 1.0D15

! Life and Death element
Key_EKILL  = 0
Ele_Killed_Each_Load_Step(1:Max_Step,1:1000)  =0
EKILL_Weaken_Factor = 1.0D-6

! Sparse Matrix Storage K
Key_K_Sparse =0

! Thermal stress
Key_Thermal_Stress =0
Key_Scheme_Thermal_Stress = 1
Key_Initial_Temperature   = 1
Thermal_Str_Temper(1:100) = ZR

! Define a matrix to store the node numbers of the 12 edges Ele_3D_Edges_Node(12,2), 2022-04-14
! Local number
Ele_3D_Edges_Node(1,1:2)  = [1,2]
Ele_3D_Edges_Node(2,1:2)  = [2,3]
Ele_3D_Edges_Node(3,1:2)  = [3,4]
Ele_3D_Edges_Node(4,1:2)  = [4,1]
Ele_3D_Edges_Node(5,1:2)  = [1,5]
Ele_3D_Edges_Node(6,1:2)  = [2,6]
Ele_3D_Edges_Node(7,1:2)  = [3,7]
Ele_3D_Edges_Node(8,1:2)  = [4,8]
Ele_3D_Edges_Node(9,1:2)  = [5,6]
Ele_3D_Edges_Node(10,1:2) = [6,7]
Ele_3D_Edges_Node(11,1:2) = [7,8]
Ele_3D_Edges_Node(12,1:2) = [8,5]

! Others
file_Sparse_K_Location_COO_bin = .False.
file_Sparse_K_Mask_COO_bin     = .False.
Key_Schollmann_Twist = 0
Key_Ave_Stress    = 1
S1_Point_Factor   = 0.2D0
num_Gauss_S1      = 216
Key_Adj_Prp_Step_3D = 0
Prp_3D_Factor_m     = 0.3333D0
Adj_Prp_Step_3D_Max = 1.5D0
Prp_Bisec_Factor_3D = 2.0D0
Key_Smooth_Front    = 0
Key_Smooth_Front_Twice =0
Key_3D_FluEle_Triang= 1
Key_Crack_Inner_Pressure = 0
Key_Multi_TipEnrNode=0
Key_Junction_Enrich         =  0
Key_TipEnrich_Radius_Factor =  2.0D0
Num_Check_T_Matrix   = 180
Key_Pre_Conditioner  =0
Lis_Number_Iteration = 5000
MDOF_2D              = 200
MDOF_3D              = 156
Key_Initiation       = 0
Key_Ini_Rule         = 1
Key_FS_Seque_Coupling= 0
Key_TS_Seque_Coupling= 0
Key_SIFs_DIM_Points  = 2
Key_SIFs_DIM_Method  = 1
Factor_Propagation   = 1.5
Propagation_Length   = -99.0D0
Key_Local_Mesh_Refine= 0
Key_Large_Deform     = 0
Key_Front_Segmentation = 0
Number_front_Segments  =2
Key_Adj_Ini_Crack_3D = 0
Key_Check_and_Adjust_Cracks_3D =  0
Adjust_Cracks_Delta_Angle   =  45.0D0
Adjust_Cracks_Resolution    =  6
num_Suspended_Point  = 0
key_tip_fluid_element_3d = 1
allocate(Suspended_Points(1000,3))
Suspended_Points(1:1000,1:3)=-TEN_15
Flag_Local_Refined   = .false.
Ave_Elem_L_Enrich_Unlocalrefined = ZR
Key_Ini_Cr_3D_Type   =1
Picard_Alpha         = 0.25D0
Picard_Tol           = 0.01
Max_Picard_Iter      = 50
Key_Block_Model      = 0
Num_Elem_Block_Bou   = 0
Key_InPlane_Growth   = 0
Key_Stop_Outside_Crack =0
! ---Wellbore Related (2022-04-19)----
Key_HF_Multistage_3D = 0
num_Wellbore         = 0
num_Points_WB(1:20)  = 2
Wellbore_Coors(1:20,1:20,1:3) = ZR
num_Stages_Wellbores(1:20)    = 1
num_Crs_Stages_Wellbores(1:20,1:20) = 1
Key_Gen_Ini_Crack_Wellbores   = 0
Key_Gen_Ini_Crack_Rule_Wellbores =  1
                                          ! =1, the generated initial fractures are perpendicular to the wellbore (default)
                                          ! =2, generate the initial cracks perpendicular to the direction of minimum principal stress.
                                          ! =3, user-defined.
Normal_Vector_Gen_Ini_Crack_Wellbores(1:3) = [ONE,ZR,ZR]
Size_Ini_Crack_Wellbores      = 5.0D0
Num_Poly_Edges_PolyCr_WB      = 5

Wellbores_Start_Point(1:20,1:3) = -TEN_15
Wellbores_End_Point(1:20,1:3)   = -TEN_15
Injection_Q_Stages_Wellbores(1:20,1:20)= -TEN_15
Injection_T_Stages_Wellbores(1:20,1:20)= -TEN_15
! ----3D Crack Leading Edge S1 or K Smooth Treatment - (2022-04-25)----
Key_Denoise_Vertex_Value =0
Key_Smooth_Vertex_Value  =0
Smooth_Vertex_n          =4
Key_Smooth_Vertex_Value2 =0
Smooth_Vertex_n2         =3
! -------3D crack front Theta smoothing treatment (available when CFCP=3)------- Added on
! 2022-07-14. (Later abandoned)
Key_Denoise_Theta_Value= 0
Key_Smooth_Theta_Value = 0
Smooth_Theta_n         = 2
Key_Smooth_Theta_Value2= 0
Smooth_Theta_n2        = 2
! -------3D crack front edge GrowthFactor smoothing treatment (available when CFCP=3)-------Added on
! 2022-07-14.
Key_Denoise_GF_Value= 0
Key_Smooth_GF_Value = 0
Smooth_GF_n         = 2
Key_Smooth_GF_Value2= 0
Smooth_GF_n2        = 2
! ----Related to Compressive-Shear Cracks----
Key_CS_Crack(1:Max_Num_Cr) = 0
allocate(Cracks_Initial_Meshed(Max_Num_Cr_3D))
Cracks_Initial_Meshed(1:Max_Num_Cr_3D) =0
allocate(Cracks_Initial_Adjusted(Max_Num_Cr_3D))
Cracks_Initial_Adjusted(1:Max_Num_Cr_3D)=0
allocate(Cracks_Checked_and_Adjusted(Max_Num_Cr_3D))
Cracks_Checked_and_Adjusted(1:Max_Num_Cr_3D)=0

Key_EBE_Precondition      = 1
Key_EBE_Condition_Number  = 0
Key_EBE_Sym_Storage_K     = 0
! Used to mark the type and condition of cracks. 2022-06-02.
allocate(Crack_Type_Status_3D(Max_Num_Cr_3D,10))
Crack_Type_Status_3D(1:Max_Num_Cr_3D,1)=1
Crack_Type_Status_3D(1:Max_Num_Cr_3D,2)=1
Crack_Type_Status_3D(1:Max_Num_Cr_3D,3)=1
Crack_Type_Status_3D(1:Max_Num_Cr_3D,4)=0
Crack_Type_Status_3D(1:Max_Num_Cr_3D,5)=0
Flag_HF_3D = 0
Key_Cpp_Call_Fortran_Lib = 0

! The following variables are used to quickly define the initial stress field (without regional
! distinction). 2022-06-03.
InSitu_S1_3D = ZR
InSitu_S2_3D = ZR
InSitu_S3_3D = ZR
InSitu_S1_nv_3D = [1.0D0,0.0D0,0.0D0]
InSitu_S2_nv_3D = [0.0D0,1.0D0,0.0D0]
InSitu_S3_nv_3D = [0.0D0,0.0D0,1.0D0]
! Non-uniform initial stress field (in x, y, z directions). 2022-07-06. NEWFTU2022070601.
Key_Nonuniform_InSitu_X_with_Z          = 0
InSitu_Sx_3D_Seg_Strs_X_with_Z(1:100)   = ZR
InSitu_Sx_3D_Seg_Loca_X_with_Z(1:100)   = ZR
Key_Nonuniform_InSitu_X_with_Y          = 0
InSitu_Sx_3D_Seg_Strs_X_with_Y(1:100)   = ZR
InSitu_Sx_3D_Seg_Loca_X_with_Y(1:100)   = ZR
Key_Nonuniform_InSitu_Y_with_Z          = 0
InSitu_Sy_3D_Seg_Strs_Y_with_Z(1:100)   = ZR
InSitu_Sy_3D_Seg_Loca_Y_with_Z(1:100)   = ZR
Key_Nonuniform_InSitu_Y_with_X          = 0
InSitu_Sy_3D_Seg_Strs_Y_with_X(1:100)   = ZR
InSitu_Sy_3D_Seg_Loca_Y_with_X(1:100)   = ZR
Key_Nonuniform_InSitu_Z_with_X          = 0
InSitu_Sz_3D_Seg_Strs_Z_with_X(1:100)   = ZR
InSitu_Sz_3D_Seg_Loca_Z_with_X(1:100)   = ZR
Key_Nonuniform_InSitu_Z_with_Y          = 0
InSitu_Sz_3D_Seg_Strs_Z_with_Y(1:100)   = ZR
InSitu_Sz_3D_Seg_Loca_Z_with_Y(1:100)   = ZR
allocate(Cracks_FluidEle_CalP_Glo_Insitu(Max_Num_Cr_3D*Max_Max_N_CalP_3D))
Cracks_FluidEle_CalP_Glo_Insitu(:)  = ZR
! 3D Natural Fracture Related
Key_NaCr_Type_3D       = 1
Num_Poly_Edges_NaCr    = 6
Key_NaCr_Cross         = 0
Key_NaCr_Growth        = 0
NaCr_3D_n_Vector(1:3)  = [ONE,ZR,ZR]
NaCr_3D_n_Vector_Delta = ZR
NaCr_3D_Size           = -999.0D0
NaCr_3D_Sz_Delta       = ZR
NaCr_3D_Rect_Longside_Vector = ZR
NaCr_3D_Rect_L         = -999.0D0
NaCr_3D_Rect_W         = -999.0D0
NaCr_3D_Rect_L_Delta   = ZR
NaCr_3D_Rect_W_Delta   = ZR
NaCr_3D_Rect_Longside_Vector_Delta = ZR
num_XA_Input_Cracks    = 0
XA_Min_Frac_Radius     = ZR
! -------3D Natural Fracture Activation Algorithm (for 3D)-----------2023-01-07
Key_NaCr_Active_Scheme_3D  =  1
Size_Factor_of_Active_NaCr = 6.0D0
KIc_NaCr                   = ZR
St_NaCr                    = ZR
Key_CFCP_3_Type   = 1
                                           ! 2: Can be expanded if greater than zero (Ref: Tang_2019_Analysis of stress interference
                                           ! among_Eq.16)
Key_3D_Cr_Update_Scheme = 2
                                           ! 2: Non-extended crack tips do not extend; crack tip coordinates are updated in very small
                                           ! increments (default).
Schollm_Max_Theta = 55.0D0
! Surface load related.
Num_Surface_Loads = 0
Surface_Pressure(1:100) = ZR

! Other parameters allocate memory. 2022-09-04.
allocate(Cracks_FluidEle_num_3D(Max_Num_Cr_3D))
Cracks_FluidEle_num_3D(1:Max_Num_Cr_3D) = 0
allocate(Cracks_Real_CalP_Num_3D(Max_Num_Cr_3D))
Cracks_Real_CalP_Num_3D(1:Max_Num_Cr_3D)= 0
allocate(Cracks_Volume(Max_Num_Cr_3D))
Cracks_Volume     = ZR
allocate(Cracks_Volume_Old(Max_Num_Cr_3D))
Cracks_Volume_Old = ZR
allocate(Cracks_FluidEle_CalP_Glo_Info(Max_Num_Cr_3D*Max_Max_N_CalP_3D,3))
Cracks_FluidEle_CalP_Glo_Info(1:Max_Num_Cr_3D*Max_Max_N_CalP_3D,1:3)  = 0
allocate(Crack3D_Centroid(Max_Num_Cr_3D,3))
Crack3D_Centroid(1:Max_Num_Cr_3D,1:3) = ZR
allocate(Crack3D_Meshed_Outline_num(Max_Num_Cr_3D))
Crack3D_Meshed_Outline_num(1:Max_Num_Cr_3D) = 0
allocate(KI_3D(Max_Num_Cr_3D))
allocate(KII_3D(Max_Num_Cr_3D))
allocate(KIII_3D(Max_Num_Cr_3D))
allocate(KI_eq_3D(Max_Num_Cr_3D))

!--------------------------------------------------------------------------------------------------
! 2022-10-09. IMPROV2022100901. The following are the variables originally in PhiPsi_Read_Input.f.
!--------------------------------------------------------------------------------------------------
Key_Junction_Check= 1
Key_Propa_Type    = 1
Max_MNR_Iter      = 30
Max_Num_Lnsrch    = 20
Prop_Angle_Alowed = 180.0D0
Key_Unit_System = 1
k_tt = ZR
k_nn = 1.0D16
a_Ave_Shi         = 5.0
Factor_Ave_R      = 0.7
Key_HF_Multistage = 0
Delta_Factor_SIFs = ZP1
Key_Tip_Pres_Zero = 1
Key_HF_Del_Ne_Pres= 0
Coff_a_AMPT       = ZP1
Water_Pressure    = 1.0D6
Key_Leakoff       = 0
Coeff_Leak        = 1.0D-5
Key_IniPre_PassOn = 0
                          ! 0: Does not inherit, total_time reset to zero, crack opening iterates from 0 (default)
                          ! 1: Inheritance, see my notes for details, V3_P78 (does not support PNR iteration method #4)
Key_Cal_deltaTime = 1
                          ! 1: Rectangle method (default), delta_V = delta_L * delta_w1;
                          ! The calculated delta_V is too small, so delta_Time is also too small;
                          ! Picard iteration must use the rectangle method; otherwise, the calculated water pressure decreases
                          ! linearly to zero.
                          ! 2: Trapezoidal method, delta_L * 0.5 * (w1 + w2);
                          ! The calculated delta_V is too small, so delta_Time is too large.
SOR_factor        = 0.75D0
Key_Visco_Type    = 1
                          ! 2: Viscosity changes with the concentration of the support (dynamic viscosity)
Viscosity_Par_m   = TWO
Max_c             = ZP6
Factor_Check_Dis  = 3.0
Factor_L_Jun_NewCr= 5.23
Key_Kink_Point    = 1
First_XFEM_Step   = 0
!------------------
Desired_KIc  = 2.0D6


!------------------
Key_Ele_Max_Related_Cracks = 10

St_KIc_Conversion = 6.88D0
!------------------
XA_Step_Count =  0

! --------Hydraulic Fracturing Experiment Simulation Related----------
HFE_Initial_Try_Pressure = ZR
HFE_Initial_Pressure_Step_Size = 0.1D6

!------------
Key_Allow_3D_Outside_Crack = 0

Key_3D_HF_Time_Step_Method = 1


Key_3D_HF_SlipWater_fk_Type =1

!IMPROV2023081001.
SIFs_DIM_3D_Offset_Delta_Factor = 0.5D0
!SIFs_DIM_3D_r_1_Factor          = 0.1D0
!SIFs_DIM_3D_r_2_Factor          = 0.2D0
!2023-08-22.
SIFs_DIM_3D_r_1_Factor          = 1.0D0
SIFs_DIM_3D_r_2_Factor          = 2.0D0
SIFs_DIM_3D_r_k_Factor          = 1.5D0

Key_Crack_Aperture_Method =2

!2023-08-22.
Key_Print_SIFs_to_Screen = 0

!2023-08-27.
Crack_Max_Min_Aperture(1:Max_Num_Cr,1:3) = -999.0D0

!2024-02-15.
Key_Save_Crack_Radius=0

! 3D circle equivalent to polygon resolution. 2024-02-26.
Circle_3D_Eqv_Polygon_Resolution = 21

! Maximum number of iterations for 3D SlipWater NR time steps and pressure steps. IMPROV2024022801.
SlipWater_Max_Time_Steps_3D   = 20  
SlipWater_Max_Pres_Steps_3D   = 10 
SlipWater_Time_Step_Conv_Check= 1
SlipWater_Pres_Step_Conv_Check= 1
!2024-03-09.
Key_Random          = 1
Inject_Crack_Num    = 1
Key_TipEnrich       = 1
Key_CFCP_3_Type     = 1
Key_FD_Tipenrich    = 1

!2024-03-13.
MAT_ALLOW_CRACK_Initiation(1:Max_Materials) = 1
Key_Max_Num_Initiation_Cracks               = 1
Num_Initiation_Cracks                   = 0

!2024-03-17.
Name_of_Keywords = ''

!2024-03-28.
Key_Dimension  = 2

!2024-04-11.
key_min_aperture  = 0
Key_Non_Negtive_Aperture_3D = 0

!2024-04-24.
Key_Print_EBEPCG_Solution_Time = 0

!2024-05-02. IMPROV202450201.
Key_Scheme_Signed_Dis_InPlane = 1

! 2024-06-23. Related to the initial crack generation area. 2024-06-23
Key_Ini_Crack_Zone = 0
Ini_Crack_Zone_MinX = ZR
Ini_Crack_Zone_MaxX = ZR          
Ini_Crack_Zone_MinY = ZR
Ini_Crack_Zone_MaxY = ZR   
Ini_Crack_Zone_MinZ = ZR
Ini_Crack_Zone_MaxZ = ZR    

! Elastic modulus Weibull processing related parameters. 2024-06-24. NEWFTU2024062402.
Flag_Weibull_E = .False.
allocate(Key_Weibull_E(Max_Materials))
Key_Weibull_E(1:Max_Materials) = 0
allocate(Weibull_Parameters_E(Max_Materials,3))
Weibull_Parameters_E(1:Max_Materials,1) = TWO
Weibull_Parameters_E(1:Max_Materials,2) = ZR
Weibull_Parameters_E(1:Max_Materials,3) = Con_Big_20

! Jagged crack pattern. 2024-06-25. For demonstration purposes only.
Key_Sawtooth_Crack_Path = 0

! User-defined 2D crack path related. 2024-11-17.
Key_User_Defined_2D_Crack_Path = 0
Num_User_Defined_2D_Crack_Path = 0

! Initialization of other variables.
Key_Static_Fatigue = 0
Key_Proppant_Creep = 0
Key_Rand_Circ_Incl = 0
Key_Element_Break  = 0
num_CP_set_x       = 0
num_CP_set_y       = 0
Num_Cross          = 0
Enrich_Freedom     = 0
key_hf_secant_ts   = 0
KI  = ZR
KII = ZR
Yes_XFEM = .false.
cg_tol = 1.0D-6


print *," "

RETURN
END SUBROUTINE Initialize
