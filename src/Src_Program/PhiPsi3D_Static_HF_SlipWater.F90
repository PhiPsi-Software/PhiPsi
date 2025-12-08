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
 
SUBROUTINE PhiPsi3D_Static_HF_SlipWater
! Static analysis of three-dimensional problems.
! Constant water pressure. Double iterative loop. Newton-Raphson method to determine fracturing and
! time steps.
! 2022-06-28.
!
! Updated on 2023-05-23: When Key_3D_HF_Time_Step_Method=1, the Newton-Raphson method is used to
! solve the time step.
! When Key_3D_HF_Time_Step_Method=2, the time step is solved using the bisection method.
!

! The global variable Crack_Type_Status_3D(i_C,10) is used to indicate the type and status of
! cracks.
! Column 1 (Fracture Type): =1, HF Fracture; =2, Natural Fracture; =3, Post-Fracturing Hydraulic
! Fracture
! Note: Natural fractures and after-pressurized hydraulic fractures may potentially turn into HF
! fractures.
! Column 2 (Fracture Status): =1, HF fracturing not completed; =2, HF fracturing completed
! Column 3 (Can the crack continue to propagate): =1, yes; =0, no
! Column 4 (Whether the fracture has obtained a fluid node): =1, Yes; =0, No
! Column 5 (Did the crack propagate in the previous step?): =1, Yes; =0, No
! The global variable Cracks_Stages_Wellbores(i_WB, i_Stage, i_C) is used for the crack numbers
! corresponding to each segment of each wellbore.
      
!----------------------------------
! Read the public variable module.
!----------------------------------
use Global_Float_Type
use Global_Common   
use Global_Filename
use Global_Model
use Global_XFEM_Elements
use Global_Elem_Area_Vol
use Global_Crack_Common
use Global_Crack_3D
use Global_DISP
use Global_HF
use Global_Stress
use Global_Strain
use Global_POST
      
! For testing
use Global_Ragged_Array_Real_Classs

!---------------------------------------------------------------------------------------------
! Read the subroutine interface module (activate the compiler's parameter consistency check).
!---------------------------------------------------------------------------------------------
use Global_D3_HF_Get_Pres_by_NR_Method
use Global_EBE_XFEM_PCG_3D_Get_K
      
!----------------------------
! Variable type declaration.
!----------------------------
implicit none
integer isub
logical Yes_Last_Growth
real(kind=FT) Lambda
integer,ALLOCATABLE::freeDOF(:),fixedDOF(:)
real(kind=FT),ALLOCATABLE:: globalK(:,:)
real(kind=FT),ALLOCATABLE:: ori_globalK(:,:)
real(kind=FT),ALLOCATABLE::tem_DISP(:)
integer num_FreeD
real(kind=FT) Max_Disp_x,Max_Disp_y,Min_Disp_x,Min_Disp_y,Max_Disp_z,Min_Disp_z
integer i_C
logical Yes_Growth(Max_Num_Cr_3D)
integer ifra
integer num_FixedD
real(kind=FT),ALLOCATABLE:: F_U(:),F(:)
real(kind=FT),ALLOCATABLE:: delta_x(:)
integer date_time(8)
integer(LIT) F_time
character*10  current_data
real(kind=FT),ALLOCATABLE:: Contact_DISP(:)
real(kind=FT),ALLOCATABLE:: DISP2(:)
real(kind=FT),ALLOCATABLE::diag_precon_no_invert(:)
integer i_WB,i_Stage,i_Prop
real(kind=FT) c_Stage_Q,c_Stage_Time,c_Time,c_Pres
real(kind=FT) Time_Initial_Guess,Pres_Initial_Guess
real(kind=FT) f_Vol_Tol
real(kind=FT) Max_KI_eq_3D,Min_KI_eq_3D,Ave_KI_eq_3D,KIc
integer Max_Prop_Steps,i_Time
real(kind=FT) Old_Time
real(kind=FT) NR_delta_Pres,NR_delta_Time
logical Tem_Condition_1,Tem_Condition_2
logical Yes_T_Convergent
real(kind=FT) f_K(2),f_K_Tol
real(kind=FT) Output_Pres
real(kind=FT) Applied_Time,d_f_K
real(kind=FT) Conv_Time
real(kind=FT) Time_Tol,Pres_Tol
integer,ALLOCATABLE::tem_all_local(:,:)
integer c_Crack_Max_Num,c_Vertex_Max_Num
integer c_Elem_num,c_mat_num
real(kind=FT) c_NR_delta_Time
integer i_Bisection,Max_Bisection
real(kind=FT) Bisec_Start_Time,Bisec_End_Time,Bisec_Tol_Time,Bisec_Mid_Time
real(kind=FT) c_Bisec_Gap,Bisec_F_KIc_Start,Bisec_F_KIc_End
real(kind=FT) Bisec_F_KIc_Mid
real(kind=FT) Bisec_Start_Time_Initial,Bisec_End_Time_Initial
real(kind=FT),ALLOCATABLE:: K_CSR_aa(:)
integer,ALLOCATABLE:: K_CSR_ja(:)
integer,ALLOCATABLE:: K_CSR_ia(:)
integer(kind=LIT) K_CSR_NNZ_Max
integer K_CSR_NNZ
integer Total_Num_G_P
    
!------------------------------------------
! Initial value of the temporary variable.
!------------------------------------------
Yes_Last_Growth = .False.
      
!-----------------------------
! Formatted output statement.
!-----------------------------
1002 FORMAT('     Force factor is ',F5.3)   
1021 FORMAT(5X,'Range of displacement x:  ',E14.7,' m to ',E14.7,' m')  
1022 FORMAT(5X,'Range of displacement y:  ',E14.7,' m to ',E14.7,' m')  
1023 FORMAT(5X,'Range of displacement z:  ',E14.7,' m to ',E14.7,' m')     
3001 FORMAT(5X,'Elapsed CPU time: ',I7,' s, about ',F7.2,' mins')     
4001 FORMAT(5X,'Number of crack is ',I3,' / ',I5)   
1171 FORMAT('  >> Prop step ',I3,' of stage ',I3,' /',I3,' of WB ',I3,' /',I3,' started:')   
1172 FORMAT('     Prop step count:',I3)        
4007 FORMAT(5X,'......................................................')   
4008 FORMAT(5X,'N-R Time iteration failed after ',I3,' tries!')  
4018 FORMAT(5X,'Bisection Time iteration failed after ',I3,' tries!') 
4009 FORMAT('     >> Time step ',I5,' of ',I5,T39,' started (Newton-Raphson)<<')    
4019 FORMAT('     >> Bisection step ',I5,' of ',I5,T39,' started <<')   
5001 FORMAT(5X,'Converged pressure (MPa)      :',F11.3)   
5000 FORMAT(5X,'Converged time (s)            :',F11.3)    
5002 FORMAT(5X,'Converged time (mins)         :',F11.3)   
5003 FORMAT(5X,'Converged Max_KIeq (MPa*m^1/2):',F11.3) 
5004 FORMAT(5X,'Material Para KI_c (MPa*m^1/2):',F11.3)  
5005 FORMAT(5x,'>> Time step iteration done after ',I3,' tries <<')
5006 FORMAT(5X,'Converged Ave_KIeq (MPa*m^1/2):',F11.3)
5007 FORMAT(5X,'Converged Min_KIeq (MPa*m^1/2):',F11.3)

    
!-------------------------------
! Calculate the material matrix
! D for each material.
!-------------------------------
call Get_Material_Matrix_3D
!If has composite material.
if (any(Material_Type == 5))then   
  !Get the rotation matrix for the composite material for each element.
  call Get_Ele_Composite_Mat_Rot_Matrix 
endif

!-------------------------------------
! Generate initial natural fractures. 
! NEWFTU2022061001.
!-------------------------------------
if (Key_Random_NaCr==1) then 
  call Tool_Generate_Natural_Fractures_3D
endif

!-----------------------------------------------
! Save wellbore coordinate data (if available). 
! NEWFTU2022110901.
!-----------------------------------------------
if(num_Wellbore>=1)then
    print *, "    Saving *.wbfp file of wellbore..." 
    call Save_HF_wbfp_File
endif  
    
!-------------------------------------------------------------
! Read initial natural fractures from the FracMan *.fab file. 
! NEWFTU2022072401.
!-------------------------------------------------------------
if (Key_Random_NaCr==2) then 
  call Tool_Read_fab_Natural_Fractures_3D
endif 

!----------------------------------------------------------
! Set natural fractures based on the data in the kpp file.
! 2023-08-24. NEWFTU2023082401.
!----------------------------------------------------------
if (Key_Random_NaCr==3) then 
  call Tool_Set_Natural_Fractures_by_kpp_3D
endif  

!----------------------------------------------------------------------------        
! Segment each initial crack (generate triangular geometric crack elements).
!---------------------------------------------------------------------------- 
if(num_Crack.ne.0)then
  print *,'    Meshing initial 3D cracks...'  
  call D3_Mesh_Initial_Crack
endif
      
!----------------------------------------------------
! Inspect and adjust each initial crack. 2022-08-03.
!----------------------------------------------------
if(num_Crack.ne.0 .and. Key_Check_and_Adjust_Cracks_3D>=1)then
  print *,'    Check and adjust initial 3D cracks...'  
  call D3_Check_and_Adjust_Cracks
endif 

!-----------------------------------------------------------------
! Calculate the traditional degree-of-freedom displacement field 
! under in-situ stress (far-field stress), and calculate the 
! in-situ stress level, then save it.
! Nodal stress and Gauss point stress.
! Note: There are no enhancement nodes in the model at this time.
!-----------------------------------------------------------------
if(Key_InSitu_Strategy /=0)then
  ! Key_InSitu_Strategy=4 does not require calculating the initial stress field, but is directly
  ! defined or read from a file, 2022-06-03.
  if(Key_InSitu_Strategy /=4)then
      call Cal_InSitu_Stress_3D
  endif
endif

!---------------------------------------------------------------------
! Initial ground stress treatment. Key_InSitu_Strategy=4. 2022-06-03.
! Ref: Cruz_2019_An XFEM implementation in Abaqus to model 
! intersections between fractures in porous rocks.pdf, Section 3.5
! NEWFTU2022060301.
!---------------------------------------------------------------------
if(Key_InSitu_Strategy ==4)then
      ! Obtain internal force (stress at each conventional Gauss point): obtained according to user custom
      ! settings or input file.
      ! Store variables such as InSitu_Strs_Gaus_xx(num_Elem, Num_Gauss_P_FEM_3D).
      print *,'    Getting internal stress/strain of Gauss points.'  
      call D3_Get_Internal_Stress_and_Strain
      ! Obtain the constraint reactions for all nodes and store them in Reactions_Nodes_3D(num_Node, 3).
      ! call D3_Get_InSitu_Node_Reactions ! Abandon this plan
endif
      
!---------------------------------------------------------------------------------
! Make minor adjustments to each initial crack. Hydraulic fracturing analysis 
! does not support this because the negative effects are significant. 2022-07-08.
!---------------------------------------------------------------------------------
!if(num_Crack.ne.0 .and. Key_Adj_Ini_Crack_3D==1)then
!    call D3_Adjust_Initial_Crack
!endif


Flag_HF_3D = 1
      
!-------------------------------
! Cycle between each load step.
!-------------------------------
isub = 0
!@@@@@@@@@@@@@@@@@@@@@@@
! Wellbore circulation.
!@@@@@@@@@@@@@@@@@@@@@@@
do i_WB = 1,num_Wellbore
   !@@@@@@@@@@@@@@@@@
   ! Segmented loop.
   !@@@@@@@@@@@@@@@@@
   do i_Stage = 1,num_Stages_Wellbores(i_WB)
    ! Current fracturing fluid flow rate for the current wellbore segment, unit: m^3/s
    c_Stage_Q = Injection_Q_Stages_Wellbores(i_WB,i_Stage)           
    ! Fracturing fluid injection time for the current section of the current wellbore, unit: s
    c_Stage_Time = Injection_T_Stages_Wellbores(i_WB,i_Stage)  
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Key variable definitions (for solving time steps using the Newton-Raphson method).
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Key_3D_HF_Time_Step_Method = 1 ! Time step solved using the Newton-Raphson method.
    if(Key_3D_HF_Time_Step_Method ==  1)then 
        ! Suitable for small models within 1 meter. 2024-02-21. IMPROV2024022101.
        if(Max_Model_Range < 1.0D0) then
            Time_Initial_Guess =  1.0D0
            Pres_Initial_Guess =  9.87654321D6
                                                 !which would cause F=0 during the initial step calculation, making it impossible to compute.
                                                 !BUGFIX2022091802.
            Max_Prop_Steps     =  1000
            ! SlipWater_Max_Time_Steps_3D     =  20      ! Maximum time step iterations, global variable.
            ! SlipWater_Max_Pres_Steps_3D     =  10       !Maximum pressure step iterations, global variable.
            
            ! NR_delta_Pres      =  0.01D6    !Small pressure increment for N-R iteration (default: 0.01D6)
            ! NR_delta_Time       = 0.01D0    ! Tiny time increment for N-R iteration (default: 0.01D0)
            ! f_Vol_Tol          =  1.0D-6    !Convergence tolerance for pressure iteration (meaning volume)
            ! f_K_Tol            =  0.01D6    ! Convergence tolerance for time step iteration (meaning is fracture toughness)
            ! Pres_Tol           =  1.0D-5    !Pressure change ratio. abs((tem_Pres-Old_Pres)/Old_Pres) < Pres_Tol
            ! Time_Tol           =  1.0D-5    !Time change ratio. abs((c_Time-Old_Time)/Old_Time) < Time_Tol

            NR_delta_Pres      =  0.01D6
            NR_delta_Time      =  0.01D0
            f_Vol_Tol          =  1.0D-6
            f_K_Tol            =  0.01D6
            Pres_Tol           =  1.0D-5
            Time_Tol           =  1.0D-4
        ! Suitable for small models within the range of 1 meter to 100 meters. 2022-09-10.
        elseif(Max_Model_Range < 100.0D0 .and. Max_Model_Range >= 1.0D0) then
            Time_Initial_Guess =  1.0D0
            Pres_Initial_Guess =  9.87654321D6
                                                 !which would cause F=0 during the initial step calculation, making it impossible to compute.
                                                 !BUGFIX2022091802.
            Max_Prop_Steps     =  1000
            ! SlipWater_Max_Time_Steps_3D     =  20      ! Maximum time step iterations, global variable.
            ! SlipWater_Max_Pres_Steps_3D     =  10       !Maximum pressure step iterations, global variable.
            NR_delta_Pres      =  0.01D6
            NR_delta_Time      =  0.01D0
            f_Vol_Tol          =  1.0D-6
            f_K_Tol            =  0.01D6
            Pres_Tol           =  1.0D-4
            Time_Tol           =  1.0D-4
        ! Applicable to small models, kilometer-scale models, and large models. 2022-09-10.
        ! IMPROV2022091001.
        elseif(Max_Model_Range>=100.0D0) then
            Time_Initial_Guess =  1.0D0
            Pres_Initial_Guess =  9.87654321D6
            Max_Prop_Steps     =  1000
            ! SlipWater_Max_Time_Steps_3D     =  20      ! Maximum time step iterations, global variable.
            ! SlipWater_Max_Pres_Steps_3D     =  10       !Maximum pressure step iterations, global variable.
            NR_delta_Pres      =  0.01D6
            NR_delta_Time      =  0.01D0
            f_Vol_Tol          =  1.0D-5
            f_K_Tol            =  0.2D6
            ! Pres_Tol           =  1.0D-2   ! Pressure change ratio. abs((tem_Pres-Old_Pres)/Old_Pres) < Pres_Tol
            ! Time_Tol           =  2.0D-2   ! Time change ratio. abs((c_Time-Old_Time)/Old_Time) < Time_Tol
            
            !2023-09-21. IMPROV2023092101.
            Pres_Tol           =  5.0D-2
            Time_Tol           =  5.0D-2
        endif
        
    endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Key variable definitions (for solving time steps using the bisection method). 2023-05-23.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !2023-05-18.
    if(Key_3D_HF_Time_Step_Method ==  2)then 
        ! Used for pressure step Newton-Raphson iteration.
        Pres_Initial_Guess =  9.87654321D6
        Max_Prop_Steps     =  1000
        ! SlipWater_Max_Pres_Steps_3D     =  10       !Maximum pressure step cycles
        NR_delta_Pres      =  0.01D6
        f_Vol_Tol          =  1.0D-5
        Pres_Tol           =  1.0D-2
        ! Used for time-step bisection iteration.
        Max_Bisection = 40
        Bisec_Start_Time_Initial   = 0.001D0
        ! Bisec_Start_Time_Initial = -3600.0D0    !Start time of the bisection method time step.
        ! Bisec_Start_Time_Initial = 0.1D0    !Starting point of the bisection time step.
        ! Bisec_Start_Time_Initial = 1.0D0    !Starting point of the bisection time step.
        ! Bisec_End_Time_Initial   = 3600.0D0      ! Bisection method time step endpoint (1 hour).
        Bisec_End_Time_Initial   = c_Stage_Time
        ! Bisec_End_Time_Initial = 1.0D7 !Bisection method time step endpoint (take a relatively large
        ! value, approximately 115 days).
        ! Bisec_End_Time_Initial = 1.0D8 !Bisection method time step endpoint (take a relatively large
        ! value, approximately 3.2 years).
        Bisec_Tol_Time     = 1.0D0
    endif
    
    
    !@@@@@@@@@@@@@@@@@@@@@@
    ! Fracture step cycle.
    !@@@@@@@@@@@@@@@@@@@@@@
    do i_Prop = 1,Max_Prop_Steps
      ! If it is the first rupture step.
      if(i_Prop==1) then
          c_Time = Time_Initial_Guess
          c_Pres = Pres_Initial_Guess
      endif

      print *, "  " 
          
      isub = isub + 1
      
      ! Information output.
      WRITE(*,1171) i_Prop,i_Stage,num_Stages_Wellbores(i_WB),i_WB,num_Wellbore
      WRITE(*,1172) isub          

      !********************************************************
      ! Stage fracturing generates initial cracks, 2022-05-31.
      !********************************************************
      ! Generate initial fractures along the wellbore (WB).
      call D3_HF_Generate_Initial_Cracks_of_WB(i_Prop,i_WB,i_Stage)
      
      !***************************************************************************************
      ! Perform mesh segmentation on each initial crack (generating triangular 
      ! geometric crack elements), only for cracks that have not been discretely initialized.
      !***************************************************************************************
      if(num_Crack.ne.0)then
          print *,'    Meshing initial 3D cracks...'  
          call D3_Mesh_Initial_Crack
      endif
          
      !****************************************************
      ! Inspect and adjust each initial crack. 2022-08-03.
      !****************************************************
      if(num_Crack.ne.0 .and. Key_Check_and_Adjust_Cracks_3D>=1)then
          print *,'    Check and adjust initial 3D cracks...'  
          call D3_Check_and_Adjust_Cracks 
      endif   
      
      !******************************************
      ! Mesh the newly added cracks. 2023-01-09.
      !******************************************
      print *,'    Meshing for 3D new cracks...'  
      call D3_Mesh_Initial_Crack
      
      !****************************************************
      ! Inspect and adjust each initial crack. 2023-01-09.
      !****************************************************
      if(num_Crack.ne.0 .and. Key_Check_and_Adjust_Cracks_3D>=1)then
          print *,'    Check and adjust initial 3D cracks...'  
          call D3_Check_and_Adjust_Cracks 
      endif   
      
      !****************************
      ! Obtain the payload factor.
      !****************************
      call Force_Factor(Lambda,isub,Yes_Last_Growth)
      WRITE(*,1002) Lambda
      
      !************************************
      ! If there are no cracks, then stop.
      !************************************
      if (num_Crack.eq.0)then 
          ! The number of cracks must be greater than 0, otherwise an error will pop up.
          print *,'    Warning :: num_Crack=0!'
          call Warning_Message('S',Keywords_Blank) 
      end if
          

      !******************************************************
      ! Check whether the maximum allowable number of cracks 
      ! has been exceeded (2021-08-21).
      !******************************************************
      write (*,4001) num_Crack,Max_Num_Cr_3D
      if(num_Crack > Max_Num_Cr_3D)then
          print *, '    Error :: num_Crack > Max_Num_Cr_3D!'
          print *, '             Error in PhiPsi3D_Static.f.'
          call Warning_Message('S',Keywords_Blank)  
      endif
      
      !***********************************************
      ! Adjust the crack edge nodes (to prevent them 
      ! from being on the surface of solid elements).
      !***********************************************
      print *,'    Adjusting crack vertexs...'   
      call D3_Crack_Vertex_Adjust(isub)
      
      !****************************************************************
      ! Calculate the local coordinate system of the crack edge nodes.
      !****************************************************************
      print *,'    Calculating local coordinate system of crack vertexs...'   
      call D3_Crack_Vertex_Local_Coor_Sys(isub)

      !***************************
      ! Confirm enhancement node.
      !***************************
      ifra = 1
      print *,'    Determine enriched node...'
      call Determine_Enriched_Nodes_3D(ifra,isub)
  
      
      !******************************************
      ! Assign a number to the enhancement node.
      !******************************************
      call Number_Enriched_Nodes_3D(isub)     
      print *,'    Total_FD:',Total_FD  
      print *,'    Enrich_Freedom:',Enrich_Freedom  

      !***********************************
      ! Consider the boundary conditions.
      !***********************************
      ALLOCATE(freeDOF(Total_FD))
      ALLOCATE(fixedDOF(Total_FD))
      
      call Boundary_Cond_3D(Total_FD,isub,freeDOF,num_FreeD,fixedDOF,num_FixedD)    

      !*************************************
      ! Show degrees of freedom information
      !*************************************
1101  FORMAT(5X,'Total DOFs of solid:       ',I8)
1102  FORMAT(5X,'Total Enriched DOFs:       ',I8) 
1103  FORMAT(5X,'Total active DOFs of solid:',I8)   
      WRITE(*,1101) Total_FD 
      WRITE(*,1102) Enrich_Freedom  
      WRITE(*,1103) num_FreeD


      !*****************************************************************
      ! Obtain the list of newly added XFEM elements and the 
      ! list of XFEM elements that need to update the stiffness matrix.
      ! 2022-06-24. NEWFTU2022062403.
      !*****************************************************************
      call D3_Get_New_and_Updated_XFEM_Elements(isub) 

      
      !******************************************************************************
      ! Calculate related data at fracture calculation points (fluid element nodes).
      !******************************************************************************
1104  FORMAT(5X,'Total DOFs of fluid:',I8)   

      call D3_Cal_HF_Crack_Points_Info_Linear(isub)
      
      num_Tol_CalP_Water = 0
      do i_C=1,Num_Crack
        ! If it is an HF fracture
        if(Crack_Type_Status_3D(i_C,1) == 1) then
            num_Tol_CalP_Water= num_Tol_CalP_Water +Cracks_Real_CalP_Num_3D(i_C)        
        endif
      end do 
      
      WRITE(*,1104) num_Tol_CalP_Water                                  

      !************************************************
      ! Allocate space for other related dynamic data.
      !************************************************
      if (.not. allocated(F_U))ALLOCATE(F_U(Total_FD))
      if (.not. allocated(Coupled_Q_3D)) then
          ! ALLOCATE(Coupled_Q_3D(Total_FD, num_Tol_CalP_Water))! Fluid-structure coupling matrix
          ! Using an uneven array. 2022-09-16
          ALLOCATE(Coupled_Q_3D(Total_FD))          
      endif
      if (.not. allocated(Coupled_Q_3D_Index)) then
          ALLOCATE(Coupled_Q_3D_Index(Total_FD))          
      endif      
      
      if (.not. allocated(F))ALLOCATE(F(Total_FD))
      !ALLOCATE(CalP_Pres(num_Tol_CalP_Water))        
      if (.not. allocated(delta_x))ALLOCATE(delta_x(num_FreeD))
      
      !****************************
      ! Initialize some variables.
      !****************************
      F(1:Total_FD)             = ZR
      F_U(1:Total_FD)           = ZR

      !**********************************
      ! Calculate the coupling matrix Q.
      !**********************************
      if(Flag_HF_3D==1) then
#ifndef Silverfrost
          call Cal_HF_Matrix_Q_Linear_3D(isub,Total_FD,num_Tol_CalP_Water)  
#endif
#ifdef Silverfrost
          print *,'    ERROR :: Silverfrost compiler failed to compile Cal_HF_Matrix_Q_Linear_3D.f90!'
          print *,'             In PhiPsi3D_Static_HF_SlipWater.F90.'
          call Warning_Message('S',Keywords_Blank)
#endif              
      endif     

      !************************************************
      ! Assembly load vector (without fluid pressure).
      !************************************************
      call Force_Vector_3D(Total_FD,isub,Lambda,F_U)
      
      !*****************************************************
      ! Initialization of variables for the PCG-EBE solver.
      !*****************************************************
      if (.not. allocated(DISP))ALLOCATE(DISP(Total_FD))
      DISP  = ZR
      if (.not. allocated(DISP2))ALLOCATE(DISP2(Total_FD))
      DISP2 = ZR   
      
      if(Key_SLOE==11)then  
          !Element by Element, diagonallypreconditioned conjugate gradient solver 
          print *,'    PCG-EBE: Preparing...' 
          if (.not. allocated(storK_XFEM)) then
              ALLOCATE(storK_XFEM(num_XFEM_Elem))
          endif
          ! storK_FEM is calculated only once because it will not change afterward and will only decrease.
          if(Key_EBE_Sym_Storage_K==0)then
              if (allocated(storK_FEM).eqv. .False.) then
                  ALLOCATE(storK_FEM(24,24,num_FEM_Elem))  
              endif
          endif
          if(Key_EBE_Sym_Storage_K==1)then
              if (allocated(storK_FEM_Sym).eqv. .False.) then
                  ALLOCATE(storK_FEM_Sym(300,num_FEM_Elem))  
              endif
          endif      
          if (.not. allocated(Size_Local_3D)) then
              ALLOCATE(Size_Local_3D(Num_Elem))
          endif
          if (.not. allocated(All_Local_3D))  then
              ALLOCATE(All_Local_3D(MDOF_3D,Num_Elem))   
          endif

          ! Expand the all_local array. !IMPROV2023011101. 2023-01-11. Avoid the problem of insufficient
          ! MDOF_3D.
          if(MDOF_3D > Old_MDOF_3D)then
#ifndef Silverfrost  
              print *,'    Extending matrix All_Local_3D...'  
              allocate(tem_all_local(MDOF_3D,Num_Elem))  
              tem_all_local(1:MDOF_3D,1:Num_Elem)  = 0
              tem_all_local(1:Old_MDOF_3D,1:Num_Elem) = All_Local_3D(1:Old_MDOF_3D,1:Num_Elem)
              deallocate(All_Local_3D)
              call move_alloc(tem_all_local,All_Local_3D)
#endif
#ifdef Silverfrost  
              print *, '    Error :: Extending matrix All_Local_3D is not valid for Silverfrost compiler!'
              print *, '             move_alloc function in Silverfrost is faulty, thus code is commented out!'
              print *, '             In PhiPsi3D_Static_HF_SlipWater.F90.'
              call Warning_Message('S',Keywords_Blank)
#endif
          endif
      endif
      
      if (allocated(diag_precon_no_invert).eqv. .False.) then
          ALLOCATE(diag_precon_no_invert(0:num_FreeD))
      endif
          
      !********************************
      ! Assemble the stiffness matrix.
      !********************************
      if(Key_SLOE==11)then  
          call EBE_XFEM_PCG_3D_Get_K(isub,num_FreeD,freeDOF(1:num_FreeD),      &
                                   diag_precon_no_invert)
          if (allocated(K_CSR_aa)) DEALLOCATE(K_CSR_aa)
          if (allocated(K_CSR_ja)) DEALLOCATE(K_CSR_ja)
          if (allocated(K_CSR_ia)) DEALLOCATE(K_CSR_ia)
          K_CSR_NNZ = 1
          ALLOCATE(K_CSR_aa(K_CSR_NNZ))
          ALLOCATE(K_CSR_ja(K_CSR_NNZ))
          ALLOCATE(K_CSR_ia(num_FreeD+1))
      !NEWFTU-2025120201.
      else
#ifdef gfortran
!#ifndef github
          print *,'    Get max NNZ before assembling K......'
          call Assemble_Stiffness_Matrix_SPARS_XFEM_3D_Get_MaxNNZ(isub,&
                      freeDOF(1:num_FreeD),num_FreeD,Total_FD,K_CSR_NNZ_Max)
          print *,'    K_CSR_NNZ_Max:',K_CSR_NNZ_Max    
          
          if (allocated(K_CSR_aa)) DEALLOCATE(K_CSR_aa)
          if (allocated(K_CSR_ja)) DEALLOCATE(K_CSR_ja)
          if (allocated(K_CSR_ia)) DEALLOCATE(K_CSR_ia)
          ALLOCATE(K_CSR_aa(K_CSR_NNZ_Max))
          ALLOCATE(K_CSR_ja(K_CSR_NNZ_Max))
          ALLOCATE(K_CSR_ia(num_FreeD+1))
          print *,'    Assemble the sparse K...'
          ! Assemble the sparse stiffness matrix (automatically considers boundary conditions and removes
          ! sparse matrix elements corresponding to constrained degrees of freedom)
          call Assemble_Stiffness_Matrix_SPARS_XFEM_3D(isub, &
               freeDOF(1:num_FreeD),num_FreeD,&
               fixedDOF(1:num_FixedD),num_FixedD,     &
               K_CSR_aa,K_CSR_ja,K_CSR_ia,&
               K_CSR_NNZ_Max,K_CSR_NNZ,&
               Total_FD,Total_Num_G_P)
          diag_precon_no_invert =  ONE
!#endif
#endif
      endif
      
      ! Save the number of Gaussian integration points for each element to a *.elgn file for
      ! post-processing. NEWFTU2022071601.
      if (Key_Post_Elements_Gauss_Num== 1) then
          call Save_Gauss_Num_of_Elements(isub)
      endif
      
      
      !*****************
      ! Time step loop.
      !*****************
      print *,'    Start time-step loop...' 
      Yes_T_Convergent = .False.
      
      !////////////////////////////////////////////////////////////////////////
      !                                                                      !
      !                                                                      !
      ! Solve the time step using the Newton-Raphson method.                 !
      ! Note: Currently, convergence is poor under non-uniform stress fields.!
      !                                                                      !
      !////////////////////////////////////////////////////////////////////////
      if(Key_3D_HF_Time_Step_Method==1) then
        do i_Time = 1,SlipWater_Max_Time_Steps_3D
            ! Screen output.
            print *,' '
            write(*,4009) i_Time,SlipWater_Max_Time_Steps_3D
            
            ! Current time.
            Applied_Time = c_Time
            
            call D3_HF_Get_Pres_by_NR_Method(i_WB,i_Stage,i_Prop,                             &
                    isub,i_Time,c_Stage_Q,Applied_Time,SlipWater_Max_Pres_Steps_3D,           &
                    num_FreeD,Total_FD,num_Tol_CalP_Water,freeDOF(1:num_FreeD),F_U,F,Lambda,  &
                    diag_precon_no_invert(0:num_FreeD),                                       &
                    DISP,NR_delta_Pres,f_Vol_Tol,Pres_Tol,c_Pres,Output_Pres,                 &
                    K_CSR_NNZ,K_CSR_aa(1:K_CSR_NNZ),K_CSR_ja(1:K_CSR_NNZ),K_CSR_ia(1:num_FreeD+1))   
                    
            c_Pres = Output_Pres

            print *,' '
          
            ! Calculate the stress intensity factor (including the equivalent stress intensity factor).
            print *,'       Calculating stress intensity factors...'   
            call Cal_SIFs_DIM_3D(isub,DISP,8)
            ! Extract the maximum, minimum, and average equivalent stress intensity factors of the HF fracture.
            call D3_Get_Max_Min_Ave_KIeq_of_HF_Cracks(Max_KI_eq_3D,Min_KI_eq_3D,Ave_KI_eq_3D,c_Crack_Max_Num,c_Vertex_Max_Num)
            print *,'       Max_KI_eq_3D (MPa*m^1/2):',Max_KI_eq_3D/1.0D6
            print *,'       Min_KI_eq_3D (MPa*m^1/2):',Min_KI_eq_3D/1.0D6
            print *,'       Ave_KI_eq_3D (MPa*m^1/2):',Ave_KI_eq_3D/1.0D6
            ! KIc
            ! Obtain the element number of the vertex at the leading edge of the crack with the maximum
            ! equivalent stress. 2023-05-16. NEWFTU2023051601.
            c_Elem_num = Cr3D_Meshed_Node_in_Ele_Num(c_Crack_Max_Num)%row(c_Vertex_Max_Num)
            c_mat_num = Elem_Mat(c_Elem_num)
            KIc   = Material_Para(c_mat_num,6)
            
            ! Calculate f_K(1).
            if(Key_3D_HF_SlipWater_fk_Type==1) then
                f_K(1) = Max_KI_eq_3D - KIc
            elseif(Key_3D_HF_SlipWater_fk_Type==2)then
                f_K(1) = Ave_KI_eq_3D - KIc
            elseif(Key_3D_HF_SlipWater_fk_Type==3)then
                f_K(1) = Min_KI_eq_3D - KIc
            endif
                
            c_NR_delta_Time = NR_delta_Time
            Applied_Time = c_Time + c_NR_delta_Time
            
            call D3_HF_Get_Pres_by_NR_Method(i_WB,i_Stage,i_Prop,                               &
                    isub,i_Time,c_Stage_Q,Applied_Time,SlipWater_Max_Pres_Steps_3D,             &
                    num_FreeD,Total_FD,num_Tol_CalP_Water,freeDOF(1:num_FreeD),F_U,F,Lambda,    &
                    diag_precon_no_invert(0:num_FreeD),                                         &
                    DISP2,NR_delta_Pres,f_Vol_Tol,Pres_Tol,c_Pres,Output_Pres,                  &   
                    K_CSR_NNZ,K_CSR_aa(1:K_CSR_NNZ),K_CSR_ja(1:K_CSR_NNZ),K_CSR_ia(1:num_FreeD+1))  
            
            print *,' '
            ! Calculate the stress intensity factor (including the equivalent stress intensity factor).
            print *,'       Calculating stress intensity factors...'   
            call Cal_SIFs_DIM_3D(isub,DISP2,8)
            ! Extract the maximum, minimum, and average equivalent stress intensity factors of the HF fracture.
            call D3_Get_Max_Min_Ave_KIeq_of_HF_Cracks(Max_KI_eq_3D,Min_KI_eq_3D,Ave_KI_eq_3D,c_Crack_Max_Num,c_Vertex_Max_Num)
            print *,'       Max_KI_eq_3D (MPa*m^1/2):',Max_KI_eq_3D/1.0D6
            print *,'       Min_KI_eq_3D (MPa*m^1/2):',Min_KI_eq_3D/1.0D6
            print *,'       Ave_KI_eq_3D (MPa*m^1/2):',Ave_KI_eq_3D/1.0D6
            
            ! Calculate f_K(2).
            if(Key_3D_HF_SlipWater_fk_Type==1) then
                f_K(2) = Max_KI_eq_3D - KIc
            elseif(Key_3D_HF_SlipWater_fk_Type==2)then
                f_K(2) = Ave_KI_eq_3D - KIc
            elseif(Key_3D_HF_SlipWater_fk_Type==3)then
                f_K(2) = Min_KI_eq_3D - KIc
            endif
            
            ! Calculate the derivative d_f_K of f_K.
            print *,'       Calculate the derivative of f_K...'
            d_f_K = (f_K(2)-f_K(1))/c_NR_delta_Time
            print *,'       The derivative of f_K:',d_f_K 
            
            ! Update the time step according to the Newton-Raphson method.
            Old_Time = c_Time
            c_Time   = c_Time - f_K(1)/d_f_K
            print *,'       Updated time (s):',c_Time
            
            ! Convergence check.
            Tem_Condition_1 = .False.
            Tem_Condition_2 = .False.
            
            if(abs(f_K(1))<= f_K_Tol) then                
               Tem_Condition_1 =  .True.
            endif
            print *, '       Abs(f_K(1)),f_K_Tol:', abs(f_K(1))/1.0D6,f_K_Tol/1.0D6
            
            Conv_Time = abs((c_Time-Old_Time)/Old_Time)

            
            if(Conv_Time <= Time_Tol) then                         
               Tem_Condition_2 =  .True.
            endif  
            print *, '       Conv_Time,Time_Tol:', Conv_Time,Time_Tol 
            
            
            if(Tem_Condition_1 .and. Tem_Condition_2) then
               ! Update DISP. 2023-05-16. BUGFIX2023051602.
               DISP = DISP2
               
               print *,' '
               write(*,5005) i_Time
               write(*,*) '    .........................................'
               Yes_T_Convergent = .True.
               write(*,5001) c_Pres/1.0D6
               write(*,5000) c_Time
               write(*,5002) c_Time/60.0D0
               write(*,5003) Max_KI_eq_3D/1.0D6
               write(*,5006) Ave_KI_eq_3D/1.0D6
               write(*,5007) Min_KI_eq_3D/1.0D6
               write(*,5004) KIc/1.0D6
               write(*,*) '    .........................................'
               print *,' '
               ! Save pressure and time data to a *.wbpt file
               call Save_HF_wbpt_File(i_WB,i_Stage,i_Prop,c_Pres,c_Time)
               exit
            endif
             ! If convergence is not achieved after all steps are completed, an error message will be displayed.
            if(i_Time ==SlipWater_Max_Time_Steps_3D .and. SlipWater_Time_Step_Conv_Check==1) then
               write(*,4007) 
               write(*,4008) SlipWater_Max_Time_Steps_3D
               write(*,4007) 
               call Warning_Message('S',Keywords_Blank) 
            endif
        enddo
      
      !////////////////////////////////////////////////////
      !                                                  !
      !                                                  !
      ! Bisection method to solve time step. 2023-05-18. !
      ! NEWFTU2023051801.                                !
      ! https://en.wikipedia.org/wiki/Bisection_method   !
      !                                                  !
      !////////////////////////////////////////////////////
      elseif(Key_3D_HF_Time_Step_Method==2) then 
        ! Binary loop.
        do i_Bisection=1,Max_Bisection
        
            ! Screen output.
            print *,' '
            write(*,4019) i_Bisection,Max_Bisection
            
            !=======================================================
            ! The first bisection step requires determining 
            ! whether the start and end points have opposite signs.
            !=======================================================
            if(i_Bisection==1)then
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                ! Initialize the start and end time, without
                ! inheriting the results from the previous step.
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                Bisec_Start_Time = Bisec_Start_Time_Initial
                Bisec_End_Time = Bisec_End_Time_Initial
                print *,'       Bisec_Start_Time:',Bisec_Start_Time
                print *,'       Bisec_End_Time:',Bisec_End_Time
                
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                ! Perform the pressure iteration corresponding
                ! to the Bisec_Start_Time time point.
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                Applied_Time = Bisec_Start_Time
                
                c_Bisec_Gap = (Bisec_End_Time-Bisec_Start_Time)/TWO
                
                call D3_HF_Get_Pres_by_NR_Method(i_WB,i_Stage,i_Prop,                              &
                        isub,i_Bisection,c_Stage_Q,Applied_Time,SlipWater_Max_Pres_Steps_3D,       &
                        num_FreeD,Total_FD,num_Tol_CalP_Water,freeDOF(1:num_FreeD),F_U,F,Lambda,   &
                        diag_precon_no_invert(0:num_FreeD),                                        &
                        DISP,NR_delta_Pres,f_Vol_Tol,Pres_Tol,c_Pres,Output_Pres,                  &   
                        K_CSR_NNZ,K_CSR_aa(1:K_CSR_NNZ),K_CSR_ja(1:K_CSR_NNZ),K_CSR_ia(1:num_FreeD+1))  
                        
                c_Pres = Output_Pres
                
                print *,'       Applied pressure at start time (MPa):',c_Pres/1.0D6   
                
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                ! Calculate the stress intensity factor (including the equivalent 
                ! stress intensity factor) corresponding to the Bisec_Start_Time time point.
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                call Cal_SIFs_DIM_3D(isub,DISP,8)
                ! Extract the maximum, minimum, and average equivalent stress intensity factors of the HF fracture.
                call D3_Get_Max_Min_Ave_KIeq_of_HF_Cracks(Max_KI_eq_3D,Min_KI_eq_3D,Ave_KI_eq_3D,c_Crack_Max_Num,c_Vertex_Max_Num)
                print *,'       Max_KI_eq_3D at start time (MPa*m^1/2):',Max_KI_eq_3D/1.0D6
                print *,'       Min_KI_eq_3D at start time (MPa*m^1/2):',Min_KI_eq_3D/1.0D6
                print *,'       Ave_KI_eq_3D at start time (MPa*m^1/2):',Ave_KI_eq_3D/1.0D6
                !KIc
                ! Obtain the element number of the vertex at the leading edge of the crack with the maximum
                ! equivalent stress. 2023-05-16. NEWFTU2023051601.
                c_Elem_num = Cr3D_Meshed_Node_in_Ele_Num(c_Crack_Max_Num)%row(c_Vertex_Max_Num)
                c_mat_num = Elem_Mat(c_Elem_num)
                KIc   = Material_Para(c_mat_num,6)
                
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                ! Calculate the value of function F at the 
                ! time point Bisec_Start_Time.
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                !Bisec_F_KIc_Start = Max_KI_eq_3D - KIc
                if(Key_3D_HF_SlipWater_fk_Type==1) then
                    Bisec_F_KIc_Start = Max_KI_eq_3D - KIc
                elseif(Key_3D_HF_SlipWater_fk_Type==2)then
                    Bisec_F_KIc_Start = Ave_KI_eq_3D - KIc  
                elseif(Key_3D_HF_SlipWater_fk_Type==3)then
                    Bisec_F_KIc_Start = Min_KI_eq_3D - KIc  
                endif
                print *,'       Bisec_F_KIc_Start (MPa*m^(1/2)):',Bisec_F_KIc_Start/1.0D6
                
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                ! Perform the pressure iteration corresponding
                ! to the Bisec_End_Time time point.
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                Applied_Time = Bisec_End_Time
                !Applied_Time = 200.0D0
                call D3_HF_Get_Pres_by_NR_Method(i_WB,i_Stage,i_Prop,                               &
                        isub,i_Bisection,c_Stage_Q,Applied_Time,SlipWater_Max_Pres_Steps_3D,        &
                        num_FreeD,Total_FD,num_Tol_CalP_Water,freeDOF(1:num_FreeD),F_U,F,Lambda,    &
                        diag_precon_no_invert(0:num_FreeD),                                         &
                        DISP2,NR_delta_Pres,f_Vol_Tol,Pres_Tol,c_Pres,Output_Pres,                  &   
                        K_CSR_NNZ,K_CSR_aa(1:K_CSR_NNZ),K_CSR_ja(1:K_CSR_NNZ),K_CSR_ia(1:num_FreeD+1))    
                !c_Pres = Output_Pres
                print *,'       Applied pressure at end time (MPa):',Output_Pres/1.0D6 
                
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                ! Calculate the stress intensity factor (including the equivalent
                ! stress intensity factor) corresponding to the Bisec_End_Time point.
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                call Cal_SIFs_DIM_3D(isub,DISP2,8)
                ! Extract the maximum, minimum, and average equivalent stress intensity factors of the HF fracture.
                call D3_Get_Max_Min_Ave_KIeq_of_HF_Cracks(Max_KI_eq_3D,Min_KI_eq_3D,Ave_KI_eq_3D,c_Crack_Max_Num,c_Vertex_Max_Num)
                print *,'       Max_KI_eq_3D at end time (MPa*m^1/2):',Max_KI_eq_3D/1.0D6
                print *,'       Min_KI_eq_3D at end time (MPa*m^1/2):',Min_KI_eq_3D/1.0D6
                print *,'       Ave_KI_eq_3D at end time (MPa*m^1/2):',Ave_KI_eq_3D/1.0D6
                
                
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                ! Calculate the value of function F at the
                ! time point Bisec_End_Time.
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                !Bisec_F_KIc_End = Max_KI_eq_3D - KIc
                if(Key_3D_HF_SlipWater_fk_Type==1) then
                    Bisec_F_KIc_End = Max_KI_eq_3D - KIc
                elseif(Key_3D_HF_SlipWater_fk_Type==2)then
                    Bisec_F_KIc_End = Ave_KI_eq_3D - KIc  
                elseif(Key_3D_HF_SlipWater_fk_Type==3)then
                    Bisec_F_KIc_End = Min_KI_eq_3D - KIc  
                endif
                print *,'       Bisec_F_KIc_End (MPa*m^(1/2)):',Bisec_F_KIc_End/1.0D6
                
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                ! Determine whether the signs of the starting point 
                ! Bisec_F_KIc_Start and the ending point Bisec_F_KIc_End are opposite.
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                if(Bisec_F_KIc_Start*Bisec_F_KIc_End < ZR) then
                
                else
                   print *,'       Error-2023052301 :: Bisec_F_KIc_Start*Bisec_F_KIc_End >0!'
                   print *,'                           Bisec_F_KIc_Start (MPa*m^(1/2)):',Bisec_F_KIc_Start/1.0D6
                   print *,'                           Bisec_F_KIc_End (MPa*m^(1/2)):',Bisec_F_KIc_End/1.0D6
                   call Warning_Message('S',Keywords_Blank) 
                endif
            !=============
            ! Next steps.
            !=============
            else
                !BBBBBBBBBBBBBBBBBBB
                !  Update Mid time.
                !BBBBBBBBBBBBBBBBBBB
                Bisec_Mid_Time = (Bisec_Start_Time + Bisec_End_Time)/TWO
                
                c_Bisec_Gap = (Bisec_End_Time-Bisec_Start_Time)/TWO
                
                print *,'       Bisec_Start_Time:',Bisec_Start_Time
                print *,'       Bisec_End_Time:',Bisec_End_Time
                
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                ! Perform the pressure iteration corresponding
                ! to the Bisec_Mid_Time time point.
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                Applied_Time = Bisec_Mid_Time
                
                call D3_HF_Get_Pres_by_NR_Method(i_WB,i_Stage,i_Prop,                             &
                        isub,i_Bisection,c_Stage_Q,Applied_Time,SlipWater_Max_Pres_Steps_3D,      &
                        num_FreeD,Total_FD,num_Tol_CalP_Water,freeDOF(1:num_FreeD),F_U,F,Lambda,  &
                        diag_precon_no_invert(0:num_FreeD),                                       &
                        DISP,NR_delta_Pres,f_Vol_Tol,Pres_Tol,c_Pres,Output_Pres,                  &   
                        K_CSR_NNZ,K_CSR_aa(1:K_CSR_NNZ),K_CSR_ja(1:K_CSR_NNZ),K_CSR_ia(1:num_FreeD+1))    
                        
                c_Pres = Output_Pres

                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                ! Calculate the stress intensity factor (including the equivalent 
                ! stress intensity factor) corresponding to the Bisec_Mid_Time time point.
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                call Cal_SIFs_DIM_3D(isub,DISP,8)
                ! Extract the maximum, minimum, and average equivalent stress intensity factors of the HF fracture.
                call D3_Get_Max_Min_Ave_KIeq_of_HF_Cracks(Max_KI_eq_3D,Min_KI_eq_3D,Ave_KI_eq_3D,c_Crack_Max_Num,c_Vertex_Max_Num)
                !KIc
                ! Obtain the element number of the vertex at the leading edge of the crack with the maximum
                ! equivalent stress. 2023-05-16. NEWFTU2023051601.
                c_Elem_num = Cr3D_Meshed_Node_in_Ele_Num(c_Crack_Max_Num)%row(c_Vertex_Max_Num)
                c_mat_num = Elem_Mat(c_Elem_num)
                KIc   = Material_Para(c_mat_num,6)
                
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                ! Calculate the value of the F function at 
                ! the Bisec_Mid_Time point.
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                !Bisec_F_KIc_Mid = Max_KI_eq_3D - KIc
                if(Key_3D_HF_SlipWater_fk_Type==1) then
                    Bisec_F_KIc_Mid = Max_KI_eq_3D - KIc
                elseif(Key_3D_HF_SlipWater_fk_Type==2)then
                    Bisec_F_KIc_Mid = Ave_KI_eq_3D - KIc  
                elseif(Key_3D_HF_SlipWater_fk_Type==3)then
                    Bisec_F_KIc_Mid = Min_KI_eq_3D - KIc  
                endif
                print *,'       Bisec_F_KIc_Mid:',Bisec_F_KIc_Mid
                
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                ! Symbol judgment.
                !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                print *,'       Bisec_Gap:',c_Bisec_Gap
                if(Bisec_F_KIc_Mid*Bisec_F_KIc_Start > ZR) then
                    Bisec_Start_Time  = Bisec_Mid_Time
                    Bisec_F_KIc_Start = Bisec_F_KIc_Mid
                else
                    Bisec_End_Time  = Bisec_Mid_Time
                    Bisec_F_KIc_End = Bisec_F_KIc_Mid
                endif
            endif
            
            
            !BBBBBBBBBBBBBBBBBBBBB
            ! Convergence check.
            !BBBBBBBBBBBBBBBBBBBBB
            if(c_Bisec_Gap < Bisec_Tol_Time) then
               c_Time = Bisec_Mid_Time
               
               ! Update DISP. 2023-05-16. BUGFIX2023051602.
               DISP = DISP2
               
               print *,' '
               write(*,5005) i_Bisection
               write(*,*) '    .........................................'
               Yes_T_Convergent = .True.
               write(*,5001) c_Pres/1.0D6
               write(*,5000) c_Time
               write(*,5002) c_Time/60.0D0
               write(*,5003) Max_KI_eq_3D/1.0D6
               write(*,5006) Ave_KI_eq_3D/1.0D6
               write(*,5007) Min_KI_eq_3D/1.0D6
               write(*,5004) KIc/1.0D6
               write(*,*) '    .........................................'
               print *,' '
               ! Save pressure and time data to a *.wbpt file
               call Save_HF_wbpt_File(i_WB,i_Stage,i_Prop,c_Pres,c_Time)
               exit
            endif
            
            !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
            ! If convergence is still not achieved after 
            ! all steps are completed, an error message will be
            ! displayed.
            !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
            if(i_Bisection == Max_Bisection .and. SlipWater_Time_Step_Conv_Check==1) then
               write(*,4007) 
               write(*,4018) SlipWater_Max_Time_Steps_3D
               write(*,4007) 
               call Warning_Message('S',Keywords_Blank) 
            endif
        enddo
      endif
      
      
      !***********************************************************
      ! Save node payload (including enhanced nodes)
      ! For Matlab post-processing, 2022-04-23, NEWFTU2022042301.
      !***********************************************************
      call Save_Dof_Force_3D(isub,F,Total_FD)
      
      
      !********************
      ! Save displacement.
      !********************
      Max_Disp_x = maxval(DISP(1:3*Num_Node:3))
      Min_Disp_x = minval(DISP(1:3*Num_Node:3))
      Max_Disp_y = maxval(DISP(2:3*Num_Node:3))
      Min_Disp_y = minval(DISP(2:3*Num_Node:3))
      Max_Disp_z = maxval(DISP(3:3*Num_Node:3))
      Min_Disp_z = minval(DISP(3:3*Num_Node:3))
      WRITE(*,1021) Min_Disp_x,Max_Disp_x
      WRITE(*,1022) Min_Disp_y,Max_Disp_y
      WRITE(*,1023) Min_Disp_z,Max_Disp_z   
      call Save_Disp(isub,1)
      
      !**********************************************
      ! Calculate fracture permeability. 2022-11-26.
      !**********************************************
      print *,'    Calculating crack and element permeability......'   
#ifndef github
      call Cal_Crack_Permeability_3D(isub)
#endif
    
      !*********************************************************************
      ! Save crack-related files (including enhanced node numbering c_POS).
      !*********************************************************************
      call Save_Files_Crack_3D(isub)    
      
      !***************************************
      ! Save crack width related information.
      !***************************************
      call Save_Files_Cr_CalP_3D(isub)   
      
      !****************************************************
      ! Save VTK CRACK file, 2023-07-28. NEWFTU2023072801.
      !****************************************************
      call Save_vtk_file_for_Crack(isub)
          
      !*****************************************************************************
      ! If the time exceeds the scheduled fracturing time for that stage, 
      ! no further fracture propagation calculations will be performed. 2022-07-11.
      !NEWFTU2022071101.
      !*****************************************************************************
      if(c_Time >= c_Stage_Time)then
          goto 880
      endif
      
      !***************************************************************************
      ! If the crack is allowed to extend, then determine whether it will extend.
      ! If an extension occurs, calculate the new crack tip coordinates, 
      ! add new crack segments, and update.
      ! Crack coordinates.
      !***************************************************************************
      ! Determine the direction and size of crack propagation
      call Check_Crack_Grows_3D(1,isub,Yes_Growth)
      ! Smooth treatment of crack front edge, 2021-09-08
      if(Key_Smooth_Front>=1 .and. Key_Smooth_Front/=6)then
           call Crack_Front_Smooth_3D(1,isub)
      endif
      ! Save the maximum principal stress direction vectors of the discrete crack edge nodes
      call Save_Files_Crack_3D_Extra(isub)
      ! If a crack has propagated, update the crack propagation flag: Yes_Last_Growth
      if (any(Yes_Growth).eqv..True.)then
           Yes_Last_Growth = .True.
      ! If no crack has propagated, exit the program.
      else
           Yes_Last_Growth = .False.
           if(isub < SlipWater_Max_Time_Steps_3D)then
               print *,'    ---$---$---$---$---$---$---$---$---$---S---S---S---'
               print *,'    Warning :: No crack propapated, program was ended! '
               print *,'    ---$---$---$---$---$---$---$---$---$---S---S---S---'
           elseif(isub == SlipWater_Max_Time_Steps_3D)then
               print *,'    |<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|'
               print *,'    |    All propagation steps done!  |'
               print *,'    |<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|'
           endif
      endif
      
880   continue
        
      !*************************************************************
      ! Determine the new number of cracks and update the next step
      !num_Crack_Log(isub+1).
      !*************************************************************
      num_Crack_Log(isub+1)  = num_Crack
          
      
      !****************************
      ! Computational Node Stress.
      !****************************
      if(Key_Post_CS_N_Strs==1)then
          ALLOCATE(Stress_xx_Node(num_Node),Stress_yy_Node(num_Node),Stress_zz_Node(num_Node))
          ALLOCATE(Stress_xy_Node(num_Node),Stress_yz_Node(num_Node),Stress_xz_Node(num_Node),Stress_vm_Node(num_Node))           
          ! The precise calculation takes a very long time (issue fixed on 2020-12-27). 
          !1 indicates calculating stress, and the second 1 indicates the Cartesian coordinate system.
          ! The first True, Yes_Add_Insitu, indicates whether in-situ stress is superimposed.
          call Get_Node_Str_XFEM_IN_OUT_3D(.True., isub,DISP,1,1,Stress_xx_Node,Stress_yy_Node,Stress_zz_Node,   &   
                                           Stress_xy_Node,Stress_yz_Node,Stress_xz_Node,Stress_vm_Node)
          ! Save nodal stress.
          call Save_Stress_Node(isub,1)
      endif    

      !**************************************
      ! Calculation Node Strain (2021-09-10)
      !**************************************
      if(Key_Post_CS_N_Stra==1)then
          ALLOCATE(Strain_xx_Node(num_Node),Strain_yy_Node(num_Node),Strain_zz_Node(num_Node))
          ALLOCATE(Strain_xy_Node(num_Node),Strain_yz_Node(num_Node),Strain_xz_Node(num_Node),Strain_vm_Node(num_Node))      
          ! Accurate calculation, extremely time-consuming (issue fixed on 2020-12-27). 
          !2 indicates strain calculation, 1 indicates Cartesian coordinate results
          ! The first True, Yes_Add_Insitu, indicates whether in-situ stress is superimposed.
          call Get_Node_Str_XFEM_IN_OUT_3D(.True., isub, DISP,2,1,Strain_xx_Node,Strain_yy_Node,Strain_zz_Node,  &  
           Strain_xy_Node,Strain_yz_Node,Strain_xz_Node,Strain_vm_Node)
          !Save nodal strain.
          call Save_Strain_Node(isub,1)
      endif    
          
      !****************************
      ! Save VTK file, 2021-07-16.
      !****************************
      call Save_vtk_file(isub)
      
      !*********************************************************************
      ! Check the connectivity of the cracks. 2022-06-10. NEWFTU2022061001.
      ! And update the crack type based on the connectivity relationship.
      !IMPROV2024022202.
      !*********************************************************************
      call D3_Check_Crack_Connections(i_WB,i_Stage,i_Prop,isub)      
      
      !*******************************
      ! Clearing dynamic arrays, etc.
      !*******************************
      if (allocated(globalK)) DEALLOCATE(globalK)
      if(Key_Contact/=0 .and. Key_Contact/=5)then
          if (allocated(Contact_DISP)) DEALLOCATE(Contact_DISP)          
          if (allocated(ori_globalK)) DEALLOCATE(ori_globalK)
      endif          
      if (allocated(freeDOF)) DEALLOCATE(freeDOF)
      if (allocated(fixedDOF)) DEALLOCATE(fixedDOF)
      if (allocated(DISP)) DEALLOCATE(DISP)
      if (allocated(DISP2)) DEALLOCATE(DISP2)
      if (allocated(tem_DISP)) DEALLOCATE(tem_DISP)
      if (allocated(F)) DEALLOCATE(F)
      if (allocated(F_U)) DEALLOCATE(F_U)
      if (allocated(Coupled_Q_3D)) DEALLOCATE(Coupled_Q_3D)
      if (allocated(Coupled_Q_3D_Index)) DEALLOCATE(Coupled_Q_3D_Index)
      if (allocated(delta_x)) DEALLOCATE(delta_x)   
      
      if(Key_Post_CS_N_Strs==1)then
          if (allocated(Stress_xx_Node))DEALLOCATE(Stress_xx_Node)
          if (allocated(Stress_yy_Node))DEALLOCATE(Stress_yy_Node)
          if (allocated(Stress_zz_Node))DEALLOCATE(Stress_zz_Node)
          if (allocated(Stress_xy_Node))DEALLOCATE(Stress_xy_Node)
          if (allocated(Stress_yz_Node))DEALLOCATE(Stress_yz_Node)
          if (allocated(Stress_xz_Node))DEALLOCATE(Stress_xz_Node)
          if (allocated(Stress_vm_Node))DEALLOCATE(Stress_vm_Node)              
      endif
      
      if(Key_Post_CS_N_Stra==1)then
          if (allocated(Strain_xx_Node))DEALLOCATE(Strain_xx_Node)
          if (allocated(Strain_yy_Node))DEALLOCATE(Strain_yy_Node)
          if (allocated(Strain_zz_Node))DEALLOCATE(Strain_zz_Node)
          if (allocated(Strain_xy_Node))DEALLOCATE(Strain_xy_Node)
          if (allocated(Strain_yz_Node))DEALLOCATE(Strain_yz_Node)
          if (allocated(Strain_xz_Node))DEALLOCATE(Strain_xz_Node)
          if (allocated(Strain_vm_Node))DEALLOCATE(Strain_vm_Node)              
      endif
          
      !Clear solver 11 (EBE) related temporary variables.
      if (allocated(storK_XFEM)) DEALLOCATE(storK_XFEM)      
      if (allocated(diag_precon_no_invert)) DEALLOCATE(diag_precon_no_invert)    
      
      !*********************************************************************
      ! If all cracks have stopped expanding and no new cracks have formed, 
      ! terminate the program.
      !*********************************************************************
      if(num_Crack.ne.0)then
          if(Yes_Last_Growth .eqv. .False.)then
              print *,' '
              print *,"    All crack stop propagating!"
              print *,' '
              goto 800
          endif
      endif
      
      !********************************************************
      ! If the time exceeds the fracturing time of the current
      ! stage, exit the current stage fracture step loop.
      !********************************************************
      if(c_Time >= c_Stage_Time)then
          print *,' '
          print *,"    Current stage done!"
          print *,' '
          exit
      endif
      
      !***********************
      ! CPU time consumption.
      !***********************
      call Tool_Get_Current_Time(current_data,date_time,F_time)
      WRITE(*,3001) F_time-S_time,(dble(F_time)-dble(S_time))/Con_60

     
    enddo
   enddo
enddo
      
800 continue
      

RETURN
END SUBROUTINE PhiPsi3D_Static_HF_SlipWater
