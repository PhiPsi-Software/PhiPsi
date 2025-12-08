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
 
SUBROUTINE PhiPsi2D_Static

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
use Global_DISP
use Global_HF
use Global_Stress
use Global_POST
use Global_Contact
use Global_Field_Problem
use Global_Inclusion
use Global_Cross
use omp_lib

!----------------------------------------------------------------------------------
! Read subroutine interface module (activate compiler parameter consistency check)
!----------------------------------------------------------------------------------
use Global_Inter_Matrix_Solve_LSOE
use Global_Inter_Determine_Contact_State_by_Iteration
use Global_Inter_Assemble_Stiffness_Matrix_SPARS_FEM
use Global_Inter_Assemble_Stiffness_Matrix_SPARS_XFEM
use Global_Inter_Matrix_Solve_LSOE_Sparse

!---------------------------
! Variable Type Declaration
!---------------------------
implicit none
integer isub
logical Yes_Last_Growth
real(kind=FT) Lambda
integer,ALLOCATABLE::freeDOF(:)
real(kind=FT),ALLOCATABLE:: globalK(:,:),globalF(:)
real(kind=FT),ALLOCATABLE:: ori_globalK(:,:)
real(kind=FT),ALLOCATABLE:: tem_DISP(:)
real(kind=FT),ALLOCATABLE:: Contact_DISP(:)
! Sparse Matrix K-related (storage format: CSR)
real(kind=FT),ALLOCATABLE:: K_CSR_aa(:)
integer,ALLOCATABLE:: K_CSR_ja(:)
integer,ALLOCATABLE:: K_CSR_ia(:)
integer K_CSR_NNZ_Max
integer K_CSR_NNZ
integer num_FreeD
real(kind=FT) Max_Disp_x,Max_Disp_y,Min_Disp_x,Min_Disp_y
integer i_C
logical Yes_Growth(Max_Num_Cr,2)
integer ifra,Total_Num_G_P
real(kind=FT) Max_vm_Node,Max_vm_Gauss
real(kind=FT) Ger_Min_X_Coor,Ger_Min_Y_Coor
real(kind=FT) Ger_Max_X_Coor,Ger_Max_Y_Coor
integer Result_File_Num
real(kind=FT),ALLOCATABLE::Cracks_CalP_Tan_Aper(:,:)
logical Yes_Generated
character(200) c_File_name_1
integer,ALLOCATABLE::fixedDOF(:)
integer num_FixedD
integer num_Ele_Killed
real(kind=FT),ALLOCATABLE:: Coupled_Q(:,:)
real(kind=FT),ALLOCATABLE:: CalP_Pres(:)
integer num_Count,i_CalP
real(kind=FT) Max_StrX_Node,Min_StrX_Node,Max_StrY_Node,Min_StrY_Node
integer i_E,max_El
real(kind=FT) St,Max_S1
real(kind=FT),ALLOCATABLE::storK(:,:,:),diag_precon(:)
integer,ALLOCATABLE::size_local(:)
integer,ALLOCATABLE::all_local(:,:)
logical Inactive_Key_Multi_TipEnrNode
logical New_Crack_Flag
integer i_User_Crack,i_User_Seg,c_User_Crack,c_User_Seg
integer User_Max_num_Points,tem_isub
integer,ALLOCATABLE::User_Defined_Crack_and_Point_Count(:,:)
#ifdef Silverfrost
real(kind=FT), allocatable:: tem_vector(:)
#endif

!-------------------------------------
! Initial value of temporary variable
!-------------------------------------
Yes_Last_Growth = .False.
New_Crack_Flag = .False.

!----------------------------
! Formatted output statement
!----------------------------
1001 FORMAT('  >> Load step ',I5,' of ',I5,' started:')
1002 FORMAT('     Force factor is ',F5.3)
1021 FORMAT(5X,'Range of displacement x:   ',E14.6,' m to ',E14.6,' m')
1022 FORMAT(5X,'Range of displacement y:   ',E14.6,' m to ',E14.6,' m')
1121 FORMAT(5X,'Range of displacement x:   ',E14.6,' mm to ',E14.5,' mm')
1122 FORMAT(5X,'Range of displacement y:   ',E14.6,' mm to ',E14.5,' mm')
1131 FORMAT(5X,'KI and KII of crack ',I5,' tip 1 are ', E16.7,' and ',E16.7,' MPa.m^(1/2)')
1132 FORMAT(5X,'KI and KII of crack ',I5,' tip 2 are ', E16.7,' and ',E16.7,' MPa.m^(1/2)')
1041 FORMAT(5X,'Max Von-mises stress of all nodes:   ',E15.5,' MPa')
1042 FORMAT(5X,'Max Von-mises stress of Gauss points:   ',E15.5,' MPa')
1051 FORMAT(5X,'Max stress x of all nodes:   ',E15.5,' MPa')
1052 FORMAT(5X,'Min stress x of all nodes:   ',E15.5,' MPa')
1053 FORMAT(5X,'Max stress y of all nodes:   ',E15.5,' MPa')
1054 FORMAT(5X,'Min stress y of all nodes:   ',E15.5,' MPa')
2036 FORMAT(5X,'Total conductivity:   ',E14.5,' μm^2·cm')
2037 FORMAT(5X,'Average conductivity: ',E14.5,' μm^2·cm')
3042 FORMAT(5X,'Max damage factor of Gauss points:   ',F8.4)
4001 FORMAT(5X,'Number of elements to be killed in this step:',I5)
909 FORMAT(5X,'WARNING :: element ',I7,' is broken due to S1>St!')
Num_Frac = Num_Substeps
! If the fracture zone is defined, initial inclusions and cracks will only be generated within the
! fracture zone.
if(Key_Fracture_Zone==1)then
  Ger_Min_X_Coor = Frac_Zone_MinX
  Ger_Min_Y_Coor = Frac_Zone_MinY
  Ger_Max_X_Coor = Frac_Zone_MaxX
  Ger_Max_Y_Coor = Frac_Zone_MaxY
elseif (Key_Fracture_Zone==0)then
  Ger_Min_X_Coor = Min_X_Coor
  Ger_Min_Y_Coor = Min_Y_Coor
  Ger_Max_X_Coor = Max_X_Coor
  Ger_Max_Y_Coor = Max_Y_Coor
endif

!---------------------------------------------------
! Calculate the material matrix D for each material
!---------------------------------------------------
call Get_Material_Matrix

!---------------------------------------------------------------------------------
! If it is to calculate the support crack width with Key_Propped_Width = 1, then:
! Read related result files of hydraulic fracturing analysis
!---------------------------------------------------------------------------------
#ifndef github
if(Key_Propped_Width==1)then
  ! Obtain the file number from hydraulic fracturing analysis
  call Tool_Get_Last_Result_Number(Result_File_Num)
  ! Read the fracture coordinate files (crax files and cray files), and determine the number of
  ! fractures, the number of fracture coordinate points, and the fracture coordinates from them
  call Read_crxo_and_cryo_file_to_Crack_Coor(Result_File_Num)
endif
#endif

!------------------------------------
! Generate initial natural fractures
!------------------------------------
if (Key_Random_NaCr==1) then
  ! The total number of natural fractures is the number of fractures generated randomly.
  num_Na_Crack = num_Rand_Na_Crack
  call Tool_Generate_Natural_Fractures(  &
             Ger_Min_X_Coor,Ger_Max_X_Coor,&
             Ger_Min_Y_Coor,Ger_Max_Y_Coor,&
             NaCr_Length,NaCr_Orientation,   &
             NaCr_Len_Delta,NaCr_Ori_Delta,   &
             num_Na_Crack,NaCr_Length*1.7D0, &
             Crack_Coor(num_Crack+1:num_Crack+num_Na_Crack,1:2,1:2))
  num_Crack = num_Crack + num_Rand_Na_Crack
  Each_Cr_Poi_Num(2:num_Na_Crack+1) =2
endif

!--------------------------
! Generate initial mixture
!--------------------------
! Round inclusion
if (Key_Rand_Circ_Incl==1) then
  ! Randomly generate a circle, then generate an inscribed polygon (only within the rectangular area)
  call Tool_Generate_Circles( &
               Ger_Min_X_Coor,Ger_Max_X_Coor,&
               Ger_Min_Y_Coor,Ger_Max_Y_Coor,&
               Rand_Circ_Incl_R,Rand_Circ_Inc_R_Delta,&
               num_Rand_Circ_Incl,&
               (Rand_Circ_Incl_R+Rand_Circ_Inc_R_Delta)*2.5D0,&
               (Rand_Circ_Incl_R+Rand_Circ_Inc_R_Delta)*0.5D0,&
               Circ_Inclu_Coor(1:num_Rand_Circ_Incl,1:2),&
               Circ_Inclu_Coor(1:num_Rand_Circ_Incl,3))
  num_Circ_Incl = num_Rand_Circ_Incl
  Circ_Inclu_Mat_Num(1:num_Rand_Circ_Incl) =2
endif
! Regular polygon hybrid
if (Key_Rand_Poly_Incl==1) then
  ! Randomly generate a circle, then generate an inscribed polygon (only within the rectangular area)
  call Tool_Generate_Inscribed_Poly(&
                 Ger_Min_X_Coor,Ger_Max_X_Coor,&
                 Ger_Min_Y_Coor,Ger_Max_Y_Coor,&
                 Rand_Poly_Incl_R,Rand_Poly_Inc_R_Delta,&
                 num_Rand_Poly_Incl,num_Vert_Poly_Incl,&
                 (Rand_Poly_Incl_R+Rand_Poly_Inc_R_Delta)*2.5D0,  &
                 (Rand_Poly_Incl_R+Rand_Poly_Inc_R_Delta)*0.5D0, &
                 Poly_Incl_Coor_x(1:num_Rand_Poly_Incl,1:num_Vert_Poly_Incl),&
                 Poly_Incl_Coor_y(1:num_Rand_Poly_Incl,1:num_Vert_Poly_Incl))
  num_Poly_Incl = num_Rand_Poly_Incl
  Poly_Inclu_Edges_Num(1:num_Rand_Poly_Incl)= num_Vert_Poly_Incl
  Poly_Inclu_Mat_Num(1:num_Rand_Poly_Incl) =2
endif
! Irregular polygon inclusion, added on 2017-10-28
if (Key_Rand_Poly_Incl_Irregular ==1) then
  ! Randomly generate a circle, then generate an inscribed polygon (only within the rectangular area)
  num_Rand_Poly_Incl=sum(num_Rand_Poly_Incl_for_Each_Type(1:10))
  call Tool_Generate_Inscribed_Poly_Irregular(&
                 Ger_Min_X_Coor,Ger_Max_X_Coor,&
                 Ger_Min_Y_Coor,Ger_Max_Y_Coor,&
                 num_Rand_Poly_Incl,     &
                 num_Vert_Poly_Incl,      &
                 num_Rand_Poly_Incl_for_Each_Type(1:10),   &
                 Rand_Poly_Incl_R_Min_and_Max(1:10,1:2),   &
                 Rand_Poly_Incl_Irregular_R_Delta_Factor,   &
                 Rand_Poly_Incl_Irregular_Angle_Factor,     &
                 Rand_Poly_Incl_Irregular_Extension_Factor,   &
                 Rand_Poly_Incl_Irregular_Inclination,      &
                 1.2D0,    &
                 0.5D0,    &
                 Poly_Incl_Coor_x(1:num_Rand_Poly_Incl,        &
                                  1:num_Vert_Poly_Incl),   &
                 Poly_Incl_Coor_y(1:num_Rand_Poly_Incl,    &
                                  1:num_Vert_Poly_Incl))
  num_Poly_Incl = num_Rand_Poly_Incl
  Poly_Inclu_Edges_Num(1:num_Rand_Poly_Incl)= num_Vert_Poly_Incl
  Poly_Inclu_Mat_Num(1:num_Rand_Poly_Incl) =2
endif

! Total number of inclusions
num_Inclusion =  num_Circ_Incl + num_Poly_Incl  + num_Ellip_Incl


!------------------------
! Generate initial holes
!------------------------
! Circular hole
if (Key_Random_Hole==1) then
  call Tool_Generate_Circles(&
               Ger_Min_X_Coor,Ger_Max_X_Coor,&
               Ger_Min_Y_Coor,Ger_Max_Y_Coor,&
               Rand_Hole_R,Rand_Hole_R_Delta,&
               num_Rand_Hole,&
              (Rand_Hole_R+Rand_Hole_R_Delta)*2.5D0,&
               (Rand_Hole_R+Rand_Hole_R_Delta)*0.5D0,&
               Hole_Coor(1:num_Rand_Hole,1:2),&
               Hole_Coor(1:num_Rand_Hole,3))
  num_Circ_Hole = num_Rand_Hole
  num_Hole = num_Circ_Hole
else
  num_Hole =  num_Circ_Hole + num_Ellip_Hole
endif



!-------------------------------------------------------------------------------------
! Calculate the Biot consolidation matrix c (Biot_c_MAT), see Smith 5th edition, p.55
!-------------------------------------------------------------------------------------
if(Key_PoreP==1 .or. Key_FS_Seque_Coupling==1)then
  ALLOCATE(Biot_c_MAT(Num_Node*2,Num_Node))
  ALLOCATE(Porous_P(Num_Node))
  call Cal_Biot_c_MAT
endif

!------------------------------------------------------------------------------------------------
! Calculate the traditional degree-of-freedom displacement field under in-situ stress (far-field
! stress), and calculate the in-situ stress level, then save it.
! Node stress and Gauss point stress
! Note: There are no enhanced nodes in the model at this time
!------------------------------------------------------------------------------------------------
if(Key_InSitu_Strategy /=0)then
  if( num_Hole.ne.0 .or. num_Inclusion.ne.0)then
      call Cal_InSitu_Stress_for_Hole_and_Inclusion
  else
      call Cal_InSitu_Stress    
  endif
endif

!---------------------------
! If it is fatigue analysis
!---------------------------
if (Key_Static_Fatigue==1)then
  Num_Substeps = Num_Fatigue_N
endif


!---------------------------------------
! If it is the type of damaged material
!---------------------------------------
if (Key_Damage==1) then
  allocate(Ele_Damage(Num_Elem,Num_Gauss_Points))
  Ele_Damage(1:Num_Elem,1:Num_Gauss_Points) =  ZR
endif

!-------------------------------------------------------------------------------
! If the user customizes the crack propagation path analysis. NEWFTU2024111701.
!-------------------------------------------------------------------------------
if (Key_User_Defined_2D_Crack_Path==1)then
  ! Gradually expand from a crack-free state to all points along the entire predefined crack path.
  Num_Substeps = sum(User_Defined_2D_Crack_Path_Num_Points(1:Num_User_Defined_2D_Crack_Path)-1)
  User_Max_num_Points = maxval(User_Defined_2D_Crack_Path_Num_Points(1:Num_User_Defined_2D_Crack_Path))
  ! The crack number and point number corresponding to each step.
  allocate(User_Defined_Crack_and_Point_Count(Num_Substeps,2))
  !isub_User_Defined_2D_Crack_Path = 0
  tem_isub = 0
  do i_User_Crack =1,Num_User_Defined_2D_Crack_Path
      do i_User_Seg=1,User_Defined_2D_Crack_Path_Num_Points(i_User_Crack)-1
          tem_isub  = tem_isub  + 1
          User_Defined_Crack_and_Point_Count(tem_isub,1) = i_User_Crack
          User_Defined_Crack_and_Point_Count(tem_isub,2) = i_User_Seg
      enddo
  enddo
  ! Do not check whether the crack is expanding.
  Key_Propagation =0
endif

!---------------------------------------
! Cycle of each crack propagation step.
!---------------------------------------
do isub = 1,Num_Substeps
  print *, "  "
  ! Initially, mark it as a non-XFEM analysis
  Yes_XFEM = .False.
  ! Load Factor
  call Force_Factor(Lambda,isub,Yes_Last_Growth)
  WRITE(*,1001) isub,Num_Substeps
  WRITE(*,1002) Lambda
  ! Elements to kill at each step
  if(Key_EKILL==1)then
      num_Ele_Killed = count(Ele_Killed_Each_Load_Step(isub,:) > 0)
      write(*,4001) num_Ele_Killed
      ! Save killed elements
      call Save_Killed_Elements(isub)
  endif
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! If the user customizes the crack propagation path analysis. NEWFTU2024111701.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (Key_User_Defined_2D_Crack_Path==1)then
      ! Determine which crack you are currently at and the coordinate point of the crack based on isub.
      c_User_Crack = User_Defined_Crack_and_Point_Count(isub,1)
      c_User_Seg   = User_Defined_Crack_and_Point_Count(isub,2)
      ! Assign coordinates to the cracks involved in the calculation.
      num_Crack = c_User_Crack
      do i_User_Crack = 1,c_User_Crack
        if(i_User_Crack < c_User_Crack) then
            Each_Cr_Poi_Num(i_User_Crack) = User_Defined_2D_Crack_Path_Num_Points(i_User_Crack)
            Crack_Coor(i_User_Crack,1:Each_Cr_Poi_Num(i_User_Crack),1) = &
                User_Defined_2D_Crack_Path(i_User_Crack,1:2*Each_Cr_Poi_Num(i_User_Crack):2)
            Crack_Coor(i_User_Crack,1:Each_Cr_Poi_Num(i_User_Crack),2) = &
                User_Defined_2D_Crack_Path(i_User_Crack,2:2*Each_Cr_Poi_Num(i_User_Crack):2)
        else
            Each_Cr_Poi_Num(i_User_Crack)  = c_User_Seg + 1
            Crack_Coor(i_User_Crack,1:Each_Cr_Poi_Num(i_User_Crack),1) = &
                User_Defined_2D_Crack_Path(i_User_Crack,1:2*Each_Cr_Poi_Num(i_User_Crack):2)
            Crack_Coor(i_User_Crack,1:Each_Cr_Poi_Num(i_User_Crack),2) = &
                User_Defined_2D_Crack_Path(i_User_Crack,2:2*Each_Cr_Poi_Num(i_User_Crack):2)
        endif
      enddo
  endif

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! If it is an indirect fluid-structure interaction calculation
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(Key_FS_Seque_Coupling==1)then
#ifndef github
      print *,'    '
      print *,'    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      print *,'    STEP 1: SOLVING THE FIELD PROBLEM      '
      print *,'    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call PhiPsi2D_Static_Field_Problem_Step(isub,Num_Substeps)   
      print *,'    '
      print *,'    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      print *,'    STEP 2: SOLVING THE SOLID PROBLEM      '
      print *,'    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
#endif  
  endif

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! If it is an indirect thermosetting coupling calculation       
  ! Date: 2021-11-03               
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(Key_TS_Seque_Coupling==1)then
      print *,'    '
      print *,'    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      print *,'    STEP 1: SOLVING THE THERMO FIELD PROBLEM  '
      print *,'    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      !call PhiPsi2D_Static_Field_Problem_Step(isub,Num_Substeps)
      print *,'    '
      print *,'    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      print *,'    STEP 2: SOLVING THE SOLID PROBLEM      '
      print *,'    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
  endif
  

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! If there are cracks, holes, or inclusions (XFEM) 
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(num_Crack.ne.0 .or. num_Hole.ne.0 .or. num_Inclusion.ne.0 )then
      ! Marked as XFEM analysis, not FEM analysis
      Yes_XFEM = .True.
      !***************************
      ! Confirm enhancement node.
      !***************************
      ifra = isub
      if (Key_Local_Mesh_Refine == 0)then
          call Determine_Enriched_Nodes(ifra,isub)
      endif

      !***************************************************
      ! Enhanced element Local Optimization (2021-06-27).
      ! Attention: This operation will change the total
      ! number of elements, nodes, and coordinates of
      ! added nodes. The enriched scheme is referred to
      ! my paper.
      !***************************************************
#ifndef github
      if (Key_Local_Mesh_Refine > 0)then
          Flag_Local_Refined =.False.
          if (isub>=2)then
              ! Reload the original geometric model file, otherwise the encryption will be applied on top of the
              ! previously encrypted mesh, resulting in the mesh becoming increasingly dense.
              call Read_Geo
          endif
          !Determine enriched nodes and elements for the first time.
          if (Key_Multi_TipEnrNode   ==  1) then
              Key_Multi_TipEnrNode   =  0
              Inactive_Key_Multi_TipEnrNode = .True.
          endif
          call Determine_Enriched_Nodes(ifra,isub)
          print *,'    Local mesh refinement...'
          call Refine_Enriched_Elements(ifra,isub)
          Flag_Local_Refined = .True.
          !Determine enriched nodes and elements again.
          if(Inactive_Key_Multi_TipEnrNode)then
              Key_Multi_TipEnrNode          =  1
              Inactive_Key_Multi_TipEnrNode = .False.
          endif
          call Determine_Enriched_Nodes(ifra,isub)
      endif
#endif
      !***************************************
      ! Assign a number to the enhanced node.
      !***************************************
      call Number_Enriched_Nodes(isub)
      print *,'    Total_FD:          ',Total_FD
      print *,'    Enrich_Freedom:    ',Enrich_Freedom
      print *,'    H enriched nodes:  ',n_h_Node
      print *,'    Tip enriched nodes:',n_t_Node
      
      !********************************************************************
      ! Save crack-related files (including enhanced node numbering c_POS)
      !********************************************************************
      call Save_Files_Crack(isub)
      ! If there are holes, save the hole coordinates and other related files.
      if(num_Hole.ne.0)then
          call Save_Files_Holes(isub)
      endif
      ! If it contains Cross, save the related files.
      if(num_Cross.ne.0)then
          call Save_Files_Cross(isub)
      endif
      ! If inclusions are present, save the coordinates of the inclusions
      if(num_Inclusion.ne.0)then
          call Save_Files_Inclusions(isub)
      endif
      
      !**********************
      ! Save VTK CRACK file.
      !**********************
      call Save_vtk_file_for_Crack(isub)
          
      !****************************************************
      ! Calculate data related to crack calculation points
      !****************************************************
      !call Cal_Crack_Points_Info(isub)
      ! The following processing (related to water pressure) is aimed at obtaining 
      ! calculation point coordinates consistent with hydraulic fracturing analysis, 
      ! in order to examine the distribution of fracture openings.
      if(num_Crack.ne.0)then
          Cracks_HF_State(1:Max_Num_Cr) = 1
          !call Stat_Crack_Connection(isub)
          
          call Stat_Crack_Connection(isub)
          
          call Cal_HF_Crack_Points_Info_Linear(isub)
1104      FORMAT(5X,'Total DOFs of fluid:',I8)
          if(Key_Crack_Inner_Pressure==1) then
              WRITE(*,1104) num_Tol_CalP_Water
          endif
      endif

      !******************************
      ! Consider boundary conditions
      !******************************
      ALLOCATE(freeDOF(Total_FD))
      ALLOCATE(fixedDOF(Total_FD))
      call Boundary_Cond_output_Fixed_Freedom(Total_FD,isub,freeDOF,num_FreeD,fixedDOF,num_FixedD)

      !********************************************
      ! Allocate other relevant dynamic data space
      !********************************************
      if(num_Crack.ne.0)then
          ! ALLOCATE(DISP(Total_FD)) ! The displacement obtained in the current step, including active and
          ! fixed degrees of freedom
          ALLOCATE(Coupled_Q(Total_FD,num_Tol_CalP_Water))
          ! ALLOCATE(F(Total_FD))                     !External load vector considering fluid pressure
          ALLOCATE(CalP_Pres(num_Tol_CalP_Water))
          ! ALLOCATE(delta_x(num_FreeD)) !Newton Raphson iteration displacement and pressure increment, note
          ! the dimensions: num_freeDOF num_free_CalP
      endif
      !*********************************************************************************
      ! Calculate the coupling matrix Q (if considering fluid pressure within the seam)
      !*********************************************************************************
      if(Key_Crack_Inner_Pressure==1) then
          call Cal_HF_Matrix_Q_Linear(isub,Coupled_Q,Total_FD)
      endif

      !************************************************************************
      ! Convert water pressure into the pressure vector at all
      ! calculation points (if considering fluid pressure within the fracture)
      !************************************************************************
      if(Key_Crack_Inner_Pressure==1) then
          num_Count = 0
          do i_C = 1,num_Crack
              do i_CalP=1,Cracks_CalP_Num(i_C)
                  num_Count = num_Count + 1
                  CalP_Pres(num_Count) = Crack_Pressure(i_C)
                  Cracks_CalP_Pres(i_C,i_CalP)=Crack_Pressure(i_C)
              end do
          end do
      endif

      !*************************************************************
      ! If it is to calculate the crack opening width with support,
      ! Key_Propped_Width=1, then:
      !*************************************************************
      ! Read the wpnp (Width of Proppant No Pressure) file, which contains
      ! Crack opening of the support crack under an uncompressed state (obtained from fluid-solid coupling
      ! analysis)
#ifndef github
      if(Key_Propped_Width==1)then
          ! Read the corresponding wpnp file into the Cracks_CalP_wpnp variable
          call Read_Cr_CalP_Wpnp_to(Result_File_Num,Cracks_CalP_wpnp,Max_Num_Cr,Max_Num_Cr_CalP)
          ! Read the crack opening calculated under water pressure
          call Read_Cr_CalP_Aper_to(Result_File_Num,Cracks_CalP_wpor,Max_Num_Cr,Max_Num_Cr_CalP)
      endif
#endif

      !*************
      ! Load vector
      !*************
      ALLOCATE(globalF(Total_FD))
      call Force_Vector(Total_FD,isub,.False.,Lambda,globalF)
      ! If there is water pressure within the joint
      if(Key_Crack_Inner_Pressure==1)then
          globalF(freeDOF(1:num_FreeD))=globalF(freeDOF(1:num_FreeD))+ &
                          MATMUL(Coupled_Q(freeDOF(1:num_FreeD),&
                                 1:num_Tol_CalP_Water),CalP_Pres(1:num_Tol_CalP_Water))
      endif
      print *,'    Sum of globalF:',sum(globalF)

      !****************************************************
      !       If solver 11: Element by Element, diagonally
      !       preconditioned conjugate gradient solver
      !****************************************************
      if(Key_SLOE==11)then
          ALLOCATE(DISP(Total_FD))
          ALLOCATE(storK(MDOF_2D,MDOF_2D,Num_Elem))
          ALLOCATE(diag_precon(0:num_FreeD))
          ALLOCATE(size_local(Num_Elem))
          ALLOCATE(all_local(MDOF_2D,Num_Elem))
          call Ele_by_Ele_XFEM_PCG(isub,Lambda,cg_tol,max_cg, &
                      num_FreeD,freeDOF(1:num_FreeD),&
                      globalF(freeDOF(1:num_FreeD)),DISP,&
                      Total_Num_G_P,&
                      storK,size_local,all_local,&
                      diag_precon(0:num_FreeD))
          goto 555
      endif

      !*****************************************************************************
      ! Assemble the stiffness matrix and solve for displacements
      !---------------------------------------------------------
      ! Depending on whether the embedding of the crack surface is considered, there 
      ! are two situations. If embedding of the crack surface is allowed,
      ! Then it can be calculated directly; if not allowed, iterative calculation 
      ! is required.
      !*****************************************************************************
      !++++++++++++++++++++++++++++++++++++++
      ! Allow embedding of the crack surface
      ! (ignoring crack surface contact).
      !++++++++++++++++++++++++++++++++++++++
      if(Key_Contact==0) then
          ! Assembly stiffness matrix
          if(Key_K_Sparse==0)then
              ALLOCATE(globalK(Total_FD,Total_FD))
              print *,'    Assembling K...'
              call Assemble_Stiffness_Matrix_XFEM(isub,globalK,Total_FD,Total_Num_G_P)
              print *,'    Sum of globalK: ',  sum(globalK)
          elseif(Key_K_Sparse==1)then
#ifdef gfortran
#ifdef notcbfortran
#ifndef github
              ! K_CSR_NNZ_Max = CEILING(dble(Total_FD)**2*Sparse_Ratio)  ! Maximum number of non-zero elements
              print *,'    Get max NNZ before assembling K......'
              call Assemble_Stiffness_Matrix_SPARS_XFEM_Get_MaxNNZ(isub,&
                          freeDOF(1:num_FreeD),num_FreeD,Total_FD,K_CSR_NNZ_Max)
              print *,'    K_CSR_NNZ_Max:',K_CSR_NNZ_Max    
          
              if (allocated(K_CSR_aa)) DEALLOCATE(K_CSR_aa)
              if (allocated(K_CSR_ja)) DEALLOCATE(K_CSR_ja)
              if (allocated(K_CSR_ia)) DEALLOCATE(K_CSR_ia)
              ALLOCATE(K_CSR_aa(K_CSR_NNZ_Max))
              ALLOCATE(K_CSR_ja(K_CSR_NNZ_Max))
              ALLOCATE(K_CSR_ia(num_FreeD+1))
              print *,'    Assemble the sparse K...'
              call Assemble_Stiffness_Matrix_SPARS_XFEM(isub, &
                   freeDOF(1:num_FreeD),num_FreeD,&
                   fixedDOF(1:num_FixedD),num_FixedD,     &
                   K_CSR_aa,K_CSR_ja,K_CSR_ia,&
                   K_CSR_NNZ_Max,K_CSR_NNZ,&
                   Total_FD,Total_Num_G_P)
#endif
#endif
#endif
          endif
          ! Solve for displacement
          ALLOCATE(DISP(Total_FD))
          DISP(1:Total_FD) = ZR
          ALLOCATE(tem_DISP(num_FreeD))
          print *,'    Solving displacements...'
          if(Key_K_Sparse==0)then
              call Matrix_Solve_LSOE(0,1,Key_SLOE,&
                         globalK(freeDOF(1:num_FreeD),freeDOF(1:num_FreeD)),&
                         globalF(freeDOF(1:num_FreeD)), &
                         tem_DISP,num_FreeD)
          elseif(Key_K_Sparse==1)then
#ifdef gfortran
#ifdef notcbfortran
#ifndef github
              call Matrix_Solve_LSOE_Sparse(0,1,Key_SLOE,&
                        K_CSR_NNZ,&
                        K_CSR_aa(1:K_CSR_NNZ),&
                        K_CSR_ja(1:K_CSR_NNZ),&
                        K_CSR_ia(1:num_FreeD+1),&
                        globalF(freeDOF(1:num_FreeD)),&
                        tem_DISP,num_FreeD)
#endif
#endif
#endif
          endif

          DISP(freeDOF(1:num_FreeD)) = tem_DISP
      !+++++++++++++++++++++++++++++++++++++++++++++++++
      ! Embedding of cracked surfaces is not allowed!!!
      !+++++++++++++++++++++++++++++++++++++++++++++++++
      elseif(Key_Contact /= 0) then
          if(Key_K_Sparse==0)then
            ALLOCATE(globalK(Total_FD,Total_FD))
            call Assemble_Stiffness_Matrix_XFEM(isub,globalK,Total_FD,Total_Num_G_P)
            ! ------Used for contacting iterative ori_globalK
            ALLOCATE(ori_globalK(Total_FD,Total_FD))
            ori_globalK = globalK
          elseif(Key_K_Sparse==1)then
#ifdef gfortran
#ifdef notcbfortran
#ifndef github
            ! K_CSR_NNZ_Max = CEILING(dble(Total_FD)**2*Sparse_Ratio)  ! Maximum number of non-zero elements
            print *,'    Get max NNZ before assembling K......'
            call Assemble_Stiffness_Matrix_SPARS_XFEM_Get_MaxNNZ(isub,&
                          freeDOF(1:num_FreeD),num_FreeD,Total_FD,K_CSR_NNZ_Max)
            print *,'    K_CSR_NNZ_Max:',K_CSR_NNZ_Max   
              
            if (allocated(K_CSR_aa)) DEALLOCATE(K_CSR_aa)
            if (allocated(K_CSR_ja)) DEALLOCATE(K_CSR_ja)
            if (allocated(K_CSR_ia)) DEALLOCATE(K_CSR_ia)
            ALLOCATE(K_CSR_aa(K_CSR_NNZ_Max))
            ALLOCATE(K_CSR_ja(K_CSR_NNZ_Max))
            ALLOCATE(K_CSR_ia(num_FreeD+1))
            print *,'    Assembling the sparse K......'
            call Assemble_Stiffness_Matrix_SPARS_XFEM(isub,&
                   freeDOF(1:num_FreeD),num_FreeD,&
                   fixedDOF(1:num_FixedD),num_FixedD,  &
                   K_CSR_aa,K_CSR_ja,K_CSR_ia,&
                   K_CSR_NNZ_Max,&
                   K_CSR_NNZ,Total_FD,Total_Num_G_P)
#endif
#endif
#endif
          endif

          ALLOCATE(DISP(Total_FD))
          ALLOCATE(Contact_DISP(Total_FD))
          ALLOCATE(tem_DISP(num_FreeD))
          print *,'    Solving displacements......'
          if(Key_K_Sparse==0)then
              call Matrix_Solve_LSOE(0,1,Key_SLOE,&
                         globalK(freeDOF(1:num_FreeD),&
                         freeDOF(1:num_FreeD)),&
                         globalF(freeDOF(1:num_FreeD)),&
                         tem_DISP,num_FreeD)
          elseif(Key_K_Sparse==1)then
#ifdef gfortran
#ifdef notcbfortran
#ifndef github
              call Matrix_Solve_LSOE_Sparse(0,1,Key_SLOE,&
                        K_CSR_NNZ,&
                        K_CSR_aa(1:K_CSR_NNZ),&
                        K_CSR_ja(1:K_CSR_NNZ),&
                        K_CSR_ia(1:num_FreeD+1),&
                        globalF(freeDOF(1:num_FreeD)),&
                        tem_DISP,num_FreeD)
#endif
#endif
#endif
          endif
          !************************************************************************
          ! Contact iteration, calculate the displacement of the contact iteration
          !************************************************************************
          Contact_DISP(1:Total_FD) = ZR
          Contact_DISP(freeDOF(1:num_FreeD))=tem_DISP
          if(Key_K_Sparse==0)then
              call Determine_Contact_State_by_Iteration(&
                         isub,isub,isub,&
                  Contact_DISP,Total_FD, &
                  Usual_Freedom,Enrich_Freedom,  &
                  freeDOF,num_FreeD,globalF,ori_globalK,globalK)
          elseif(Key_K_Sparse==1)then
              !to be done.
          endif
          ! Return the displacement calculated after the new contact consideration
          DISP = Contact_DISP
      endif

555   continue

      !====================================================
      !Element by Element Contact Calculation, 2020-03-31.
      !====================================================
      if(Key_SLOE==11 .and. Key_Contact /= 0) then
          call EBE_Determine_Contact_State_by_Iteration(isub,&
                  cg_tol,max_cg,     &
                  num_FreeD,freeDOF(1:num_FreeD),&
                  globalF,DISP,    &
                  storK,size_local,all_local,&
                  diag_precon(0:num_FreeD))
      endif
      
#ifndef Silverfrost
        Max_Disp_x = maxval(DISP(1:2*Num_Node:2))
        Min_Disp_x = minval(DISP(1:2*Num_Node:2))
        Max_Disp_y = maxval(DISP(2:2*Num_Node:2))
        Min_Disp_y = minval(DISP(2:2*Num_Node:2))
#endif   
     
#ifdef Silverfrost
        allocate(tem_vector(Num_Node))
        tem_vector = DISP(1:2*Num_Node:2)
        Max_Disp_x = maxval(tem_vector)
        tem_vector = DISP(1:2*Num_Node:2)
        Min_Disp_x = minval(tem_vector)
        
        tem_vector = DISP(2:2*Num_Node:2)
        Max_Disp_y = maxval(tem_vector)
        tem_vector = DISP(2:2*Num_Node:2)
        Min_Disp_y = minval(tem_vector)
        deallocate(tem_vector)
#endif

      WRITE(*,1021) Min_Disp_x,Max_Disp_x
      WRITE(*,1022) Min_Disp_y,Max_Disp_y
      ! Calculate crack width.
      if(num_Crack.ne.0)then
          call Cal_Crack_Aperture(isub,DISP)
#ifndef Silverfrost
          print *,'    Max aperture of crack 1 (mm):',&
            maxval(Cracks_CalP_Aper(1,1:Cracks_CalP_Num(1)))*1000.0D0
          print *,'    Min aperture of crack 1 (mm):',&
            minval(Cracks_CalP_Aper(1,1:Cracks_CalP_Num(1)))*1000.0D0
#endif 
      endif
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%                                                                 %
  !%                                                                 %
  !%                                                                 %
  !%      If there are no cracks, holes, or inclusions (FEM)         %
  !%                                                                 %
  !%                                                                 %
  !%                                                                 %
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif (num_Crack.eq.0 .and. num_Hole.eq.0 .and. num_Inclusion.eq.0 )then
      !**************************
      ! Total degrees of freedom
      !**************************
      Total_FD = 2*Num_Node
      print *,'    Total_FD:',Total_FD
      !*****************
      ! Cluster payload
      !*****************
      ALLOCATE(globalF(Total_FD))
      call Force_Vector(Total_FD,isub,.False.,Lambda,globalF)
      print *,'    Sum of globalF:',sum(globalF)

      !******************************
      ! Consider boundary conditions
      !******************************
      ALLOCATE(freeDOF(Total_FD))
      ALLOCATE(fixedDOF(Total_FD))
      call Boundary_Cond(Total_FD,isub,freeDOF,num_FreeD)

      !****************************************************
      !       If solver 11: Element by Element, diagonally
      !       preconditioned conjugate gradient solver
      !****************************************************
      if(Key_SLOE==11)then
          ALLOCATE(DISP(Total_FD))
          call Ele_by_Ele_FEM_PCG(isub,Lambda,cg_tol,max_cg,   &
                              num_FreeD,freeDOF(1:num_FreeD),&
                              globalF(freeDOF(1:num_FreeD)),DISP,&
                              Total_Num_G_P)
          goto 666
      endif

      !***************************
      ! Assembly stiffness matrix
      !***************************
      if(Key_K_Sparse==0)then
          print *,'    Assembling K......'
          ALLOCATE(globalK(Total_FD,Total_FD))
          call Assemble_Stiffness_Matrix_FEM(isub,globalK,Total_FD,Total_Num_G_P)
          print *,'    Sum of globalK: ',  sum(globalK)
      elseif(Key_K_Sparse==1)then
#ifdef gfortran
#ifdef notcbfortran
#ifndef github
          ! K_CSR_NNZ_Max = CEILING(Total_FD**2*Sparse_Ratio)      !Maximum number of non-zero elements
          print *,'    Get max NNZ before assembling K......'
          call Assemble_Stiffness_Matrix_SPARS_FEM_Get_MaxNNZ(isub,&
                          freeDOF(1:num_FreeD),num_FreeD,Total_FD,K_CSR_NNZ_Max)
          print *,'    K_CSR_NNZ_Max:',K_CSR_NNZ_Max    
              
          if (allocated(K_CSR_aa)) DEALLOCATE(K_CSR_aa)
          if (allocated(K_CSR_ja)) DEALLOCATE(K_CSR_ja)
          if (allocated(K_CSR_ia)) DEALLOCATE(K_CSR_ia)
          ALLOCATE(K_CSR_aa(K_CSR_NNZ_Max))
          ALLOCATE(K_CSR_ja(K_CSR_NNZ_Max))
          ALLOCATE(K_CSR_ia(num_FreeD+1))
          ! Assemble the sparse stiffness matrix (automatically considers 
          ! boundary conditions and removes sparse matrix elements corresponding
          ! to constrained degrees of freedom)
          print *,'    Assembling the sparse K......'
          call Assemble_Stiffness_Matrix_SPARS_FEM(isub,&
                  freeDOF(1:num_FreeD),num_FreeD,&
                  K_CSR_aa,K_CSR_ja,K_CSR_ia,&
                  K_CSR_NNZ_Max,K_CSR_NNZ,Total_FD,Total_Num_G_P)
#endif
#endif
#endif
      endif

      !************************
      ! Solve for displacement
      !************************
      ALLOCATE(DISP(Total_FD))
      DISP(1:Total_FD) = ZR
      ALLOCATE(tem_DISP(num_FreeD))
      print *,'    Solving displacements......'
      if(Key_K_Sparse==0)then
          call Matrix_Solve_LSOE(0,1,Key_SLOE,&
                     globalK(freeDOF(1:num_FreeD),&
                     freeDOF(1:num_FreeD)),&
                     globalF(freeDOF(1:num_FreeD)),&
                     tem_DISP,num_FreeD)
      elseif(Key_K_Sparse==1)then
#ifdef gfortran
#ifdef notcbfortran
#ifndef github
          call Matrix_Solve_LSOE_Sparse(0,1,Key_SLOE,&
                    K_CSR_NNZ,&
                    K_CSR_aa(1:K_CSR_NNZ),&
                    K_CSR_ja(1:K_CSR_NNZ),&
                    K_CSR_ia(1:num_FreeD+1),&
                    globalF(freeDOF(1:num_FreeD)),&
                    tem_DISP,num_FreeD)
#endif
#endif
#endif
      endif
      DISP(freeDOF(1:num_FreeD)) = tem_DISP

666   continue

      Max_Disp_x = maxval(DISP(1:Total_FD:2))
      Min_Disp_x = minval(DISP(1:Total_FD:2))
      Max_Disp_y = maxval(DISP(2:Total_FD:2))
      Min_Disp_y = minval(DISP(2:Total_FD:2))
      if(Key_Unit_System==1)then
          WRITE(*,1021) Min_Disp_x,Max_Disp_x
          WRITE(*,1022) Min_Disp_y,Max_Disp_y
      elseif(Key_Unit_System==2)then
          WRITE(*,1121) Min_Disp_x,Max_Disp_x
          WRITE(*,1122) Min_Disp_y,Max_Disp_y
      endif
  end if  
  

  !*******************
  ! Save displacement
  !*******************
  call Save_Disp(isub,1)
          
  !*******************************************************
  ! If the contact condition of the crack surface is
  ! considered, save the element contact state elcs file.
  !*******************************************************
  if((num_Crack.ne.0) .and. (Key_Contact/=0)) then
      call Save_Ele_Contac_State(isub)
  endif

  !******************************************************************************************
  ! Save information related to crack width, including the coordinates of calculation points
  !******************************************************************************************
  if(num_Crack.ne.0) then
      call Save_Files_Cr_CalP(isub)
  endif

  !*********************************************************************************************
  ! If it is for calculating the support crack width, Key_Propped_Width = 1, then:
  !-----------
  !2016-08-27
  !*********************************************************************************************
  ! Calculate the flow-conductivity coefficient C_f at each calculation point of each fracture and
  ! save it to a *.cond file, and also calculate the system's total flow-conductivity
  ! Total_Conductivity
#ifndef github
  if(Key_Propped_Width==1)then
      call Cal_Propped_Cr_CalP_Conductivity(ifra,isub)
      call Save_Files_Cr_CalP_Propp_Conductivity(ifra,isub)
      print *,'    ++++++++++++++++++++++++++++++++++++++++++++'
      write(*,2036) Total_Conductivity*1.0D14
      write(*,2037) Ave_Conductivity*1.0D14
      print *,'    ++++++++++++++++++++++++++++++++++++++++++++'
  endif
#endif  

  !................................................................................................
  ! The following post-processing calculations are executed on different cores because they do not
  ! affect each other.
  ! Note: Later gave up
  !................................................................................................
  !-------------------------------------------
  !  If needed, calculate and save the stress
  !-------------------------------------------
  if (Key_Post_CS_N_Strs ==1) then
      ! Allocate memory space for node stress
      ALLOCATE(Stress_xx_Node(num_Node))
      ALLOCATE(Stress_yy_Node(num_Node))
      ALLOCATE(Stress_xy_Node(num_Node))
      ALLOCATE(Stress_vm_Node(num_Node))
      ! If thermal stress is considered, allocate memory for nodal thermal stress.
      if(Key_Thermal_Stress==1)then
          ALLOCATE(TStress_xx_Node(num_Node))
          ALLOCATE(TStress_yy_Node(num_Node))
          ALLOCATE(TStress_xy_Node(num_Node))
          ALLOCATE(TStress_vm_Node(num_Node))
      endif
      ! Calculate the nodal stress and store it in the global variables Stress_xx_Node, Stress_yy_Node,
      ! Stress_xy_Node, Stress_vm_Node
      if(Yes_XFEM.eqv..True.)then
          call Get_Node_Stress_XFEM(isub,DISP)
      else
          call Get_Node_Stress_FEM(isub)
      endif
      ! Maximum Mises stress at the screen output node
      Max_vm_Node   = maxval(Stress_vm_Node(1:num_Node))/1.0D6
      Max_StrX_Node = maxval(Stress_xx_Node(1:num_Node))/1.0D6
      Min_StrX_Node = minval(Stress_xx_Node(1:num_Node))/1.0D6
      Max_StrY_Node = maxval(Stress_yy_Node(1:num_Node))/1.0D6
      Min_StrY_Node = minval(Stress_yy_Node(1:num_Node))/1.0D6
      WRITE(*,1041) Max_vm_Node
      WRITE(*,1051) Max_StrX_Node
      WRITE(*,1052) Min_StrX_Node
      WRITE(*,1053) Max_StrY_Node
      WRITE(*,1054) Min_StrY_Node
      ! Save Node Stress
      call Save_Stress_Node(isub,1)
  end if

  !--------------------------------------------------------------------------------
  ! If needed, calculate and save the Gauss point coordinates for post-processing.
  !--------------------------------------------------------------------------------
  if (Key_Post_CS_G_Coor==1) then
      ALLOCATE(Gauss_CoorX(Total_Num_G_P))
      ALLOCATE(Gauss_CoorY(Total_Num_G_P))
      ! Obtain Gauss point coordinates (Note: In fact, the Gauss point coordinates have already been
      ! indirectly calculated during the assembly of the stiffness matrix, this calculation is done again
      ! here)
      call Cal_Gauss_Coors(isub,Total_Num_G_P,Gauss_CoorX,Gauss_CoorY)
      ! Save Gauss point coordinates
      call Save_Gauss_Coors(isub,Total_Num_G_P,Gauss_CoorX,Gauss_CoorY)
  end if
  

  

  !**********************************************************************************
  ! If needed, calculate and save the Gauss point displacements for post-processing.
  !**********************************************************************************
  if (Key_Post_CS_G_Disp==1) then
      ALLOCATE(DISP_x_Gauss(Total_Num_G_P))
      ALLOCATE(DISP_y_Gauss(Total_Num_G_P))
      ! Obtain Gauss point displacement
      call Get_Gauss_Disps(isub,DISP,Total_Num_G_P)
      ! Save Gauss point displacement
      call Save_Gauss_Disps(isub,Total_Num_G_P)
  end if


  !***********************************************************************************
  ! If needed, calculate and save the Gauss point stress for post-processing display.
  !***********************************************************************************
  if (Key_Post_CS_G_Strs==1) then
      ALLOCATE(Stress_xx_Gauss(Total_Num_G_P))
      ALLOCATE(Stress_yy_Gauss(Total_Num_G_P))
      ALLOCATE(Stress_xy_Gauss(Total_Num_G_P))
      ALLOCATE(Stress_vm_Gauss(Total_Num_G_P))
      ! If thermal stress is considered, allocate memory for thermal stress at Gauss points
      if(Key_Thermal_Stress==1)then
          ALLOCATE(TStress_xx_Gauss(Total_Num_G_P))
          ALLOCATE(TStress_yy_Gauss(Total_Num_G_P))
          ALLOCATE(TStress_xy_Gauss(Total_Num_G_P))
          ALLOCATE(TStress_vm_Gauss(Total_Num_G_P))
      endif
      ! Calculate Gauss point stress
      if(Yes_XFEM.eqv..True.)then
          call Get_Gauss_Stress_XFEM(isub,DISP)
      else
          call Get_Gauss_Stress_FEM(isub)
      endif
      ! Display the maximum Mises stress at Gauss points on the screen
      Max_vm_Gauss=maxval(Stress_vm_Gauss(1:Total_Num_G_P))/1.0D6
      WRITE(*,1042) Max_vm_Gauss
      ! Save Gauss point stress
      call Save_Gauss_Stress(isub,Total_Num_G_P)
  end if

  !*********************************************************
  ! If necessary, look for the damaged elements, 2020-02-28
  ! Only one element can be destroyed at a time
  !*********************************************************
  if(Key_Element_Break==1) then
      if(Key_Element_Break_Rule ==1)then
          Max_S1 =1.0D-20
          St = Material_Para(Break_Mat_Num,5)
          do i_E = 1,Num_Elem
              if (Elem_Mat(i_E)==Break_Mat_Num) then
                  if(Elem_Ave_Gauss_S1(i_E)>=Max_S1)then
                      Max_S1 = Elem_Ave_Gauss_S1(i_E)
                      Max_El = i_E
                  endif
              endif
          enddo
          print *,'    Max s1 (MPa) of elements:',Max_S1/1.0D6
          if(Max_S1>=St)then
              !Elem_Break(Max_El) =.True.
              write(*,909) i_E
          endif
      endif
      ! Save killed (destroyed) elements
      call Save_Killed_Elements(isub)
  endif
  

  !----------------------------------------------------
  ! If it contains damaged material, material type = 3
  !----------------------------------------------------
  if(Key_Damage ==1)then
#ifndef github
      ! Obtain the damage variable values and D matrix at Gauss points, and 
      ! store them in Ele_Damage(Num_Elem,4) and Ele_Damage_D(Num_Elem,3,3)
      allocate(Damage_Gauss(Total_Num_G_P))
      if(Yes_XFEM.eqv..True.)then
          call Get_Gauss_Damage_XFEM(isub,DISP)
      else
          call Get_Gauss_Damage_FEM(isub)
      endif
      ! Screen output of the maximum damage factor at Gauss points
      WRITE(*,3042) maxval(Damage_Gauss(1:Total_Num_G_P))
      ! Save Gauss point damage factor
      call Save_Gauss_Damage(isub,Total_Num_G_P)
#endif
  endif


  !*******************************************************
  ! If the contact condition of the crack surface is 
  ! considered, save the element contact state elcs file.
  !*******************************************************
  if((num_Crack.ne.0) .and. (Key_Contact/=0)) then
      call Save_Ele_Contac_State(isub)
  end if

  !***************************************************************************
  ! If the crack is allowed to propagate and is under maximum circumferential 
  ! tensile stress.
  ! According to the guidelines, calculate the stress intensity factor.
  !***************************************************************************
  if(num_Crack.ne.0)then
      if ((Key_Propagation ==1) .and. (CFCP == 1)) then
          if (Key_SIFs_Method==1) then
              call Cal_SIFs_DIM(isub,DISP)
          elseif (Key_SIFs_Method==2) then
              call Cal_SIFs_IIM(isub,.True.,DISP)
          end if
          do i_C = 1,num_Crack
              WRITE(*,1131)i_C,KI(i_C,1)/1.0D6,KII(i_C,1)/1.0D6
              WRITE(*,1132)i_C,KI(i_C,2)/1.0D6,KII(i_C,2)/1.0D6
          end do
      end if
      call Save_SIFs_KI_and_KII(isub)
  endif

  !***************************
  ! Save VTK file, 2021-07-16
  !***************************
  call Save_vtk_file(isub)

  ! Save the number of iterations for each step
  call Save_HF_time(1,isub,isub,0.0D0)
  !--------------------------------------------------------------------------------
  ! If the crack is allowed to extend, then determine whether it extends.
  ! If an extension occurs, calculate the new crack tip coordinates, add new crack
  ! segments, and update.
  ! Crack coordinates
  !--------------------------------------------------------------------------------
  if(num_Crack.ne.0)then
    if (Key_Propagation ==1) then
      ! If it is a fatigue analysis, and it follows the Paris criterion, and it is not the maximum
      ! circumferential tensile stress criterion, then
      ! Calculate the stress intensity factor, because the Paris law requires knowing the stress intensity
      ! factor.
      if(Key_Static_Fatigue==1 .and. CFCP == 2 .and. Key_Fatigue_Cri == 1)then
          if (Key_SIFs_Method==1) then
              call Cal_SIFs_DIM(isub,DISP)
          elseif (Key_SIFs_Method==2) then
              call Cal_SIFs_IIM(isub,.True.,DISP)
          end if
          do i_C = 1,num_Crack
              WRITE(*,1131)i_C,KI(i_C,1)/1.0D6,KII(i_C,1)/1.0D6
              WRITE(*,1132)i_C,KI(i_C,2)/1.0D6,KII(i_C,2)/1.0D6
          end do
      endif
      ! Determine the direction and extent of crack propagation
      call Check_Crack_Grows(1,isub,Yes_Growth)
      ! If it is a life-and-Death element analysis, the program is terminated based on whether the crack
      ! propagates or not.
      if(Key_Ekill==1 ) then
          ! If there are elements to be excavated (killed) in the next step, proceed to the next load step.
          num_Ele_Killed = count(Ele_Killed_Each_Load_Step(isub+1,:) > 0)
          !if(num_Ele_Killed>=1) then
          if(num_Ele_Killed>=1 .or. isub==1) then
              Yes_Last_Growth =.True.
              goto 199
          endif
      endif

      ! If a crack has propagated, update the crack propagation flag: Yes_Last_Growth
      if (any(Yes_Growth).eqv..True.)then
          Yes_Last_Growth = .True.
      ! If no cracks have propagated, exit the program.
      else
          Yes_Last_Growth = .False.
          if(isub < Num_Substeps)then
              print *,'    ---$---$---$---$---$---$---$---$---$---S---S---S---'
              print *,'    Warning :: No crack propapated, program was ended!'
              print *,'    ---$---$---$---$---$---$---$---$---$---S---S---S---'
              goto 200
          elseif(isub == Num_Substeps)then
              print *,'    |<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|'
              print *,'    |    All propagation steps done!  |'
              print *,'    |<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|'
          endif
          goto 199
        endif
    end if
  endif

  199  continue
  
  !---------------------------------------
  ! Check for any new cracks, 2024-06-22.
  !---------------------------------------
  New_Crack_Flag =.false.
  if (Key_Initiation ==1 .and. (Num_Initiation_Cracks < Key_Max_Num_Initiation_Cracks)) then
      call Check_Crack_initialize_2D(1,isub,New_Crack_Flag) 
  endif

  !------------------------------------------------------------------------------------------
  ! If it contains damaged material and may generate cracks (Material_Dam2Frac_Value < ONE),
  ! Then perform the inspection.
  !------------------------------------------------------------------------------------------
  Yes_Generated = .False.
  if(Key_Damage ==1 .and. Material_Dam2Frac_Value < ONE)then 
      call Check_Damage_Crack_Generate(1,isub,Total_Num_G_P,Yes_Generated)
  endif

  !---------------------------------------------------------------
  ! If hole generation is allowed, then only detection is needed.
  !---------------------------------------------------------------
  if(num_Hole.ne.0 .and. Key_Hole_Crack_Generate==1)then
      call Check_Hole_Crack_Generate(1,isub,Yes_Generated)
  endif

  !-------------------------------------------------------------------------------------------------
  ! If necessary, calculate and save the tangential relative displacement at each crack calculation
  ! point.
  !-------------------------------------------------------------------------------------------------
  if (Key_Post_S_TanDisp==1)then
     print *,'    Cal and save tangent relative disp of cracks.'
     allocate(Cracks_CalP_Tan_Aper(num_Crack,Max_Num_Cr_CalP))
     call Cal_Crack_Tan_Relative_Disp(isub,DISP,Cracks_CalP_Tan_Aper)
     call Save_Crack_Tan_Relative_Disp(isub,Cracks_CalP_Tan_Aper)
     deallocate(Cracks_CalP_Tan_Aper)
  endif

  !---------------------
  ! Clear dynamic array
  !---------------------

  if(Key_K_Sparse==0)then
      if (allocated(globalK)) DEALLOCATE(globalK)
  elseif(Key_K_Sparse==1)then
#ifdef gfortran
#ifdef notcbfortran
      if (allocated(K_CSR_aa)) DEALLOCATE(K_CSR_aa)
      if (allocated(K_CSR_ja)) DEALLOCATE(K_CSR_ja)
      if (allocated(K_CSR_ia)) DEALLOCATE(K_CSR_ia)
#endif
#endif
  endif

  if (allocated(freeDOF)) DEALLOCATE(freeDOF)
  if (allocated(fixedDOF)) DEALLOCATE(fixedDOF)
  if (allocated(globalF)) DEALLOCATE(globalF)
  if (allocated(DISP)) DEALLOCATE(DISP)
  if (allocated(tem_disp)) DEALLOCATE(tem_disp)
  if(num_Crack.ne.0)then
      if (allocated(Coupled_Q)) DEALLOCATE(Coupled_Q)
      if (allocated(CalP_Pres)) DEALLOCATE(CalP_Pres)
  endif
  if((num_Crack.ne.0) .and. (Key_Contact/=0)) then
      if (allocated(Contact_DISP)) DEALLOCATE(Contact_DISP)
  endif
  ! Clear the memory space for node stress
  if (Key_Post_CS_N_Strs ==1) then
      if(allocated(Stress_xx_Node)) DEALLOCATE(Stress_xx_Node)
      if(allocated(Stress_yy_Node)) DEALLOCATE(Stress_yy_Node)
      if(allocated(Stress_xy_Node)) DEALLOCATE(Stress_xy_Node)
      if(allocated(Stress_vm_Node)) DEALLOCATE(Stress_vm_Node)
      if(allocated(TStress_xx_Node)) DEALLOCATE(TStress_xx_Node)
      if(allocated(TStress_yy_Node)) DEALLOCATE(TStress_yy_Node)
      if(allocated(TStress_xy_Node)) DEALLOCATE(TStress_xy_Node)
      if(allocated(TStress_vm_Node)) DEALLOCATE(TStress_vm_Node)
  end if
  ! Clear Gauss point coordinate memory
  if (Key_Post_CS_G_Coor==1) then
      if (allocated(Gauss_CoorX)) DEALLOCATE(Gauss_CoorX)
      if (allocated(Gauss_CoorY)) DEALLOCATE(Gauss_CoorY)
  endif
  ! Clear Gauss point displacement memory
  if (Key_Post_CS_G_Disp==1) then
      if (allocated(DISP_x_Gauss)) DEALLOCATE(DISP_x_Gauss)
      if (allocated(DISP_y_Gauss)) DEALLOCATE(DISP_y_Gauss)
  endif
  ! Clear Gauss point stress memory space
  if (Key_Post_CS_G_Strs==1) then
    if(allocated(Stress_xx_Gauss)) DEALLOCATE(Stress_xx_Gauss)
    if(allocated(Stress_yy_Gauss)) DEALLOCATE(Stress_yy_Gauss)
    if(allocated(Stress_xy_Gauss)) DEALLOCATE(Stress_xy_Gauss)
    if(allocated(Stress_vm_Gauss)) DEALLOCATE(Stress_vm_Gauss)
    if(allocated(TStress_xx_Gauss)) DEALLOCATE(TStress_xx_Gauss)
    if(allocated(TStress_yy_Gauss)) DEALLOCATE(TStress_yy_Gauss)
    if(allocated(TStress_xy_Gauss)) DEALLOCATE(TStress_xy_Gauss)
    if(allocated(TStress_vm_Gauss)) DEALLOCATE(TStress_vm_Gauss)
  endif
  if(Key_Damage ==1)then
      if(allocated(Damage_Gauss)) deallocate(Damage_Gauss)
  endif
  if(Key_Contact/=0 .and. num_Crack.ne.0)then
      if(allocated(ori_globalK)) DEALLOCATE(ori_globalK)
  endif
  !Clear solver 11 (EBE) related temporary variables.
  if (allocated(storK)) DEALLOCATE(storK)
  if (allocated(diag_precon)) DEALLOCATE(diag_precon)
  if (allocated(size_local)) DEALLOCATE(size_local)
  if (allocated(all_local)) DEALLOCATE(all_local)

  ! If all cracks have stopped propagating, no new cracks have formed, and they are not part of the
  ! user-defined path analysis, then terminate the program.
  if(num_Crack.ne.0)then
      if((Yes_Last_Growth .eqv. .False.) .and. &
         (Yes_Generated .eqv. .False.)   .and. &
         (New_Crack_Flag .eqv. .False.)  .and. &
         (Key_User_Defined_2D_Crack_Path/=1)  )then
          goto 200
      endif
  endif

end do

200 continue

!-----------------------------------------------------------
! Delete temporary binary files for the group sparse matrix
!-----------------------------------------------------------
if (file_Sparse_K_Location_COO_bin .eqv. .True.)then
c_File_name_1=trim(Full_Pathname)//'_Sparse_K_Location_COO.bin'
print *,'    Deleting temporary binary files......'
Open(1234,File =c_File_name_1,Access='Direct',Form = 'Unformatted', RecL = 4)
close(1234, status='delete')
endif
if (file_Sparse_K_Mask_COO_bin .eqv. .True.)then
c_File_name_1=trim(Full_Pathname)//'_Sparse_K_Mask_COO.bin'
print *,'    Deleting temporary binary files......'
Open(1234,File =c_File_name_1,Access='Direct',Form = 'Unformatted', RecL = 4)
close(1234, status='delete')
endif

!----------------------------------------------------------------
! Release memory after calculation to prevent unexpected issues.
!----------------------------------------------------------------
if (Key_Contact /= 0) then
  if (allocated(Elem_Conta_Sta)) DEALLOCATE(Elem_Conta_Sta)
  if (allocated(Elem_Conta_Sta_Last)) then
      DEALLOCATE(Elem_Conta_Sta_Last)
  endif
  if (allocated(Ele_NumCalP)) DEALLOCATE(Ele_NumCalP)
  if (allocated(Ele_CalPNum)) DEALLOCATE(Ele_CalPNum)
end if
if(Key_PoreP==1)then
  if (allocated(Biot_c_MAT)) DEALLOCATE(Biot_c_MAT)
  if (allocated(Porous_P)) DEALLOCATE(Porous_P)
endif
if (allocated(storK)) DEALLOCATE(storK)
if (allocated(diag_precon)) DEALLOCATE(diag_precon)
if (allocated(size_local)) DEALLOCATE(size_local)
if (allocated(all_local)) DEALLOCATE(all_local)


RETURN
END SUBROUTINE PhiPsi2D_Static
