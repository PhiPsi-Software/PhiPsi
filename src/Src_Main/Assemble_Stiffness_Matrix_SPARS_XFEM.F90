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
 
SUBROUTINE Assemble_Stiffness_Matrix_SPARS_XFEM(isub, &
                     freeDOF,num_FreeD,&
                     fixedDOF,num_FixedD,&
                     K_CSR_aa,K_CSR_ja,K_CSR_ia,&
                     K_CSR_NNZ_Max,K_CSR_NNZ,&
                     T_Freedom,Total_Num_G_P)
! Assemble the stiffness matrix, including the problem with cracks (store K in COO sparse matrix
! format).

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
use Global_Contact
use omp_lib
use CMD_Progress
use ISO_FORTRAN_ENV
use module_INTERFACE_coocsr
use Global_Inclusion,only:num_Inclusion,num_Circ_Incl,c_POS_Incl,Enriched_Node_Type_Incl,num_Poly_Incl,&
                          Poly_Inclu_Mat_Num,Circ_Inclu_Mat_Num
use Global_Cross,only:num_cross,c_POS_Cross,Enriched_Node_Type_Cross


!**********************
! Variable Declaration
!**********************
implicit none
integer,intent(in)::isub,num_FreeD,T_Freedom
integer,intent(in)::K_CSR_NNZ_Max
integer,intent(in)::freeDOF(1:num_FreeD)
integer,intent(in)::num_FixedD
integer,intent(in)::fixedDOF(1:num_FixedD)
integer,intent(out)::Total_Num_G_P
real(kind=FT),intent(out)::K_CSR_aa(K_CSR_NNZ_Max)
integer,intent(out)::K_CSR_ja(K_CSR_NNZ_Max)
integer,intent(out)::K_CSR_ia(num_FreeD+1)
integer,intent(out)::K_CSR_NNZ
!---------------------
logical Yes_Gauss_in_Incl
integer c_Incl_Num,c_MatNum
integer i_C,i_E,i_G
integer i_H,i_Cross,i_Incl
real(kind=FT) c_thick,c_D(3,3)
real(kind=FT) c_X_NODES(4),c_Y_NODES(4)
integer c_NN(4)
! Gauss point data of enhanced elements, integration scheme 2
real(kind=FT) kesi_Enr(Num_Gauss_Points),yita_Enr(Num_Gauss_Points),&
                 weight_Enr(Num_Gauss_Points)
real(kind=FT) kesi_N_Enr(4),           &
                 yita_N_Enr(4),&
                 weight_N_Enr(4)
real(kind=FT) kesi(Num_Gauss_Points),yita(Num_Gauss_Points),&
                 weight(Num_Gauss_Points)
integer:: Location_ESM(MDOF_2D)
! integer:: Location_ESM_NoFEM(160) !Does not include FEM degrees of freedom, 32*5, each element has
! at most 32 enhanced degrees of freedom, each element is assumed to be related to at most 5 cracks
integer::Location_ESM_C_Crack(80),Location_ESM_C_Cr_NoFEM(60)
! integer::Loc_ESM_C_Crack(50) ! Numbering of element enhanced degrees of freedom in the total
! stiffness (including conventional finite element degrees of freedom)
! integer::Loc_ESM_C_Cr_NoFEM(32) ! The numbering of element enhanced degrees of freedom in the
! global stiffness matrix (including only extended finite element numbers)
integer num_Loc_ESM_C_Cr_NoFEM
real(kind=FT) B(3,80),tem_B(3,80)
integer num_B,num_tem_B
integer num_Loc_ESM
integer num_Loc_ESM_C_Crack
integer c_Num_Gauss_Point
real(kind=FT) detJ
! The following variables are related to the set of sparse matrices
integer i_Check,K_COO_NNZ_try,i_delete,i_row,i_col
real(kind=FT),ALLOCATABLE::K_COO_aa_try(:)
integer,ALLOCATABLE::K_COO_ia_try(:)
integer,ALLOCATABLE::K_COO_ja_try(:)
! Used for option1 (fast, memory-intensive sparse matrix set algorithm)
! integer*4, ALLOCATABLE:: K_num_NNZ(:,:) ! This matrix stores the non-zero element number
! corresponding to each element of the stiffness matrix
                                              ! This variable consumes a lot of memory, but it can significantly speed up the algorithm.
! Used for option2 (slow, memory-efficient sparse matrix grouping algorithm)
! integer, ALLOCATABLE :: tem_K_COO_ia_try(:) ! Row numbers in the overall K matrix corresponding to
! the K_COO_NNZ_try-th non-zero element
! integer, ALLOCATABLE :: tem_K_COO_ja_try(:) ! The column number of the overall K matrix
! corresponding to the K_COO_NNZ_try-th non-zero element
!**********************************************************************************************************
logical(KIND=1),ALLOCATABLE::Mask_COO(:,:)
integer(kind=4),ALLOCATABLE::Location_COO(:,:)
!**********************************************************************************************************
integer c_row,c_col,c_num_NNZ
real(kind=FT) c_G_x,c_G_y
real(kind=FT) c_N(2,8)
integer real_row,real_col
real(kind=FT) localK(MDOF_2D,MDOF_2D)
! For conversion to CSR format
real(kind=FT),ALLOCATABLE::  K_COO_aa_try2(:)
integer,ALLOCATABLE::  K_COO_ia_try2(:)
integer,ALLOCATABLE::  K_COO_ja_try2(:)
logical c_Yes_In
logical values
integer,ALLOCATABLE:: iwork(:)
integer job,value2
integer indu(num_FreeD),iwk(num_FreeD+1)
logical logi_Yes_Mask_COO
integer c_Location
integer(KIND=1) int_0,int_1
integer(KIND=1) c_int_0_or_1
integer(KIND=1),ALLOCATABLE:: int_0_vector(:)
integer(KIND=4) location_0
integer(KIND=4),ALLOCATABLE:: location_vector(:)
type( CLS_CMD_Progress ) ::Progress
type( CLS_CMD_Progress ) ::Progress_Mask_COO
type( CLS_CMD_Progress ) ::Progress_Location_COO
character(200) c_File_name_1,c_File_name_2
integer i,j
integer(kind=int64) Rec_num
integer num_Ele_Killed
integer, parameter :: integer_1 = selected_int_kind(1)

integer Scheme 
! 2024-11-07. Scheme = 1 New algorithm related.
real(kind=FT),allocatable::Temp_K_Each(:,:)
real(kind=FT),allocatable::K_Each(:,:)
integer K_Each_NNZ_Max  
real(kind=FT),allocatable::Combined_K_Each(:,:)
integer Temp_K_Each_Count,K_Each_Count,num_rows_Combined_K
integer,allocatable:: K_rows(:), K_cols(:)
real(kind=FT),allocatable:: K_vals(:)
integer :: progress_counter = 0
real(kind=FT) :: progress_percent
integer :: bar_width = 50
integer :: filled
character(len=1) :: spin(4) = ['|','/','-','\']
integer :: spin_idx = 1



1112 FORMAT(5X,'Real density ratio of K is ',F12.4,'%')


! Stiffness Matrix Assembly Progress Bar Settings
#ifndef Silverfrost
call Progress % Set( N = Num_Elem , L = 25 )
Progress % Prefix = "     K assembling process:  "
Progress % M = ">"
Progress % O = "."
! Progress_Mask_COO Progress Bar Settings
call Progress_Mask_COO % Set( N = T_Freedom/2 , L = 20 )
Progress_Mask_COO % Prefix ="     Mask_COO generating process:      "
Progress_Mask_COO % M = ">"
Progress_Mask_COO % O = "."
! Progress_Location_COO Progress Bar Settings
call Progress_Location_COO % Set( N = T_Freedom/2 , L = 20 )
Progress_Location_COO % Prefix ="     Location_COO generating process:  "
Progress_Location_COO % M = ">"
Progress_Location_COO % O = "."
#endif




Scheme = 1
            ! =2, old algorithm. Deprecated.

!////////////////////////
!                      /
! New Algorithm /
!       2024-11-07     /
!   NEWFTU2024110703.  /
!                      /
!////////////////////////
if (Scheme ==1) then
    !***************************************************************
    ! Step 0: Initialization of some variables and data preparation
    !***************************************************************
    K_CSR_aa(1:K_CSR_NNZ_Max)   = ZR
    K_CSR_ja(1:K_CSR_NNZ_Max)   = 0
    K_CSR_ia(1:num_FreeD+1)   = 0
    K_CSR_NNZ = 0
    Total_Num_G_P  = 0
    
    ! Life-and-death unit preparation
    if(Key_EKILL==1)then
        num_Ele_Killed = count(Ele_Killed_Each_Load_Step(1:isub,:)>0)
    endif

    ! Standard 64 Gauss points
    if (Key_Integral_Sol  == 2)then
        call Cal_Gauss_Points_QUAD(Num_Gauss_Points,kesi_Enr,yita_Enr,weight_Enr)
    ! Subdivision of quadrilaterals to calculate Gauss points (the Gauss points are unrelated to cracks,
    ! and are for several regular quadrilaterals)
    elseif (Key_Integral_Sol  == 3)then
        call Cal_Gauss_Points_QUAD_for_SUBQUAD(Num_Sub_Quads,kesi_Enr,yita_Enr,weight_Enr)
        Num_Gauss_Points = Num_Sub_Quads*4
    endif
    call Cal_Gauss_Points_QUAD(4,kesi_N_Enr,yita_N_Enr,weight_N_Enr)
    
    !********************************************************
    ! Step 1: Obtain the maximum number of non-zero elements
    !********************************************************
    EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gauss_Points)= .False.

    
    Temp_K_Each_Count = 0
    
    do i_E = 1,Num_Elem
        ! Progress bar.
        !bar_progress = dble(i_E)/Num_Elem
        !bar = '['//repeat('=', int(bar_progress*50))//repeat(' ',50-int(bar_progress*50))//']'
        !write(*, '(a,a,f6.1,a,a,a)', advance='no') char(13), '    ',bar_progress*100, '%  ', trim(bar)

          c_NN    = G_NN(:,i_E)
         
          ! ------------------------------
          ! Point Scheme 1: Triangulation
          ! ------------------------------
          if(Key_Integral_Sol.eq.1)then
              !call Cal_Gauss_Points_Subtriangle(kesi,yita,weight)
              !c_Num_Gauss_Point  = size(kesi,2);
          ! -------------------------------------
          ! Point Plan 2 or 3: Even Distribution
          ! -------------------------------------
          elseif(Key_Integral_Sol.eq.2 .or. Key_Integral_Sol.eq.3)then
              ! Starting number of Gauss points for each unit
              Ele_GP_Start_Num(i_E) = Total_Num_G_P + 1
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Determine the Gauss points of the current element
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !if the current element are enriched element, then 8x8 gauss points is suggested:
              if(num_Crack/= 0 .and. (maxval(Enriched_Node_Type(c_NN,1:num_Crack)).ne.0))then
                  c_Num_Gauss_Point = Num_Gauss_Points
                  Total_Num_G_P     = Total_Num_G_P +c_Num_Gauss_Point 
              ! If it is a Hole enhanced node, then 8x8 Gauss points are suggested:
              elseif(num_Hole/= 0 .and.(maxval(Enriched_Node_Type_Hl(c_NN,1:num_Hole)).ne.0))then
                  c_Num_Gauss_Point = Num_Gauss_Points
                  Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
              ! If it is a Cross-enhanced node, then 8x8 Gauss points are suggested:
              elseif(num_Cross/= 0 .and. (maxval(Enriched_Node_Type_Cross(c_NN,1:num_Cross)).ne.0))then
                  c_Num_Gauss_Point = Num_Gauss_Points
                  Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
              ! If there are embedded reinforcement nodes, then 8x8 Gauss points are suggested:
              elseif(num_Inclusion/= 0 .and.(maxval(Enriched_Node_Type_Incl(c_NN,1:num_Inclusion)).ne.0))then
                  c_Num_Gauss_Point = Num_Gauss_Points
                  Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point  
              !if the current element are not enriched element, then 2x2 gauss points:
              else
                  c_Num_Gauss_Point = 4
                  Total_Num_G_P     = Total_Num_G_P +c_Num_Gauss_Point
              end if
              
              ! Save the number of Gauss points for each element to a global variable
              num_GP_Elem(i_E) = c_Num_Gauss_Point

              Location_ESM(1:MDOF_2D)       = 0
              !Location_ESM_NoFEM(1:160) = 0
              num_Loc_ESM       = 0
              !num_Loc_ESM_NoFEM = 0
              
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            !Decide the location of each element stiffness 
            !matrix in the global stiffness matrix for cracks.   
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if(num_Crack/=0)then
                  do i_C =1,num_Crack
                      call Location_Element_Stiff_Matrix(i_E,i_C,&
                                                c_POS(:,i_C),&
                                                Location_ESM_C_Crack,&
                                                num_Loc_ESM_C_Crack,&
                                                Location_ESM_C_Cr_NoFEM,&
                                                num_Loc_ESM_C_Cr_NoFEM)
                      ! Includes FEM degrees of freedom
                      Location_ESM(num_Loc_ESM+1:&
                              num_Loc_ESM+num_Loc_ESM_C_Crack) =&
                              Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
                      num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
                  end do
              endif
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            !Decide the location of each element stiffness 
            !matrix in the global stiffness matrix for hole.   
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if(num_Hole/=0)then
                do i_H =1,num_Hole
                    call Location_Element_Stiff_Matrix_Hl(i_E,i_H, &
                                         c_POS_Hl(1:Num_Node,i_H), &
                                         Location_ESM_C_Crack, &
                                         num_Loc_ESM_C_Crack, &
                                         Location_ESM_C_Cr_NoFEM,&
                                         num_Loc_ESM_C_Cr_NoFEM)
                    ! Includes FEM degrees of freedom
                    Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = &
                               Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
                    num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
                end do
            endif
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            !Decide the location of each element stiffness 
            !matrix in the global stiffness matrix for cross.   
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if(num_Cross/=0)then
                do i_Cross =1,num_Cross
                    call Location_Element_Stiff_Matrix_Cross(i_E,i_Cross,&
                                         c_POS_Cross(1:Num_Node,i_Cross),&
                                         Location_ESM_C_Crack,&
                                         num_Loc_ESM_C_Crack,&
                                         Location_ESM_C_Cr_NoFEM,&
                                         num_Loc_ESM_C_Cr_NoFEM)
                    ! Includes FEM degrees of freedom
                    Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = &
                               Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
                    num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
                end do
            endif
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            !Decide the location of each element stiffness 
            !matrix in the global stiffness matrix for inclusion.   
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if(num_Inclusion/=0)then
                do i_Incl =1,num_Inclusion
                    call Location_Element_Stiff_Matrix_Incl(i_E,i_Incl, &
                                            c_POS_Incl(1:Num_Node,i_Incl),&
                                            Location_ESM_C_Crack,&
                                            num_Loc_ESM_C_Crack,&
                                            Location_ESM_C_Cr_NoFEM,&
                                            num_Loc_ESM_C_Cr_NoFEM)
                    ! Includes FEM degrees of freedom
                    Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = &
                               Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
                    num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
                end do
            endif

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Extract the non-zero elements of the stiffness matrix.
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do i_row = 1,num_Loc_ESM
              do i_col = 1,num_Loc_ESM
                  c_row = Location_ESM(i_row)
                  c_col = Location_ESM(i_col)
                  ! The corresponding degree of freedom needs to be a functional degree of freedom.
                  if(Flag_FreeDOF(c_row)==1  .and. Flag_FreeDOF(c_col)==1) then
                    Temp_K_Each_Count =  Temp_K_Each_Count + 1
                  endif
              end do
            end do
          endif
    end do
    
    
    !**************************************************************
    ! Step 2: Obtain the non-zero elements of the stiffness matrix
    !**************************************************************
    K_Each_NNZ_Max = Temp_K_Each_Count
    allocate(Temp_K_Each(K_Each_NNZ_Max,3))
    EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gauss_Points)= .False.
    Total_Num_G_P = 0
    Temp_K_Each_Count = 0
    progress_counter = 0
    do i_E = 1,Num_Elem
        ! Progress bar.
        !bar_progress = dble(i_E)/Num_Elem
        !bar = '['//repeat('=', int(bar_progress*50))//repeat(' ',50-int(bar_progress*50))//']'
        !write(*, '(a,a,f6.1,a,a,a)', advance='no') char(13), '    ',bar_progress*100, '%  ', trim(bar)
        
        c_thick = thick(Elem_Mat(i_E))
        c_D     = D(Elem_Mat(i_E),:,:)
          
        ! Elastic modulus of the Weibull distribution. 2024-06-25. NEWFTU2024062402.
        if(Flag_Weibull_E)then
            if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
                c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
            endif
        endif
          
        ! Life-and-death unit processing: weakening D matrix (elastic modulus)
        if(Key_EKILL==1 .and. num_Ele_Killed>=1) then
            if(any(Ele_Killed_Each_Load_Step(1:isub,1:num_Ele_Killed)==i_E))then
                c_D= c_D*EKILL_Weaken_Factor
            endif
        endif
        
        ! Damage unit processing: weakened D matrix (elastic modulus), 2020-02-28
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
              !call Cal_Gauss_Points_Subtriangle(kesi,yita,weight)
              !c_Num_Gauss_Point  = size(kesi,2);
          ! -------------------------------------
          ! Point Plan 2 or 3: Even Distribution
          ! -------------------------------------
          elseif(Key_Integral_Sol.eq.2 .or. Key_Integral_Sol.eq.3)then
              ! Starting number of Gauss points for each unit
              Ele_GP_Start_Num(i_E) = Total_Num_G_P + 1
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Determine the Gauss points of the current element
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !if the current element are enriched element, then 8x8 gauss points is suggested:
              if(num_Crack/= 0 .and. maxval(Enriched_Node_Type(c_NN,1:num_Crack)).ne.0)then
                  kesi(1:Num_Gauss_Points)    = kesi_Enr
                  yita(1:Num_Gauss_Points)    = yita_Enr
                  weight(1:Num_Gauss_Points)  = weight_Enr
                  c_Num_Gauss_Point = Num_Gauss_Points
                  Total_Num_G_P     = Total_Num_G_P +c_Num_Gauss_Point 
              ! If it is a Hole enhanced node, then 8x8 Gauss points are suggested:
              elseif(num_Hole/= 0 .and.(maxval(Enriched_Node_Type_Hl(c_NN,1:num_Hole)).ne.0))then
                  kesi(1:Num_Gauss_Points)    = kesi_Enr
                  yita(1:Num_Gauss_Points)    = yita_Enr
                  weight(1:Num_Gauss_Points)  = weight_Enr
                  c_Num_Gauss_Point = Num_Gauss_Points
                  Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
              ! If it is a Cross-enhanced node, then 8x8 Gauss points are suggested:
              elseif(num_Cross/= 0 .and. (maxval(Enriched_Node_Type_Cross(c_NN,1:num_Cross)).ne.0))then
                  kesi(1:Num_Gauss_Points)    = kesi_Enr
                  yita(1:Num_Gauss_Points)    = yita_Enr
                  weight(1:Num_Gauss_Points)  = weight_Enr
                  c_Num_Gauss_Point = Num_Gauss_Points
                  Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
              ! If there are embedded reinforcement nodes, then 8x8 Gauss points are suggested:
              elseif(num_Inclusion/= 0 .and.(maxval(Enriched_Node_Type_Incl(c_NN,1:num_Inclusion)).ne.0))then
                  kesi(1:Num_Gauss_Points)    = kesi_Enr
                  yita(1:Num_Gauss_Points)    = yita_Enr
                  weight(1:Num_Gauss_Points)  = weight_Enr
                  c_Num_Gauss_Point = Num_Gauss_Points
                  Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point  
              !if the current element are not enriched element, then 2x2 gauss points:
              else
                  kesi(1:4)    = kesi_N_Enr
                  yita(1:4)    = yita_N_Enr
                  weight(1:4)  = weight_N_Enr
                  c_Num_Gauss_Point = 4
                  Total_Num_G_P     = Total_Num_G_P +c_Num_Gauss_Point
              end if

              
              ! Save the number of Gauss points for each element to a global variable
              num_GP_Elem(i_E) = c_Num_Gauss_Point

              Location_ESM(1:MDOF_2D)       = 0
              !Location_ESM_NoFEM(1:160) = 0
              num_Loc_ESM       = 0
              
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !Decide the location of each element stiffness
              !matrix in the global stiffness matrix.
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !num_Loc_ESM_NoFEM = 0
              if(num_Crack/=0)then
                  do i_C =1,num_Crack
                      call Location_Element_Stiff_Matrix(i_E,i_C,&
                                                c_POS(:,i_C),&
                                                Location_ESM_C_Crack,&
                                                num_Loc_ESM_C_Crack,&
                                                Location_ESM_C_Cr_NoFEM,&
                                                num_Loc_ESM_C_Cr_NoFEM)
                      ! Includes FEM degrees of freedom
                      Location_ESM(num_Loc_ESM+1:&
                              num_Loc_ESM+num_Loc_ESM_C_Crack) =&
                              Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
                      num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
                  end do
              endif
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            !Decide the location of each element stiffness 
            !matrix in the global stiffness matrix for hole.   
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if(num_Hole/=0)then
                do i_H =1,num_Hole
                    call Location_Element_Stiff_Matrix_Hl(i_E,i_H, &
                                         c_POS_Hl(1:Num_Node,i_H), &
                                         Location_ESM_C_Crack, &
                                         num_Loc_ESM_C_Crack, &
                                         Location_ESM_C_Cr_NoFEM,&
                                         num_Loc_ESM_C_Cr_NoFEM)
                    ! Includes FEM degrees of freedom
                    Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = &
                               Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
                    num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
                end do
            endif
            
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            !Decide the location of each element stiffness 
            !matrix in the global stiffness matrix for cross.   
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if(num_Cross/=0)then
                do i_Cross =1,num_Cross
                    call Location_Element_Stiff_Matrix_Cross(i_E,i_Cross,&
                                         c_POS_Cross(1:Num_Node,i_Cross),&
                                         Location_ESM_C_Crack,&
                                         num_Loc_ESM_C_Crack,&
                                         Location_ESM_C_Cr_NoFEM,&
                                         num_Loc_ESM_C_Cr_NoFEM)
                    ! Includes FEM degrees of freedom
                    Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = &
                               Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
                    num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
                end do
            endif
            
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            !Decide the location of each element stiffness 
            !matrix in the global stiffness matrix for inclusion.   
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if(num_Inclusion/=0)then
                do i_Incl =1,num_Inclusion
                    call Location_Element_Stiff_Matrix_Incl(i_E,i_Incl, &
                                            c_POS_Incl(1:Num_Node,i_Incl),&
                                            Location_ESM_C_Crack,&
                                            num_Loc_ESM_C_Crack,&
                                            Location_ESM_C_Cr_NoFEM,&
                                            num_Loc_ESM_C_Cr_NoFEM)
                    ! Includes FEM degrees of freedom
                    Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = &
                               Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
                    num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
                end do
            endif

              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Loop over each Gauss point to assemble the element stiffness matrix.
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              localK(1:MDOF_2D,1:MDOF_2D) = ZR
              do i_G = 1,c_Num_Gauss_Point
                  !B(1:3,1:50) = ZR
                  B(1:3,1:80) = ZR
                  num_B = 0
                  call Cal_N(kesi(i_G),yita(i_G),c_N)
                  call Cal_detJ(kesi(i_G),yita(i_G),c_X_NODES,c_Y_NODES,detJ)
                  !Global coordinates of the gauss point.
                  c_G_x = DOT_PRODUCT(c_N(1,1:7:2),c_X_NODES(1:4))
                  c_G_y = DOT_PRODUCT(c_N(1,1:7:2),c_Y_NODES(1:4))  
                  
                  !Calculate the B Matrix, Loop through each crack.
                  if(num_Crack/=0)then
                      do i_C =1,num_Crack
                          call Cal_B_Matrix_Crack(kesi(i_G),yita(i_G),&
                                            i_C,i_E,i_G,&
                                            c_NN,c_X_NODES,c_Y_NODES,&
                                            tem_B,num_tem_B)
                          B(1:3,num_B+1:num_B+num_tem_B) =tem_B(1:3,1:num_tem_B)
                          num_B = num_B + num_tem_B
                      end do
                  endif   
                  
                 !Calculate the B Matrix, Loop through each hole.
                  if(num_Hole/=0)then
                      do i_H =1,num_Hole 
                          call Cal_B_Matrix_Hl(kesi(i_G),yita(i_G), &
                                             i_H,i_E,i_G, &
                                             c_NN,c_X_NODES,c_Y_NODES, &
                                             tem_B,num_tem_B)
                          B(1:3,num_B+1:num_B+num_tem_B) =  tem_B(1:3,1:num_tem_B)
                          num_B = num_B + num_tem_B
                      end do
                  endif
                  !Calculate the B Matrix, Loop through each cross.
                  if(num_Cross/=0)then
                      do i_Cross =1,num_Cross
                          call Cal_B_Matrix_Cross(kesi(i_G),yita(i_G), &
                                             i_Cross,i_E,i_G, &
                                             c_NN,c_X_NODES,c_Y_NODES, &
                                             tem_B,num_tem_B)
                          B(1:3,num_B+1:num_B+num_tem_B) = tem_B(1:3,1:num_tem_B)
                          num_B = num_B + num_tem_B
                      end do
                  endif
                  !Calculate the B Matrix, Loop through each circle inclusion.
                  if(num_Circ_Incl/=0)then
                      do i_Incl =1,num_Circ_Incl 
                          call Cal_B_Matrix_Circ_Incl(kesi(i_G), &
                                            yita(i_G),i_Incl,i_E,i_G, &
                                             c_NN,c_X_NODES,c_Y_NODES, &
                                             tem_B,num_tem_B)
                          B(1:3,num_B+1:num_B+num_tem_B) =  tem_B(1:3,1:num_tem_B)
                          num_B = num_B + num_tem_B
                      end do
                      ! Determine whether the current Gauss point is within a certain inclusion
                      call Tool_Yes_Point_in_Inclusions(c_G_x,c_G_y,Yes_Gauss_in_Incl,c_Incl_Num)
                      if(Yes_Gauss_in_Incl)then
                          ! Obtain the current composite material matrix D
                          c_MatNum = Circ_Inclu_Mat_Num(c_Incl_Num)
                          c_D     = D(c_MatNum,:,:)   
                      endif
                  endif
                  !Calculate the B Matrix, Loop through each polygon inclusion.
                  if(num_Poly_Incl/=0)then
                      do i_Incl =1,num_Poly_Incl
                          call Cal_B_Matrix_Poly_Incl(kesi(i_G), &
                                            yita(i_G),i_Incl,i_E,i_G, &
                                             c_NN,c_X_NODES,c_Y_NODES, &
                                             tem_B,num_tem_B)
                          B(1:3,num_B+1:num_B+num_tem_B) =  tem_B(1:3,1:num_tem_B)
                          num_B = num_B + num_tem_B
                      end do
                      ! Determine whether the current Gauss point is inside an inclusion (including circular and polygonal
                      ! inclusions)
                      call Tool_Yes_Point_in_Inclusions(c_G_x,c_G_y,Yes_Gauss_in_Incl,c_Incl_Num)
                      if(Yes_Gauss_in_Incl)then
                          ! Obtain the current composite material matrix D
                          c_MatNum = Poly_Inclu_Mat_Num(c_Incl_Num)
                          c_D     = D(c_MatNum,:,:)   
                      endif
                  endif
     
                  localK(1:num_B,1:num_B) = localK(1:num_B,1:num_B) +&
                               c_thick*weight(i_G)*detJ*&
                               MATMUL(MATMUL(transpose(B(1:3,1:num_B)),c_D),B(1:3,1:num_B))
              enddo

                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ! Extract the non-zero elements of the stiffness matrix.
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do i_row = 1,num_Loc_ESM
                  do i_col = 1,num_Loc_ESM
                      c_row = Location_ESM(i_row)
                      c_col = Location_ESM(i_col)

                      ! The corresponding degree of freedom needs to be a functional degree of freedom.
                      if(Flag_FreeDOF(c_row)==1  .and. Flag_FreeDOF(c_col)==1) then
                         ! Row numbers corresponding to the degrees of freedom after removing constraints
                         ! Real_FreeDOF_Index is the stiffness matrix degree of freedom index after removing the constrained
                         ! degrees of freedom corresponding to all degrees of freedom. 2024-11-07.
                         real_row =  Real_FreeDOF_Index(c_row) 
                         real_col =  Real_FreeDOF_Index(c_col) 
                        Temp_K_Each_Count =  Temp_K_Each_Count + 1
                        
                        Temp_K_Each(Temp_K_Each_Count,1) = real(real_row)
                        Temp_K_Each(Temp_K_Each_Count,2) = real(real_col)
                        Temp_K_Each(Temp_K_Each_Count,3) = localK(i_row,i_col)
                      endif
                  end do
                end do
          endif
          
        ! Progress bar.
!        progress_counter = progress_counter + 1
! if(mod(progress_counter, max(1,Num_Elem/200)) == 0 .or. progress_counter==Num_Elem) then ! Update
! approximately every 0.5% completion
!            progress_percent = real(progress_counter)/real(Num_Elem)*100
!            filled = nint(bar_width * progress_counter / real(Num_Elem))
! ! Update rotation animation
!            spin_idx = mod(spin_idx, 4) + 1
! ! Output progress bar
!            write(*,'(A,A8,A1)',advance='no') char(13), '     Now', ' '
!            write(*,'(A1)',advance='no') '['
! ! Output completed part
!            write(*,'(A)',advance='no') repeat('=', filled)
! ! Output unfinished part
!            write(*,'(A)',advance='no') repeat(' ', bar_width-filled)
! ! Output percentages and counts
!            if(progress_counter==Num_Elem) spin_idx=1
!            write(*,'(A1,F6.2,A2,A1,A2,I0,A1,I0,A1)',advance='no') &
!                ']', progress_percent, '% ', spin(spin_idx), ' (', &
!                progress_counter, '/', Num_Elem, ')'
!            call flush(6)
!        endif

    end do
    
! write (*,'(a)') ' '    !You need to activate this line to use the progress bar.
    
    !*******************************
    !               Step 3
    ! K merge duplicate data (sum).
    !*******************************
    K_Each_Count = Temp_K_Each_Count
    allocate(K_Each(K_Each_Count,3))
    K_Each(1:K_Each_Count,1:3) = Temp_K_Each(1:K_Each_Count,1:3)
    deallocate(Temp_K_Each)
    allocate(Combined_K_Each(K_Each_Count,3))
    !write ( * ,'(a)') 'Combining K duplicate lines...' 
    call Matrix_n_x_3_Combine_Same_i_j_Lines(K_Each_Count,K_Each(1:K_Each_Count,1:3),&
                                             Combined_K_Each(1:K_Each_Count,1:3), &
                                             num_rows_Combined_K)
                                             
    !**************************
    !           Step 4
    ! Regenerate the merged K.
    !**************************
    deallocate(K_Each)
    allocate(K_Each(num_rows_Combined_K,3))
    K_Each(1:num_rows_Combined_K,1:3)  = Combined_K_Each(1:num_rows_Combined_K,1:3)
    deallocate(Combined_K_Each)
    
    
    !***************************************************
    !           Step 5
    ! Convert matrix K into a COO sparse format matrix.
    !***************************************************
    K_CSR_NNZ = 0
    do i=1,num_rows_Combined_K
        if(K_Each(i,3)/=ZR) then
            K_CSR_NNZ = K_CSR_NNZ +1 
        endif      
    enddo 
    
    if(K_CSR_NNZ > K_CSR_NNZ_Max) then
        print *,'    ERROR :: K_CSR_NNZ > K_CSR_NNZ_Max in Assemble_Stiffness_Matrix_SPARS_XFEM.f90!'
        print *,"         K_CSR_NNZ:     ",K_CSR_NNZ
        print *,"         K_CSR_NNZ_Max: ",K_CSR_NNZ_Max
        call Warning_Message('S',Keywords_Blank)
    endif
    
    
    ! Allocate row indices, column indices, and non-zero value arrays
    allocate(K_rows(K_CSR_NNZ), K_cols(K_CSR_NNZ), K_vals(K_CSR_NNZ))
    
    ! Fill the rows, columns, and values of a sparse matrix
    K_CSR_NNZ = 0
    do i=1,num_rows_Combined_K
        !if(abs(MXS(i,3))>Tol_15) then
        !if(abs(MXS(i,3))>Tol_7) then  
        if(K_Each(i,3)/=ZR) then
            K_CSR_NNZ = K_CSR_NNZ +1 
            K_rows(K_CSR_NNZ) = int(K_Each(i, 1))
            K_cols(K_CSR_NNZ) = int(K_Each(i, 2))
            K_vals(K_CSR_NNZ) = K_Each(i, 3)
        endif     
    enddo 
    
    write(*,1112) DBLE(K_CSR_NNZ)/T_Freedom/T_Freedom*100.0D0
    
    !************************
    !           Step 6
    ! Convert to CSR format.
    !************************
    call coocsr(num_FreeD,K_CSR_NNZ,K_vals(1:K_CSR_NNZ),K_rows(1:K_CSR_NNZ),K_cols(1:K_CSR_NNZ),&
                       K_CSR_aa(1:K_CSR_NNZ),K_CSR_ja(1:K_CSR_NNZ),K_CSR_ia(1:num_FreeD+1))
                       
    deallocate(K_rows,K_cols,K_vals)    
    

endif

!////////////////////////
!                      /
! Old algorithm /
!                      /
!////////////////////////
if (Scheme ==2) then
    !***************************************************************
    ! Step 0: Initialization of some variables and data preparation
    !***************************************************************
    K_CSR_aa(1:K_CSR_NNZ_Max)   = ZR
    K_CSR_ja(1:K_CSR_NNZ_Max)   = 0
    K_CSR_ia(1:num_FreeD+1)   = 0
    
    allocate(K_COO_aa_try(K_CSR_NNZ_Max))
    allocate(K_COO_ia_try(K_CSR_NNZ_Max))
    allocate(K_COO_ja_try(K_CSR_NNZ_Max))

    K_COO_aa_try(1:K_CSR_NNZ_Max)   = ZR
    K_COO_ia_try(1:K_CSR_NNZ_Max)   = 0
    K_COO_ja_try(1:K_CSR_NNZ_Max)   = 0
    K_CSR_NNZ = 0
    K_COO_NNZ_try = 0
    Total_Num_G_P = 0
    int_0 = 0_integer_1
    int_1 = 1_integer_1
    location_0 = 0
    
    

    ! Life-and-death unit preparation
    if(Key_EKILL==1)then
          num_Ele_Killed = count(Ele_Killed_Each_Load_Step(1:isub,:)>0)
    endif

    ! Standard 64 Gauss points
    if (Key_Integral_Sol  == 2)then
        call Cal_Gauss_Points_QUAD(Num_Gauss_Points,kesi_Enr,yita_Enr,weight_Enr)
    ! Subdivision of quadrilaterals to calculate Gauss points (Gauss points are unrelated to cracks,
    ! consisting of several regular quadrilaterals)
    elseif (Key_Integral_Sol  == 3)then
        call Cal_Gauss_Points_QUAD_for_SUBQUAD(Num_Sub_Quads,kesi_Enr,yita_Enr,weight_Enr)
          Num_Gauss_Points = Num_Sub_Quads*4
    endif
    call Cal_Gauss_Points_QUAD(4,kesi_N_Enr,yita_N_Enr,weight_N_Enr)


    !------------------------------------------------------------
    ! OPTION-1: Store Location_COO and Mask_COO in a binary file
    !------------------------------------------------------------
    ! Open the Sparse_K_Location_COO.bin file, which stores the indices of the non-zero elements
    ! corresponding to the stiffness matrix elements.
    if(Sparse_Store_Method==1)then
          !###############
          ! Mask_COO File
          !###############
          if(isub<=1)then
              print *,"    Generating Mask_COO binary file..."
          else
              print *,"    Refreshing Mask_COO binary file..."
          endif
          !***********************************************************************
          ! option 1: Mask_COO file zeroing algorithm, one by one, low efficiency
          !***********************************************************************
    !C         c_File_name_2=
    !C    &          trim(Full_Pathname)//'_Sparse_K_Mask_COO.bin'
    !C         Open(1001,File = c_File_name_2,Access = 'Direct',
    !C    &             Form = 'Unformatted', RecL = 1)
    !C         do i=1,T_Freedom
    !C           do j=1,T_Freedom
    !C             Write(1001, Rec=(i-1)*T_Freedom+j) int_0
    !C           enddo
    !C         enddo
          !***********************************************************************************
          ! Option 2: Mask_COO file zeroing algorithm, assemble the vector when writing, high
          ! efficiency
          !***********************************************************************************
          c_File_name_2=trim(Full_Pathname)//'_Sparse_K_Mask_COO.bin'

    !C         Open(1001,File = c_File_name_2,Access = 'Direct',
    ! C & Form = 'Unformatted', RecL = T_Freedom) !Note that the record length has been changed to
    ! T_Freedom
    !C         allocate(int_0_vector(T_Freedom))
    !C         int_0_vector(1:T_Freedom) = int_0
    !C
    !C         do j=1,T_Freedom
    !C             Write(1001,Rec=j) int_0_vector(1:T_Freedom)
    !C         enddo
    !C         close(1001)
    !C         deallocate(int_0_vector)

          ! ~~~~~~~~~~~~~~~The above code segment can be replaced with the following code to improve
          ! efficiency~~~~~~~~~~~~~~~
          if(isub<=1)then
              Open(1001,File = c_File_name_2,Access = 'Direct',&
                 Form = 'Unformatted', status='unknown',&
                 RecL = T_Freedom*2)
          else
              Open(1001,File = c_File_name_2,Access = 'Direct',&
                 Form = 'Unformatted', status='old',&
                 RecL = T_Freedom*2)
          endif
          allocate(int_0_vector(T_Freedom*2))
          int_0_vector(1:T_Freedom*2) = int_0
          do j=1,T_Freedom/2
              Write(1001,Rec=j) int_0_vector
              ! Progress bar
#ifndef Silverfrost
              if(mod(j,100)==0 .and. j/=T_Freedom/2)then
                  call Progress_Mask_COO %Put(j, CMD_PROGRESS_ABSOLUTE)
              endif
              if(j==T_Freedom/2)then
                  call Progress_Mask_COO %Put(j, CMD_PROGRESS_ABSOLUTE)
              endif
#endif
          enddo
          close(1001)
          deallocate(int_0_vector)
          ! Open again
          Open(1001,File = c_File_name_2,Access = 'Direct',Form = 'Unformatted', RecL = 1)
          !###################
          ! Location_COO File
          !###################
          if(isub<=1)then
              print *,"    Generating Location_COO binary file..."
          else
              print *,"    Refreshing Location_COO binary file..."
          endif
          c_File_name_1=trim(Full_Pathname)//'_Sparse_K_Location_COO.bin'

    !C         Open(1002,File =c_File_name_1,Access='Direct',
    !C    &             Form = 'Unformatted', RecL = 4*T_Freedom)
    !C         allocate(location_vector(T_Freedom))
    !C         location_vector(1:T_Freedom) = location_0
    !C         do j=1,T_Freedom
    ! C Write(1002, Rec=j) location_vector(1:T_Freedom) ! Assign zero value to Location_COO in the
    ! binary file
    !C         enddo
    !C         close(1002)
    !C         deallocate(location_vector)
          ! ~~~~~~~~~~~~~~~The above code segment can be replaced with the following code to improve
          ! efficiency~~~~~~~~~~~~~~~
          if(isub<=1)then
              Open(1002,File =c_File_name_1,Access='Direct',&
                 Form = 'Unformatted',status='unknown',&
                 RecL = 4*T_Freedom*2)
          else
              Open(1002,File =c_File_name_1,Access='Direct',&
                 Form = 'Unformatted',status='old',&
                 RecL = 4*T_Freedom*2)
          endif
          allocate(location_vector(T_Freedom*2))
          location_vector(1:T_Freedom*2) = location_0
          do j=1,T_Freedom/2
              Write(1002, Rec=j)  location_vector
              ! Progress bar
#ifndef Silverfrost
              if(mod(j,100)==0 .and. j/=T_Freedom/2)then
                  call Progress_Location_COO %Put(j, CMD_PROGRESS_ABSOLUTE)
              endif
              if(j==T_Freedom/2)then
                  call Progress_Location_COO %Put(j, CMD_PROGRESS_ABSOLUTE)
              endif
#endif
          enddo
          close(1002)
          deallocate(location_vector)
          ! Open again
          Open(1002,File =c_File_name_1,Access='Direct',Form = 'Unformatted', RecL = 4)

          file_Sparse_K_Location_COO_bin = .True.
          file_Sparse_K_Mask_COO_bin = .True.
    endif

    !-----------------------------------------------------
    ! OPTION-2: Store Location_COO and Mask_COO in memory
    !-----------------------------------------------------
    if(Sparse_Store_Method==2)then
          ! Allocate memory for intermediate variables
          allocate(Location_COO(T_Freedom,T_Freedom))
          Location_COO(1:T_Freedom,1:T_Freedom)= 0
          allocate(Mask_COO(T_Freedom,T_Freedom))
          Mask_COO(1:T_Freedom,1:T_Freedom) = .False.
    endif

    !**************************************************************
    ! Step 1: Obtain the non-zero elements before division by zero
    !**************************************************************
    EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gauss_Points)= .False.
    do i_E = 1,Num_Elem
          c_thick = thick(Elem_Mat(i_E))
          c_D     = D(Elem_Mat(i_E),:,:)
          
          ! Elastic modulus of the Weibull distribution. 2024-06-25. NEWFTU2024062402.
          if(Flag_Weibull_E)then
              if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
                  c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
              endif
          endif
          
          ! Life-and-death unit processing: weakening D matrix (elastic modulus)
          if(Key_EKILL==1 .and. num_Ele_Killed>=1) then
            if(any(Ele_Killed_Each_Load_Step(1:isub,1:num_Ele_Killed)==i_E))then
                c_D= c_D*EKILL_Weaken_Factor
            endif
          endif
          c_NN    = G_NN(:,i_E)
          c_X_NODES = G_X_NODES(:,i_E)
          c_Y_NODES = G_Y_NODES(:,i_E)           
          ! ------------------------------
          ! Point Scheme 1: Triangulation
          ! ------------------------------
          if(Key_Integral_Sol.eq.1)then
              !call Cal_Gauss_Points_Subtriangle(kesi,yita,weight)
              !c_Num_Gauss_Point  = size(kesi,2);
          ! -------------------------------------
          ! Point Plan 2 or 3: Even Distribution
          ! -------------------------------------
          elseif(Key_Integral_Sol.eq.2 .or. Key_Integral_Sol.eq.3)then
              ! Starting number of Gauss points for each unit
              Ele_GP_Start_Num(i_E) = Total_Num_G_P + 1
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Determine the Gauss points of the current element
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !if the current element are enriched element, then 8x8 gauss points is suggested:
              if(maxval(Enriched_Node_Type(c_NN,1:num_Crack)).ne.0)then
                  kesi(1:Num_Gauss_Points)    = kesi_Enr
                  yita(1:Num_Gauss_Points)    = yita_Enr
                  weight(1:Num_Gauss_Points)  = weight_Enr
                  c_Num_Gauss_Point = Num_Gauss_Points
                  Total_Num_G_P     = Total_Num_G_P +c_Num_Gauss_Point
              !if the current element are not enriched element, then 2x2 gauss points:
              else
                  kesi(1:4)    = kesi_N_Enr
                  yita(1:4)    = yita_N_Enr
                  weight(1:4)  = weight_N_Enr
                  c_Num_Gauss_Point = 4
                  Total_Num_G_P     = Total_Num_G_P +c_Num_Gauss_Point
              end if
              ! Save the number of Gauss points for each element to a global variable
              num_GP_Elem(i_E) = c_Num_Gauss_Point
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !Decide the location of each element stiffness
              !matrix in the global stiffness matrix.
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              Location_ESM(1:MDOF_2D)       = 0
              !Location_ESM_NoFEM(1:160) = 0
              num_Loc_ESM       = 0
              !num_Loc_ESM_NoFEM = 0
              do i_C =1,num_Crack
                  call Location_Element_Stiff_Matrix(i_E,i_C,&
                                            c_POS(:,i_C),&
                                            Location_ESM_C_Crack,&
                                            num_Loc_ESM_C_Crack,&
                                            Location_ESM_C_Cr_NoFEM,&
                                            num_Loc_ESM_C_Cr_NoFEM)
                  ! Includes FEM degrees of freedom
                  Location_ESM(num_Loc_ESM+1:&
                          num_Loc_ESM+num_Loc_ESM_C_Crack) =&
                          Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
                  num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
              end do

              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Loop over each Gauss point to assemble the element stiffness matrix.
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              localK(1:MDOF_2D,1:MDOF_2D) = ZR
              do i_G = 1,c_Num_Gauss_Point
                  !B(1:3,1:50) = ZR
                  B(1:3,1:80) = ZR
                  num_B = 0
                  call Cal_detJ(kesi(i_G),yita(i_G),c_X_NODES,c_Y_NODES,detJ)
                  !Calculate the B Matrix, Loop through each crack.
                  do i_C =1,num_Crack
                      call Cal_B_Matrix_Crack(kesi(i_G),yita(i_G),&
                                        i_C,i_E,i_G,&
                                        c_NN,c_X_NODES,c_Y_NODES,&
                                        tem_B,num_tem_B)
                      B(1:3,num_B+1:num_B+num_tem_B) =tem_B(1:3,1:num_tem_B)
                      num_B = num_B + num_tem_B
                  end do
                  localK(1:num_B,1:num_B) = localK(1:num_B,1:num_B) +&
                               c_thick*weight(i_G)*detJ*&
                               MATMUL(MATMUL(transpose(B(1:3,1:num_B)),c_D),B(1:3,1:num_B))
              enddo

              !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              ! OPTION-1: Binary file (Sparse_K_Location_COO.bin) stores Location_COO
              ! and Mask_COO
              !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              if(Sparse_Store_Method==1)then
                do i_row = 1,num_Loc_ESM
                  do i_col = 1,num_Loc_ESM
                      c_row = Location_ESM(i_row)
                      c_col = Location_ESM(i_col)
                      ! If it is in the lower triangle, skip this iteration (for the upper triangle, including the
                      ! diagonal, the condition is: c_col >= c_row)
                      if(c_col<c_row)then
                          cycle
                      endif
                      ! The corresponding degree of freedom needs to be a functional degree of freedom.
    !C                     if(any(freeDOF(1:num_FreeD).eq.c_row)  .and.
    !C    &                   any(freeDOF(1:num_FreeD).eq.c_col)) then

                      !if(all(fixedDOF(1:num_FixedD).ne.c_row)  .and.all(fixedDOF(1:num_FixedD).ne.c_col)) then
                      
                      !IMPROV2024110702.
                      if(Flag_FreeDOF(c_row)==1  .and. Flag_FreeDOF(c_col)==1) then

                         ! Row numbers corresponding to the degrees of freedom after removing constraints
                         !--------------------------
                         ! option 1: Low efficiency
                         !--------------------------
    !C                        real_row = minloc(freeDOF(1:num_FreeD),1,MASK = (freeDOF(1:num_FreeD).eq.c_row))
    !C                        real_col = minloc(freeDOF(1:num_FreeD),1,MASK = (freeDOF(1:num_FreeD).eq.c_col))
                         !---------------------------
                         ! option 2: high efficiency
                         !---------------------------
    ! C call Vector_Location_Int_Ascending_Order( !Replacing Vector_Location_Int with
    ! Vector_Location_Int_Ascending_Order can improve computation efficiency by 38%
    !C    &                            num_FreeD,freeDOF(1:num_FreeD),
    !C    &                                  c_row,real_row,c_Yes_In)
    ! C call Vector_Location_Int_Ascending_Order( !Replacing Vector_Location_Int with
    ! Vector_Location_Int_Ascending_Order can improve computation efficiency by 38%
    !C    &                            num_FreeD,freeDOF(1:num_FreeD),
    !C    &                                  c_col,real_col,c_Yes_In)
                         !------------------------------
                         ! option 3: Most efficient
                         !2024-11-07. IMPROV2024110703.
                         !------------------------------
                         ! Real_FreeDOF_Index is the stiffness matrix degree of freedom index after removing the constrained
                         ! degrees of freedom corresponding to all degrees of freedom. 2024-11-07.
                         real_row =  Real_FreeDOF_Index(c_row) 
                         real_col =  Real_FreeDOF_Index(c_col) 
                        
                         Rec_num = (c_row-1)*T_Freedom+c_col
                         Read(1001, Rec =Rec_num) c_int_0_or_1

                         if (c_int_0_or_1 ==int_0)then
                             ! Upper triangle
                             K_COO_NNZ_try = K_COO_NNZ_try + 1
                             K_COO_ia_try(K_COO_NNZ_try) = real_row
                             K_COO_ja_try(K_COO_NNZ_try) = real_col
                             K_COO_aa_try(K_COO_NNZ_try) = localK(i_row,i_col)
                             Rec_num = (c_row-1)*T_Freedom+c_col
                             Write(1001, Rec=Rec_num) int_1
                             Write(1002, Rec=Rec_num) K_COO_NNZ_try
                             ! Lower triangle
                             if(c_row /= c_col)then
                                 K_COO_NNZ_try = K_COO_NNZ_try + 1
                                 K_COO_ia_try(K_COO_NNZ_try) = real_col
                                 K_COO_ja_try(K_COO_NNZ_try) = real_row
                                 K_COO_aa_try(K_COO_NNZ_try) = localK(i_row,i_col)
                                 Rec_num = (c_col-1)*T_Freedom + c_row
                                 Write(1001, Rec=Rec_num) int_1
                                 Write(1002, Rec=Rec_num) K_COO_NNZ_try
                             endif
                          ! If there is already a contribution in that position of the global stiffness matrix, it is
                          ! accumulated without increasing the number of non-zero elements.
                          elseif(c_int_0_or_1 ==int_1)then
                             ! Upper triangle
                             Rec_num = (c_row-1)*T_Freedom+c_col
                             Read(1002, Rec =Rec_num) c_Location
                             K_COO_aa_try(c_Location)=K_COO_aa_try(c_Location)+ localK(i_row,i_col)
                             ! Lower triangle
                             if(c_row /= c_col)then
                                 Rec_num = (c_col-1)*T_Freedom+c_row
                                 Read(1002, Rec =Rec_num) c_Location
                                 K_COO_aa_try(c_Location)=K_COO_aa_try(c_Location)+ localK(i_row,i_col)
                             endif
                          endif
                      endif
                  end do
                end do
              endif

              !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              ! OPTION-2: Store temporary variables Location_COO and Mask_COO in memory
              !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              if(Sparse_Store_Method==2)then
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ! Extract the non-zero elements of the stiffness matrix.
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do i_row = 1,num_Loc_ESM
                  do i_col = 1,num_Loc_ESM
                      c_row = Location_ESM(i_row)
                      c_col = Location_ESM(i_col)
                      ! If it is in the lower triangle, skip this iteration (for the upper triangle, including the
                      ! diagonal, the condition is: c_col >= c_row)
                      if(c_col<c_row)then
                          cycle
                      endif
                      ! The corresponding degree of freedom needs to be a functional degree of freedom.

                      !if(all(fixedDOF(1:num_FixedD).ne.c_row)  .and.all(fixedDOF(1:num_FixedD).ne.c_col)) then

                      !IMPROV2024110702.
                      if(Flag_FreeDOF(c_row)==1  .and. Flag_FreeDOF(c_col)==1) then

                         ! Row numbers corresponding to the degrees of freedom after removing constraints
                         !--------------------------
                         ! Option 1: Low efficiency
                         !--------------------------
    !C                        real_row = minloc(freeDOF(1:num_FreeD),1,
    !C    &                         MASK = (freeDOF(1:num_FreeD).eq.c_row))
    !C                        real_col = minloc(freeDOF(1:num_FreeD),1,
    !C    &                         MASK = (freeDOF(1:num_FreeD).eq.c_col))
                         !---------------------------
                         ! option 2: high efficiency
                         !---------------------------
    !C                        call Vector_Location_Int(num_FreeD,
    !C    &                                  freeDOF(1:num_FreeD),
    !C    &                                  c_row,real_row,c_Yes_In)
    !C                        call Vector_Location_Int(num_FreeD,
    !C    &                                  freeDOF(1:num_FreeD),
    !C    &                                  c_col,real_col,c_Yes_In)
                         !------------------------------
                         ! option 3: highest efficiency
                         !2024-11-07. IMPROV2024110703.
                         !------------------------------
                         ! Real_FreeDOF_Index is the stiffness matrix degree of freedom index after removing the constrained
                         ! degrees of freedom corresponding to all degrees of freedom. 2024-11-07.
                         real_row =  Real_FreeDOF_Index(c_row) 
                         real_col =  Real_FreeDOF_Index(c_col) 
                         
                         logi_Yes_Mask_COO = Mask_COO(c_row,c_col)
                         if (logi_Yes_Mask_COO .eqv. .False.)then
                             ! Upper triangle
                             K_COO_NNZ_try = K_COO_NNZ_try + 1
                             K_COO_ia_try(K_COO_NNZ_try) = real_row
                             K_COO_ja_try(K_COO_NNZ_try) = real_col
                             K_COO_aa_try(K_COO_NNZ_try) = localK(i_row,i_col)
                             Mask_COO(c_row,c_col) = .True.
                             Location_COO(c_row,c_col) =K_COO_NNZ_try
                             ! Lower triangle
                             if(c_row /= c_col)then
                                 K_COO_NNZ_try = K_COO_NNZ_try + 1
                                 K_COO_ia_try(K_COO_NNZ_try) = real_col
                                 K_COO_ja_try(K_COO_NNZ_try) = real_row
                                 K_COO_aa_try(K_COO_NNZ_try) = localK(i_row,i_col)
                                 Mask_COO(c_col,c_row) = .True.
                                 Location_COO(c_col,c_row)=K_COO_NNZ_try
                             endif
                          ! If there is already a contribution in that position of the global stiffness matrix, it is
                          ! accumulated without increasing the number of non-zero elements.
                          elseif(logi_Yes_Mask_COO .eqv. .True.)then
                             ! Upper triangle
                             c_Location = Location_COO(c_row,c_col)
                             K_COO_aa_try(c_Location)=K_COO_aa_try(c_Location)+ localK(i_row,i_col)
                             ! Lower triangle
                             if(c_row /= c_col)then
                                 c_Location = Location_COO(c_col,c_row)
                                 K_COO_aa_try(c_Location)=K_COO_aa_try(c_Location)+ localK(i_row,i_col)
                             endif
                          endif
                      endif
                  end do
                end do
              endif
          endif
          ! Percentage of Feedback Completed
#ifndef Silverfrost
          if(i_E==1.or.i_E==5.or.i_E==10.or.i_E==20.or.i_E==30)then
              call Progress % Put(i_E, CMD_PROGRESS_ABSOLUTE)
          endif
          if(Mod(i_E,50)==0 .and. (i_E/=Num_Elem))then
              !write(*,2111) DBLE(i_E)/DBLE(Num_Elem)*100.0
              call Progress % Put(i_E, CMD_PROGRESS_ABSOLUTE)
          endif
          if(i_E==Num_Elem)then
              !write(*,2111) DBLE(i_E)/DBLE(Num_Elem)*100.0
              call Progress % Put(i_E, CMD_PROGRESS_ABSOLUTE)
          endif
#endif
    end do

    if(Sparse_Store_Method==1)then
        close( 1001)
        close( 1002)
    endif
    !Clear memory variables
    if(Sparse_Store_Method==2)then
          deallocate(Location_COO)
          deallocate(Mask_COO)
    endif

    !*********************************************************************************************
    ! Step 2: Remove zero elements (including very small values) from the sparse matrix, unsorted
    !*********************************************************************************************

    allocate(K_COO_aa_try2(K_COO_NNZ_try))
    allocate(K_COO_ia_try2(K_COO_NNZ_try))
    allocate(K_COO_ja_try2(K_COO_NNZ_try))
    K_COO_aa_try2(1:K_COO_NNZ_try)   = ZR
    K_COO_ia_try2(1:K_COO_NNZ_try)   = 0
    K_COO_ja_try2(1:K_COO_NNZ_try)   = 0
    K_CSR_NNZ = 0
    do i_delete = 1,K_COO_NNZ_try
          if(abs(K_COO_aa_try(i_delete))> 1.0D-8)then
              K_CSR_NNZ = K_CSR_NNZ + 1
              K_COO_ia_try2(K_CSR_NNZ) = K_COO_ia_try(i_delete)
              K_COO_ja_try2(K_CSR_NNZ) = K_COO_ja_try(i_delete)
              K_COO_aa_try2(K_CSR_NNZ) = K_COO_aa_try(i_delete)
          endif
    enddo

    write(*,1112) DBLE(K_COO_NNZ_try)/T_Freedom/T_Freedom*100.0D0


    !************************************************************
    ! Step 3: Convert to CSR format (SPARSEKIT function package)
    !************************************************************
    call coocsr(num_FreeD,K_CSR_NNZ,K_COO_aa_try2(1:K_CSR_NNZ), &
                                    K_COO_ia_try2(1:K_CSR_NNZ), &
                                    K_COO_ja_try2(1:K_CSR_NNZ),&
                                    K_CSR_aa(1:K_CSR_NNZ),&
                                    K_CSR_ja(1:K_CSR_NNZ),&
                                    K_CSR_ia(1:num_FreeD+1))
    deallocate(K_COO_aa_try2)
    deallocate(K_COO_ia_try2)
    deallocate(K_COO_ja_try2)

    !*************************************************************
    ! Step 4: Sort the CSR sparse matrix (in ascending row order)
    !*************************************************************
    job = 3
    value2 = 1
    call clncsr(job,value2,num_FreeD,K_CSR_aa(1:K_CSR_NNZ), &
                                   K_CSR_ja(1:K_CSR_NNZ), &
                                   K_CSR_ia(1:num_FreeD+1), &
                                   indu,iwk)

    !**************************************************************************
    ! Step 5: Sort the CSR sparse matrix (in ascending order by column number)
    !**************************************************************************
    ALLOCATE(iwork(2*K_CSR_NNZ))
    call csort(num_FreeD,K_CSR_aa(1:K_CSR_NNZ),K_CSR_ja(1:K_CSR_NNZ),&
               K_CSR_ia(1:num_FreeD+1),iwork,values)

endif

!.................................................................
! Handling of non-zero displacement boundary conditions:
! Finite Element Methods in Engineering, Zeng Pan, Fourth Edition
!.................................................................
if(Num_Boux_nonzero >0)then
    print *, '    Error :: Num_Boux_nonzero >0 not avaliable yet when *Key_K_Sparse=1!'
    print *, '    Error :: In Assemble_Stiffness_Matrix_SPARS_XFEM.F90!'
    call Warning_Message('S',Keywords_Blank)
endif
if(Num_Bouy_nonzero >0)then
    print *, '    Error :: Num_Bouy_nonzero >0 not avaliable yet when *Key_K_Sparse=1!'
    print *, '    Error :: In Assemble_Stiffness_Matrix_SPARS_XFEM.F90!'
    call Warning_Message('S',Keywords_Blank)
endif

!...................................................................................
! Penalty stiffness method for node coupling, theoretical basis:
!http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch09.d/IFEM.Ch09.pdf
!...................................................................................
if (num_CP_x_nodes >0)then
    print *, '    Error :: num_CP_x_nodes >0 not avaliable yet when *Key_K_Sparse=1!'
    print *, '    Error :: In Assemble_Stiffness_Matrix_SPARS_XFEM.F90!'
    call Warning_Message('S',Keywords_Blank)
endif
if (num_CP_y_nodes > 0)then
    print *, '    Error :: num_CP_y_nodes >0 not avaliable yet when *Key_K_Sparse=1!'
    print *, '    Error :: In Assemble_Stiffness_Matrix_SPARS_XFEM.F90!'
    call Warning_Message('S',Keywords_Blank)
endif

RETURN
END SUBROUTINE Assemble_Stiffness_Matrix_SPARS_XFEM
