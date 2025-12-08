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
 
SUBROUTINE Assemble_Stiffness_Matrix_SPARS_XFEM_3D(isub, &
                     freeDOF,num_FreeD,&
                     fixedDOF,num_FixedD,&
                     K_CSR_aa,K_CSR_ja,K_CSR_ia,&
                     K_CSR_NNZ_Max,K_CSR_NNZ,&
                     T_Freedom,Total_Num_G_P)
! Assemble the stiffness matrix for 3D problems, including problems with cracks (store K in COO
! sparse matrix format).

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
use Global_POST,only:Key_Save_Nothing,Key_Post_Elements_Gauss_Num
use Global_Crack_3D,only:c_POS_3D,Enriched_Node_Type_3D,Elem_num_Related_Cracks

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
integer i_C,i_E,i_G
real(kind=FT) c_D(6,6)
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
integer c_NN(8)  
real(kind=FT) kesi_Enr(Num_Gau_Points_3D),yita_Enr(Num_Gau_Points_3D),&
              zeta_Enr(Num_Gau_Points_3D), weight_Enr(Num_Gau_Points_3D)
real(kind=FT) kesi_No_Enr(Num_Gauss_P_FEM_3D),yita_No_Enr(Num_Gauss_P_FEM_3D),&
              zeta_No_Enr(Num_Gauss_P_FEM_3D),weight_No_Enr(Num_Gauss_P_FEM_3D)
real(kind=FT),ALLOCATABLE::kesi_Enr_Cubes(:),yita_Enr_Cubes(:),zeta_Enr_Cubes(:),weight_Enr_Cubes(:)
real(kind=FT)    kesi_Enr_MC(Num_Gau_Points_3D_MC), yita_Enr_MC(Num_Gau_Points_3D_MC),& 
                 zeta_Enr_MC(Num_Gau_Points_3D_MC), weight_Enr_MC(Num_Gau_Points_3D_MC)
real(kind=FT) ,ALLOCATABLE::kesi(:),yita(:),zeta(:),weight(:)                 
integer:: Location_ESM(MDOF_3D)
integer::Location_ESM_C_Crack(MDOF_3D),Location_ESM_C_Cr_NoFEM(MDOF_3D)
integer num_Loc_ESM_C_Cr_NoFEM
real(kind=FT) B(6,MDOF_3D),tem_B(6,MDOF_3D)
integer num_B,num_tem_B
integer num_Loc_ESM,num_Loc_ESM_C_Crack
integer c_Num_Gauss_Point
real(kind=FT) detJ
real(kind=FT) Rot_c_D_Comp(6,6),c_D_Comp(6,6),Volume_Ratio
real(kind=FT) T_Matrix(6,6),TT_Matrix(6,6)  
integer mat_num
real(kind=FT) localK(MDOF_3D,MDOF_3D)
integer c_POS_3D_c_Ele(8)
integer c_row,c_col
integer real_row,real_col
integer i
integer k_each_nnz_max
integer, parameter :: integer_1 = selected_int_kind(1)
integer i_row,i_col
real(kind=FT),allocatable::Temp_K_Each(:,:)
real(kind=FT),allocatable::K_Each(:,:)
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


!***************************************************************
! Step 0: Initialization of some variables and data preparation
!***************************************************************
K_CSR_aa(1:K_CSR_NNZ_Max)   = ZR
K_CSR_ja(1:K_CSR_NNZ_Max)   = 0
K_CSR_ia(1:num_FreeD+1)   = 0
K_CSR_NNZ = 0
Total_Num_G_P  = 0

call Cal_Gauss_Points_3D_8nodes(Num_Gauss_P_FEM_3D,kesi_No_Enr,yita_No_Enr,zeta_No_Enr,weight_No_Enr)

! If it is a fixed-point calculation.
if (Key_Integral_Sol  == 2) then    
    call Cal_Gauss_Points_3D_8nodes(Num_Gau_Points_3D,kesi_Enr,yita_Enr,zeta_Enr,weight_Enr)
    call Cal_Gauss_Points_3D_8nodes(Num_Gau_Points_3D_MC,kesi_Enr_MC,yita_Enr_MC,zeta_Enr_MC,weight_Enr_MC)
! If it is 3D fixed block integration. 2022-07-29. NEWFTU2022072901.
elseif (Key_Integral_Sol  == 3) then     
    allocate(kesi_Enr_Cubes(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
    allocate(yita_Enr_Cubes(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
    allocate(zeta_Enr_Cubes(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
    allocate(weight_Enr_Cubes(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
    call Cal_Gauss_Points_3D_for_SUBCUBES(Num_Sub_3D_Cubes,Num_Gau_Points_3D_Cube,kesi_Enr_Cubes,&
                                      yita_Enr_Cubes,zeta_Enr_Cubes,weight_Enr_Cubes)                                
endif

!********************************************************
! Step 1: Obtain the maximum number of non-zero elements
!********************************************************
EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gau_Points_3D)= .False.

Temp_K_Each_Count = 0

do i_E = 1,Num_Elem
    c_NN    = G_NN(:,i_E)
    select case(Key_Integral_Sol)
    !////////////////
    ! Fixed Points /
    !////////////////
    case(2)
    ! For the enhancement element
    if (num_Crack/= 0 .and. &
        (maxval(Enriched_Node_Type_3D(c_NN,1:num_Crack)).ne.0))then
          c_Num_Gauss_Point            = Num_Gau_Points_3D
          Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
    else 
          c_Num_Gauss_Point            = Num_Gauss_P_FEM_3D
          Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
    end if
                
    ! For Integration Scheme 2: If the element is related to multiple cracks, more Gauss points are
    ! used. 2022-07-16. NEWFTU2022071604.
    if(Elem_num_Related_Cracks(i_E)>=2) then  
          c_Num_Gauss_Point               = Num_Gau_Points_3D_MC
          Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
    endif
    !///////////////////////////
    ! Fixed-block integration /
    !///////////////////////////
    !2022-07-29. NEWFTU2022072901.
    case(3)
    ! For enhancement elements
    if (num_Crack/= 0 .and. (maxval(Enriched_Node_Type_3D(c_NN,1:num_Crack)).ne.0))then
          c_Num_Gauss_Point            = Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube
          Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
    else 
          c_Num_Gauss_Point            = Num_Gauss_P_FEM_3D
          Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
    end if

    !////////////////////////////////////////////////
    ! Integration by parts (temporarily abandoned) /
    !////////////////////////////////////////////////
    !2022-07-27.
    case(4)
    !To Be Done.
    end select  

    ! Save the number of Gauss integration points for each element to a global variable for
    ! post-processing. 2022-07-16.
    if(Key_Save_Nothing==0) then
        if(Key_Post_Elements_Gauss_Num==1) then
            Elements_Gauss_Num(i_E)      =  c_Num_Gauss_Point
        endif
    endif

    ! Save the number of Gauss points for each element to a global variable
    num_GP_Elem(i_E) = c_Num_Gauss_Point


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !Decide the location of each element stiffness 
    !matrix in the global stiffness matrix.   
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Location_ESM(1:MDOF_3D)   = 0
    num_Loc_ESM           = 0
    if(num_Crack/=0)then
        do i_C =1,num_Crack
            c_POS_3D_c_Ele(1:8) = c_POS_3D(c_NN,i_C)
            if(i_C >1 .and. sum(c_POS_3D_c_Ele(1:8))==0) cycle
            call Location_Element_Stiff_Matrix_3D(i_E,i_C,c_POS_3D_c_Ele(1:8),Location_ESM_C_Crack, &
                                                    num_Loc_ESM_C_Crack,Location_ESM_C_Cr_NoFEM,num_Loc_ESM_C_Cr_NoFEM)
            ! Includes FEM degrees of freedom
            Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
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
end do

!**************************************************************
! Step 2: Obtain the non-zero elements of the stiffness matrix
!**************************************************************
K_Each_NNZ_Max = Temp_K_Each_Count
allocate(Temp_K_Each(K_Each_NNZ_Max,3))
EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gau_Points_3D)= .False.

Temp_K_Each_Count = 0
progress_counter = 0
do i_E = 1,Num_Elem
    mat_num = Elem_Mat(i_E)
    localK(1:MDOF_3D,1:MDOF_3D) = ZR
    c_D(1:6,1:6)     = D(Elem_Mat(i_E),1:6,1:6)     
    !If the mat is composite material.
    if (Material_Type(mat_num)==5)then
        Volume_Ratio = Material_Para_Added(mat_num,10)
        c_D_Comp = D_Comp(mat_num,1:6,1:6)
        T_Matrix = Ele_ComMat_RotMatrix(i_E,1:6,1:6)
        TT_Matrix= TRANSPOSE(T_Matrix)
        Rot_c_D_Comp = MATMUL(TT_Matrix,c_D_Comp)
        Rot_c_D_Comp = MATMUL(Rot_c_D_Comp,T_Matrix)
        c_D =(ONE-Volume_Ratio)*c_D + Volume_Ratio*Rot_c_D_Comp
    endif           
    c_NN    = G_NN(1:8,i_E)
    c_X_NODES = G_X_NODES(1:8,i_E)
    c_Y_NODES = G_Y_NODES(1:8,i_E)    
    c_Z_NODES = G_Z_NODES(1:8,i_E) 
      
    select case(Key_Integral_Sol)
    !////////////////
    ! Fixed Points /
    !////////////////
    case(2)
        ! For enhancement elements
        if (num_Crack/= 0 .and. (maxval(Enriched_Node_Type_3D(c_NN,1:num_Crack)).ne.0))then
            if(allocated(kesi)) deallocate(kesi)
            if(allocated(yita)) deallocate(yita)
            if(allocated(zeta)) deallocate(zeta)
            if(allocated(weight)) deallocate(weight)              
            allocate(kesi(Num_Gau_Points_3D))
            allocate(yita(Num_Gau_Points_3D))
            allocate(zeta(Num_Gau_Points_3D))
            allocate(weight(Num_Gau_Points_3D))
            kesi(1:Num_Gau_Points_3D)    = kesi_Enr
            yita(1:Num_Gau_Points_3D)    = yita_Enr
            zeta(1:Num_Gau_Points_3D)    = zeta_Enr
            weight(1:Num_Gau_Points_3D)  = weight_Enr
            c_Num_Gauss_Point            = Num_Gau_Points_3D
            Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
            else 
            if(allocated(kesi)) deallocate(kesi)
            if(allocated(yita)) deallocate(yita)
            if(allocated(zeta)) deallocate(zeta)
            if(allocated(weight)) deallocate(weight)
            allocate(kesi(Num_Gauss_P_FEM_3D))
            allocate(yita(Num_Gauss_P_FEM_3D))
            allocate(zeta(Num_Gauss_P_FEM_3D))
            allocate(weight(Num_Gauss_P_FEM_3D))              
            kesi(1:Num_Gauss_P_FEM_3D)   = kesi_No_Enr
            yita(1:Num_Gauss_P_FEM_3D)   = yita_No_Enr
            zeta(1:Num_Gauss_P_FEM_3D)   = zeta_No_Enr
            weight(1:Num_Gauss_P_FEM_3D) = weight_No_Enr
            c_Num_Gauss_Point            = Num_Gauss_P_FEM_3D
            Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
        end if
            
        ! For Integration Scheme 2: If the element is related to multiple cracks, more Gauss points are
        ! used. 2022-07-16. NEWFTU2022071604.
        if(Elem_num_Related_Cracks(i_E)>=2) then  
            if(allocated(kesi)) deallocate(kesi)
            if(allocated(yita)) deallocate(yita)
            if(allocated(zeta)) deallocate(zeta)
            if(allocated(weight)) deallocate(weight)
            allocate(kesi(Num_Gau_Points_3D_MC))
            allocate(yita(Num_Gau_Points_3D_MC))
            allocate(zeta(Num_Gau_Points_3D_MC))
            allocate(weight(Num_Gau_Points_3D_MC))    
            kesi(1:Num_Gau_Points_3D_MC)    = kesi_Enr_MC
            yita(1:Num_Gau_Points_3D_MC)    = yita_Enr_MC
            zeta(1:Num_Gau_Points_3D_MC)    = zeta_Enr_MC
            weight(1:Num_Gau_Points_3D_MC)  = weight_Enr_MC
            c_Num_Gauss_Point               = Num_Gau_Points_3D_MC
            Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
        endif
    !///////////////////////////
    ! Fixed-block integration /
    !///////////////////////////
    !2022-07-29. NEWFTU2022072901.
    case(3)
        ! For enhancement elements
        if (num_Crack/= 0 .and. (maxval(Enriched_Node_Type_3D(c_NN,1:num_Crack)).ne.0))then
            if(allocated(kesi)) deallocate(kesi) 
            if(allocated(yita)) deallocate(yita)
            if(allocated(zeta)) deallocate(zeta)
            if(allocated(weight)) deallocate(weight)              
            allocate(kesi(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
            allocate(yita(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
            allocate(zeta(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
            allocate(weight(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
            kesi(1:Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube)    = kesi_Enr_Cubes
            yita(1:Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube)    = yita_Enr_Cubes
            zeta(1:Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube)    = zeta_Enr_Cubes
            weight(1:Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube)  = weight_Enr_Cubes
            c_Num_Gauss_Point            = Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube
            Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
        else 
            if(allocated(kesi)) deallocate(kesi)  
            if(allocated(yita)) deallocate(yita)
            if(allocated(zeta)) deallocate(zeta)
            if(allocated(weight)) deallocate(weight)
            allocate(kesi(Num_Gauss_P_FEM_3D))
            allocate(yita(Num_Gauss_P_FEM_3D))
            allocate(zeta(Num_Gauss_P_FEM_3D))
            allocate(weight(Num_Gauss_P_FEM_3D))              
            kesi(1:Num_Gauss_P_FEM_3D)   = kesi_No_Enr
            yita(1:Num_Gauss_P_FEM_3D)   = yita_No_Enr
            zeta(1:Num_Gauss_P_FEM_3D)   = zeta_No_Enr
            weight(1:Num_Gauss_P_FEM_3D) = weight_No_Enr
            c_Num_Gauss_Point            = Num_Gauss_P_FEM_3D
            Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
        end if
    !////////////////////////////////////////////////
    ! Integration by parts (temporarily abandoned) /
    !////////////////////////////////////////////////
    !2022-07-27.
    case(4)
        !To Be Done.
    end select      

    ! Save the number of Gauss points for each element to a global variable
    num_GP_Elem(i_E) = c_Num_Gauss_Point

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !Decide the location of each element stiffness 
    !matrix in the global stiffness matrix.   
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Location_ESM(1:MDOF_3D)   = 0
    num_Loc_ESM           = 0
    if(num_Crack/=0)then
        do i_C =1,num_Crack
            c_POS_3D_c_Ele(1:8) = c_POS_3D(c_NN,i_C)
            if(i_C >1 .and. sum(c_POS_3D_c_Ele(1:8))==0) cycle
            call Location_Element_Stiff_Matrix_3D(i_E,i_C,c_POS_3D_c_Ele(1:8),Location_ESM_C_Crack, &
                                                    num_Loc_ESM_C_Crack,Location_ESM_C_Cr_NoFEM,num_Loc_ESM_C_Cr_NoFEM)
            ! Includes FEM degrees of freedom
            Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
            num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
        end do
    endif
      
      
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Loop through each Gauss point to assemble the stiffness matrix of the current element
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    do i_G = 1,c_Num_Gauss_Point
        B(1:6,1:MDOF_3D) = ZR
        num_B = 0
        call Cal_detJ_3D(kesi(i_G),yita(i_G),zeta(i_G),c_X_NODES,c_Y_NODES,c_Z_NODES,detJ)  
        !Calculate the B Matrix, Loop through each crack.
        if(num_Crack/=0)then
            do i_C =1,num_Crack 
                c_POS_3D_c_Ele(1:8) = c_POS_3D(c_NN,i_C)
                if(i_C >1 .and. sum(c_POS_3D_c_Ele(1:8))==0) cycle
                call Cal_B_Matrix_Crack_3D(kesi(i_G),yita(i_G),zeta(i_G),i_C,i_E,i_G,  & 
                         c_NN,c_X_NODES,c_Y_NODES,c_Z_NODES,tem_B,num_tem_B)
                B(1:6,num_B+1:num_B+num_tem_B) = tem_B(1:6,1:num_tem_B)
                num_B = num_B + num_tem_B
            end do
        endif
        localK(1:num_B,1:num_B) = localK(1:num_B,1:num_B) + weight(i_G)*detJ*   &
               MATMUL(MATMUL(transpose(B(1:6,1:num_B)),c_D),B(1:6,1:num_B))
    end do
      
      
      
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Extract the non-zero elements of the stiffness matrix.
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    do i_row = 1,num_Loc_ESM
        do i_col = 1,num_Loc_ESM
        c_row = Location_ESM(i_row)
        c_col = Location_ESM(i_col)
        !IMPROV2024110702.
        if(Flag_FreeDOF(c_row)==1  .and. Flag_FreeDOF(c_col)==1) then
            ! Row numbers corresponding to the degrees of freedom after removing constraints
            ! Location_FreeDOF is the stiffness matrix degree of freedom index corresponding to all degrees of
            ! freedom after removing the constrained degrees of freedom.
            ! The function of Location_FreeDOF is the same as Real_FreeDOF_Index in 2D problems. 2024-11-15.
            real_row =  Location_FreeDOF(c_row) 
            real_col =  Location_FreeDOF(c_col) 
            Temp_K_Each_Count =  Temp_K_Each_Count + 1

            Temp_K_Each(Temp_K_Each_Count,1) = real(real_row)
            Temp_K_Each(Temp_K_Each_Count,2) = real(real_col)
            Temp_K_Each(Temp_K_Each_Count,3) = localK(i_row,i_col)
            endif
        end do
    end do
end do


!*******************************
!               Step 3
! K merge duplicate data (sum).
!*******************************
K_Each_Count = Temp_K_Each_Count
allocate(K_Each(K_Each_Count,3))
K_Each(1:K_Each_Count,1:3) = Temp_K_Each(1:K_Each_Count,1:3)
deallocate(Temp_K_Each)
allocate(Combined_K_Each(K_Each_Count,3))
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


!*******************************************************
!           Step 5
! Convert the matrix K into a COO sparse matrix format.
!*******************************************************
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



!----------------------------------------------------------------------------------------------------
! Penalty function method for handling compression-shear fracture surfaces, 2022-05-12.
! NEWFTU2022051201.
! Ref: \theory_documents\034 Efficient Algorithm for 3D Compression-Shear Fracture Processing Based
! on Penalty Function Method - Faculty Visit - 2022-05-11.pdf
! Ref: \theory_documents\033 Three-Point Constraint Penalty Function Method - Equations 9.8-9.9-9.10
! - 2022-05-11.pdf
!
! Note, the types of enhanced nodes are not distinguished here, which may cause issues. For example,
! there may be multiple tip-enriched elements. Refer to EBE_XFEM_PCG_3D_Get_K.f90. 2023-08-27.
!
!----------------------------------------------------------------------------------------------------
if (Key_CS_Natural_Crack==1) then
    print *, '    Error :: Key_CS_Natural_Crack=1 not avaliable yet when *Key_K_Sparse=1!'
    print *, '    Error :: In Assemble_Stiffness_Matrix_SPARS_XFEM_3D.F90!'
    call Warning_Message('S',Keywords_Blank)
endif

RETURN
END SUBROUTINE Assemble_Stiffness_Matrix_SPARS_XFEM_3D
