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
 
SUBROUTINE Assemble_Stiffness_Matrix_SPARS_FEM(isub, freeDOF,num_FreeD, &
               K_CSR_aa,K_CSR_ja,K_CSR_ia,K_CSR_NNZ_Max,K_CSR_NNZ,   &
               T_Freedom,Total_Num_G_P)
! Assemble the stiffness matrix.
! Stored in compressed row format.

!***************************
! Global variable reference
!***************************
use Global_Float_Type
use Global_Model
use Global_Filename
use Global_Common
use Global_Material
use module_INTERFACE_coocsr

!**********************
! Variable Declaration
!**********************
implicit none
!include 'omp_lib.h'
integer,intent(in)::isub,num_FreeD,T_Freedom,K_CSR_NNZ_Max
integer,intent(in)::freeDOF(1:num_FreeD)
real(kind=FT),intent(out)::K_CSR_aa(K_CSR_NNZ_Max)
integer,intent(out)::K_CSR_ja(K_CSR_NNZ_Max)
integer,intent(out)::K_CSR_ia(num_FreeD+1)
integer,intent(out)::K_CSR_NNZ
integer,intent(out)::Total_Num_G_P
!-------
integer i_E
real(kind=FT) c_thick,c_D(3,3)
real(kind=FT) c_X_NODES(4),c_Y_NODES(4)
integer c_NN(4) 
real(kind=FT) kesi(Num_Gauss_P_FEM),yita(Num_Gauss_P_FEM),weight(Num_Gauss_P_FEM)             
real(kind=FT) localK(8,8)
integer local(8),i_row,i_col,nIndex
! The following variables are related to the set of sparse matrices
integer i_Check,K_COO_NNZ_try,i_delete 
! integer K_COO_ia_try(K_CSR_NNZ_Max) !The row number in the K matrix corresponding to the active
! degree of freedom of the K_COO_NNZ_try-th nonzero element
! integer K_COO_ja_try(K_CSR_NNZ_Max) ! Column number in K matrix corresponding to the active degree
! of freedom for the K_COO_NNZ_try-th nonzero element

real(kind=FT),ALLOCATABLE::K_COO_aa_try(:)
integer,ALLOCATABLE::K_COO_ia_try(:)
integer,ALLOCATABLE::K_COO_ja_try(:)

logical(KIND=1),ALLOCATABLE::Mask_COO(:,:)
integer c_row,c_col,c_num_NNZ
integer real_row,real_col
real(kind=FT),ALLOCATABLE::  K_COO_aa_try2(:)
integer,ALLOCATABLE::  K_COO_ia_try2(:)
integer,ALLOCATABLE::  K_COO_ja_try2(:)

logical values
integer,ALLOCATABLE:: iwork(:)
integer job,value2
integer indu(num_FreeD),iwk(num_FreeD+1)
integer num_Ele_Killed
integer Scheme 

! 2024-11-07. Scheme = 1 New algorithm related.
real(kind=FT),allocatable::Temp_K_Each(:,:)
real(kind=FT),allocatable::K_Each(:,:)
real(kind=FT),allocatable::Combined_K_Each(:,:)
integer Temp_K_Each_Count,K_Each_Count,num_rows_Combined_K
integer K_Each_NNZ_Max
integer,allocatable:: K_rows(:), K_cols(:)
real(kind=FT),allocatable:: K_vals(:)
integer i
real(FT) :: progress
character(len=52) :: bar
integer :: progress_counter = 0
real(FT) :: progress_percent
integer :: bar_width = 50
integer :: filled
character(len=1) :: spin(4) = ['|','/','-','\']
integer :: spin_idx = 1

!---------------------------------------------------------------------------------------------
1112 FORMAT(5X,'Real density ratio of K is ',F12.4,'%')    

Scheme = 1
            ! =2, old algorithm

!////////////////////////
!                      /
! New Algorithm /
!       2024-11-07     /
!   NEWFTU2024110702.  /
!                      /
!////////////////////////
if (Scheme ==1) then
    !***************************************************************
    ! Step 1: Initialization of some variables and data preparation
    !***************************************************************
    K_CSR_aa(1:K_CSR_NNZ_Max)   = ZR
    K_CSR_ja(1:K_CSR_NNZ_Max)   = 0
    K_CSR_ia(1:num_FreeD+1)     = 0
    K_CSR_NNZ = 0
    
    !K_COO_NNZ_try = 0
    !Mask_COO(1:T_Freedom,1:T_Freedom) = .False.
    Total_Num_G_P = 0
    ! Life-and-death unit preparation
    if(Key_EKILL==1)then
        num_Ele_Killed = count(Ele_Killed_Each_Load_Step(1:isub,:)>0)
    endif

    
    !**************************************
    !               Step 2
    ! Get the number of non-zero elements.
    !**************************************
    Temp_K_Each_Count = 0
    do i_E = 1,Num_Elem
        ! Progress bar.
        !progress = dble(i_E)/Num_Elem
        !bar = '['//repeat('=', int(progress*50))//repeat(' ',50-int(progress*50))//']'
        !write(*, '(a,a,f6.1,a,a,a)', advance='no') char(13), '    ',progress*100, '%  ', trim(bar)
        
        ! Starting number of Gauss points for each unit
        Ele_GP_Start_Num(i_E) = Total_Num_G_P + 1
        Total_Num_G_P = Total_Num_G_P +4
        c_NN    = G_NN(:,i_E)
        !Traditional index locations
        local=[c_NN(1)*2-1,c_NN(1)*2,c_NN(2)*2-1,c_NN(2)*2, &
               c_NN(3)*2-1,c_NN(3)*2,c_NN(4)*2-1,c_NN(4)*2]

        ! Extract the non-zero elements of the stiffness matrix (without summing for now)
        do i_row = 1,8
            do i_col = 1,8
                c_row = local(i_row)
                c_col = local(i_col)
                ! The corresponding degree of freedom needs to be a functional degree of freedom.
                if(Flag_FreeDOF(c_row)==1  .and. Flag_FreeDOF(c_col)==1) then
                    Temp_K_Each_Count =  Temp_K_Each_Count + 1
                endif
            end do
        end do
    end do
    
    
    !********************************************************************************
    !               Step 3
    ! Obtain an unordered COO format sparse matrix (considering boundary conditions)
    ! Save to K_Each(K_CSR_NNZ_Max, 3)
    !********************************************************************************
    call Cal_Gauss_Points_QUAD(4,kesi,yita,weight)
    
    !2024-11-07.
    K_Each_NNZ_Max = Temp_K_Each_Count
    allocate(Temp_K_Each(K_Each_NNZ_Max,3))
    
    Temp_K_Each_Count = 0
    !nIndex = 0
    
    progress_counter = 0
    do i_E = 1,Num_Elem
        ! Progress bar.
        !progress = dble(i_E)/Num_Elem
        !bar = '['//repeat('=', int(progress*50))//repeat(' ',50-int(progress*50))//']'
        !write(*, '(a,a,f6.1,a,a,a)', advance='no') char(13), '    ',progress*100, '%  ', trim(bar)
        
        ! Starting number of Gauss points for each unit
        c_thick = thick(Elem_Mat(i_E))
        c_D     = D(Elem_Mat(i_E),:,:)   

        ! Elastic modulus of the Weibull distribution. 2024-06-25. NEWFTU2024062402.
        if(Flag_Weibull_E)then
            if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
            c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
            endif
        endif

        ! Life and Death Unit Processing: Weakening the D Matrix
        if(Key_EKILL==1 .and. num_Ele_Killed>=1) then
            if(any(Ele_Killed_Each_Load_Step(1:isub,1:num_Ele_Killed)==i_E))then
                c_D= c_D*EKILL_Weaken_Factor
            endif
        endif
        c_NN    = G_NN(:,i_E)
        c_X_NODES = G_X_NODES(:,i_E)
        c_Y_NODES = G_Y_NODES(:,i_E)           
        !Traditional index locations
        local=[c_NN(1)*2-1,c_NN(1)*2,c_NN(2)*2-1,c_NN(2)*2, &
        c_NN(3)*2-1,c_NN(3)*2,c_NN(4)*2-1,c_NN(4)*2]
        !Get the element stiffness matrix of the current element	
        call Cal_Ele_Stiffness_Matrix_N4(c_X_NODES,c_Y_NODES,  &
             c_thick,c_D,kesi,yita,weight,localK)
        ! Extract the non-zero elements of the stiffness matrix (without summing for now)
        do i_row = 1,8
            do i_col = 1,8
                c_row = local(i_row)
                c_col = local(i_col)
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
        
! !Progress bar.
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
    
    !*********************************************************************************************
    !               Step 4
    ! K merges duplicate data (summing them), and this function is very time-consuming. It can be
    ! optimized.
    !*********************************************************************************************
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
    !         Step 5
    ! Regenerate the merged K.
    !**************************
    deallocate(K_Each)
    allocate(K_Each(num_rows_Combined_K,3))
    K_Each(1:num_rows_Combined_K,1:3)  = Combined_K_Each(1:num_rows_Combined_K,1:3)
    deallocate(Combined_K_Each)
    
    
    !***************************************************
    !           Step 6
    ! Convert matrix K into a COO sparse format matrix.
    !***************************************************
    K_CSR_NNZ = 0
    do i=1,num_rows_Combined_K
        if(K_Each(i,3)/=ZR) then
            K_CSR_NNZ = K_CSR_NNZ +1 
        endif      
    enddo 
    
    ! Check whether the number of non-zero elements exceeds the limit.
    if(K_CSR_NNZ > K_CSR_NNZ_Max) then
        print *,'    ERROR :: K_CSR_NNZ > K_CSR_NNZ_Max in Assemble_Stiffness_Matrix_SPARS_FEM.f90!'
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
    !           Step 7
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
    ! Step 1: Initialization of some variables and data preparation
    !***************************************************************
    K_CSR_aa(1:K_CSR_NNZ_Max)   = ZR
    K_CSR_ja(1:K_CSR_NNZ_Max)   = 0
    K_CSR_ia(1:num_FreeD+1)     = 0
    K_CSR_NNZ = 0
    K_COO_NNZ_try = 0
    allocate(Mask_COO(T_Freedom,T_Freedom))
    Mask_COO(1:T_Freedom,1:T_Freedom) = .False.
    
    allocate(K_COO_aa_try(K_CSR_NNZ_Max))
    allocate(K_COO_ia_try(K_CSR_NNZ_Max))
    allocate(K_COO_ja_try(K_CSR_NNZ_Max))
    
    Total_Num_G_P = 0
    ! Life-and-death unit preparation
    if(Key_EKILL==1)then
        num_Ele_Killed = count(Ele_Killed_Each_Load_Step(1:isub,:)>0)
    endif
    
    !****************************************************************************************
    ! Step 2: Obtain the unsorted COO-format sparse matrix (considering boundary conditions)
    !****************************************************************************************
    call Cal_Gauss_Points_QUAD(4,kesi,yita,weight)
    nIndex = 0
    do i_E = 1,Num_Elem
        ! Starting number of Gauss points for each unit
        Ele_GP_Start_Num(i_E) = Total_Num_G_P + 1
        Total_Num_G_P = Total_Num_G_P +4
        c_thick = thick(Elem_Mat(i_E))
        c_D     = D(Elem_Mat(i_E),:,:)   

        ! Elastic modulus of the Weibull distribution. 2024-06-25. NEWFTU2024062402.
        if(Flag_Weibull_E)then
            if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
            c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
            endif
        endif

        ! Life and Death Unit Processing: Weakening the D Matrix
        if(Key_EKILL==1 .and. num_Ele_Killed>=1) then
            if(any(Ele_Killed_Each_Load_Step(1:isub,1:num_Ele_Killed)==i_E))then
                c_D= c_D*EKILL_Weaken_Factor
            endif
        endif
        !C         c_NN    = G_NN(i_E,:)
        !C         c_X_NODES = G_X_NODES(i_E,:)
        !C         c_Y_NODES = G_Y_NODES(i_E,:)  
        c_NN    = G_NN(:,i_E)
        c_X_NODES = G_X_NODES(:,i_E)
        c_Y_NODES = G_Y_NODES(:,i_E)           
        !Traditional index locations
        local=[c_NN(1)*2-1,c_NN(1)*2,c_NN(2)*2-1,c_NN(2)*2, &
        c_NN(3)*2-1,c_NN(3)*2,c_NN(4)*2-1,c_NN(4)*2]
        !Get the element stiffness matrix of the current element	
        call Cal_Ele_Stiffness_Matrix_N4(c_X_NODES,c_Y_NODES,  &
             c_thick,c_D,kesi,yita,weight,localK)
        ! Extract the non-zero elements of the stiffness matrix (without summing for now)
        do i_row = 1,8
            do i_col = 1,8
                c_row = local(i_row)
                c_col = local(i_col)
                ! The corresponding degree of freedom needs to be a functional degree of freedom.
                !if(any(freeDOF(1:num_FreeD).eq.c_row)  .and. 
                !   any(freeDOF(1:num_FreeD).eq.c_col)) then

                !IMPROV2024110701.
                if(Flag_FreeDOF(c_row)==1  .and. Flag_FreeDOF(c_col)==1) then
                    ! The corresponding row number after removing the constrained degrees of freedom.
                    ! OPTION-1: Low efficiency.
                    !real_row = minloc(freeDOF(1:num_FreeD),1,MASK = (freeDOF(1:num_FreeD).eq.c_row))
                    !real_col = minloc(freeDOF(1:num_FreeD),1,MASK = (freeDOF(1:num_FreeD).eq.c_col))
                    
                    !OPTION 2: 2024-11-07. IMPROV2024110703.
                    ! Real_FreeDOF_Index is the stiffness matrix degree of freedom index after removing the constrained
                    ! degrees of freedom corresponding to all degrees of freedom. 2024-11-07.
                    real_row =  Real_FreeDOF_Index(c_row) 
                    real_col =  Real_FreeDOF_Index(c_col) 
                    ! If there has been no contribution at this position in the overall stiffness matrix, then add a
                    ! non-zero element.
                    if (.not. Mask_COO(c_row,c_col))then
                        K_COO_NNZ_try = K_COO_NNZ_try + 1
                        K_COO_ia_try(K_COO_NNZ_try) = real_row
                        K_COO_ja_try(K_COO_NNZ_try) = real_col 
                        K_COO_aa_try(K_COO_NNZ_try) = localK(i_row,i_col)
                        Mask_COO(c_row,c_col) = .True.
                    ! If there is already a contribution in that position of the global stiffness matrix, it is
                    ! accumulated without increasing the number of non-zero elements.
                    elseif(Mask_COO(c_row,c_col))then
                        ! Extract non-zero element number c_num_NNZ
                        do i_Check=1,K_COO_NNZ_try
                            if (K_COO_ia_try(i_Check)==real_row)then
                                if(K_COO_ja_try(i_Check)==real_col)then
                                    c_num_NNZ = i_Check    
                                    exit
                                endif
                            endif
                        enddo  
                        K_COO_aa_try(c_num_NNZ)=K_COO_aa_try(c_num_NNZ) + localK(i_row,i_col)
                    endif
                endif
            end do
        end do
    end do

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
    if(abs(K_COO_aa_try(i_delete))> 1.0D-5)then
          K_CSR_NNZ = K_CSR_NNZ + 1
          K_COO_ia_try2(K_CSR_NNZ) = K_COO_ia_try(i_delete)
          K_COO_ja_try2(K_CSR_NNZ) = K_COO_ja_try(i_delete)
          K_COO_aa_try2(K_CSR_NNZ) = K_COO_aa_try(i_delete)
    endif
    enddo
    
    write(*,1112) DBLE(K_CSR_NNZ)/T_Freedom/T_Freedom*100.0D0

    !************************************************************
    ! Step 3: Convert to CSR format (SPARSEKIT function package)
    !************************************************************
    call coocsr(num_FreeD,K_CSR_NNZ,K_COO_aa_try2(1:K_CSR_NNZ), &
                                  K_COO_ia_try2(1:K_CSR_NNZ), &
                                  K_COO_ja_try2(1:K_CSR_NNZ), &
                                  K_CSR_aa(1:K_CSR_NNZ),      &
                                  K_CSR_ja(1:K_CSR_NNZ),      &
                                  K_CSR_ia(1:num_FreeD+1))

    !*************************************************************
    ! Step 4: Sort the CSR sparse matrix (in ascending row order)
    !*************************************************************
    job = 3
    value2 = 1   
    call clncsr(job,value2,num_FreeD,K_CSR_aa(1:K_CSR_NNZ), &
                                 K_CSR_ja(1:K_CSR_NNZ),   &
                                 K_CSR_ia(1:num_FreeD+1), &
                                 indu,iwk)

    !**************************************************************************
    ! Step 5: Sort the CSR sparse matrix (in ascending order by column number)
    !**************************************************************************
    ALLOCATE(iwork(2*K_CSR_NNZ))
    call csort(num_FreeD,K_CSR_aa(1:K_CSR_NNZ),                &
                       K_CSR_ja(1:K_CSR_NNZ),                &
                       K_CSR_ia(1:num_FreeD+1),iwork,values)

    !****************
    ! Step 6: Others
    !****************
    ! Save the number of Gauss points for each element to a global variable
    num_GP_Elem(1:num_Elem) = Num_Gauss_P_FEM

endif

!.................................................................
! Handling of non-zero displacement boundary conditions:
! Finite Element Methods in Engineering, Zeng Pan, Fourth Edition
!.................................................................
if(Num_Boux_nonzero >0)then
    print *, '    Error :: Num_Boux_nonzero >0 not avaliable yet when *Key_K_Sparse=1!'
    print *, '    Error :: In Assemble_Stiffness_Matrix_SPARS_FEM.F90!'
    call Warning_Message('S',Keywords_Blank)
endif
if(Num_Bouy_nonzero >0)then
    print *, '    Error :: Num_Bouy_nonzero >0 not avaliable yet when *Key_K_Sparse=1!'
    print *, '    Error :: In Assemble_Stiffness_Matrix_SPARS_FEM.F90!'
    call Warning_Message('S',Keywords_Blank)
endif

!...................................................................................
! Penalty stiffness method for node coupling, theoretical basis:
!http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch09.d/IFEM.Ch09.pdf
!...................................................................................
if (num_CP_x_nodes >0)then
    print *, '    Error :: num_CP_x_nodes >0 not avaliable yet when *Key_K_Sparse=1!'
    print *, '    Error :: In Assemble_Stiffness_Matrix_SPARS_FEM.F90!'
    call Warning_Message('S',Keywords_Blank)
endif
if (num_CP_y_nodes > 0)then
    print *, '    Error :: num_CP_y_nodes >0 not avaliable yet when *Key_K_Sparse=1!'
    print *, '    Error :: In Assemble_Stiffness_Matrix_SPARS_FEM.F90!'
    call Warning_Message('S',Keywords_Blank)
endif

RETURN
END SUBROUTINE Assemble_Stiffness_Matrix_SPARS_FEM
