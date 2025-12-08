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
 
SUBROUTINE Assemble_Stiffness_Matrix_SPARS_FEM_3D(isub,freeDOF,num_FreeD, &
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
integer,intent(in)::isub,num_FreeD,T_Freedom
integer(kind=LIT),intent(in)::K_CSR_NNZ_Max
integer,intent(in)::freeDOF(1:num_FreeD)
real(kind=FT),intent(out)::K_CSR_aa(K_CSR_NNZ_Max)
integer,intent(out)::K_CSR_ja(K_CSR_NNZ_Max)
integer,intent(out)::K_CSR_ia(num_FreeD+1)
integer,intent(out)::K_CSR_NNZ
integer,intent(out)::Total_Num_G_P
integer i_E   
real(kind=FT) c_D(6,6)
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
integer c_NN(8)  
real(kind=FT) kesi(Num_Gauss_P_FEM_3D),yita(Num_Gauss_P_FEM_3D)
real(kind=FT) zeta(Num_Gauss_P_FEM_3D),weight(Num_Gauss_P_FEM_3D)     
real(kind=FT) localK(24,24)     
integer local(24),i_row,i_col
integer mat_num
integer c_row,c_col
integer real_row,real_col
real(kind=FT),allocatable::Temp_K_Each(:,:)
real(kind=FT),allocatable::K_Each(:,:)
real(kind=FT),allocatable::Combined_K_Each(:,:)
integer Temp_K_Each_Count,K_Each_Count,num_rows_Combined_K
integer K_Each_NNZ_Max
integer,allocatable:: K_rows(:), K_cols(:)
real(kind=FT),allocatable:: K_vals(:)
real(kind=FT) Rot_c_D_Comp(6,6),c_D_Comp(6,6),Volume_Ratio
real(kind=FT) T_Matrix(6,6),TT_Matrix(6,6)
integer i
integer :: progress_counter = 0
real(8) :: progress_percent
integer :: bar_width = 50
integer :: filled
character(len=1) :: spin(4) = ['|','/','-','\']
integer :: spin_idx = 1

1112 FORMAT(5X,'Real density ratio of K is ',F12.4,'%')    


!***************************************************************
! Step 1: Initialization of some variables and data preparation
!***************************************************************
K_CSR_aa(1:K_CSR_NNZ_Max)   = ZR
K_CSR_ja(1:K_CSR_NNZ_Max)   = 0
K_CSR_ia(1:num_FreeD+1)     = 0
K_CSR_NNZ = 0

Total_Num_G_P = 0

!*************************************************
! Step 2: Obtain the number of non-zero elements.
!*************************************************
Temp_K_Each_Count = 0
do i_E = 1,Num_Elem
    ! Starting number of Gauss points for each element
    Ele_GP_Start_Num(i_E) = Total_Num_G_P + 1
    Total_Num_G_P = Total_Num_G_P + Num_Gauss_P_FEM_3D
    c_NN    = G_NN(:,i_E)
    !Traditional index locations
    local=[c_NN(1)*3-2,c_NN(1)*3-1,c_NN(1)*3, &
           c_NN(2)*3-2,c_NN(2)*3-1,c_NN(2)*3, &
           c_NN(3)*3-2,c_NN(3)*3-1,c_NN(3)*3, &
           c_NN(4)*3-2,c_NN(4)*3-1,c_NN(4)*3, &
           c_NN(5)*3-2,c_NN(5)*3-1,c_NN(5)*3, &
           c_NN(6)*3-2,c_NN(6)*3-1,c_NN(6)*3, &
           c_NN(7)*3-2,c_NN(7)*3-1,c_NN(7)*3, &
           c_NN(8)*3-2,c_NN(8)*3-1,c_NN(8)*3]   
    ! Extract the non-zero elements of the stiffness matrix (without summing for now)
    do i_row = 1,24
        do i_col = 1,24
            c_row = local(i_row)
            c_col = local(i_col)
            ! The corresponding degree of freedom needs to be a functional degree of freedom.
            if(Flag_FreeDOF(c_row)==1  .and. Flag_FreeDOF(c_col)==1) then
                ! Location_FreeDOF is the stiffness matrix degree of freedom index after removing the constrained
                ! degrees of freedom corresponding to all degrees of freedom.
                ! The function of Location_FreeDOF is the same as Real_FreeDOF_Index in 2D problems. 2024-11-15.
                real_row =  Location_FreeDOF(c_row) 
                real_col =  Location_FreeDOF(c_col) 
                Temp_K_Each_Count =  Temp_K_Each_Count + 1
            endif
        end do
    end do
end do


!********************************************************************************
!               Step 3:
! Obtain an unordered COO format sparse matrix (considering boundary conditions)
! Save to K_Each(K_CSR_NNZ_Max, 3)
!********************************************************************************
call Cal_Gauss_Points_3D_8nodes(Num_Gauss_P_FEM_3D,kesi,yita,zeta,weight)

!2024-11-07.
K_Each_NNZ_Max = Temp_K_Each_Count
allocate(Temp_K_Each(K_Each_NNZ_Max,3))

Temp_K_Each_Count = 0

progress_counter = 0
do i_E = 1,Num_Elem
    mat_num = Elem_Mat(i_E)
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
    c_NN    = G_NN(:,i_E)
    c_X_NODES = G_X_NODES(:,i_E)
    c_Y_NODES = G_Y_NODES(:,i_E)    
    c_Z_NODES = G_Z_NODES(:,i_E)        
    !Traditional index locations
    local=[c_NN(1)*3-2,c_NN(1)*3-1,c_NN(1)*3, &
           c_NN(2)*3-2,c_NN(2)*3-1,c_NN(2)*3, &
           c_NN(3)*3-2,c_NN(3)*3-1,c_NN(3)*3, &
           c_NN(4)*3-2,c_NN(4)*3-1,c_NN(4)*3, &
           c_NN(5)*3-2,c_NN(5)*3-1,c_NN(5)*3, &
           c_NN(6)*3-2,c_NN(6)*3-1,c_NN(6)*3, &
           c_NN(7)*3-2,c_NN(7)*3-1,c_NN(7)*3, &
           c_NN(8)*3-2,c_NN(8)*3-1,c_NN(8)*3]   
    !Get the element stiffness matrix of the current element	
    call Cal_Ele_Stiffness_Matrix_3D_8nodes(i_E, &
              Num_Gauss_P_FEM_3D,c_X_NODES,c_Y_NODES,c_Z_NODES, &
              c_D,kesi,yita,zeta,weight,localK)  
    ! Extract the non-zero elements of the stiffness matrix (without summing for now)
    do i_row = 1,24
        do i_col = 1,24
            c_row = local(i_row)
            c_col = local(i_col)
            ! The corresponding degree of freedom needs to be a functional degree of freedom.
            if(Flag_FreeDOF(c_row)==1  .and. Flag_FreeDOF(c_col)==1) then
                ! Row numbers corresponding to the degrees of freedom after removing constraints
                ! Location_FreeDOF is the stiffness matrix degree of freedom index after removing the constrained
                ! degrees of freedom corresponding to all degrees of freedom.
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

!*********************************************************************************************
!               Step 4:
! K merges duplicate data (summing them), and this function is very time-consuming. It can be
! optimized.
!*********************************************************************************************
K_Each_Count = Temp_K_Each_Count
allocate(K_Each(K_Each_Count,3))
K_Each(1:K_Each_Count,1:3) = Temp_K_Each(1:K_Each_Count,1:3)
deallocate(Temp_K_Each)
allocate(Combined_K_Each(K_Each_Count,3))
write ( * ,'(a)') '     Combining duplicate lines of K...' 
call Matrix_n_x_3_Combine_Same_i_j_Lines(K_Each_Count,K_Each(1:K_Each_Count,1:3),&
                                         Combined_K_Each(1:K_Each_Count,1:3), &
                                         num_rows_Combined_K)
                                       
!**************************
!           Step 5:
! Regenerate the merged K.
!**************************
deallocate(K_Each)
allocate(K_Each(num_rows_Combined_K,3))
K_Each(1:num_rows_Combined_K,1:3)  = Combined_K_Each(1:num_rows_Combined_K,1:3)
deallocate(Combined_K_Each)


!***************************************************
!           Step 6:
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
    print *,'    ERROR :: K_CSR_NNZ > K_CSR_NNZ_Max in Assemble_Stiffness_Matrix_SPARS_FEM_3D.f90!'
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
!           Step 7:
! Convert to CSR format.
!************************
write ( * ,'(a)') '     Converting K to CSR format...' 
call coocsr(num_FreeD,K_CSR_NNZ,K_vals(1:K_CSR_NNZ),K_rows(1:K_CSR_NNZ),K_cols(1:K_CSR_NNZ),&
                   K_CSR_aa(1:K_CSR_NNZ),K_CSR_ja(1:K_CSR_NNZ),K_CSR_ia(1:num_FreeD+1))                   
deallocate(K_rows,K_cols,K_vals)                   


RETURN
END SUBROUTINE Assemble_Stiffness_Matrix_SPARS_FEM_3D
