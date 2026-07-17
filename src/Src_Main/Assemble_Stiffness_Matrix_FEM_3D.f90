!-----------------------------------------------------------
! Brief: Assemble the full 3D FEM global stiffness matrix
!        for 8-node hexahedral elements.
!
! Parameters:
!   Input:  isub            - current load step index
!           c_Total_Freedom - total DOF count
!   Output: c_globalK       - assembled stiffness matrix
!           Total_Num_G_P   - total number of Gauss points
!
! Notes:   Supports composite material rotation, element
!          damage softening, and OpenMP parallelization.
!-----------------------------------------------------------

subroutine Assemble_Stiffness_Matrix_FEM_3D(isub,c_globalK, c_Total_Freedom,Total_Num_G_P)
!     Assemble the stiffness matrix (3D problem).
use Global_Float_Type
use Global_Model
use Global_Filename
use Global_Common
use Global_Material
use omp_lib
! Program Interface Testing. 2022-09-06.
use Global_Cal_Ele_Stiffness_Matrix_3D_8nodes

implicit none
!include 'omp_lib.h'
integer,intent(in)::isub,c_Total_Freedom
integer,intent(out)::Total_Num_G_P
real(kind=FT) ,intent(out)::c_globalK(c_Total_Freedom, c_Total_Freedom)
integer :: i_E
integer R_num,C_num
real(kind=FT) c_D(6,6)
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
integer :: c_NN(8)
real(kind=FT) kesi(Num_Gauss_P_FEM_3D),yita(Num_Gauss_P_FEM_3D), zeta(Num_Gauss_P_FEM_3D),weight(Num_Gauss_P_FEM_3D)
real(kind=FT) localK(24,24)
integer local(24),i_row,i_col,nIndex 
integer :: mat_num
real(kind=FT) Rot_c_D_Comp(6,6),c_D_Comp(6,6),Volume_Ratio
real(kind=FT) T_Matrix(6,6),TT_Matrix(6,6)

c_globalK(1:Total_FD,1:Total_FD) = ZR


call Cal_Gauss_Points_3D_8nodes(Num_Gauss_P_FEM_3D, kesi,yita,zeta,weight)
nIndex = 0
Total_Num_G_P = 0

!.............................
! OpenMP multi-core computing
!.............................
!$OMP PARALLEL do DEFAULT(SHARED) private(i_E,i_row,i_col,localK, &
!$OMP& c_D,c_NN,c_X_NODES,c_Y_NODES,c_Z_NODES,local, &
!$OMP& Volume_Ratio,c_D_Comp,T_Matrix,TT_Matrix, &
!$OMP& Rot_c_D_Comp,mat_num) SCHEDULE(static)

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

    ! Damage element processing: weakening D matrix (elastic modulus). 2026-02-02. NEWFTU-2026020201.
    if(Key_Element_Break==1 .and. Elem_Break(i_E))then
        !c_D= c_D*Tol_3
        c_D= c_D*Tol_6
    endif  

    c_NN    = G_NN(1:8,i_E)
    c_X_NODES = G_X_NODES(1:8,i_E)
    c_Y_NODES = G_Y_NODES(1:8,i_E)    
    c_Z_NODES = G_Z_NODES(1:8,i_E) 
    !$omp critical      
    Total_Num_G_P = Total_Num_G_P + Num_Gauss_P_FEM_3D
    !$omp end critical  
    !Traditional index locations
    local=[c_NN(1)*3-2,c_NN(1)*3-1,c_NN(1)*3, c_NN(2)*3-2,c_NN(2)*3-1,c_NN(2)*3, c_NN(3)*3-2,c_NN(3)*3-1,c_NN(3)*3, &
    c_NN(4)*3-2,c_NN(4)*3-1,c_NN(4)*3, c_NN(5)*3-2,c_NN(5)*3-1,c_NN(5)*3, c_NN(6)*3-2,c_NN(6)*3-1,c_NN(6)*3, &
    c_NN(7)*3-2,c_NN(7)*3-1,c_NN(7)*3, c_NN(8)*3-2,c_NN(8)*3-1,c_NN(8)*3]
    ! Calculate the stiffness matrix of the current element
    call Cal_Ele_Stiffness_Matrix_3D_8nodes(i_E, Num_Gauss_P_FEM_3D,c_X_NODES,c_Y_NODES,c_Z_NODES, &
    c_D,kesi,yita,zeta,weight, localK)
    !$omp critical          
    ! Global stiffness matrix of the assembly
    do i_row = 1,24
        do i_col = 1,24
            c_globalK(local(i_row),local(i_col)) = c_globalK(local(i_row),local(i_col)) + localK(i_row,i_col)
        end do
    end do  
    !$omp end critical          
end do
!$omp end parallel do 



return
END subroutine Assemble_Stiffness_Matrix_FEM_3D
