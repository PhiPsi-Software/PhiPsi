!-----------------------------------------------------------
! Brief: Compute the maximum number of nonzero entries
!        needed to store the 3D FEM stiffness matrix.
!
! Parameters:
!   Input:  isub          - current load step index
!           freeDOF       - mapping of free DOFs
!           num_FreeD     - number of free DOFs
!           T_Freedom     - total DOF count
!   Output: K_CSR_NNZ_Max - maximum nonzero count
!
! Notes:   Used as a pre-pass to allocate CSR arrays
!          before the actual stiffness assembly.
!-----------------------------------------------------------

SUBROUTINE Assemble_Stiffness_Matrix_SPARS_FEM_3D_Get_MaxNNZ(isub,freeDOF,num_FreeD,T_Freedom,K_CSR_NNZ_Max)
! Obtain the maximum number of non-zero elements K_CSR_NNZ_Max before assembling the stiffness
! matrix.
! 2024-11-15.
!

!***************************
! Global variable reference
!***************************
use Global_Float_Type
use Global_Model,only:G_NN,Num_Elem,Location_FreeDOF,Flag_FreeDOF
use Global_Filename
use Global_Common

!**********************
! Variable Declaration
!**********************
implicit none
integer,intent(in)::isub,num_FreeD,T_Freedom
integer(kind=LIT),intent(out)::K_CSR_NNZ_Max
integer,intent(in)::freeDOF(1:num_FreeD)
integer i_E   
integer c_NN(8)  
integer local(24),i_row,i_col
integer c_row,c_col
integer real_row,real_col
integer(kind=LIT) Temp_K_Each_Count

Temp_K_Each_Count = 0
do i_E = 1,Num_Elem
    c_NN    = G_NN(:,i_E)
    !Traditional index locations
local=[c_NN(1)*3-2,c_NN(1)*3-1,c_NN(1)*3, c_NN(2)*3-2,c_NN(2)*3-1,c_NN(2)*3, c_NN(3)*3-2,c_NN(3)*3-1,c_NN(3)*3, &
c_NN(4)*3-2,c_NN(4)*3-1,c_NN(4)*3, c_NN(5)*3-2,c_NN(5)*3-1,c_NN(5)*3, c_NN(6)*3-2,c_NN(6)*3-1,c_NN(6)*3, &
c_NN(7)*3-2,c_NN(7)*3-1,c_NN(7)*3, c_NN(8)*3-2,c_NN(8)*3-1,c_NN(8)*3]
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

K_CSR_NNZ_Max = Temp_K_Each_Count
RETURN
END SUBROUTINE Assemble_Stiffness_Matrix_SPARS_FEM_3D_Get_MaxNNZ
