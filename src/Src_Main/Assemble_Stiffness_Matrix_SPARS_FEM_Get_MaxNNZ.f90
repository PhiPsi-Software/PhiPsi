!-----------------------------------------------------------
! Brief: Compute the maximum nonzero count for the 2D
!        FEM stiffness matrix in CSR sparse storage.
!
! Parameters:
!   Input:  isub          - current load step index
!           freeDOF       - mapping of free DOFs
!           num_FreeD     - number of free DOFs
!           T_Freedom     - total DOF count
!   Output: K_CSR_NNZ_Max - maximum nonzero count
!
! Notes:   Pre-pass to size the sparse storage buffers
!          before assembly of the global stiffness.
!-----------------------------------------------------------

SUBROUTINE Assemble_Stiffness_Matrix_SPARS_FEM_Get_MaxNNZ(isub, freeDOF,num_FreeD,T_Freedom,K_CSR_NNZ_Max)
! Obtain the maximum number of non-zero elements K_CSR_NNZ_Max before assembling the global
! stiffness matrix.
! 2024-11-16.
!

!***************************
! Global variable reference
!***************************
use Global_Model,only:flag_freedof,G_NN,Num_Elem

!**********************
! Variable Declaration
!**********************
implicit none
integer,intent(in)::isub,num_FreeD,T_Freedom
integer,intent(out)::K_CSR_NNZ_Max
integer,intent(in)::freeDOF(1:num_FreeD)
!-------
integer i_E
integer c_NN(4) 
integer local(8),i_row,i_col
integer c_row,c_col
integer Temp_K_Each_Count

Temp_K_Each_Count = 0
do i_E = 1,Num_Elem
    c_NN    = G_NN(:,i_E)
    !Traditional index locations
local=[c_NN(1)*2-1,c_NN(1)*2,c_NN(2)*2-1,c_NN(2)*2, c_NN(3)*2-1,c_NN(3)*2,c_NN(4)*2-1,c_NN(4)*2]

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

K_CSR_NNZ_Max  =  Temp_K_Each_Count               

RETURN
END SUBROUTINE Assemble_Stiffness_Matrix_SPARS_FEM_Get_MaxNNZ
