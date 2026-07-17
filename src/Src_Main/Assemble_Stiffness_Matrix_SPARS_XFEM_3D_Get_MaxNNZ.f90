!-----------------------------------------------------------
! Brief: Compute the maximum nonzero count for the 3D
!        XFEM stiffness matrix including crack DOFs.
!
! Parameters:
!   Input:  isub      - current load step index
!           freeDOF   - mapping of free DOFs
!           num_FreeD - number of free DOFs
!           T_Freedom - total DOF count
!   Output: K_CSR_NNZ_Max - maximum nonzero count
!
! Notes:   Pre-pass that walks each element, locates the
!          stiffness entries (FEM + enriched), and counts
!          only those belonging to free DOFs.
!-----------------------------------------------------------

SUBROUTINE Assemble_Stiffness_Matrix_SPARS_XFEM_3D_Get_MaxNNZ(isub,freeDOF,num_FreeD,T_Freedom,K_CSR_NNZ_Max)
! Obtain the maximum number of non-zero elements K_CSR_NNZ_Max before assembling the global
! stiffness matrix.
! 2024-11-15.
!

!***************************
! Global variable reference
!***************************
use Global_Float_Type
use Global_Model,only:G_NN,Num_Elem,Flag_FreeDOF
use Global_Common,only:MDOF_3D
use Global_Crack_Common,only:num_Crack
use Global_Crack_3D,only:c_POS_3D

!**********************
! Variable Declaration
!**********************
implicit none
!include 'omp_lib.h'
integer,intent(in)::isub,num_FreeD,T_Freedom
integer(kind=LIT),intent(out)::K_CSR_NNZ_Max
integer,intent(in)::freeDOF(1:num_FreeD)
!-----------------
integer i_E   
integer c_NN(8)  
integer i_row,i_col
integer c_row,c_col
integer(kind=LIT) Temp_K_Each_Count
integer:: Location_ESM(MDOF_3D)
integer::Location_ESM_C_Crack(MDOF_3D),Location_ESM_C_Cr_NoFEM(MDOF_3D)
integer num_Loc_ESM_C_Cr_NoFEM
integer num_Loc_ESM,num_Loc_ESM_C_Crack
integer i_C
integer c_POS_3D_c_Ele(8)

Temp_K_Each_Count = 0
do i_E = 1,Num_Elem
    c_NN    = G_NN(:,i_E)
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !Decide the location of each element stiffness 
    !matrix in the global stiffness matrix.   
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Location_ESM(1:MDOF_3D)   = 0
    num_Loc_ESM               = 0
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

K_CSR_NNZ_Max = Temp_K_Each_Count
RETURN
END SUBROUTINE Assemble_Stiffness_Matrix_SPARS_XFEM_3D_Get_MaxNNZ
