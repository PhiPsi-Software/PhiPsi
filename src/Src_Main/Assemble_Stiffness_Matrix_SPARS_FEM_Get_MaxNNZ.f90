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
 
SUBROUTINE Assemble_Stiffness_Matrix_SPARS_FEM_Get_MaxNNZ(isub,&
                          freeDOF,num_FreeD,T_Freedom,K_CSR_NNZ_Max)
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

K_CSR_NNZ_Max  =  Temp_K_Each_Count               

RETURN
END SUBROUTINE Assemble_Stiffness_Matrix_SPARS_FEM_Get_MaxNNZ
