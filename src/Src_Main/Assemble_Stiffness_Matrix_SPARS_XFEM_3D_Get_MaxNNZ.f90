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
