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
 
SUBROUTINE Assemble_Stiffness_Matrix_SPARS_XFEM_Get_MaxNNZ(isub,&
                          freeDOF,num_FreeD,T_Freedom,K_CSR_NNZ_Max)
! Obtain the maximum number of non-zero elements K_CSR_NNZ_Max before assembling the stiffness
! matrix.
! 2024-11-16.
!

!**********************
! Read public variable
!**********************
use Global_Float_Type
use Global_Crack
use Global_Crack_Common
use Global_Model
use Global_Filename
use Global_Common
use Global_Model,only:flag_freedof,G_NN,Num_Elem
use Global_Inclusion,only:num_Inclusion,c_POS_Incl
use Global_Cross,only:num_cross,c_POS_Cross

!**********************
! Variable Declaration
!**********************
implicit none
integer,intent(in)::isub,num_FreeD,T_Freedom
integer,intent(out)::K_CSR_NNZ_Max
integer,intent(in)::freeDOF(1:num_FreeD)
!---------------------
integer i_C,i_E
integer c_NN(4)
integer:: Location_ESM(MDOF_2D)
integer::Location_ESM_C_Crack(80),Location_ESM_C_Cr_NoFEM(60)
integer num_Loc_ESM_C_Cr_NoFEM
integer num_Loc_ESM
integer num_Loc_ESM_C_Crack
! The following variables are related to the set of sparse matrices
integer i_row,i_col
integer c_row,c_col
integer Temp_K_Each_Count
integer i_H,i_Cross,i_Incl

EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gauss_Points)= .False.
    
Temp_K_Each_Count = 0
    
do i_E = 1,Num_Elem
    c_NN    = G_NN(:,i_E)

    ! ------------------------------
    ! Point Scheme 1: Triangulation
    ! ------------------------------
    if(Key_Integral_Sol.eq.1)then
        !call Cal_Gauss_Points_Subtriangle(kesi,yita,weight)
        !c_Num_Gauss_Point  = size(kesi,2);
    ! ----------------------------------------
    ! Point scheme 2 or 3: evenly distributed
    ! ----------------------------------------
    elseif(Key_Integral_Sol.eq.2 .or. Key_Integral_Sol.eq.3)then
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !Decide the location of each element stiffness
        !matrix in the global stiffness matrix.
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Location_ESM(1:MDOF_2D)       = 0
        !Location_ESM_NoFEM(1:160) = 0
        num_Loc_ESM       = 0
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
    
K_CSR_NNZ_Max = Temp_K_Each_Count

RETURN
END SUBROUTINE Assemble_Stiffness_Matrix_SPARS_XFEM_Get_MaxNNZ
