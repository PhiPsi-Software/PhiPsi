!-----------------------------------------------------------
! Brief: Compute the global DOF layout for a 3D XFEM 8-node
!        hexahedral element, with and without the standard FEM
!        contribution.
!
! Parameters:
!   Input:  i_E, i_C              - element and crack indices
!           c_POS_3D_c_Ele(8)     - per-node DOF positions
!   Output: Location_ESM_C_Crack(MDOF_3D)   - global positions incl. FEM
!           num_Loc_ESM_C_Crack              - count of FEM+XFEM entries
!           Location_ESM_C_Cr_NoFEM(MDOF_3D) - XFEM-only positions
!           num_Loc_ESM_C_Cr_NoFEM           - XFEM-only count
!
! Notes:   Counts tip / Heaviside / junction enrichments and uses
!          Num_F_Functions for tip DOFs per node.
!-----------------------------------------------------------

subroutine Location_Element_Stiff_Matrix_3D(i_E,i_C, c_POS_3D_c_Ele, Location_ESM_C_Crack, num_Loc_ESM_C_Crack, &
Location_ESM_C_Cr_NoFEM, num_Loc_ESM_C_Cr_NoFEM)
!     Position of the element stiffness matrix in the global stiffness (3D problem).
use Global_Float_Type
use Global_Crack_3D
use Global_Model
use Global_Common

implicit none

integer,intent(in)::i_E,i_C
integer,intent(in)::c_POS_3D_c_Ele(8)
integer,intent(out)::Location_ESM_C_Crack(MDOF_3D)
integer,intent(out)::Location_ESM_C_Cr_NoFEM(MDOF_3D)
integer,intent(out)::num_Loc_ESM_C_Crack,num_Loc_ESM_C_Cr_NoFEM

integer::Location_FEM(24)     
integer c_NN(8),i_f,i_N      
integer num_t,num_h,num_j,cnt
integer:: Location_XFEM(MDOF_3D)  
integer :: tem

!Enriched_Node_c_Crack = Enriched_Node_Type_3D(1:Num_Node,i_C)
c_NN    = G_NN(1:8,i_E)

! Vector Initialization
Location_FEM(1:24)                 = 0
Location_XFEM(1:MDOF_3D)           = 0
Location_ESM_C_Crack(1:MDOF_3D)    = 0
Location_ESM_C_Cr_NoFEM(1:MDOF_3D) = 0
num_Loc_ESM_C_Crack                = 0
num_Loc_ESM_C_Cr_NoFEM             = 0

if (i_C .eq. 1)then
    do i_N =1,8
        Location_FEM(3*i_N-2) = 3*c_NN(i_N) - 2
        Location_FEM(3*i_N-1) = 3*c_NN(i_N) - 1
        Location_FEM(3*i_N)   = 3*c_NN(i_N)
    end do
end if

!////////////////////////////////////////////
!If the element has no enriched nodes, then:
!////////////////////////////////////////////
if(sum(Enriched_Node_Type_3D(c_NN,i_C)) == 0) then
    if (i_C .eq. 1)then
        do i_N =1,8
            Location_FEM(3*i_N-2) = 3*c_NN(i_N) - 2
            Location_FEM(3*i_N-1) = 3*c_NN(i_N) - 1
            Location_FEM(3*i_N)   = 3*c_NN(i_N)
        end do

        Location_ESM_C_Crack(1:24) = Location_FEM
        num_Loc_ESM_C_Crack =24
        num_Loc_ESM_C_Cr_NoFEM = 0
    end if

    return

endif

!/////////////////////////////////////////
!If the element has enriched nodes, then:      
!/////////////////////////////////////////
if(sum(Enriched_Node_Type_3D(c_NN,i_C)) > 0) then
    !Get the number of the tip enriched nodes of the element.
    num_t = count(Enriched_Node_Type_3D(c_NN,i_C).eq.1)
    !Get the number of the Heaviside enriched nodes of the element.
    num_h = count(Enriched_Node_Type_3D(c_NN,i_C).eq.2)
    !Get the number of the tip junction nodes of the element.     
    num_j = count(Enriched_Node_Type_3D(c_NN,i_C).eq.3)

    tem = 3*(num_h*1 + num_t*Num_F_Functions + num_j*1) 
    cnt = 0

    do i_N = 1,8
        if (Enriched_Node_Type_3D(c_NN(i_N),i_C) .eq. 2  )then
            cnt = cnt + 1
            Location_XFEM(3*cnt-2) = 3*c_POS_3D_c_Ele(i_N) - 2
            Location_XFEM(3*cnt-1) = 3*c_POS_3D_c_Ele(i_N) - 1
            Location_XFEM(3*cnt  ) = 3*c_POS_3D_c_Ele(i_N)
        elseif(Enriched_Node_Type_3D(c_NN(i_N),i_C)  .eq.  1)then
            do i_f=1,Num_F_Functions
                cnt = cnt + 1
                Location_XFEM(3*cnt-2) = 3*(c_POS_3D_c_Ele(i_N)+i_f-1) - 2
                Location_XFEM(3*cnt-1) = 3*(c_POS_3D_c_Ele(i_N)+i_f-1) - 1
                Location_XFEM(3*cnt  ) = 3*(c_POS_3D_c_Ele(i_N)+i_f-1)
            end do               
        elseif (Enriched_Node_Type_3D(c_NN(i_N),i_C) .eq.  3)then
            cnt = cnt + 1
            Location_XFEM(3*cnt-2) = 3*c_POS_3D_c_Ele(i_N) - 2
            Location_XFEM(3*cnt-1) = 3*c_POS_3D_c_Ele(i_N) - 1
            Location_XFEM(3*cnt )  = 3*c_POS_3D_c_Ele(i_N)
        end if
    end do   

    if (i_C .eq. 1)then
        ! Includes FEM degrees of freedom
        Location_ESM_C_Crack(1:24)      = Location_FEM
        Location_ESM_C_Crack(25:24+tem) = Location_XFEM(1:tem)
        num_Loc_ESM_C_Crack             = 24+tem    
        ! Does not include FEM degrees of freedom
        Location_ESM_C_Cr_NoFEM(1:tem)= Location_XFEM(1:tem)
        num_Loc_ESM_C_Cr_NoFEM        = tem
    else
        ! Includes FEM degrees of freedom
        Location_ESM_C_Crack(1:tem)   = Location_XFEM(1:tem)
        num_Loc_ESM_C_Crack           = tem 
        ! Does not include FEM degrees of freedom
        Location_ESM_C_Cr_NoFEM(1:tem)= Location_XFEM(1:tem)
        num_Loc_ESM_C_Cr_NoFEM        = tem              
    end if
end if

return
END subroutine Location_Element_Stiff_Matrix_3D
