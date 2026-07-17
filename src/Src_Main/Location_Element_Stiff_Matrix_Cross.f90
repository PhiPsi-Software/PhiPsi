!-----------------------------------------------------------
! Brief: Compute the global DOF layout for a 2D element enriched at
!        a cross (intersection of two cracks) point.
!
! Parameters:
!   Input:  i_E, i_Cross      - element and cross indices
!           POS_c_Cross(Num_Node) - per-node DOF positions
!   Output: Location_ESM_C_Crack(80)   - FEM + cross DOF positions
!           num_Loc_ESM_C_Crack        - count
!           Location_ESM_C_Cr_NoFEM(60) - cross-only positions
!           num_Loc_ESM_C_Cr_NoFEM     - count
!
! Notes:   Non-Heaviside nodes are excluded from the cross DOF list.
!-----------------------------------------------------------

subroutine Location_Element_Stiff_Matrix_Cross(i_E,i_Cross, POS_c_Cross, Location_ESM_C_Crack, num_Loc_ESM_C_Crack, &
Location_ESM_C_Cr_NoFEM, num_Loc_ESM_C_Cr_NoFEM)
!     The position of the element stiffness matrix in the overall stiffness (corresponding to the intersection point of the cross, excluding the two cracks)
!     Nodes other than the Heaviside-enhanced nodes.

use Global_Float_Type
use Global_Model
use Global_Common
use Global_Cross

implicit none

integer,intent(in)::i_E,i_Cross,POS_c_Cross(Num_Node)
integer,intent(out)::Location_ESM_C_Crack(80)
integer,intent(out)::Location_ESM_C_Cr_NoFEM(60)
integer,intent(out)::num_Loc_ESM_C_Crack,num_Loc_ESM_C_Cr_NoFEM      
integer::Location_FEM(8)     
real(kind=FT) c_X_NODES(4),c_Y_NODES(4)
integer c_NN(4),i_N,cnt      
integer :: Enriched_Node_c_Cross(Num_Node)
integer :: num_crs
integer:: Location_XFEM(60)
integer :: tem

Enriched_Node_c_Cross = Enriched_Node_Type_Cross(:,i_Cross)

c_NN    = G_NN(:,i_E)
c_X_NODES = G_X_NODES(:,i_E)
c_Y_NODES = G_Y_NODES(:,i_E)
! Vector Initialization
Location_FEM(1:8)=0
Location_XFEM(1:60)=0
Location_ESM_C_Crack(1:80)=0
Location_ESM_C_Cr_NoFEM(1:60) = 0

num_Loc_ESM_C_Crack = 0
num_Loc_ESM_C_Cr_NoFEM = 0

if (i_Cross .eq. 1)then
    do i_N =1,4
        Location_FEM(2*i_N-1) = 2*c_NN(i_N) - 1
        Location_FEM(2*i_N)   = 2*c_NN(i_N)
    end do
end if

!If the element has no enriched nodes, then:
if(maxval(Enriched_Node_c_Cross(c_NN)).eq.0 .and. minval(Enriched_Node_c_Cross(c_NN)).eq.0) then
    if (i_Cross .eq. 1)then

        Location_ESM_C_Crack(1:8) = Location_FEM
        num_Loc_ESM_C_Crack =8
        num_Loc_ESM_C_Cr_NoFEM = 0
    end if
    !If the element has enriched nodes, then:      
else

    !Get the number of the cross enriched nodes of the element.
    num_crs = count(Enriched_Node_c_Cross(c_NN).eq.1)
    tem = 2*(num_crs)   
    cnt = 0

    do i_N = 1,4
        if (Enriched_Node_c_Cross(c_NN(i_N)) .eq. 1)then
            cnt = cnt + 1
            Location_XFEM(2*cnt-1) = 2*POS_c_Cross(c_NN(i_N)) - 1
            Location_XFEM(2*cnt  ) = 2*POS_c_Cross(c_NN(i_N))
        end if
    end do   
    !-----------------------
    if (i_Cross .eq. 1)then
        ! Includes FEM degrees of freedom
        Location_ESM_C_Crack(1:8)     = Location_FEM
        Location_ESM_C_Crack(9:8+tem) = Location_XFEM(1:tem)
        num_Loc_ESM_C_Crack           = 8+tem    
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
END subroutine Location_Element_Stiff_Matrix_Cross
