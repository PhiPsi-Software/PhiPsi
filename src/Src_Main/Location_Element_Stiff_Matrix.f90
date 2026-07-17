!-----------------------------------------------------------
! Brief: Compute the global-stiffness DOF layout (with and without
!        the FEM part) for a 2D XFEM element relative to a crack.
!
! Parameters:
!   Input:  i_E, i_C            - element and crack indices
!           POS_c_Crack(Num_Node) - per-node DOF positions for crack
!   Output: Location_ESM_C_Crack(80)   - global positions including FEM
!           num_Loc_ESM_C_Crack        - number of FEM+XFEM entries
!           Location_ESM_C_Cr_NoFEM(60) - XFEM-only positions
!           num_Loc_ESM_C_Cr_NoFEM     - XFEM-only count
!
! Notes:   Selects the number of tip branch functions from
!          Key_TipEnrich (1, 4, or 1) and skips FEM entries when
!          i_C is not the first crack.
!-----------------------------------------------------------

subroutine Location_Element_Stiff_Matrix(i_E,i_C,POS_c_Crack, Location_ESM_C_Crack, num_Loc_ESM_C_Crack, &
Location_ESM_C_Cr_NoFEM, num_Loc_ESM_C_Cr_NoFEM)
!     The position of the element stiffness matrix in the global stiffness matrix.
use Global_Float_Type
use Global_Crack
use Global_Model
use Global_Common

implicit none

integer,intent(in)::i_E,i_C,POS_c_Crack(Num_Node)
integer,intent(out)::Location_ESM_C_Crack(80)
integer,intent(out)::Location_ESM_C_Cr_NoFEM(60)
integer,intent(out)::num_Loc_ESM_C_Crack,num_Loc_ESM_C_Cr_NoFEM

integer::Location_FEM(8)     
real(kind=FT) c_X_NODES(4),c_Y_NODES(4)
integer c_NN(4),i_f,i_N      
integer :: Enriched_Node_c_Crack(Num_Node)
integer num_t,num_h,num_j,cnt
integer:: Location_XFEM(60)
integer :: tem

Enriched_Node_c_Crack = Enriched_Node_Type(:,i_C)

c_NN    = G_NN(:,i_E)
c_X_NODES = G_X_NODES(:,i_E)
c_Y_NODES = G_Y_NODES(:,i_E)
! Vector initialization
Location_FEM(1:8)=0
Location_XFEM(1:60)=0
Location_ESM_C_Crack(1:80)=0
Location_ESM_C_Cr_NoFEM(1:60) = 0

num_Loc_ESM_C_Crack = 0
num_Loc_ESM_C_Cr_NoFEM = 0

if (i_C .eq. 1)then
    do i_N =1,4
        Location_FEM(2*i_N-1) = 2*c_NN(i_N) - 1
        Location_FEM(2*i_N)   = 2*c_NN(i_N)
    end do
end if

!If the element has no enriched nodes, then:
if(maxval(Enriched_Node_c_Crack(c_NN)).eq.0 .and. minval(Enriched_Node_c_Crack(c_NN)).eq.0) then
    if (i_C .eq. 1)then

        Location_ESM_C_Crack(1:8) = Location_FEM
        num_Loc_ESM_C_Crack =8
        num_Loc_ESM_C_Cr_NoFEM = 0

    end if
    !If the element has enriched nodes, then:      
else

    !Get the number of the tip enriched nodes of the element.
    num_t = count(Enriched_Node_c_Crack(c_NN).eq.1)
    !Get the number of the Heaviside enriched nodes of the element.
    num_h = count(Enriched_Node_c_Crack(c_NN).eq.2)
    !Get the number of the tip junction nodes of the element.     
    num_j = count(Enriched_Node_c_Crack(c_NN).eq.3) + &
    count(Enriched_Node_c_Crack(c_NN).eq.6)

    if (Key_TipEnrich==0) then   
        tem = 2*(num_h*1 + num_j*1)   
    elseif (Key_TipEnrich==1) then   
        tem = 2*(num_h*1 + num_t*4 + num_j*1)   
    elseif((Key_TipEnrich==2).or.(Key_TipEnrich==3) .or. (Key_TipEnrich==4).or.(Key_TipEnrich==5)) then
        tem = 2*(num_h*1 + num_t + num_j*1)   
    end if
    cnt = 0

    do i_N = 1,4
        if (Enriched_Node_c_Crack(c_NN(i_N)) .eq. 2  )then
            cnt = cnt + 1
            Location_XFEM(2*cnt-1) = 2*POS_c_Crack(c_NN(i_N)) - 1
            Location_XFEM(2*cnt  ) = 2*POS_c_Crack(c_NN(i_N))
        elseif( Enriched_Node_c_Crack(c_NN(i_N))  .eq.  1)then
            !------------------------
            ! Tip Enhancement Plan 1
            !------------------------
            if (Key_TipEnrich==1) then   
                do i_f=1,4
                    cnt = cnt + 1
                    Location_XFEM(2*cnt-1) = 2*(POS_c_Crack(c_NN(i_N))+i_f-1) - 1
                    Location_XFEM(2*cnt  ) = 2*(POS_c_Crack(c_NN(i_N))+i_f-1)
                end do
                !-------------------------------
                ! Tip Enhancement Plans 2 and 3
                !-------------------------------
                ! 2: Keep only the first item F, 3: Heaviside smooth transition scheme (see my PhD thesis for details), (4) cohesive crack tip
            elseif((Key_TipEnrich==2).or.(Key_TipEnrich==3) .or. (Key_TipEnrich==4).or.(Key_TipEnrich==5) ) then
                cnt = cnt + 1
                Location_XFEM(2*cnt-1) = 2*(POS_c_Crack(c_NN(i_N))) - 1
                Location_XFEM(2*cnt  ) = 2*(POS_c_Crack(c_NN(i_N)))
            endif
        elseif (Enriched_Node_c_Crack(c_NN(i_N)) .eq.  3)then
            cnt = cnt + 1
            Location_XFEM(2*cnt-1) = 2*POS_c_Crack(c_NN(i_N)) - 1
            Location_XFEM(2*cnt )  = 2*POS_c_Crack(c_NN(i_N))
        elseif (Enriched_Node_c_Crack(c_NN(i_N)) .eq.  6)then
            cnt = cnt + 1
            Location_XFEM(2*cnt-1) = 2*POS_c_Crack(c_NN(i_N)) -1
            Location_XFEM(2*cnt )  = 2*POS_c_Crack(c_NN(i_N))
        elseif (Enriched_Node_c_Crack(c_NN(i_N)) .eq.  4)then
            cnt = cnt + 1
            Location_XFEM(2*cnt-1) = 2*POS_c_Crack(c_NN(i_N)) -1
            Location_XFEM(2*cnt )  = 2*POS_c_Crack(c_NN(i_N))
        end if
    end do   
    if (i_C .eq. 1)then
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
END subroutine Location_Element_Stiff_Matrix
