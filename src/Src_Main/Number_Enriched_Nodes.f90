!-----------------------------------------------------------
! Brief: Assign global DOF indices to the 2D XFEM enriched nodes
!        (Heaviside, tip, junction, cross, hole, inclusion) for the
!        current load step.
!
! Parameters:
!   Input:  isub - sub-step index (used for output labeling only)
!   Output: (none - fills c_POS, c_POS_Hl, c_POS_Cross,
!           c_POS_Incl, n_*_Node, Total_FD, Usual_Freedom,
!           Enrich_Freedom)
!
! Notes:   Branches on Key_TipEnrich to pick the number of tip
!          branch functions (0, 1, or 4).
!-----------------------------------------------------------

subroutine Number_Enriched_Nodes(isub)
!     Number the enriched nodes.
use Global_Float_Type
use Global_Common
use Global_Crack
use Global_Crack_Common
use Global_Model
use Global_Filename
use Global_Inclusion
use Global_Cross

implicit none

integer,intent(in)::isub
integer i_C,i_N,i_H,i_Incl
integer :: max_Cr_Node
character(5) temp
integer i_Cross,i_HC

print *,'    Numbering enriched nodes...'

n_h_Node =0
n_t_Node =0
n_j_Node =0
n_c_Node =0     
n_hl_Node=0
n_cross_Node=0
n_Incl_Node    = 0
Total_FD       = 0
Usual_Freedom  = 0
Enrich_Freedom = 0

Usual_Freedom  = 2*Num_Node

! Initialize c_POS
c_POS(1:Num_Node,1:num_Crack) = 0
! If there are holes, initialize c_POS_Hl
if(num_Hole /=0)then
    c_POS_Hl(1:Num_Node,1:num_Hole) = 0
endif
! If it contains Cross, initialize c_POS_Cross
if(num_Cross /=0)then
    c_POS_Cross(1:Num_Node,1:num_Cross) = 0
endif
! If there are inclusions, initialize c_POS_Incl
if(num_Inclusion /=0)then
    c_POS_Incl(1:Num_Node,1:num_Inclusion) = 0
endif

!Loop through each crack.
do i_C = 1,num_Crack
    !Loop through each node.
    do i_N = 1,Num_Node
        !----------------------------
        ! No crack tip reinforcement
        !----------------------------
        if (Key_TipEnrich==0) then
            if (Enriched_Node_Type(i_N,i_C).eq.2) then
                c_POS(i_N,i_C) = (Num_Node+n_h_Node+n_j_Node+n_c_Node) + 1
                n_h_Node = n_h_Node + 1
            elseif (Enriched_Node_Type(i_N,i_C) .eq.3)then
                c_POS(i_N,i_C) = (Num_Node+n_h_Node+n_j_Node+n_c_Node) + 1
                n_j_Node = n_j_Node + 1
            elseif (Enriched_Node_Type(i_N,i_C) .eq.6)then
                c_POS(i_N,i_C) = (Num_Node+n_h_Node+n_j_Node+n_c_Node) + 1
                n_j_Node = n_j_Node + 1
            elseif (Enriched_Node_Type(i_N,i_C) .eq.4)then
                c_POS(i_N,i_C) = (Num_Node+n_h_Node+n_j_Node+n_c_Node) + 1
                n_c_Node = n_c_Node + 1
            end if
            !------------------------------------
            ! Standard Tip Enhancement (4 items)
            !------------------------------------
        elseif (Key_TipEnrich==1) then
            if (Enriched_Node_Type(i_N,i_C).eq.2) then
                c_POS(i_N,i_C) = (Num_Node + n_h_Node + n_t_Node*4 + n_j_Node+ n_c_Node) + 1
                n_h_Node = n_h_Node + 1
            elseif (Enriched_Node_Type(i_N,i_C) .eq.1)then
                c_POS(i_N,i_C) = (Num_Node + n_h_Node + n_t_Node*4 + n_j_Node+ n_c_Node) + 1
                n_t_Node = n_t_Node + 1
            elseif (Enriched_Node_Type(i_N,i_C) .eq.3)then
                c_POS(i_N,i_C) = (Num_Node + n_h_Node + n_t_Node*4 + n_j_Node+ n_c_Node) + 1
                n_j_Node = n_j_Node + 1
            elseif (Enriched_Node_Type(i_N,i_C) .eq.6)then
                c_POS(i_N,i_C) = (Num_Node + n_h_Node + n_t_Node*4 + n_j_Node+ n_c_Node) + 1
                n_j_Node = n_j_Node + 1
            elseif (Enriched_Node_Type(i_N,i_C) .eq.4)then
                c_POS(i_N,i_C) = (Num_Node + n_h_Node + n_t_Node*4 + n_j_Node+ n_c_Node) + 1
                n_c_Node = n_c_Node + 1
            end if
            !-------------------------------
            ! Tip Enhancement Plans 2 and 3
            !-------------------------------
        elseif(Key_TipEnrich==2 .or. Key_TipEnrich==3 .or. Key_TipEnrich==4 .or. Key_TipEnrich==5) then
            ! 2: Keep only the first item F,
            ! 3: Heaviside smooth transition scheme (see my PhD thesis for details), 
            ! 4: Viscous crack tip
            ! 5: Asymptotic H function, 2026-05-26.
            if (Enriched_Node_Type(i_N,i_C).eq.2) then
                c_POS(i_N,i_C) = (Num_Node + n_h_Node + n_t_Node + n_j_Node+ n_c_Node) + 1
                n_h_Node = n_h_Node + 1
            elseif (Enriched_Node_Type(i_N,i_C) .eq.1)then
                c_POS(i_N,i_C) = (Num_Node + n_h_Node + n_t_Node + n_j_Node+ n_c_Node) + 1
                n_t_Node = n_t_Node + 1
            elseif (Enriched_Node_Type(i_N,i_C) .eq.3)then
                c_POS(i_N,i_C) = (Num_Node + n_h_Node + n_t_Node + n_j_Node+ n_c_Node) + 1
                n_j_Node = n_j_Node + 1
            elseif (Enriched_Node_Type(i_N,i_C) .eq.6)then
                c_POS(i_N,i_C) = (Num_Node + n_h_Node + n_t_Node + n_j_Node+ n_c_Node) + 1
                n_j_Node = n_j_Node + 1
            elseif (Enriched_Node_Type(i_N,i_C) .eq.4)then
                c_POS(i_N,i_C) = (Num_Node + n_h_Node + n_t_Node + n_j_Node+ n_c_Node) + 1
                n_c_Node = n_c_Node + 1
            end if
        end if
    end do
end do     

! Maximum enhancement node number (of the crack)
if(num_Crack >=1)then
    max_Cr_Node = maxval(c_POS(1:Num_Node,1:num_Crack))
else
    max_Cr_Node = Num_Node
endif

!Loop through each hole.
if(num_Hole /=0)then
    do i_H = 1,num_Hole
        !Loop through each node.
        do i_N = 1,Num_Node
            if (Enriched_Node_Type_Hl(i_N,i_H) /= 0) then
                c_POS_Hl(i_N,i_H) = max_Cr_Node+n_hl_Node +1
                n_hl_Node = n_hl_Node + 1
            end if
        end do
    end do  
endif

!Loop through each cross.
if(num_Cross /=0)then
    do i_Cross = 1,num_Cross
        !Loop through each node.
        do i_N = 1,Num_Node
            if(Enriched_Node_Type_Cross(i_N,i_Cross)/= 0)then
                c_POS_Cross(i_N,i_Cross) = max_Cr_Node+n_cross_Node +1
                n_cross_Node = n_cross_Node + 1
            end if
        end do
    end do  
endif


! Various mixed loops.
if(num_Inclusion /=0)then
    do i_Incl = 1,num_Inclusion
        !Loop through each node.
        do i_N = 1,Num_Node
            if (Enriched_Node_Type_Incl(i_N,i_Incl) /= 0)then
                c_POS_Incl(i_N,i_Incl) = max_Cr_Node+n_Incl_Node +1
                n_Incl_Node = n_Incl_Node + 1
            end if
        end do
    end do  
endif

!Total degrees of freedom. 
if (Key_TipEnrich==0) then    
    Total_FD=2*(Num_Node+n_h_Node+n_j_Node+ n_c_Node)
elseif (Key_TipEnrich==1) then
    Total_FD=2*(Num_Node+n_h_Node+n_t_Node*4+n_j_Node+ n_c_Node)
elseif (Key_TipEnrich==2 .or. Key_TipEnrich==3 .or. Key_TipEnrich==4 .or. Key_TipEnrich==5) then
    ! 2: Keep only the first item F.
    ! 3: Heaviside smooth transition scheme (see my PhD thesis for details).
    ! 4: Viscous crack tip.
    ! 5: Asymptotic H function, 2026-05-26.
    Total_FD=2*(Num_Node+n_h_Node+n_t_Node+n_j_Node+ n_c_Node)
endif
if(num_Hole /=0)then
    Total_FD =  Total_FD + 2*n_hl_Node
endif
if(num_Cross /=0)then
    Total_FD =  Total_FD + 2*n_cross_Node
endif 
if(num_Inclusion /=0)then
    Total_FD =  Total_FD + 2*n_Incl_Node
endif

Enrich_Freedom = Total_FD - Usual_Freedom

return
END subroutine Number_Enriched_Nodes
