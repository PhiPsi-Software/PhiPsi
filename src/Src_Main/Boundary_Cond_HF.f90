!-----------------------------------------------------------
! Brief: Apply hydraulic-fracturing pressure boundary
!        conditions, fixing the fluid pressure at crack tips.
!
! Parameters:
!   Input:  iter - current HF iteration index
!   Output: freeDOF_HF        - free pressure DOFs
!           num_free_CalP     - number of free pressure pts
!           Local_freeDOF_HF  - local-to-global index map
!
! Notes:   Crack-tip pressure is set to zero when the
!          toughness-storage parameter is small and the
!          Key_Tip_Pres_Zero keyword is enabled.
!-----------------------------------------------------------

subroutine Boundary_Cond_HF(iter, freeDOF_HF,num_free_CalP, Local_freeDOF_HF)
!     Boundary condition (the water pressure at the crack tip is 0, so the water pressure degree of freedom at the crack tip is constrained)

! num_Tol_CalP: Total number of calculation points
! num_free_CalP: Number of calculation points after removing constraints

!     ----------------------------
!     Read public variable module
!     ----------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack
use Global_Crack_Common
use Global_Material

!     --------------------------
!     Variable Type Declaration
!     --------------------------
implicit none

integer,intent(in)::iter
integer,intent(out)::freeDOF_HF(num_Tol_CalP_Water),num_free_CalP
integer,intent(out)::Local_freeDOF_HF(num_Tol_CalP_Water)

integer i_CalP,i_C,i_fd
integer fixedDOF(num_Tol_CalP_Water),num_Count

print *,'    Dealing with boundary condition of HF...'   

!Initialize vector of fixed DOFs.
freeDOF_HF(1:num_Tol_CalP_Water) = 0
fixedDOF(1:num_Tol_CalP_Water)   = 0

num_Count = 0
do i_C = 1,num_Crack
    ! If the current crack is driven by fracturing fluid, then:
    if (Cracks_HF_State(i_C) == 1) then  
        ! Current Crack Calculation Point Cycle
        do i_CalP=1,Cracks_CalP_Num(i_C)
            num_Count = num_Count + 1
            if (Cracks_CalP_Type(i_C,i_CalP,1) == 5) then
                ! See Cal_HF_Crack_Points_Info_Linear(isub) for details
                ! If viscosity dominates and the keyword for setting crack tip water pressure to 0 is enabled, then
                ! set the crack tip water pressure to 0.
                if(Global_K_m < ONE .and. Key_Tip_Pres_Zero==1) then
                    fixedDOF(num_Count) =  num_Count  
                endif
            end if

        end do
    end if
end do

num_free_CalP = 0              
do i_fd =1,num_Tol_CalP_Water
    if ( .not.(ANY( fixedDOF(1:num_Count) .eq. i_fd)) ) then
        num_free_CalP = num_free_CalP +1
        freeDOF_HF(num_free_CalP) = Total_FD+i_fd
        Local_freeDOF_HF(num_free_CalP) = i_fd
    end if
end do

return
END subroutine Boundary_Cond_HF
