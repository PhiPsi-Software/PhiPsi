!-----------------------------------------------------------
! Brief: Patch the HF Jacobian with the -delta_t*H fluid block.
!
! Parameters:
!   Input:  Counter         - current iteration counter
!           Total_FD        - total displacement DOFs
!           num_free_CalP   - number of active fluid DOFs
!           freeDOF_HF      - global fluid-DOF numbering
!           Local_freeDOF_HF - local fluid-DOF indexing into H
!           delta_Time      - time increment
!           H               - fluid coupling matrix
!   In/Out: NR_Deri         - full Newton-Raphson Jacobian
!
! Notes:   Inserts -delta_Time*H at the fluid-DOF block of the
!          full Jacobian to enforce the time-discrete flow
!          coupling in subsequent NR iterations.
!-----------------------------------------------------------

subroutine Cal_HF_Modify_Jacobian(Counter,NR_Deri,Total_FD, num_free_CalP, freeDOF_HF, Local_freeDOF_HF,delta_Time,H)
! Hydraulic fracturing Newton-Raphson iteration algorithm Jacobian matrix (subsequent iteration
! steps)

!     ----------------------------
!     Read public variable module
!     ----------------------------
use Global_Float_Type
use Global_Crack
use Global_Crack_Common

implicit none
integer, intent(in)::Counter,Total_FD, num_free_CalP
real(kind=FT),intent(in)::H(num_Tol_CalP_Water, num_Tol_CalP_Water), delta_Time
integer, intent(in)::freeDOF_HF(num_Tol_CalP_Water), Local_freeDOF_HF(num_Tol_CalP_Water)
real(kind=FT),intent(inout):: NR_Deri(Total_FD+num_Tol_CalP_Water,Total_FD+num_Tol_CalP_Water)
real(kind=FT) k,max_K(3)

integer :: i

print *,'    Modifying the Jacobian matrix......'

NR_Deri(freeDOF_HF(1:num_free_CalP),freeDOF_HF(1:num_free_CalP))= -H(Local_freeDOF_HF(1:num_free_CalP), &
Local_freeDOF_HF(1:num_free_CalP))*delta_Time

return 
end subroutine Cal_HF_Modify_Jacobian      
