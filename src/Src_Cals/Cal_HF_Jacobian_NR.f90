!-----------------------------------------------------------
! Brief: Assemble the dense Newton-Raphson Jacobian for the
!        coupled HF mechanical-fluid system.
!
! Parameters:
!   Input:  Counter, Total_FD, num_FreeD, num_free_CalP - sizes
!           globalK           - mechanical stiffness matrix
!           Coupled_Q         - coupling matrix
!           freeDOF/HF        - DOF index maps
!           Local_freeDOF_HF  - local-to-global HF map
!           delta_Time, H     - time step and H matrix
!   Output: NR_Deri           - assembled Jacobian
!
! Notes:   Layout: top-left K, off-diagonals +/-Q, bottom-right
!          -delta_Time*H.
!-----------------------------------------------------------

subroutine Cal_HF_Jacobian_NR(Counter,NR_Deri,Total_FD, num_FreeD, num_free_CalP, globalK,Coupled_Q, freeDOF,freeDOF_HF, &
Local_freeDOF_HF,delta_Time,H)
! Hydraulic Fracturing Newton-Raphson Iterative Algorithm Assembled Derivative Matrix (Jacobian
! Matrix)

!     ----------------------------
!     Read public variable module
!     ----------------------------
use Global_Float_Type
use Global_Crack
use Global_Crack_Common

implicit none
integer, intent(in)::Counter,Total_FD,num_FreeD, num_free_CalP
real(kind=FT),intent(in):: Coupled_Q(Total_FD,num_Tol_CalP_Water), globalK(Total_FD,Total_FD),delta_Time, &
H(num_Tol_CalP_Water,num_Tol_CalP_Water)
integer, intent(in)::freeDOF(Total_FD), freeDOF_HF(num_Tol_CalP_Water), Local_freeDOF_HF(num_Tol_CalP_Water)
real(kind=FT),intent(out):: NR_Deri(Total_FD+num_Tol_CalP_Water,Total_FD+num_Tol_CalP_Water)
print *,'    Calculating the Jacobian matrix of ' //'NR interation......'

NR_Deri(1:Total_FD+num_Tol_CalP_Water, 1:Total_FD+num_Tol_CalP_Water)=ZR

NR_Deri(freeDOF(1:num_FreeD),freeDOF(1:num_FreeD)) = globalK(freeDOF(1:num_FreeD),freeDOF(1:num_FreeD))
NR_Deri(freeDOF(1:num_FreeD),freeDOF_HF(1:num_free_CalP)) = -Coupled_Q(freeDOF(1:num_FreeD), &
Local_freeDOF_HF(1:num_free_CalP))
NR_Deri(freeDOF_HF(1:num_free_CalP),freeDOF(1:num_FreeD)) = -transpose(Coupled_Q(freeDOF(1:num_FreeD), &
Local_freeDOF_HF(1:num_free_CalP)))
NR_Deri(freeDOF_HF(1:num_free_CalP),freeDOF_HF(1:num_free_CalP))= -H(Local_freeDOF_HF(1:num_free_CalP), &
Local_freeDOF_HF(1:num_free_CalP))*delta_Time

return 
end subroutine Cal_HF_Jacobian_NR          
