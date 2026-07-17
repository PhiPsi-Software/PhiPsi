!-----------------------------------------------------------
! Brief: Augment the HF Jacobian with penalty terms that
!        enforce equal fluid pressure at crack junctions.
!
! Parameters:
!   Input:  Counter, Total_FD, num_FreeD, num_free_CalP - sizes
!           freeDOF, freeDOF_HF - DOF index maps
!   In/Out: NR_Deri           - Jacobian (modified in place)
!
! Notes:   Penalty stiffness k = 1.0D4 * max|NR_Deri entries|;
!          applied to every pair listed in Cracks_JunPair.
!-----------------------------------------------------------

subroutine Cal_HF_Jacobian_PenaltyFunc(Counter,NR_Deri,Total_FD, num_FreeD, num_free_CalP, freeDOF,freeDOF_HF)
! Penalty function method for handling junction points with the same water pressure
! For details on the penalty function method, see 'The Finite Element Method in Engineering'
! translated by Zeng Pan, 3rd edition, Section 3.8

!     ----------------------------
!     Read public variable module
!     ----------------------------
use Global_Float_Type
use Global_Crack
use Global_Crack_Common

implicit none
integer, intent(in)::Counter,Total_FD,num_FreeD, num_free_CalP
integer, intent(in)::freeDOF(Total_FD), freeDOF_HF(num_Tol_CalP_Water)
real(kind=FT),intent(inout):: NR_Deri(Total_FD+num_Tol_CalP_Water,Total_FD+num_Tol_CalP_Water)
real(kind=FT) k,max_K(3)
integer i_Pair,CalP_1,CalP_2
print *,'    Modifing the Jacobian matrix by the penalty function' //' method......'

max_K(1) = maxval(abs(NR_Deri(freeDOF(1:num_FreeD), freeDOF_HF(1:num_free_CalP))))
max_K(2) = maxval(abs(NR_Deri(freeDOF_HF(1:num_free_CalP), freeDOF(1:num_FreeD))))
max_K(3) = maxval(abs(NR_Deri(freeDOF_HF(1:num_free_CalP), freeDOF_HF(1:num_free_CalP))))
k = 1.0D4*maxval(max_K)
do i_Pair = 1,Num_JunPair
    CalP_1 = Cracks_JunPair(i_Pair,1)
    CalP_2 = Cracks_JunPair(i_Pair,2)


    NR_Deri(Total_FD+CalP_1,Total_FD+CalP_1) = NR_Deri(Total_FD+CalP_1,Total_FD+CalP_1) + k
    NR_Deri(Total_FD+CalP_1,Total_FD+CalP_2) = NR_Deri(Total_FD+CalP_1,Total_FD+CalP_2) - k
    NR_Deri(Total_FD+CalP_2,Total_FD+CalP_1) = NR_Deri(Total_FD+CalP_2,Total_FD+CalP_1) - k
    NR_Deri(Total_FD+CalP_2,Total_FD+CalP_2) = NR_Deri(Total_FD+CalP_2,Total_FD+CalP_2) + k
end do
return 
end subroutine Cal_HF_Jacobian_PenaltyFunc     
