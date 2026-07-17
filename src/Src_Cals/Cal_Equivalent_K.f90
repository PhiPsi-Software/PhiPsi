!-----------------------------------------------------------
! Brief: Compute the equivalent stress intensity factor from KI and KII for a crack tip.
!
! Parameters:
!   Input:  i_C     - crack index
!           i_Tip   - crack-tip index
!           c_KI    - mode-I stress intensity factor
!           c_KII   - mode-II stress intensity factor
!   Output: c_Kequ  - equivalent SIF used for propagation criterion
!
! Notes:   Returns -1e14 (a sentinel) for boundary / junction crack tips that are
!          not allowed to propagate freely. For free tips the maximum-tangential-
!          stress based formula in terms of KI and KII is applied.
!-----------------------------------------------------------

subroutine Cal_Equivalent_K(i_C,i_Tip,c_KI,c_KII,c_Kequ)

!     This function calculates equivalent stress intensity factors.
use Global_Float_Type
use Global_Crack
implicit none
integer,intent(in)::i_C,i_Tip
real(kind=FT),intent(in)::c_KI,c_KII
real(kind=FT),intent(out)::c_Kequ  
real(kind=FT) delta_angle,tem1


if ((Crack_Tip_Type(i_C,i_Tip) .ne. 1) .and. &
(Crack_Tip_Type(i_C,i_Tip) .ne. -1).and. &
(Crack_Tip_Type(i_C,i_Tip) .ne. -2)) then
delta_angle = TWO* atan(-2*c_KII/c_KI/(ONE+sqrt(ONE+EIG*(c_KII/c_KI)**2)))
tem1  = delta_angle/TWO
c_Kequ = cos(tem1)*(c_KI*cos(tem1)**2-1.5D0*sin(delta_angle))
else
    c_Kequ = -100.0D12
endif
return 
end subroutine Cal_Equivalent_K                 
