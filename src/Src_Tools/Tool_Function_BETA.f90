!-----------------------------------------------------------
! Brief: Compute the Beta function B(P,Q)
!
! Parameters:
!   Input:  P  - first parameter (> 0)
!   Input:  Q  - second parameter (> 0)
!   Output: BT - computed Beta value
!
! Notes:   B(P,Q) = Gamma(P)*Gamma(Q)/Gamma(P+Q).
!-----------------------------------------------------------

subroutine Tool_Function_BETA(P,Q,BT)
!       ==========================================
!       Purpose: Compute the beta function B(p,q)
!       Input :  p  --- Parameter  ( p > 0 )
!                q  --- Parameter  ( q > 0 )
!       Output:  BT --- B(p,q)
!       Routine called: GAMMA
!       ==========================================

use Global_Float_Type
IMPLICIT real(kind=FT) (A-H,O-Z)

CALL Tool_Function_GAMMA(P,GP)
CALL Tool_Function_GAMMA(Q,GQ)
PPQ=P+Q
CALL Tool_Function_GAMMA(PPQ,GPQ)
BT=GP*GQ/GPQ

return
end SUBROUTINE Tool_Function_BETA
