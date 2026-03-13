 
subroutine Cal_Shape_Functions_Hex8(xi, eta, zeta, N, dN_dxi, dN_deta, dN_dzeta)

use Global_Float_Type      
implicit none
real(kind=FT), intent(in) :: xi, eta, zeta
real(kind=FT), intent(out) :: N(8), dN_dxi(8), dN_deta(8), dN_dzeta(8)

real(kind=FT) :: xip, xim, etap, etam, zetap, zetam

xip = ONE + xi
xim = ONE - xi
etap = ONE + eta
etam = ONE - eta
zetap = ONE + zeta
zetam = ONE - zeta

N(1) = ZP125 * xim * etam * zetam
N(2) = ZP125 * xip * etam * zetam
N(3) = ZP125 * xip * etap * zetam
N(4) = ZP125 * xim * etap * zetam
N(5) = ZP125 * xim * etam * zetap
N(6) = ZP125 * xip * etam * zetap
N(7) = ZP125 * xip * etap * zetap
N(8) = ZP125 * xim * etap * zetap

dN_dxi(1) = -ZP125 * etam * zetam
dN_dxi(2) =  ZP125 * etam * zetam
dN_dxi(3) =  ZP125 * etap * zetam
dN_dxi(4) = -ZP125 * etap * zetam
dN_dxi(5) = -ZP125 * etam * zetap
dN_dxi(6) =  ZP125 * etam * zetap
dN_dxi(7) =  ZP125 * etap * zetap
dN_dxi(8) = -ZP125 * etap * zetap

dN_deta(1) = -ZP125 * xim * zetam
dN_deta(2) = -ZP125 * xip * zetam
dN_deta(3) =  ZP125 * xip * zetam
dN_deta(4) =  ZP125 * xim * zetam
dN_deta(5) = -ZP125 * xim * zetap
dN_deta(6) = -ZP125 * xip * zetap
dN_deta(7) =  ZP125 * xip * zetap
dN_deta(8) =  ZP125 * xim * zetap

dN_dzeta(1) = -ZP125 * xim * etam
dN_dzeta(2) = -ZP125 * xip * etam
dN_dzeta(3) = -ZP125 * xip * etap
dN_dzeta(4) = -ZP125 * xim * etap
dN_dzeta(5) =  ZP125 * xim * etam
dN_dzeta(6) =  ZP125 * xip * etam
dN_dzeta(7) =  ZP125 * xip * etap
dN_dzeta(8) =  ZP125 * xim * etap
    
end subroutine Cal_Shape_Functions_Hex8
