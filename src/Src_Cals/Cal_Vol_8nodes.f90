!-----------------------------------------------------------
! Brief: Volume of a 3D 8-node hexahedral element.
!
! Parameters:
!   Input:  i_Elem - element index (currently unused)
!           x(8), y(8), z(8) - nodal coordinates
!   Output: vol    - element volume (signed)
!
! Notes:   Reorders nodes to the (A,B,C,D,E,F,G,H) convention of
!          Tool_Volume_Hexahedron before computing the volume.
!-----------------------------------------------------------

subroutine Cal_Vol_8nodes(i_Elem,x,y,z,vol)

! Calculate the volume of a 3D 8-node hexahedral element
use Global_Float_Type
implicit none
integer ,intent(in)::i_Elem
real(kind=FT),intent(in):: x(8)
real(kind=FT),intent(in):: y(8)
real(kind=FT),intent(in):: z(8)
real(kind=FT) ,intent(out)::vol
real(kind=FT) A(3),B(3),C(3),D(3),E(3),F(3),G(3),H(3)

A(1) = x(8);A(2) = y(8);A(3) = z(8)
B(1) = x(5);B(2) = y(5);B(3) = z(5)
C(1) = x(6);C(2) = y(6);C(3) = z(6)
D(1) = x(7);D(2) = y(7);D(3) = z(7)

E(1) = x(4);E(2) = y(4);E(3) = z(4)
F(1) = x(1);F(2) = y(1);F(3) = z(1)      
G(1) = x(2);G(2) = y(2);G(3) = z(2)
H(1) = x(3);H(2) = y(3);H(3) = z(3)  

call Tool_Volume_Hexahedron(A,B,C,D,E,F,G,H,vol)

return
END subroutine Cal_Vol_8nodes
