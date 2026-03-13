subroutine Cal_JDomain_Element_of_Point_3D_Simple_HEX8(TipCoor, rJ, J_Elem,Max_J_Elem, num_J_Elem, Q_Elem_Nodes)
use Global_Float_Type
use Global_Model
implicit none
integer,intent(in) :: Max_J_Elem
real(kind=FT),intent(in) :: TipCoor(3), rJ
integer,intent(out) :: J_Elem(Max_J_Elem), num_J_Elem
real(kind=FT),intent(out) :: Q_Elem_Nodes(num_Elem,8)

integer :: i_Ele, iN
integer :: NN(8)
real(kind=FT) :: X(8),Y(8),Z(8)
real(kind=FT) :: xc(3), distc, distn, qn

Q_Elem_Nodes(:,:) = ZR
num_J_Elem = 0

do i_Ele = 1, Num_Elem
  
  
  NN = G_NN(1:8,i_Ele)
  X  = G_X_NODES(1:8,i_Ele)
  Y  = G_Y_NODES(1:8,i_Ele)
  Z  = G_Z_NODES(1:8,i_Ele)

  xc(1) = sum(X)/8.0D0
  xc(2) = sum(Y)/8.0D0
  xc(3) = sum(Z)/8.0D0

  distc = sqrt( (xc(1)-TipCoor(1))**2 + (xc(2)-TipCoor(2))**2 + (xc(3)-TipCoor(3))**2 )
  if (distc <= rJ) then
    num_J_Elem = num_J_Elem + 1
    if (num_J_Elem > size(J_Elem)) then
      num_J_Elem = num_J_Elem - 1
      exit
    endif
    
    
    J_Elem(num_J_Elem) = i_Ele

    do iN=1,8
      distn = sqrt( (X(iN)-TipCoor(1))**2 + (Y(iN)-TipCoor(2))**2 + (Z(iN)-TipCoor(3))**2 )
      qn = ONE - distn/rJ
      if (qn < ZR) qn = ZR
      if (qn > ONE) qn = ONE
      Q_Elem_Nodes(i_Ele,iN) = qn
    enddo
  endif
enddo

end subroutine Cal_JDomain_Element_of_Point_3D_Simple_HEX8