 
subroutine Cal_Ele_Stiffness_Matrix_3D_8nodes(i_E,num_Gauss,X_NODES,Y_NODES,Z_NODES, &
                             c_D,kesi,yita,zeta,weight,localK)
use Global_Float_Type
implicit none

integer,intent(in)::i_E,num_Gauss
real(kind=FT),intent(in)::c_D(6,6)
real(kind=FT),intent(in)::kesi(num_Gauss),yita(num_Gauss),zeta(num_Gauss),weight(num_Gauss)
real(kind=FT),intent(in)::X_NODES(8),Y_NODES(8),Z_NODES(8)
real(kind=FT),intent(out)::localK(24,24)          

real(kind=FT) J(3,3),detJ,dNdkesi(8,3),dNdx(8,3),B_FEM(6,24)
real(kind=FT) Inverse_J(3,3)
integer i_G,i_N
real(kind=FT) temp(3,8),Coor(3,8)
real(kind=FT) one_p_yita,one_m_yita
real(kind=FT) one_p_zeta,one_m_zeta
real(kind=FT) one_p_kesi,one_m_kesi
#ifdef Silverfrost
real(kind=FT) :: tmp1(24,6), tmp2(24,24)
#endif
localK(1:24,1:24)=ZR
    
do i_G = 1,num_Gauss      
  dNdkesi(1:8,1:3) = ZR
  one_m_yita = ONE -yita(i_G)
  one_p_yita = ONE +yita(i_G)
  one_m_zeta = ONE -zeta(i_G)
  one_p_zeta = ONE +zeta(i_G)
  one_m_kesi = ONE -kesi(i_G)
  one_p_kesi = ONE +kesi(i_G)          
  temp(1,1:8)=[-(one_m_yita)*(one_m_zeta), (one_m_yita)*(one_m_zeta), &
                (one_p_yita)*(one_m_zeta),-(one_p_yita)*(one_m_zeta), &
               -(one_m_yita)*(one_p_zeta), (one_m_yita)*(one_p_zeta), &
                (one_p_yita)*(one_p_zeta),-(one_p_yita)*(one_p_zeta)]  
  temp(2,1:8)=[-(one_m_kesi)*(one_m_zeta),-(one_p_kesi)*(one_m_zeta), &
                (one_p_kesi)*(one_m_zeta), (one_m_kesi)*(one_m_zeta), &
               -(one_m_kesi)*(one_p_zeta),-(one_p_kesi)*(one_p_zeta), &
                (one_p_kesi)*(one_p_zeta), (one_m_kesi)*(one_p_zeta)] 
  temp(3,1:8)=[-(one_m_kesi)*(one_m_yita),-(one_p_kesi)*(one_m_yita), &
               -(one_p_kesi)*(one_p_yita),-(one_m_kesi)*(one_p_yita), &
                (one_m_kesi)*(one_m_yita), (one_p_kesi)*(one_m_yita), &
                (one_p_kesi)*(one_p_yita), (one_m_kesi)*(one_p_yita)]        

  temp = temp/EIG
  dNdkesi = transpose(temp)
  Coor(1,:)=[X_NODES(1),X_NODES(2),X_NODES(3),X_NODES(4),X_NODES(5),X_NODES(6),X_NODES(7),X_NODES(8)]
  Coor(2,:)=[Y_NODES(1),Y_NODES(2),Y_NODES(3),Y_NODES(4),Y_NODES(5),Y_NODES(6),Y_NODES(7),Y_NODES(8)]
  Coor(3,:)=[Z_NODES(1),Z_NODES(2),Z_NODES(3),Z_NODES(4),Z_NODES(5),Z_NODES(6),Z_NODES(7),Z_NODES(8)]
  J = MATMUL(Coor,dNdkesi)
  detJ = J(1,1)*J(2,2)*J(3,3) +J(1,2)*J(2,3)*J(3,1) +J(2,1)*J(3,2)*J(1,3) - &
         J(1,1)*J(2,3)*J(3,2) -J(1,2)*J(2,1)*J(3,3) -J(1,3)*J(2,2)*J(3,1)
  
  call Matrix_Inverse_3x3(J,Inverse_J)  
  
  dNdx = MATMUL(dNdkesi,Inverse_J(1:3,1:3))   
  
  B_FEM(1:6,1:24) =ZR
  do i_N = 1,8
      B_FEM(1,(i_N-1)*3 + 1) = dNdx(i_N,1)
      B_FEM(2,(i_N-1)*3 + 2) = dNdx(i_N,2)
      B_FEM(3,(i_N-1)*3 + 3) = dNdx(i_N,3)
      
      B_FEM(4,(i_N-1)*3 + 1) = dNdx(i_N,2)
      B_FEM(4,(i_N-1)*3 + 2) = dNdx(i_N,1)
      
      B_FEM(5,(i_N-1)*3 + 2) = dNdx(i_N,3)
      B_FEM(5,(i_N-1)*3 + 3) = dNdx(i_N,2)
      
      B_FEM(6,(i_N-1)*3 + 1) = dNdx(i_N,3)
      B_FEM(6,(i_N-1)*3 + 3) = dNdx(i_N,1)          
  end do     
#ifndef Silverfrost
  localK = localK + weight(i_G)*detJ*MATMUL(MATMUL(transpose(B_FEM),c_D),B_FEM) 
#endif
#ifdef Silverfrost
  tmp1 = MATMUL(transpose(B_FEM), c_D)
  tmp2 = MATMUL(tmp1, B_FEM)
  localK = localK + weight(i_G)*detJ*tmp2
#endif
end do


return 
end SUBROUTINE Cal_Ele_Stiffness_Matrix_3D_8nodes  
