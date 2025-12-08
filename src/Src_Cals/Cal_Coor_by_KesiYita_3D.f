 
      subroutine Cal_Coor_by_KesiYita_3D(kesi,yita,zeta,X_NODES,Y_NODES,
     &                                   Z_NODES,Out_x,Out_y,Out_z)
      use Global_Float_Type
      implicit none
      real(kind=FT),intent(in)::kesi,yita,zeta
      real(kind=FT),intent(in)::X_NODES(8),Y_NODES(8),Z_NODES(8)
      real(kind=FT),intent(out)::Out_x,Out_y,Out_z
      real(kind=FT) N(8),N_all(3,24)
      
      call Cal_N_3D(kesi,yita,zeta,N_all)
      
      N(1:8) = N_all(1,1:24:3)
      
      Out_x  = N(1)*X_NODES(1)+N(2)*X_NODES(2)+
     &         N(3)*X_NODES(3)+N(4)*X_NODES(4)+
     &         N(5)*X_NODES(5)+N(6)*X_NODES(6)+
     &         N(7)*X_NODES(7)+N(8)*X_NODES(8)
      Out_y  = N(1)*Y_NODES(1)+N(2)*Y_NODES(2)+
     &         N(3)*Y_NODES(3)+N(4)*Y_NODES(4)+
     &         N(5)*Y_NODES(5)+N(6)*Y_NODES(6)+
     &         N(7)*Y_NODES(7)+N(8)*Y_NODES(8)
      Out_z  = N(1)*Z_NODES(1)+N(2)*Z_NODES(2)+
     &         N(3)*Z_NODES(3)+N(4)*Z_NODES(4)+
     &         N(5)*Z_NODES(5)+N(6)*Z_NODES(6)+
     &         N(7)*Z_NODES(7)+N(8)*Z_NODES(8)
     
      return 
      end SUBROUTINE Cal_Coor_by_KesiYita_3D                          
