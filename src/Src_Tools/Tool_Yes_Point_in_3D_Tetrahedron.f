 
      subroutine Tool_Yes_Point_in_3D_Tetrahedron(Point,
     &                            A,B,C,D,
     &                            Yes_in,Yes_on)

      use Global_Float_Type
      implicit none
      real(kind=FT),intent(in)::Point(3),A(3),B(3),C(3),D(3)
      logical,intent(out):: Yes_in,Yes_on
      real(kind=FT) alpha,beta,gama,delta
      real(kind=FT) c_Distance1,c_Distance2
      real(kind=FT) c_max,c_min
      real(kind=FT) coor_x_max,coor_x_min,coor_y_max,coor_y_min
      real(kind=FT) coor_z_max,coor_z_min
      real(kind=FT) Tool_Function_2Point_Dis_3D
      
      Yes_in = .False.
      Yes_on = .False.
      
      
      
      coor_x_max = max(A(1),B(1),C(1),D(1))
      if(Point(1) > (coor_x_max+Tol_8))then
          return
      endif
      coor_x_min = min(A(1),B(1),C(1),D(1))
      if(Point(1) < (coor_x_min-Tol_8))then
          return
      endif
      coor_y_max = max(A(2),B(2),C(2),D(2))
      if(Point(2) > (coor_y_max+Tol_8))then
          return
      endif
      coor_y_min = min(A(2),B(2),C(2),D(2))
      if(Point(2) < (coor_y_min-Tol_8))then
          return
      endif
      coor_z_max = max(A(3),B(3),C(3),D(3))
      if(Point(3) > (coor_z_max+Tol_8))then
          return
      endif
      coor_z_min = min(A(3),B(3),C(3),D(3))
      if(Point(3) < (coor_z_min-Tol_8))then
          return
      endif

      if(Tool_Function_2Point_Dis_3D(A,B)<Tol_20) return
      if(Tool_Function_2Point_Dis_3D(A,C)<Tol_20) return
      if(Tool_Function_2Point_Dis_3D(A,D)<Tol_20) return
      if(Tool_Function_2Point_Dis_3D(B,C)<Tol_20) return
      if(Tool_Function_2Point_Dis_3D(B,D)<Tol_20) return
      if(Tool_Function_2Point_Dis_3D(C,D)<Tol_20) return
      

      call Tool_Dis_Point_to_3D_Tri_only_Dis(Point,B,C,D,c_Distance1)
      call Tool_Dis_Point_to_3D_Tri_only_Dis(A,B,C,D,c_Distance2)
       
      
      alpha = c_Distance1/c_Distance2
      
      if(alpha > (ONE+Tol_10)) then
          return
      endif
      if(alpha < (ZR-Tol_10)) then
          return
      endif
      call Tool_Dis_Point_to_3D_Tri_only_Dis(Point,A,C,D,c_Distance1)
      call Tool_Dis_Point_to_3D_Tri_only_Dis(B,A,C,D,c_Distance2)   
      
      beta = c_Distance1/c_Distance2
      
      
      if(beta > (ONE+Tol_10)) then
          return
      endif
      if(beta < (ZR-Tol_10)) then
          return
      endif
      
      
      call Tool_Dis_Point_to_3D_Tri_only_Dis(Point,A,B,D,c_Distance1)
      call Tool_Dis_Point_to_3D_Tri_only_Dis(C,A,B,D,c_Distance2)  
      
      gama = c_Distance1/c_Distance2     
      if(gama > (ONE+Tol_10)) then
          return
      endif
      if(gama < (ZR-Tol_10)) then
          return
      endif
      
      
      call Tool_Dis_Point_to_3D_Tri_only_Dis(Point,A,B,C,c_Distance1)
      call Tool_Dis_Point_to_3D_Tri_only_Dis(D,A,B,C,c_Distance2)   
      
      delta = c_Distance1/c_Distance2
      if(delta > (ONE+Tol_10)) then
          return
      endif
      if(delta < (ZR-Tol_10)) then
          return
      endif
      
      
      c_max = max(alpha,beta,gama,delta)
      c_min = min(alpha,beta,gama,delta)
      
      
      if (c_max <= (ONE+Tol_10) .and. c_min>=(ZR-Tol_10)) then
          Yes_in = .True.
      endif

      if ((abs(c_max-ONE) <Tol_10) .and. (abs(c_min-ONE) <Tol_10))then
          Yes_on = .True.
      endif
      
      return 
      end SUBROUTINE Tool_Yes_Point_in_3D_Tetrahedron                   
