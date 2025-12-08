 
      subroutine Tool_Yes_Point_in_3D_Pyramid_with_Tol(Point,
     &                            A,B,C,D,E,Tol,
     &                            Yes_in,Yes_on)
      use Global_Float_Type
      implicit none
      real(kind=FT),intent(in)::Point(3),
     &                 A(3),B(3),C(3),D(3),E(3),Tol
      logical,intent(out):: Yes_in,Yes_on
      logical c_Yes_in_1,c_Yes_on_1,c_Yes_in_2,c_Yes_on_2
      real(kind=FT) coor_x_max,coor_x_min,coor_y_max,coor_y_min
      real(kind=FT) coor_z_max,coor_z_min
      
      Yes_in = .False.
      Yes_on = .False.
      
      coor_x_max = max(A(1),B(1),C(1),D(1),E(1))
      if(Point(1) > (coor_x_max+Tol_8))then
          return
      endif
      coor_x_min = min(A(1),B(1),C(1),D(1),E(1))
      if(Point(1) < (coor_x_min-Tol_8))then
          return
      endif
      
      coor_y_max = max(A(2),B(2),C(2),D(2),E(2))
      if(Point(2) > (coor_y_max+Tol_8))then
          return
      endif
      coor_y_min = min(A(2),B(2),C(2),D(2),E(2))
      if(Point(2) < (coor_y_min-Tol_8))then
          return
      endif
      
      coor_z_max = max(A(3),B(3),C(3),D(3),E(3))
      if(Point(3) > (coor_z_max+Tol_8))then
          return
      endif
      coor_z_min = min(A(3),B(3),C(3),D(3),E(3))
      if(Point(3) < (coor_z_min-Tol_8))then
          return
      endif
      
      call Tool_Yes_Point_in_3D_Tetrahedron_with_Tol(Point,
     &                            A,B,C,E,Tol,
     &                            c_Yes_in_1,c_Yes_on_1)
      if (c_Yes_in_1) then
          Yes_in = .True.
          return
      endif
      if (c_Yes_on_1) then
          Yes_on = .True.
          return
      endif
      
      call Tool_Yes_Point_in_3D_Tetrahedron_with_Tol(Point,
     &                            A,C,D,E,Tol,
     &                            c_Yes_in_2,c_Yes_on_2)
      if (c_Yes_in_2) then
          Yes_in = .True.
          return
      endif
      if (c_Yes_on_2) then
          Yes_on = .True.
          return
      endif

      
      return 
      end SUBROUTINE Tool_Yes_Point_in_3D_Pyramid_with_Tol                
