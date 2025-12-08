 
      subroutine Tool_Yes_Point_on_3D_Triangle(P,A,B,C,Yes_on)

      use Global_Float_Type
      implicit none
      real(kind=FT),intent(in)::P(3),A(3),B(3),C(3)
      logical,intent(out):: Yes_on
      real(kind=FT) coor_x_max,coor_x_min,coor_y_max,
     &              coor_y_min,coor_z_max,coor_z_min
      real(kind=FT) Area_All,Area_1,Area_2,Area_3
      
      Yes_on = .False.
      
      coor_x_max = max(A(1),B(1),C(1))
      coor_x_min = min(A(1),B(1),C(1))
      coor_y_max = max(A(2),B(2),C(2))
      coor_y_min = min(A(2),B(2),C(2))
      coor_z_max = max(A(3),B(3),C(3))
      coor_z_min = min(A(3),B(3),C(3))
      
      
      call Tool_Area_Tri_3D(A,B,C,Area_All)
      call Tool_Area_Tri_3D(A,B,P,Area_1)
      call Tool_Area_Tri_3D(A,C,P,Area_2)
      call Tool_Area_Tri_3D(B,C,P,Area_3)
      
      if(abs(Area_1+Area_2+Area_3-Area_All)<=Tol_11)then
          Yes_on = .True.
      endif
      
      return 
      end SUBROUTINE Tool_Yes_Point_on_3D_Triangle               
