 
      subroutine Tool_Volume_Pyramid(A,B,C,D,E,Volume)
       
      use Global_Float_Type
      IMPLICIT NONE

      real(kind=FT),intent(in)::A(3),B(3),C(3),D(3),E(3)
      real(kind=FT),intent(out)::Volume
      
      real(kind=FT) Area_BCDE,Area_BCD,Area_CDE
      real(kind=FT) Distance
      
      call Tool_Area_Tri_3D(B,C,D,Area_BCD)
      call Tool_Area_Tri_3D(C,D,E,Area_CDE)
      
      Area_BCDE = Area_BCD + Area_CDE
      
      
      
      if(Area_BCD>=Tol_20) then
          call Tool_Dis_Point_to_3D_Tri_only_Dis(A,B,C,D,Distance)
      else
          call Tool_Dis_Point_to_3D_Tri_only_Dis(A,C,D,E,Distance)
      endif
      
     
      Volume = ONE/THR*Area_BCDE*abs(Distance)
      
      RETURN
      END subroutine Tool_Volume_Pyramid
