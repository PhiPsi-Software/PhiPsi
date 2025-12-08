 
      subroutine Tool_Yes_Line_on_3D_Triangle(P1,P2,A,B,C,Yes_on,
     &                                        Point)

      use Global_Float_Type
      implicit none
      real(kind=FT),intent(in)::P1(3),P2(3),A(3),B(3),C(3)
      logical,intent(out):: Yes_on
      real(kind=FT),intent(out)::Point(3)
      logical Yes_P1_on,Yes_P2_on
      real(kind=FT) Tool_Function_Dis_Point_to_3D_Tri_v2,c_Dis
      
      Yes_on = .False.
      Point(1:3)  = -TEN_15
      
      call Tool_Yes_Point_on_3D_Triangle(P1,A,B,C,Yes_P1_on)
      call Tool_Yes_Point_on_3D_Triangle(P2,A,B,C,Yes_P2_on)
      
      if((Yes_P1_on .eqv. .False.) .and. (Yes_P2_on .eqv. .False.))then
          return
      endif
      
      if((Yes_P1_on .eqv. .True.) .and. (Yes_P2_on .eqv. .True.))then
          Yes_on = .True.
          Point  = P1
          return
      endif      
      
      if((Yes_P1_on .eqv. .True.) .and. (Yes_P2_on .eqv. .False.))then
          c_Dis = Tool_Function_Dis_Point_to_3D_Tri_v2(P2,A,B,C)
          if(abs(c_Dis)<=Tol_11) then
              Yes_on = .True.
              Point  = P1
              return
          endif
      endif       
      
      if((Yes_P2_on .eqv. .True.) .and. (Yes_P1_on .eqv. .False.))then
          c_Dis = Tool_Function_Dis_Point_to_3D_Tri_v2(P1,A,B,C)
          if(abs(c_Dis)<=Tol_11) then
              Yes_on = .True.
              Point  = P2
              return
          endif
      endif         
      return 
      end SUBROUTINE Tool_Yes_Line_on_3D_Triangle         
