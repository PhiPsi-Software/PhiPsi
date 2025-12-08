 
      subroutine Tool_2_Points_of_Cir_give_Point_L(
     &                      x,y,r,
     &                      Point_I,Length_of_Arc,Point_Status,
     &                      P1_P2)

      use Global_Float_Type
      use Global_Elem_Area_Vol
      
      implicit none
      real(kind=FT),intent(in)::x,y,r,Point_I(2),Length_of_Arc
      real(kind=FT),intent(out)::P1_P2(2,2)
      integer,intent(out)::Point_Status
      real(kind=FT) Tool_Function_2Point_Dis
      real(kind=FT) c_r
      real(kind=FT) x1,y1,r1
      real(kind=FT) Radian_Angle
      real(kind=FT) c_Inters(2,2)
      integer c_Circles_Status,c_num_Inters
      
      
      P1_P2(1:2,1:2) = ZR      
      c_r = Tool_Function_2Point_Dis(Point_I,[x,y])
      if(abs(c_r - r)<=Tol_11)then
          Point_Status = 1
          Radian_Angle = Length_of_Arc/(TWO*pi*r)*(TWO*pi)
          x1 = Point_I(1)
          y1 = Point_I(2)
          r1 = TWO*(r*sin(Radian_Angle/TWO))
          call Tool_Intersection_Circle_and_Circle(
     &                      x,y,r,
     &                      x1,y1,r1,c_Circles_Status,
     &                      c_num_Inters,c_Inters)  
          P1_P2 = c_Inters
      else
          Point_Status =-1
          P1_P2(1:2,1:2) = ZR
      endif
      
      
      return 
      end SUBROUTINE Tool_2_Points_of_Cir_give_Point_L                         
