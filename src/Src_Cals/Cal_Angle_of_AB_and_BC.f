 
      subroutine Cal_Angle_of_AB_and_BC(Line_AB,Line_BC,angle_AB_BC)
      use Global_Float_Type
      use Global_Common
      implicit none
      real(kind=FT),intent(in)::Line_AB(2,2),Line_BC(2,2)
      real(kind=FT),intent(out)::angle_AB_BC
      real(kind=FT) a_x,a_y,b_x,b_y,c_x,c_y,tem1,tem2,angle
      real(kind=FT) tt1
      real(kind=FT) a(3),b(3),Direction(3)
      
      if  ((Line_AB(2,1) .ne. Line_BC(1,1)) .or.
     &     (Line_AB(2,2) .ne. Line_BC(1,2))) then
          print *,'    ********************************'
     &                 //'*********************************'    
          print *,'    Attention should be paid here!'
          print *,'    This situation rarely' 
     &              //' appears in subroutine Cal_Angle_of_AB_and_BC!'
          print *,'    ********************************'
     &                 //'*********************************'  
      end if     
      a_x = Line_AB(1,1)
      a_y = Line_AB(1,2)
      b_x = Line_AB(2,1)
      b_y = Line_AB(2,2)
      c_x = Line_BC(2,1)
      c_y = Line_BC(2,2)

      tem1 = sqrt((a_x-b_x)**2 + (a_y-b_y)**2)
      tem2 = sqrt((c_x-b_x)**2 + (c_y-b_y)**2)
      tt1 = DOT_PRODUCT([a_x-b_x,a_y-b_y],[c_x-b_x,c_y-b_y])/(tem1*tem2)
      if (tt1 < -ONE) then
          tt1 = -ONE
      elseif(tt1 > ONE)then
          tt1 = ONE
      endif
      angle  = acos(tt1)
      
      a = [a_x-b_x,a_y-b_y,ZR]
      b = [c_x-b_x,c_y-b_y,ZR]
      
      call Vector_Cross_Product_3(a,b,Direction)
      if (Direction(3) >= ZR) then
          angle_AB_BC = angle
      elseif (Direction(3) < ZR) then
          angle_AB_BC = TWO*pi-angle
      end if
      return 
      end SUBROUTINE Cal_Angle_of_AB_and_BC              
