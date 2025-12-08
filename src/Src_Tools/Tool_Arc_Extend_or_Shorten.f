 
      subroutine Tool_Arc_Extend_or_Shorten(
     &                      delta_L_Arc,
     &                      Old_Arc_Cr,
     &                      New_Arc_Cr,Extend_Status)
      use Global_Float_Type
      use Global_Common
      
      implicit none
      real(kind=FT),intent(in)::delta_L_Arc,Old_Arc_Cr(11)
      real(kind=FT),intent(out)::New_Arc_Cr(11)
      logical,intent(out)::Extend_Status
      
      real(kind=FT) o_x,o_y
      real(kind=FT) r
      real(kind=FT) c_Arc_Center(2)
      real(kind=FT) c_Arc_Direction
      real(kind=FT) Point_I(2)
      real(kind=FT) P1_P2(2,2)
      integer       Point_Status
      real(kind=FT) c_Arc_Start_Point(2)
      logical       c_Yes_feasible
      real(kind=FT) Arc_Radian_P1,Arc_Radian_P2
      real(kind=FT) Old_Radian,Old_Arc_L
      real(kind=FT) c_Arc_Radian_S_P1,c_Arc_Radian_E_P1
      real(kind=FT) c_Arc_Radian_S_P2,c_Arc_Radian_E_P2
      
 1001 FORMAT(5X,'-- Error :: abs(-delta_L_Arc) is > Old_Arc_L!') 
 1002 FORMAT(5X,'            Error in Tool_Arc_Extend_or_Shorten.f!')  
 
      Extend_Status = .False.
      o_x               = Old_Arc_Cr(1)
      o_y               = Old_Arc_Cr(2)
      c_Arc_Center(1:2) = [o_x,o_Y]
      c_Arc_Direction   = Old_Arc_Cr(3)
      r                 = Old_Arc_Cr(4)
      Old_Radian        = Old_Arc_Cr(7)
      Old_Arc_L         = Old_Radian*r
      c_Arc_Start_Point = [Old_Arc_Cr(8),Old_Arc_Cr(9)]
      Point_I(1)        = Old_Arc_Cr(10)
      Point_I(2)        = Old_Arc_Cr(11)
      
      if(delta_L_Arc < Tol_11)then
          if(abs(delta_L_Arc) > Old_Arc_L)then
           write (*,1001)
           write (*,1002)
           call Warning_Message('S',Keywords_Blank)
          endif
      endif
      
      
      call Tool_2_Points_of_Cir_give_Point_L(
     &                      o_x,o_y,r,
     &                      Point_I,delta_L_Arc,Point_Status,
     &                      P1_P2)
      if(Point_Status ==1 .and.(sum(P1_P2(2,1:2))>Tol_11)) then       
          call Tool_Arc_r_and_Radian_Given_Coors(c_Arc_Direction,
     &                           c_Arc_Start_Point,
     &                           P1_P2(1,1:2),c_Arc_Center,
     &                           c_Yes_feasible,
     &                           r,c_Arc_Radian_S_P1,c_Arc_Radian_E_P1,
     &                           Arc_Radian_P1)
          call Tool_Arc_r_and_Radian_Given_Coors(c_Arc_Direction,
     &                           c_Arc_Start_Point,
     &                           P1_P2(2,1:2),c_Arc_Center,
     &                           c_Yes_feasible,
     &                           r,c_Arc_Radian_S_P2,c_Arc_Radian_E_P2,
     &                           Arc_Radian_P2)
          if(delta_L_Arc >= Tol_11)then
              if(Arc_Radian_P1>=Arc_Radian_P2)then
                  New_Arc_Cr     = Old_Arc_Cr
                  New_Arc_Cr(5)  = c_Arc_Radian_S_P1
                  New_Arc_Cr(6)  = c_Arc_Radian_E_P1
                  New_Arc_Cr(7)  = Arc_Radian_P1
                  New_Arc_Cr(10) = P1_P2(1,1)
                  New_Arc_Cr(11) = P1_P2(1,2)
                  Extend_Status  = .True.
              elseif(Arc_Radian_P1 < Arc_Radian_P2)then
                  New_Arc_Cr     = Old_Arc_Cr
                  New_Arc_Cr(5)  = c_Arc_Radian_S_P2
                  New_Arc_Cr(6)  = c_Arc_Radian_E_P2
                  New_Arc_Cr(7)  = Arc_Radian_P2
                  New_Arc_Cr(10) = P1_P2(2,1)
                  New_Arc_Cr(11) = P1_P2(2,2) 
                  Extend_Status  = .True.
              endif
          endif
          if(delta_L_Arc <= (-Tol_11))then
              if(Arc_Radian_P1 <= Arc_Radian_P2)then
                  New_Arc_Cr     = Old_Arc_Cr
                  New_Arc_Cr(5)  = c_Arc_Radian_S_P1
                  New_Arc_Cr(6)  = c_Arc_Radian_E_P1
                  New_Arc_Cr(7)  = Arc_Radian_P1
                  New_Arc_Cr(10) = P1_P2(1,1)
                  New_Arc_Cr(11) = P1_P2(1,2)
                  Extend_Status  = .True.
              elseif(Arc_Radian_P1 > Arc_Radian_P2)then
                  New_Arc_Cr     = Old_Arc_Cr
                  New_Arc_Cr(5)  = c_Arc_Radian_S_P2
                  New_Arc_Cr(6)  = c_Arc_Radian_E_P2
                  New_Arc_Cr(7)  = Arc_Radian_P2
                  New_Arc_Cr(10) = P1_P2(2,1)
                  New_Arc_Cr(11) = P1_P2(2,2) 
                  Extend_Status  = .True.
              endif
          endif
      endif
      return 
      end SUBROUTINE Tool_Arc_Extend_or_Shorten                     
