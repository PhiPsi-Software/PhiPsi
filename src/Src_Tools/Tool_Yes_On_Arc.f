 
      subroutine Tool_Yes_On_Arc(x,y,c_Arc_Crack_Coor,Yes_ON)

      use Global_Float_Type
      
      implicit none
      real(kind=FT),intent(in):: x,y,c_Arc_Crack_Coor(11)
      logical,intent(out):: Yes_ON
      

      real(kind=FT) Tool_Function_2Point_Dis
      real(kind=FT) o_x,o_y
      real(kind=FT) r,c_r
      real(kind=FT) c_Arc_Center(2)
      real(kind=FT) c_Angle
      real(kind=FT) c_Arc_Radian_S,c_Arc_Radian_E
      real(kind=FT) c_Arc_Direction
      
      Yes_ON = .False.
      
      o_x = c_Arc_Crack_Coor(1)
      o_y = c_Arc_Crack_Coor(2)
      c_Arc_Center(1:2) = [o_x,o_Y]
      c_Arc_Direction  = c_Arc_Crack_Coor(3)
      r   = c_Arc_Crack_Coor(4)
      
      c_r = Tool_Function_2Point_Dis(c_Arc_Center,[x,y])
      
      c_Arc_Radian_S = c_Arc_Crack_Coor(5)
      c_Arc_Radian_E = c_Arc_Crack_Coor(6)
      
      if(abs(c_r - r)<=Tol_11)then
          call Tool_Angle_of_Inclination_of_VectorAB(
     &                          c_Arc_Center,[x,y],c_Angle)
          
          if(c_Arc_Direction > HLF)then
              if (c_Arc_Radian_E >=c_Arc_Radian_S) then
                  if ((c_Angle >=c_Arc_Radian_S) .and.
     &                (c_Angle <=c_Arc_Radian_E)) then
                      Yes_ON = .True.
                  endif
              elseif(c_Arc_Radian_E < c_Arc_Radian_S)then
                  if ((c_Angle >=c_Arc_Radian_S) .and.
     &                (c_Angle <=Con_360)) then
                      Yes_ON = .True.
                  endif
                  if ((c_Angle >=ZR) .and.
     &                (c_Angle <=c_Arc_Radian_E)) then
                      Yes_ON = .True.
                  endif
              endif
          elseif(c_Arc_Direction < (-HLF))then
              if (c_Arc_Radian_E <=c_Arc_Radian_S) then
                  if ((c_Angle >=c_Arc_Radian_E) .and.
     &                (c_Angle <=c_Arc_Radian_S)) then
                      Yes_ON = .True.
                  endif
              elseif(c_Arc_Radian_E > c_Arc_Radian_S)then
                  if ((c_Angle >=c_Arc_Radian_E) .and.
     &                (c_Angle <=Con_360)) then
                      Yes_ON = .True.
                  endif
                  if ((c_Angle >=ZR) .and.
     &                (c_Angle <=c_Arc_Radian_S)) then
                      Yes_ON = .True.
                  endif
              endif  
          endif
      else
          Yes_ON = .False.
      endif
      
      return 
      end SUBROUTINE Tool_Yes_On_Arc                          
