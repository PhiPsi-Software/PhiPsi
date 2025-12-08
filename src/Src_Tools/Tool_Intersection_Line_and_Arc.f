 
      subroutine Tool_Intersection_Line_and_Arc(
     &                                   A,B,
     &                                   c_Arc_Crack_Coor,
     &                                   num_Inter,Inter,Yes_Cross) 

      use Global_Float_Type
      use Global_Elem_Area_Vol
      
      implicit none
      real(kind=FT),intent(in)::A(2),B(2),c_Arc_Crack_Coor(11)
      integer,intent(out)::num_Inter
      real(kind=FT),intent(out)::Inter(2,2)
      logical,intent(out)::Yes_Cross
      
      real(kind=FT) o_x,o_y
      real(kind=FT) r
      integer c_num_Inter,State
      real(kind=FT) c_Inter(2,2)
      integer i_Inter
      logical c_Yes_ON
      
      Yes_Cross = .False.

      o_x               = c_Arc_Crack_Coor(1)
      o_y               = c_Arc_Crack_Coor(2)
      r                 = c_Arc_Crack_Coor(4)
     
      c_num_Inter = 0
      call Tool_Intersection_Line_and_Circle(o_x,o_y,
     &                                       r,A,B,
     &                                       c_num_Inter,State,
     &                                       c_Inter)

      num_Inter = 0
      do i_Inter = 1,c_num_Inter
          call Tool_Yes_On_Arc(c_Inter(i_Inter,1),c_Inter(i_Inter,2),
     &                         c_Arc_Crack_Coor,c_Yes_ON)
          if(c_Yes_ON)then
              Yes_Cross = .True.
              num_Inter = num_Inter +1
              Inter(num_Inter,1:2) = c_Inter(i_Inter,1:2)
          endif
      enddo
      
      return 
      end SUBROUTINE Tool_Intersection_Line_and_Arc                         
