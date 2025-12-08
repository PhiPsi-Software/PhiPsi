 
      subroutine Tool_Yes_Point_on_Line_Segments_3D_Tol(Point,
     &                                     Points_Seg,Tol,
     &                                     number_Points,Yes_on)

      use Global_Float_Type 
      implicit none
      integer,intent(in) :: number_Points
      real(kind=FT),intent(in) :: Point(3),Points_Seg(number_Points,3),
     &                            Tol
      logical,intent(out) :: Yes_on
      integer i_Seg
      real(kind=FT) c_A(3),c_B(3)
      logical c_Yes_on
      
      Yes_on = .False.
      
      do i_Seg=1,number_Points-1
         c_A = Points_Seg(i_Seg,1:3)
         c_B = Points_Seg(i_Seg+1,1:3)
         call Tool_Yes_Point_on_Line_Segment_3D_Tol(c_A,c_B,Point,Tol,
     &                                              c_Yes_on)
         if(c_Yes_on)then
             Yes_on = .True.
             return
         endif
      enddo
      
      return 
      end SUBROUTINE Tool_Yes_Point_on_Line_Segments_3D_Tol                      
