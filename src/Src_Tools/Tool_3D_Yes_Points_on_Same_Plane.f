 
      subroutine Tool_3D_Yes_Points_on_Same_Plane(num_points,
     &                                            Points,Yes_on)

      use Global_Float_Type
      use Global_Common
      implicit none
      integer,intent(in)::num_points
      real(kind=FT),intent(in)::Points(num_points,3)
      logical,intent(out)::Yes_on
      real(kind=FT) Tri_P1(3),Tri_P2(3),Tri_P3(3)
      real(kind=FT) c_Distance,c_PER(3)
      integer i_Point
      logical c_Yes,c_Yes_PER_in,c_Yes_PER_on
      
      Yes_on = .True.
      
      if(num_points<=3)then
          print *,'    Error :: Number of points is less than 4!'
          print *,'             in Tool_3D_Yes_Points_on_Same_Plane.f!'
          call Warning_Message('S',Keywords_Blank)   
          return
      endif      
      
      Tri_P1(1:3) = Points(1,1:3)
      Tri_P2(1:3) = Points(2,1:3)
      Tri_P3(1:3) = Points(3,1:3)
      call Tool_Yes_Point_on_Line_Segment_3D(Tri_P1,Tri_P2,Tri_P3,c_Yes)
      if(c_Yes)then
          print *,'    Error :: Points 1-3 are collinear!'
          print *,'             in Tool_3D_Yes_Points_on_Same_Plane.f!'
          call Warning_Message('S',Keywords_Blank)   
          return
      endif
      do i_Point = 1,num_points-3
          call Tool_Dis_Point_to_3D_Tri(Points(i_Point+3,1:3),
     &                            Tri_P1,Tri_P2,Tri_P3,
     &                            c_Distance,c_PER(1:3),
     &                            c_Yes_PER_in,c_Yes_PER_on)
          if(c_Distance > Tol_10) then
              Yes_on = .False.
              return
          endif
      enddo

      
      return 
      end SUBROUTINE Tool_3D_Yes_Points_on_Same_Plane                  
