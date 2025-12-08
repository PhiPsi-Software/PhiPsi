 
      subroutine Tool_Dis_Point_to_3D_Quad(Point,
     &                            Quad_P1,Quad_P2,Quad_P3,Quad_P4,
     &                            Distance,PER,Yes_PER_in,Yes_PER_on)

      use Global_Float_Type   
      use Global_Common
      
      implicit none
      real(kind=FT),intent(in)::Point(3),
     &                 Quad_P1(3),Quad_P2(3),Quad_P3(3),Quad_P4(3)
      real(kind=FT),intent(out):: Distance,PER(3)
      logical,intent(out):: Yes_PER_in,Yes_PER_on

      
      real(kind=FT) Tri123_Distance,Tri123_PER(3)
      real(kind=FT) Tri134_Distance,Tri134_PER(3)
      logical Tri123_Yes_PER_in,Tri134_Yes_PER_in
      logical Tri123_Yes_PER_on,Tri134_Yes_PER_on

      Yes_PER_in = .False.
      Yes_PER_on = .False.
      Distance   =  ZR
      PER        =  ZR
      
      call Tool_Dis_Point_to_3D_Tri
     &               (Point,Quad_P1,Quad_P2,Quad_P3,
     &                Tri123_Distance,Tri123_PER,
     &                Tri123_Yes_PER_in,Tri123_Yes_PER_on)
      call Tool_Dis_Point_to_3D_Tri
     &               (Point,Quad_P1,Quad_P3,Quad_P4,
     &                Tri134_Distance,Tri134_PER,
     &                Tri134_Yes_PER_in,Tri134_Yes_PER_on)
      
      if(abs(Tri123_Distance-Tri134_Distance) > Tol_10)then
          print *,'    Error:: non-planar space quadrilateral!'
          print *,'            in Tool_Dis_Point_to_3D_Quad.f'
          call Warning_Message('S',Keywords_Blank)
      endif
      
      Distance = Tri123_Distance
      PER      = Tri123_PER
      if(Tri123_Yes_PER_in .or. Tri134_Yes_PER_in)then
          Yes_PER_in =.True.
      endif
      
      return 
      end SUBROUTINE Tool_Dis_Point_to_3D_Quad                  
