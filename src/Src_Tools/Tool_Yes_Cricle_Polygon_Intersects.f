 
      subroutine Tool_Yes_Cricle_Polygon_Intersects(x0,y0,r,
     &                                    Clo_Supp_Domain,num_Domain_P,
     &                                    Yes_Intersect)

      use Global_Float_Type
      use Global_Elem_Area_Vol
      use Global_Crack
      
      implicit none
      integer,intent(in)::num_Domain_P
      real(kind=FT),intent(in)::x0,y0,r
      real(kind=FT),intent(in)::Clo_Supp_Domain(num_Domain_P,2)
      logical,intent(out)::Yes_Intersect
      integer i_Edge
      real(kind=FT) edge_p1(2),edge_p2(2)
      real(kind=FT) Tool_Function_2Point_Dis
      logical Yes_p1_in,Yes_p2_in
      
      Yes_Intersect = .False.

     
      do i_Edge = 1,num_Domain_P-1
          edge_p1 = Clo_Supp_Domain(i_Edge,1:2)
          edge_p2 = Clo_Supp_Domain(i_Edge+1,1:2)
          Yes_p1_in = .False.
          if(Tool_Function_2Point_Dis([x0,y0],edge_p1)<r)then
              Yes_p1_in  = .True.
          endif
          Yes_p2_in = .False.
          if(Tool_Function_2Point_Dis([x0,y0],edge_p2)<r)then
              Yes_p2_in  = .True.
          endif
          if((Yes_p1_in .and. (Yes_p2_in.eqv..False.)) .or.
     &       (Yes_p2_in .and. (Yes_p1_in.eqv..False.)))   then
              Yes_Intersect = .True.
              return
          endif 
      enddo

      
      return 
      end SUBROUTINE Tool_Yes_Cricle_Polygon_Intersects                         
