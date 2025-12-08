 
      subroutine Tool_Yes_Crack_Polygon_Intersects(c_Crack,
     &                                    Clo_Supp_Domain,num_Domain_P,
     &                                    Yes_Intersect)
      use Global_Float_Type
      use Global_Elem_Area_Vol
      use Global_Crack
      use Global_Crack_Common
      
      implicit none
      integer,intent(in)::c_Crack,num_Domain_P
      real(kind=FT),intent(in)::Clo_Supp_Domain(num_Domain_P,2)
      logical,intent(out)::Yes_Intersect
      integer i_S,i_Edge
      real(kind=FT) crack_p1(2),crack_p2(2),edge_p1(2),edge_p2(2)
      real(kind=FT) c_X,c_Y
      logical c_Yes_Cross
      
      Yes_Intersect = .False.
      do i_S = 1,Each_Cr_Poi_Num(c_Crack)-1
          crack_p1 = Crack_Coor(c_Crack,i_S,  1:2)
          crack_p2 = Crack_Coor(c_Crack,i_S+1,1:2)
          do i_Edge = 1,num_Domain_P-1
              edge_p1 = Clo_Supp_Domain(i_Edge,1:2)
              edge_p2 = Clo_Supp_Domain(i_Edge+1,1:2)
              call Tool_Intersection(crack_p1,crack_p2,
     &                                   edge_p1,edge_p2 ,
     &                                   c_X,c_Y,c_Yes_Cross)  
              if(c_Yes_Cross)then
                  Yes_Intersect = .True.
                  return
              endif 
          enddo
      enddo
      
      return 
      end SUBROUTINE Tool_Yes_Crack_Polygon_Intersects                          
