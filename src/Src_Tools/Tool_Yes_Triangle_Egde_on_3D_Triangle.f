 
      subroutine Tool_Yes_Triangle_Egde_on_3D_Triangle(Tri_1,Tri_2,
     &                            Yes_on,Inter_Point,Side_Tri_Number,
     &                            Tri_Edge_Num)   
     

      use Global_Float_Type
      use Global_Common
      implicit none
      real(kind=FT),intent(in)::Tri_1(3,3),Tri_2(3,3)
      real(kind=FT),intent(out)::Inter_Point(3)
      logical,intent(out)::Yes_on              
      integer,intent(out)::Side_Tri_Number
      integer,intent(out)::Tri_Edge_Num
      logical Logical_Yes_x,Logical_Yes_y,Logical_Yes_z
      real(kind=FT) max_x_Tri_1,max_y_Tri_1,max_z_Tri_1
      real(kind=FT) min_x_Tri_1,min_y_Tri_1,min_z_Tri_1
      real(kind=FT) max_x_Tri_2,max_y_Tri_2,max_z_Tri_2
      real(kind=FT) min_x_Tri_2,min_y_Tri_2,min_z_Tri_2 
      integer Tri_Edges(3,3)
      integer i_Edge
      real(kind=FT) P1(3),P2(3),tem_Point(3)
      logical tem_Yes_on
      
      Yes_on  = .False.
      Inter_Point(1:3) = -TEN_15 
      Side_Tri_Number  = -999
      max_x_Tri_1 = maxval(Tri_1(1:3,1))
      max_y_Tri_1 = maxval(Tri_1(1:3,2))
      max_z_Tri_1 = maxval(Tri_1(1:3,3))
      min_x_Tri_1 = minval(Tri_1(1:3,1))
      min_y_Tri_1 = minval(Tri_1(1:3,2))
      min_z_Tri_1 = minval(Tri_1(1:3,3))      
      max_x_Tri_2 = maxval(Tri_2(1:3,1))
      max_y_Tri_2 = maxval(Tri_2(1:3,2))
      max_z_Tri_2 = maxval(Tri_2(1:3,3))
      min_x_Tri_2 = minval(Tri_2(1:3,1))
      min_y_Tri_2 = minval(Tri_2(1:3,2))
      min_z_Tri_2 = minval(Tri_2(1:3,3))         
      call Tool_Yes_Two_Ranges_Overlapped_Double(
     &                            [min_x_Tri_1,max_x_Tri_1],
     &                            [min_x_Tri_2,max_x_Tri_2],
     &                            Logical_Yes_x) 
      if(Logical_Yes_x .eqv. .False.)then
          return
      endif
      call Tool_Yes_Two_Ranges_Overlapped_Double(
     &                            [min_y_Tri_1,max_y_Tri_1],
     &                            [min_y_Tri_2,max_y_Tri_2],
     &                            Logical_Yes_y) 
      if(Logical_Yes_y .eqv. .False.)then
          return
      endif
      call Tool_Yes_Two_Ranges_Overlapped_Double(
     &                            [min_z_Tri_1,max_z_Tri_1],
     &                            [min_z_Tri_2,max_z_Tri_2],
     &                            Logical_Yes_z) 
      if(Logical_Yes_z .eqv. .False.)then
          return
      endif      
      
      Tri_Edges(1,1:2) = [1,2] 
      Tri_Edges(2,1:2) = [2,3] 
      Tri_Edges(3,1:2) = [3,1]          
      do i_Edge = 1,3
          P1(1:3) = Tri_1(Tri_Edges(i_Edge,1),1:3)
          P2(1:3) = Tri_1(Tri_Edges(i_Edge,2),1:3)
          call Tool_Yes_Line_on_3D_Triangle(P1,P2,
     &           Tri_2(1,1:3),Tri_2(2,1:3),Tri_2(3,1:3),
     &           tem_Yes_on,tem_Point)
          if(tem_Yes_on) then
              Yes_on  = .True.
              Inter_Point = tem_Point
              Side_Tri_Number = 1
              Tri_Edge_Num  = i_Edge 
              return
          endif
      enddo      

      Tri_Edges(1,1:2) = [1,2] 
      Tri_Edges(2,1:2) = [2,3] 
      Tri_Edges(3,1:2) = [3,1]          
      do i_Edge = 1,3
          P1(1:3) = Tri_2(Tri_Edges(i_Edge,1),1:3)
          P2(1:3) = Tri_2(Tri_Edges(i_Edge,2),1:3)
          call Tool_Yes_Line_on_3D_Triangle(P1,P2,
     &               Tri_1(1,1:3),Tri_1(2,1:3),Tri_1(3,1:3),
     &               tem_Yes_on,tem_Point)
          if(tem_Yes_on) then
              Yes_on  = .True.
              Inter_Point = tem_Point
              Side_Tri_Number = 2
              Tri_Edge_Num  = i_Edge 
              return
          endif
      enddo          
      
      return 
      end SUBROUTINE Tool_Yes_Triangle_Egde_on_3D_Triangle          
