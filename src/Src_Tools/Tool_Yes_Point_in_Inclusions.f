 
      subroutine Tool_Yes_Point_in_Inclusions(x,y,
     &                        Yes_Gauss_in_Incl,c_Incl_Num)

      use Global_Float_Type
      use Global_Inclusion
      
      implicit none
      integer c_Incl_Num
      real(kind=FT) x,y
      logical Yes_Gauss_in_Incl
      
      integer i_Incl,n_Vertex
      real(kind=FT) c_Incl_x,c_Incl_y,c_Incl_r,c_Dis
      real(kind=FT) Tool_Function_2Point_Dis
      
      Yes_Gauss_in_Incl = .False.
      c_Incl_Num = 0
      
      if(num_Circ_Incl/=0)then
          do i_Incl = 1,num_Circ_Incl
              c_Incl_x = Circ_Inclu_Coor(i_Incl,1)
              c_Incl_y = Circ_Inclu_Coor(i_Incl,2)
              c_Incl_r = Circ_Inclu_Coor(i_Incl,3)
              c_Dis    = Tool_Function_2Point_Dis([x,y],
     &                              [c_Incl_x,c_Incl_y])
              if(c_Dis < c_Incl_r)then
                  Yes_Gauss_in_Incl = .True.
                  c_Incl_Num = i_Incl
                  exit
              endif
          enddo
      endif
      
      if(num_Poly_Incl/=0)then
          do i_Incl = 1,num_Poly_Incl
              n_Vertex = Poly_Inclu_Edges_Num(i_Incl)
              call Tool_Yes_In_Poly(x,y,
     &                        Poly_Incl_Coor_x_Cl(i_Incl,1:n_Vertex+1),
     &                        Poly_Incl_Coor_y_Cl(i_Incl,1:n_Vertex+1),
     &                        n_Vertex+1,
     &                        Yes_Gauss_in_Incl)
     
              if(Yes_Gauss_in_Incl)then
                  c_Incl_Num = i_Incl
                  exit
              endif
          enddo
      endif
  
      return 
      end SUBROUTINE Tool_Yes_Point_in_Inclusions                          
