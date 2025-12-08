 
      subroutine Tool_Dis_Point_to_Poly(x,y,
     &                                  xpol,ypol,npol,
     &                                  Dis)

      use Global_Float_Type     
      implicit none
      integer,intent(in)::npol
      real(kind=FT),intent(in)::x,y,xpol(npol),ypol(npol)
      real(kind=FT),intent(out)::Dis
      logical Yes_in_Poly,Yes_on_Poly
      real(kind=FT) Dis_to_Edges(npol-1),Line_Edge(2,2)
      integer i_Edge
      call Tool_Yes_In_Poly(x,y,xpol,ypol,npol,Yes_in_Poly)
      
      call Tool_Yes_On_Poly(x,y,xpol,ypol,npol,Yes_ON_Poly)
      if(Yes_ON_Poly)then
          Dis = ZR
          return
      endif
      
      Dis_to_Edges(1:npol-1) = ZR
      do i_Edge = 1,npol-1
          Line_Edge(1,1) = xpol(i_Edge)
          Line_Edge(1,2) = ypol(i_Edge)
          Line_Edge(2,1) = xpol(i_Edge+1)
          Line_Edge(2,2) = ypol(i_Edge+1)
          call Tool_Dis_Point_to_Seg(Line_Edge(1,1:2),
     &                               Line_Edge(2,1:2),
     &                               [x,y],
     &                               Dis_to_Edges(i_Edge))
     
      enddo
      if(Yes_in_Poly)then
          Dis = -minval(Dis_to_Edges(1:npol-1))
      else
          Dis =  minval(Dis_to_Edges(1:npol-1))
      endif
       
      return 
      end SUBROUTINE Tool_Dis_Point_to_Poly                       
