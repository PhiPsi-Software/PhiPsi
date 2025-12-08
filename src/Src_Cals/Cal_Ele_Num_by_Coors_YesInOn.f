 
      subroutine Cal_Ele_Num_by_Coors_YesInOn(x,y,OUT_Elem,
     &                                        Yes_In,Yes_On)
      use Global_Float_Type
      use Global_Model
      
      implicit none
      
      integer i
      real(kind=FT) x,y
      integer,intent(out)::OUT_Elem
      logical,intent(out):: Yes_In, Yes_On
      
      real(kind=FT) c_x_max,c_x_min,c_y_max,c_y_min
      integer c_count,Potent_Elem(5)
      real(kind=FT) xpol(5),ypol(5)
      
      c_count = 0
      Yes_In = .False.
      Yes_On = .False.
      
      do i=1,Num_Elem
          c_x_max = maxval(G_X_NODES(:,i))
          c_x_min = minval(G_X_NODES(:,i))
          c_y_max = maxval(G_Y_NODES(:,i))
          c_y_min = minval(G_Y_NODES(:,i))
          if   ((x.le.c_x_max).and.(x.ge.c_x_min).
     &     and. (y.le.c_y_max).and.(y.ge.c_y_min))then
              c_count = c_count +1
              Potent_Elem(c_count) = i
          end if
      end do 
      
      do i =1,c_count
          xpol(1:4) = G_X_NODES(1:4,Potent_Elem(i))
          ypol(1:4) = G_Y_NODES(1:4,Potent_Elem(i))
          xpol(5) = G_X_NODES(1,Potent_Elem(i))
          ypol(5) = G_Y_NODES(1,Potent_Elem(i))
          call Tool_Yes_In_Poly(x,y,xpol,ypol,5,Yes_In)
          call Tool_Yes_On_Poly(x,y,xpol,ypol,5,Yes_On)
          if ((Yes_In.eqv..True.).or.(Yes_On.eqv..True.))then
              OUT_Elem = Potent_Elem(i)
              exit
          end if
      end do
      
      return 
      end SUBROUTINE Cal_Ele_Num_by_Coors_YesInOn                         
