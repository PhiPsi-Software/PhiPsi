 
      subroutine Cal_Equal_Division_Points(A,B,
     &                            Num_Div,Yes_Include_Endpoint,
     &                            Div_Points,Num_Div_Points)
      use Global_Float_Type
      use Global_Crack
      
      implicit none
      
      real(kind=FT),intent(in)::A(2),B(2)
      integer,intent(in)::Num_Div
      logical,intent(in)::Yes_Include_Endpoint
      real(kind=FT),intent(out)::Div_Points(Max_Num_Seg_CalP,2)
      integer,intent(out)::Num_Div_Points
      
      integer i
      real(kind=FT) a_x,a_y,b_x,b_y
      
      a_x = A(1)
      a_y = A(2)
      b_x = B(1)
      b_y = B(2)
      
      Num_Div_Points = 0
      if (Yes_Include_Endpoint.eqv..False.) then
          do i = 1,Num_Div-1
              Div_Points(i,1) = (i*b_x+(Num_Div-i)*a_x)/Num_Div
              Div_Points(i,2) = (i*b_y+(Num_Div-i)*a_y)/Num_Div
          end do
          Num_Div_Points = Num_Div-1
      else
          Div_Points(1,:)= A
          Num_Div_Points = 1
          do i = 1,Num_Div-1
              Num_Div_Points  = Num_Div_Points + 1
              Div_Points(Num_Div_Points,1)=
     &                            (i*b_x+(Num_Div-i)*a_x)/Num_Div
              Div_Points(Num_Div_Points,2)=
     &                            (i*b_y+(Num_Div-i)*a_y)/Num_Div
          end do 
          Num_Div_Points  = Num_Div_Points + 1
          Div_Points(Num_Div_Points,:) = B
      end if
      
      return 
      end SUBROUTINE Cal_Equal_Division_Points               
