 
      subroutine Tool_Fit_3D_Points_to_Circle(num_Point,
     &                                    In_Points,Out_Points)
      use Global_Float_Type     
      implicit none
      integer,intent(in)::num_Point
      real(kind=FT),intent(in)::In_Points(num_Point,3)
      real(kind=FT),intent(out)::Out_Points(num_Point,3)
      real(kind=FT) Center(3),n_Vector(3)
      real(kind=FT) Ave_Radius,Tool_Function_2Point_Dis_3D
      real(kind=FT) Nearest_Point(3)
      integer i_P
       real(kind=FT) sum_1,sum_2,sum_3
      
      
      sum_1 = ZR
      sum_2 = ZR
      sum_3 = ZR
      do i_P = 1,num_Point
          sum_1 = sum_1 + In_Points(i_P,1)
          sum_2 = sum_2 + In_Points(i_P,2)
          sum_3 = sum_3 + In_Points(i_P,3)
      enddo
      Center(1) = sum_1/dble(num_Point)
      Center(2) = sum_2/dble(num_Point)
      Center(3) = sum_3/dble(num_Point)
      
      
      Ave_Radius = ZR
      do i_P = 1,num_Point
          Ave_Radius =  Ave_Radius + 
     &      Tool_Function_2Point_Dis_3D(In_Points(i_P,1:3),Center)
      enddo
      Ave_Radius = Ave_Radius/dble(num_Point)
      

      
      call Tool_Get_Normal_Vector_of_Points_3D(
     &                 In_Points,num_Point,n_Vector)
     
      
      Out_Points = ZR
      do i_P = 1,num_Point
          call Tool_Nearest_Point_from_Point_to_Cricle_3D(
     &                 In_Points(i_P,1:3),Center,Ave_Radius,
     &                 n_Vector,Nearest_Point)     
          Out_Points(i_P,1:3) = Nearest_Point(1:3)
      enddo

      
      
      return 
      end SUBROUTINE Tool_Fit_3D_Points_to_Circle                 
