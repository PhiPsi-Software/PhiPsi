 
      subroutine Tool_Fit_3D_Points_to_Ellipse(num_Point,
     &                                    In_Points,Out_Points)

      use Global_Float_Type     
      
      implicit none
      integer,intent(in)::num_Point
      real(kind=FT),intent(in)::In_Points(num_Point,3)
      real(kind=FT),intent(out)::Out_Points(num_Point,3)
      real(kind=FT) Center(3),n_Vector(3)
      real(kind=FT) Tool_Function_2Point_Dis_3D
      real(kind=FT) Nearest_Point(3)
      real(kind=FT) all_Dis(num_Point),vector_U(3),vector_V(3)
      real(kind=FT) a,b
      integer Farthest_Point
      
      integer i_P
      
      
      Center(1) = sum(In_Points(1:num_Point,1))/dble(num_Point)
      Center(2) = sum(In_Points(1:num_Point,2))/dble(num_Point)
      Center(3) = sum(In_Points(1:num_Point,3))/dble(num_Point)
      
      all_Dis(1:num_Point) = ZR
      do i_P = 1,num_Point
          all_Dis(i_P) = Tool_Function_2Point_Dis_3D(In_Points(i_P,1:3),
     &                                               Center)
      enddo

      a = maxval(all_Dis)
      b = minval(all_Dis)
      
      Farthest_Point = maxloc(all_Dis,1)
      vector_U = In_Points(Farthest_Point,1:3)-Center(1:3)
      call Vector_Normalize(3,vector_U)   
      
      call Tool_Get_Normal_Vector_of_Points_3D(
     &                 In_Points,num_Point,n_Vector)
      
      call Vector_Normalize(3,n_Vector)   
      call Vector_Cross_Product_3(n_Vector,vector_U,vector_V)   
      call Vector_Normalize(3,vector_V)   
      
      
      do i_P = 1,num_Point
          call Tool_Nearest_Point_from_Point_to_Ellipse_3D(
     &                 In_Points(i_P,1:3),Center,a,b,vector_U,vector_V,
     &                 Nearest_Point)     
          Out_Points(i_P,1:3) = Nearest_Point(1:3)
      enddo
      
      end SUBROUTINE Tool_Fit_3D_Points_to_Ellipse             
