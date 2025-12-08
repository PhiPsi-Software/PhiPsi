 
      subroutine Tool_Get_Farthest_Two_Points_3D(num_Point,Points,
     &                                       Point1_Num,Point2_Num) 
 
      use Global_Float_Type
      use Global_Common
      
      implicit none
      integer,intent(in)::num_Point
      real(kind=FT),intent(in)::Points(num_Point,3)       
      integer,intent(out)::Point1_Num,Point2_Num 
      real(kind=FT) Dis_Matrix(num_Point,num_Point)
      real(kind=FT) Tool_Function_2Point_Dis_3D
      integer i_P,j_P
      
      if(num_Point<=1)then
        print *, '    Error :: num_Point<=1!'
        print *, '             in Tool_Get_Farthest_Points_3D.f!'
        call Warning_Message('S',Keywords_Blank)      
      endif
      Dis_Matrix = -TEN_15
      
      
      do i_P=1,num_Point
          do j_P=i_P+1,num_Point
              Dis_Matrix(i_P,j_P) = Tool_Function_2Point_Dis_3D(
     &                              Points(i_P,1:3),
     &                              Points(j_P,1:3))
              print *,i_P,j_P,Dis_Matrix(i_P,j_P) 
          enddo
      enddo
      
      call Matrix_Max_Location_Dou(num_Point,num_Point,Dis_Matrix,
     &                             Point1_Num,Point2_Num)
      if(Point1_Num== Point2_Num)then
        print *, '    Error :: Point1_Num = Point2_Num!'
        print *, '             in Tool_Get_Farthest_Points_3D.f!'
        call Warning_Message('S',Keywords_Blank)  
      endif
      
      return 
      end SUBROUTINE Tool_Get_Farthest_Two_Points_3D                    
