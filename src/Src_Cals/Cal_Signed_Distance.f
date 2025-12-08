 
      subroutine Cal_Signed_Distance(Line_AB,Point_C,S_Distance)
      use Global_Float_Type      
      implicit none
      real(kind=FT),intent(in)::Line_AB(2,2),Point_C(2)
      real(kind=FT),intent(out)::S_Distance
      
      real(kind=FT) tem_1(2,2), tem_2(2)
      real(kind=FT) tem_Det,tem_Norm
      
      
      tem_1(1,:) = Line_AB(2,:)-  Line_AB(1,:)
      tem_1(2,:) = Point_C     -  Line_AB(1,:)
      
      
      tem_2(:)   = Line_AB(2,:)-Line_AB(1,:)
      
      
      tem_Det = tem_1(1,1)*tem_1(2,2) - tem_1(2,1)*tem_1(1,2)
      
      call Vector_Norm2(2,tem_2,tem_Norm) 
      
      S_Distance = tem_Det / tem_Norm
      
      return 
      end SUBROUTINE Cal_Signed_Distance                          
