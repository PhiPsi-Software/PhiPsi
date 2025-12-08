 
      subroutine Tool_abc_of_Parabola(Point_1,Point_2,Point_3,a,b,c)
      use Global_Float_Type      
      implicit none
      real(kind=FT),intent(in):: Point_1(2),Point_2(2),Point_3(2)
      real(kind=FT),intent(out):: a,b,c
      
      real(kind=FT) D(3),K(3,3),inv_K(3,3),F(3),x1,x2,x3,y1,y2,y3
      
      x1 = Point_1(1); y1 = Point_1(2)
      x2 = Point_2(1); y2 = Point_2(2)
      x3 = Point_3(1); y3 = Point_3(2)
      
      K(1,1:3) = [x1**2, x1, ONE]
      K(2,1:3) = [x2**2, x2, ONE]
      K(3,1:3) = [x3**2, x3, ONE]

      
      F(1:3)   = [y1,y2,y3]
      
      call Matrix_Inverse(K,inv_K,3) 
      
      D = MATMUL(inv_K,F) 
      
      a = D(1)
      b = D(2)
      c = D(3)

      
      return 
      end SUBROUTINE Tool_abc_of_Parabola                          
