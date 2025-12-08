 
      subroutine Tool_Nearest_Point_from_Point_to_Cricle_3D(
     &                 In_Point,Center,r,n,Out_Point)

      use Global_Float_Type
      
      implicit none
      real(kind=FT),intent(in)::In_Point(3),Center(3),r,n(3)
      real(kind=FT),intent(out)::Out_Point(3)
      real(kind=FT) delta(3),tem(3),nor_tem
      
      delta = In_Point-Center
      tem = delta - dot_product(n,delta)*n
      call Vector_Norm2(3,tem,nor_tem)   
      Out_Point = Center + r*tem/nor_tem
      
      
      return 
      end SUBROUTINE Tool_Nearest_Point_from_Point_to_Cricle_3D                         
