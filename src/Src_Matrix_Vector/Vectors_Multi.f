 
      SUBROUTINE Vectors_Multi(Vector1,n1,Vector2,n2,Matrix_OUT)   

      use Global_Float_Type
      implicit none
      integer,intent(in):: n1,n2
      real(kind=FT),intent(in)::Vector1(n1),Vector2(n2)
      real(kind=FT),intent(out)::Matrix_OUT(n1,n2)
      
      integer i,j
      
      do i=1,n1
          do j=1,n2
              Matrix_OUT(i,j) = Vector1(i)*Vector2(j)
          end do
      end do
      
      return
      END SUBROUTINE Vectors_Multi
    


