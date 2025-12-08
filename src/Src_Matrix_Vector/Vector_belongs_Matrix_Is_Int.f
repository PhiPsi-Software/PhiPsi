 
      subroutine Vector_belongs_Matrix_Is_Int(m,n,Matrix,Vector,
     &                                              Location,Yes)   
      use Global_Float_Type     
      implicit none
      
      integer,intent(in)::m,n
      integer,intent(in)::Matrix(m,n)
      integer,intent(in)::Vector(n)
      integer,intent(out)::Location
      logical,intent(out)::Yes
      integer i                 
      integer Vector_tem(n)
      
      Location = 0
      Yes = .False.
      do i=1,m
          Vector_tem(1:n) = Matrix(i,1:n)
          if (ALL(Vector_tem .EQ. Vector)) then
              Yes = .True.
              Location = i
              exit
          end if
          
      end do
      
      return 
      end subroutine Vector_belongs_Matrix_Is_Int
    


