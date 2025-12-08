 
      subroutine Vector_belongs_Matrix_Is_Dou(m,n,Matrix,Vector,
     &                                        Location,Yes)   
      use Global_Float_Type
      implicit none
      integer,intent(in)::m,n
      real(kind=FT),intent(in)::Matrix(m,n)
      real(kind=FT),intent(in)::Vector(n)
      integer,intent(out)::Location
      logical,intent(out)::Yes
      
      integer i
      Yes = .False.

      do i=1,m
          if ((abs(MaxVal(Matrix(i,1:n)-Vector)) < Tol_11).and.
     &        (abs(MinVal(Matrix(i,1:n)-Vector)) < Tol_11)) then
              Yes = .True.
              Location = i
              exit
          end if
      end do
      
      return 
      end subroutine Vector_belongs_Matrix_Is_Dou
    


